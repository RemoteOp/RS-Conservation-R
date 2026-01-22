# ============================================================
# GLDAS NOAH Root Zone Soil Moisture (monthly) processing — AOI
# FULL FEATURE SCRIPT + master-grid fix + drought clim12 write fix
#
# Scientific setup:
# - Legacy climatology:   2000–2012
# - Recent climatology:   2013–2024
# - Drought thresholds:   2000–2024 (full record percentiles)
# - Drought frequency:    2013–2024 (recent) vs thresholds (full record)
#
# Engineering:
# - Handles mixed-geometry GeoTIFFs by aligning to a MASTER grid
# - MASTER grid is chosen as the MOST COMMON geometry among threshold files
#   (prevents "half extent" when the first file is a cropped/odd raster)
# - Ensures droughtfreq clim12 files (P10/P20/P30) are always written/updated
# ============================================================

library(terra)
library(stringr)

terraOptions(progress = 1)

# ----------------------------
# 1) USER SETTINGS
# ----------------------------

# IMPORTANT: Use forward slashes in Windows paths (R-safe)
DATA_DIR <- "C:/GLDAS_2000_2024/data"
AOI_SHP  <- "..../AOI_reprojected.shp"
OUT_R    <- "..../NOAH_analysis/R"

FILE_PATTERN <- "\\.SUB\\.tif$"

# --- PERIODS (best-practice for NOAH monthly 0.25°) ---
LEGACY_YEARS    <- 2000:2012        # non-overlap "legacy"
RECENT_YEARS    <- 2013:2024        # non-overlap "recent"
THRESHOLD_YEARS <- 2000:2024        # stable drought percentiles (full record)

# For filenames/tags
LEGACY_TAG <- sprintf("%d_%d", min(LEGACY_YEARS),    max(LEGACY_YEARS))
RECENT_TAG <- sprintf("%d_%d", min(RECENT_YEARS),    max(RECENT_YEARS))
THRESH_TAG <- sprintf("%d_%d", min(THRESHOLD_YEARS), max(THRESHOLD_YEARS))

SEASON_MONTHS <- c(1, 4, 7, 10)

# Drought percentiles (P10 severe, P20 main, etc.)
DROUGHT_PERCENTILES <- c(10, 20, 30)

EXPORT_CLIMATOLOGY      <- TRUE # set to FALSE to exclude
EXPORT_SEASONAL_TIF     <- TRUE
EXPORT_ANOMALY          <- TRUE
EXPORT_DROUGHT_FREQ     <- TRUE

EXPORT_DROUGHT_MONTHLY_TIF <- TRUE
EXPORT_DROUGHT_SEASONAL    <- TRUE
EXPORT_DROUGHT_ALLMONTHS   <- TRUE

EXPORT_DROUGHT_CSV      <- TRUE
EXPORT_DROUGHT_CSV_LONG <- TRUE

# Area thresholds for drought frequency rasters (fraction of months)
DROUGHT_AREA_THRESHOLDS <- c(0.2, 0.5, 0.8)

GDAL_OPTS <- c("COMPRESS=DEFLATE", "TILED=YES")

# ----------------------------
# 2) OUTPUT PATH CHECKS
# ----------------------------

cat("\n=== OUTPUT PATH CHECK ===\n")
cat("OUT_R:\n", OUT_R, "\n")

dir.create(OUT_R, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_R)) stop("OUT_R does not exist and could not be created: ", OUT_R)

proof_file <- file.path(OUT_R, "___script_started.txt")
writeLines(paste("Started at:", Sys.time()), proof_file)
cat("Wrote proof file:", proof_file, "\n")

# ----------------------------
# 3) HELPER FUNCTIONS
# ----------------------------

# Extract YYYYMM token from filenames like "...Ayyyymm..."
extract_yyyymm <- function(x) {
  m <- str_match(basename(x), "A(\\d{6})")[,2]
  if (is.na(m)) return(NA_character_)
  m
}

# Convert YYYYMM token to year integer.
yyyymm_to_year  <- function(yyyymm) as.integer(substr(yyyymm, 1, 4))
# Convert YYYYMM token to month integer.
yyyymm_to_month <- function(yyyymm) as.integer(substr(yyyymm, 5, 6))

month_names_short <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# Write a SpatRaster to disk with robust overwrite checks.
write_tif <- function(r, filename) {
  cat("Writing:", filename, "\n")
  if (!inherits(r, "SpatRaster")) stop("Object is not a SpatRaster: ", filename)
  if (nlyr(r) < 1) stop("Raster has 0 layers (nothing to write): ", filename)
  
  # Force overwrite even on Windows/GDAL oddities: remove first if exists
  if (file.exists(filename)) {
    ok_rm <- tryCatch({ unlink(filename); TRUE }, error = function(e) FALSE)
    if (!ok_rm && file.exists(filename)) {
      stop("Could not remove existing file before overwrite: ", filename)
    }
  }
  
  tryCatch({
    writeRaster(r, filename, overwrite = TRUE, gdal = GDAL_OPTS)
  }, error = function(e) {
    stop("writeRaster FAILED for: ", filename, "\nReason: ", conditionMessage(e))
  })
  
  if (!file.exists(filename)) stop("writeRaster reported success but file does not exist: ", filename)
}

# Build a 12-band monthly climatology (mean per calendar month).
monthly_climatology <- function(r_stack, months_vec) {
  stopifnot(nlyr(r_stack) == length(months_vec))
  out_list <- vector("list", 12)
  
  for (m in 1:12) {
    idx <- which(months_vec == m)
    if (length(idx) == 0) {
      warning(sprintf("No layers found for month %02d", m))
      next
    }
    r_m <- mean(r_stack[[idx]], na.rm = TRUE)
    names(r_m) <- month_names_short[m]
    out_list[[m]] <- r_m
  }
  
  out_list <- out_list[!vapply(out_list, is.null, logical(1))]
  rast(out_list)
}

# Select layers by month number, using layer names to avoid index shifts.
subset_month_layers <- function(r12, months) {
  keep_names <- month_names_short[months]
  keep_names <- keep_names[keep_names %in% names(r12)]
  r12[[keep_names]]
}

# Compute per-pixel quantile (used for drought thresholds).
pixel_quantile <- function(x, prob) {
  quantile(x, probs = prob, na.rm = TRUE, type = 7)
}

# ---- Robust global stats WITHOUT terra::global(fun="...") ----
# Summarize raster values (mean/median/sd/min/max) safely.
layer_stats <- function(r) {
  v <- values(r, mat = FALSE)
  v <- v[is.finite(v)]
  if (length(v) == 0) {
    return(list(mean=NA_real_, median=NA_real_, sd=NA_real_, min=NA_real_, max=NA_real_))
  }
  list(
    mean   = mean(v),
    median = stats::median(v),
    sd     = stats::sd(v),
    min    = min(v),
    max    = max(v)
  )
}

# ---- Robust area summaries WITHOUT terra::global("sum") ----
# Summarize area above thresholds using cell size in m^2.
area_above_thresholds <- function(r, thresholds) {
  a <- cellSize(r, unit = "m")  # m^2 raster
  a_vals <- values(a, mat = FALSE)
  a_vals <- a_vals[is.finite(a_vals)]
  total_area <- sum(a_vals)  # m^2
  
  out <- list(total_area_m2 = total_area)
  
  r_vals <- values(r, mat = FALSE)
  
  for (t in thresholds) {
    keep <- is.finite(r_vals) & (r_vals >= t)
    area_t <- sum(values(a, mat = FALSE)[keep], na.rm = TRUE)
    out[[sprintf("area_ge_%.2f_m2", t)]] <- area_t
    out[[sprintf("pct_area_ge_%.2f", t)]] <- if (is.finite(total_area) && total_area > 0) (area_t / total_area) * 100 else NA_real_
  }
  out
}

# ----------------------------
# 4) GEOMETRY HANDLING (mode geometry master template)
# ----------------------------

# Build a geometry signature (CRS, extent, resolution, dimensions).
geom_signature <- function(r) {
  e <- ext(r); rs <- res(r)
  paste0(
    "crs=", crs(r),
    "|ext=", paste(as.vector(e), collapse = ","),
    "|res=", paste(as.vector(rs), collapse = ","),
    "|dim=", nrow(r), "x", ncol(r)
  )
}

# Choose the most common geometry among candidate rasters as MASTER grid.
choose_master_template <- function(file_vec, max_scan = 120) {
  n <- length(file_vec)
  if (n == 0) stop("choose_master_template(): empty file list")
  
  scan_n <- min(n, max_scan)
  scan_files <- file_vec[seq_len(scan_n)]
  
  sigs <- character(scan_n)
  for (i in seq_len(scan_n)) {
    r <- rast(scan_files[i])
    sigs[i] <- geom_signature(r)
  }
  
  tab <- sort(table(sigs), decreasing = TRUE)
  top_sig <- names(tab)[1]
  
  idx <- which(sigs == top_sig)[1]
  master_file <- scan_files[idx]
  master_r <- rast(master_file)
  
  cat("\n=== MASTER GRID SELECTION ===\n")
  cat("Scanned files:", scan_n, "of", n, "\n")
  cat("Most common geometry count:", as.integer(tab[1]), "\n")
  cat("Master template file:", master_file, "\n")
  print(master_r)
  
  master_r
}

# Check if two rasters share identical grid geometry.
geom_ok <- function(a, b) {
  same.crs(a, b) &&
    isTRUE(all.equal(as.vector(ext(a)), as.vector(ext(b)), tolerance = 1e-9)) &&
    isTRUE(all.equal(as.vector(res(a)), as.vector(res(b)), tolerance = 1e-12)) &&
    nrow(a) == nrow(b) && ncol(a) == ncol(b)
}

# Reproject/resample a raster to match the MASTER grid.
align_to_template <- function(r, template, method = "bilinear") {
  if (!same.crs(template, r)) {
    r <- project(r, template, method = method)
  }
  if (!geom_ok(template, r)) {
    r <- resample(r, template, method = method)
  }
  r
}

# Load a stack of rasters aligned to the MASTER grid.
stack_aligned_to_template <- function(file_vec, template, method = "bilinear", label = "stack") {
  cat("\nBuilding stack aligned to MASTER grid:", label, "\n")
  out <- vector("list", length(file_vec))
  
  for (i in seq_along(file_vec)) {
    f <- file_vec[i]
    r <- rast(f)
    
    if (!geom_ok(template, r)) {
      cat("Aligning:", basename(f), "\n")
      r <- align_to_template(r, template, method = method)
    }
    out[[i]] <- r
  }
  rast(out)
}

# ----------------------------
# 5) Load files and parse dates
# ----------------------------

cat("\n=== INPUT CHECK ===\n")
cat("DATA_DIR:", DATA_DIR, "\n")
if (!dir.exists(DATA_DIR)) stop("DATA_DIR does not exist: ", DATA_DIR)

files <- list.files(DATA_DIR, pattern = FILE_PATTERN, full.names = TRUE)
cat("Found .SUB.tif files:", length(files), "\n")
if (length(files) == 0) stop("No .SUB.tif files found. Check DATA_DIR and FILE_PATTERN.")

yyyymm <- vapply(files, extract_yyyymm, FUN.VALUE = character(1))
ok <- !is.na(yyyymm)
files <- files[ok]
yyyymm <- yyyymm[ok]
cat("Files with Ayyyymm parsed:", length(files), "\n")
if (length(files) == 0) stop("No files with Ayyyymm pattern found in filenames.")

years  <- vapply(yyyymm, yyyymm_to_year,  FUN.VALUE = integer(1))
months <- vapply(yyyymm, yyyymm_to_month, FUN.VALUE = integer(1))

ord <- order(years, months)
files <- files[ord]
years <- years[ord]
months <- months[ord]

cat("First file:", files[1], "\n")
cat("Last  file:", files[length(files)], "\n")
cat("Year range:", min(years), "-", max(years), "\n")

# Clean: ignore files beyond THRESHOLD_YEARS max
keep <- years <= max(THRESHOLD_YEARS)
files  <- files[keep]
years  <- years[keep]
months <- months[keep]

legacy_idx    <- which(years %in% LEGACY_YEARS)
recent_idx    <- which(years %in% RECENT_YEARS)
threshold_idx <- which(years %in% THRESHOLD_YEARS)

cat("Legacy layers:   ", length(legacy_idx), "\n")
cat("Recent layers:   ", length(recent_idx), "\n")
cat("Threshold layers:", length(threshold_idx), "\n")

if (length(legacy_idx)    == 0) stop("No legacy files found for LEGACY_YEARS.")
if (length(recent_idx)    == 0) stop("No recent files found for RECENT_YEARS.")
if (length(threshold_idx) == 0) stop("No threshold files found for THRESHOLD_YEARS.")

legacy_files    <- files[legacy_idx]
recent_files    <- files[recent_idx]
threshold_files <- files[threshold_idx]

months_legacy    <- months[legacy_idx]
months_recent    <- months[recent_idx]
months_threshold <- months[threshold_idx]

# ----------------------------
# 6) Stack and clip to AOI (single master grid)
# ----------------------------

cat("\n=== AOI / STACK CHECK ===\n")
if (!file.exists(AOI_SHP)) stop("AOI shapefile not found: ", AOI_SHP)

aoi <- vect(AOI_SHP)

# Choose MASTER_TEMPLATE as the most common geometry in threshold files
MASTER_TEMPLATE <- choose_master_template(threshold_files, max_scan = 120)

# Build stacks aligned to that master grid
r_threshold <- stack_aligned_to_template(threshold_files, MASTER_TEMPLATE, method = "bilinear", label = "threshold")
r_legacy    <- stack_aligned_to_template(legacy_files,    MASTER_TEMPLATE, method = "bilinear", label = "legacy")
r_recent    <- stack_aligned_to_template(recent_files,    MASTER_TEMPLATE, method = "bilinear", label = "recent")

cat("Threshold raster layers:", nlyr(r_threshold), "\n")
cat("Legacy raster layers:   ", nlyr(r_legacy), "\n")
cat("Recent raster layers:   ", nlyr(r_recent), "\n")

# Project AOI to master CRS, then crop/mask all stacks
aoi_proj <- project(aoi, crs(MASTER_TEMPLATE))

r_threshold_clip <- mask(crop(r_threshold, aoi_proj), aoi_proj)
r_legacy_clip    <- mask(crop(r_legacy,    aoi_proj), aoi_proj)
r_recent_clip    <- mask(crop(r_recent,    aoi_proj), aoi_proj)

cat("Clipped threshold layers:", nlyr(r_threshold_clip), "\n")
cat("Clipped legacy layers:   ", nlyr(r_legacy_clip), "\n")
cat("Clipped recent layers:   ", nlyr(r_recent_clip), "\n")

# ----------------------------
# 7) Climatologies + anomalies
# ----------------------------

cat("\n=== CLIMATOLOGY / ANOMALY ===\n")

clim_legacy_12 <- monthly_climatology(r_legacy_clip, months_legacy)
clim_recent_12 <- monthly_climatology(r_recent_clip, months_recent)

if (EXPORT_CLIMATOLOGY) {
  # Outputs: 12-band monthly climatology TIFFs for legacy and recent baselines.
  write_tif(
    clim_legacy_12,
    file.path(OUT_R, sprintf("RZSM_NOAH_legacy_%s_clim12_AOI.tif", LEGACY_TAG))
  )
  write_tif(
    clim_recent_12,
    file.path(OUT_R, sprintf("RZSM_NOAH_recent_%s_clim12_AOI.tif", RECENT_TAG))
  )
}

if (EXPORT_SEASONAL_TIF) {
  # Outputs: 4-band snapshot TIFFs for representative months (Jan/Apr/Jul/Oct).
  write_tif(
    subset_month_layers(clim_legacy_12, SEASON_MONTHS),
    file.path(
      OUT_R,
      sprintf("RZSM_NOAH_legacy_%s_season_%s_AOI.tif", LEGACY_TAG, paste(SEASON_MONTHS, collapse = "_"))
    )
  )
  write_tif(
    subset_month_layers(clim_recent_12, SEASON_MONTHS),
    file.path(
      OUT_R,
      sprintf("RZSM_NOAH_recent_%s_season_%s_AOI.tif", RECENT_TAG, paste(SEASON_MONTHS, collapse = "_"))
    )
  )
}

if (EXPORT_ANOMALY) {
  # Outputs: anomaly climatology (recent minus legacy), plus seasonal snapshot subset.
  common_names <- month_names_short[
    month_names_short %in% names(clim_recent_12) &
      month_names_short %in% names(clim_legacy_12)
  ]
  anom_12 <- clim_recent_12[[common_names]] - clim_legacy_12[[common_names]]
  names(anom_12) <- common_names
  
  write_tif(
    anom_12,
    file.path(
      OUT_R,
      sprintf("RZSM_NOAH_anomaly_recent_%s_minus_legacy_%s_clim12_AOI.tif", RECENT_TAG, LEGACY_TAG)
    )
  )
  
  if (EXPORT_SEASONAL_TIF) {
    write_tif(
      subset_month_layers(anom_12, SEASON_MONTHS),
      file.path(
        OUT_R,
        sprintf(
          "RZSM_NOAH_anomaly_recent_%s_minus_legacy_%s_season_%s_AOI.tif",
          RECENT_TAG, LEGACY_TAG, paste(SEASON_MONTHS, collapse = "_")
        )
      )
    )
  }
}

# ----------------------------
# 8) Drought frequency products + CSVs
#    Thresholds: computed from THRESHOLD_YEARS (full record)
#    Frequency:  computed from RECENT_YEARS (recent) vs those thresholds
# ----------------------------

cat("\n=== DROUGHT FREQUENCY ===\n")

if (EXPORT_DROUGHT_FREQ) {
  # Outputs: monthly drought frequency climatologies per percentile, plus optional exports.
  
  long_rows_allP <- list()
  
  for (p in DROUGHT_PERCENTILES) {
    
    prob <- p / 100
    message(sprintf("Percentile P%d...", p))
    
    # 8.1 Thresholds per month (FULL RECORD)
    thr_list <- vector("list", 12)
    
    for (m in 1:12) {
      idx_t <- which(months_threshold == m)
      if (length(idx_t) == 0) {
        thr_list[[m]] <- NULL
        next
      }
      
      thr_m <- app(r_threshold_clip[[idx_t]], fun = pixel_quantile, prob = prob)
      names(thr_m) <- month_names_short[m]
      thr_list[[m]] <- thr_m
    }
    
    thr_list <- thr_list[!vapply(thr_list, is.null, logical(1))]
    if (length(thr_list) == 0) stop("No threshold layers produced for P", p)
    
    thr12 <- rast(thr_list)
    
    # 8.2 Frequency per month (RECENT)
    freq_list <- vector("list", 12)
    
    for (m in 1:12) {
      idx_r <- which(months_recent == m)
      if (length(idx_r) == 0) {
        freq_list[[m]] <- NULL
        next
      }
      
      thr_name <- month_names_short[m]
      if (!(thr_name %in% names(thr12))) {
        freq_list[[m]] <- NULL
        next
      }
      
      thr_m <- thr12[[thr_name]]
      below <- r_recent_clip[[idx_r]] < thr_m
      freq_m <- app(below, fun = mean, na.rm = TRUE)
      
      names(freq_m) <- thr_name
      freq_list[[m]] <- freq_m
    }
    
    freq_list <- freq_list[!vapply(freq_list, is.null, logical(1))]
    if (length(freq_list) == 0) stop("No recent frequency layers produced for P", p)
    
    freq12 <- rast(freq_list)
    
    # ---------------------------------------------------------
    # FIX: ALWAYS WRITE clim12 FOR EACH PERCENTILE (P10/P20/P30)
    # ---------------------------------------------------------
    write_tif(
      freq12,
      file.path(
        OUT_R,
        sprintf(
          "RZSM_NOAH_droughtfreq_recent_%s_below_threshold_%s_P%02d_clim12_AOI.tif",
          RECENT_TAG, THRESH_TAG, p
        )
      )
    )
    
    # 8.3 Monthly single-band exports
    if (EXPORT_DROUGHT_MONTHLY_TIF) {
      monthly_dir <- file.path(OUT_R, sprintf("P%02d_monthly", p))
      dir.create(monthly_dir, recursive = TRUE, showWarnings = FALSE)
      
      for (m in 1:12) {
        nm <- month_names_short[m]
        if (!(nm %in% names(freq12))) next
        
        write_tif(
          freq12[[nm]],
          file.path(
            monthly_dir,
            sprintf(
              "RZSM_NOAH_droughtfreq_recent_%s_below_threshold_%s_P%02d_%s_AOI.tif",
              RECENT_TAG, THRESH_TAG, p, nm
            )
          )
        )
      }
    }
    
    # 8.4 Seasonal export
    if (EXPORT_DROUGHT_SEASONAL && EXPORT_SEASONAL_TIF) {
      
      season_names <- month_names_short[SEASON_MONTHS]
      season_names <- season_names[season_names %in% names(freq12)]
      
      if (length(season_names) > 0) {
        write_tif(
          freq12[[season_names]],
          file.path(
            OUT_R,
            sprintf(
              "RZSM_NOAH_droughtfreq_recent_%s_below_threshold_%s_P%02d_season_%s_AOI.tif",
              RECENT_TAG, THRESH_TAG, p, paste(SEASON_MONTHS, collapse = "_")
            )
          )
        )
      } else {
        warning("Seasonal export skipped (none of season months exist in freq12) for P", p)
      }
    }
    
    # 8.5 ALLMONTHS export
    freq_all <- NULL
    
    if (EXPORT_DROUGHT_ALLMONTHS) {
      
      thr_for_recent_list <- vector("list", nlyr(r_recent_clip))
      
      for (i in 1:nlyr(r_recent_clip)) {
        m <- months_recent[i]
        nm <- month_names_short[m]
        
        if (!(nm %in% names(thr12))) {
          thr_for_recent_list[[i]] <- NULL
          next
        }
        
        thr_for_recent_list[[i]] <- thr12[[nm]]
      }
      
      keep_i <- !vapply(thr_for_recent_list, is.null, logical(1))
      
      if (!any(keep_i)) {
        warning("ALLMONTHS export skipped: no threshold layers matched recent layers for P", p)
      } else {
        
        thr_for_recent <- rast(thr_for_recent_list[keep_i])
        recent_keep    <- r_recent_clip[[which(keep_i)]]
        
        below_all <- recent_keep < thr_for_recent
        freq_all  <- app(below_all, fun = mean, na.rm = TRUE)
        
        names(freq_all) <- sprintf("DroughtFreq_below_P%02d_ALLMONTHS", p)
        
        write_tif(
          freq_all,
          file.path(
            OUT_R,
            sprintf(
              "RZSM_NOAH_droughtfreq_recent_%s_below_threshold_%s_P%02d_ALLMONTHS_AOI.tif",
              RECENT_TAG, THRESH_TAG, p
            )
          )
        )
      }
    }
    
    # 8.6 CSV exports (per percentile + long combined)
    if (EXPORT_DROUGHT_CSV || EXPORT_DROUGHT_CSV_LONG) {
      
      rows <- list()
      
      for (m in 1:12) {
        nm <- month_names_short[m]
        if (!(nm %in% names(freq12))) next
        
        st <- layer_stats(freq12[[nm]])
        ar <- area_above_thresholds(freq12[[nm]], DROUGHT_AREA_THRESHOLDS)
        
        rows[[length(rows) + 1]] <- c(
          list(
            percentile    = p,
            recent_tag    = RECENT_TAG,
            threshold_tag = THRESH_TAG,
            band_type     = "monthly",
            month_num     = m,
            month_name    = nm
          ),
          st, ar
        )
      }
      
      if (EXPORT_DROUGHT_ALLMONTHS && !is.null(freq_all)) {
        st_all <- layer_stats(freq_all)
        ar_all <- area_above_thresholds(freq_all, DROUGHT_AREA_THRESHOLDS)
        
        rows[[length(rows) + 1]] <- c(
          list(
            percentile    = p,
            recent_tag    = RECENT_TAG,
            threshold_tag = THRESH_TAG,
            band_type     = "allmonths",
            month_num     = NA,
            month_name    = "ALLMONTHS"
          ),
          st_all, ar_all
        )
      }
      
      if (length(rows) == 0) {
        warning("No CSV rows produced for P", p)
      } else {
        df <- do.call(rbind, lapply(rows, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
        
        if (EXPORT_DROUGHT_CSV) {
          out_csv <- file.path(
            OUT_R,
            sprintf(
              "RZSM_NOAH_droughtfreq_summary_recent_%s_threshold_%s_P%02d_AOI.csv",
              RECENT_TAG, THRESH_TAG, p
            )
          )
          write.csv(df, out_csv, row.names = FALSE)
          cat("Wrote CSV:", out_csv, "\n")
        }
        
        if (EXPORT_DROUGHT_CSV_LONG) {
          pct_cols <- vapply(
            DROUGHT_AREA_THRESHOLDS,
            function(t) sprintf("pct_area_ge_%.2f", t),
            FUN.VALUE = character(1)
          )
          
          keep_cols <- c(
            "percentile","recent_tag","threshold_tag","band_type","month_num","month_name",
            "mean","median","sd","min","max", pct_cols, "total_area_m2"
          )
          existing <- keep_cols[keep_cols %in% names(df)]
          df_plot <- df[, existing, drop = FALSE]
          if ("mean" %in% names(df_plot)) names(df_plot)[names(df_plot) == "mean"] <- "mean_freq"
          
          long_rows_allP[[length(long_rows_allP) + 1]] <- df_plot
        }
      }
    }
  }
  
  # 8.7 Long CSV combined across percentiles
  if (EXPORT_DROUGHT_CSV_LONG && length(long_rows_allP) > 0) {
    df_long_all <- do.call(rbind, long_rows_allP)
    suppressWarnings(df_long_all$month_num <- as.numeric(df_long_all$month_num))
    
    out_csv_long <- file.path(
      OUT_R,
      sprintf(
        "RZSM_NOAH_droughtfreq_LONG_for_plotting_recent_%s_threshold_%s_AOI.csv",
        RECENT_TAG, THRESH_TAG
      )
    )
    write.csv(df_long_all, out_csv_long, row.names = FALSE)
    cat("Wrote LONG CSV:", out_csv_long, "\n")
  }
}

cat("\n=== FINISHED SUCCESSFULLY ===\n")
