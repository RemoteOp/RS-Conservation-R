# ============================================================
# Figure 3-2 — Temporal coverage check (sanity-check only)
# - Reads filenames only (no raster reading)
# - Confirms every month is present from 2000–2024
# - Highlights missing months
# - Exports: PNG heatmap + PNG tick plot
# - NEVER overwrites existing outputs
# ============================================================

library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

# ----------------------------
# 1) USER SETTINGS
# ----------------------------
DATA_DIR <- "C:/GLDAS_2000_2024/data"
OUT_DIR  <- "C:/Users/saral/SL_Website/RootZoneSoilMoisture/NOAH_analysis/R"

FILE_PATTERN <- "\\.SUB\\.tif$"

START_YEAR <- 2000
END_YEAR   <- 2024

month_names_short <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

# If FALSE: script will error if the PNG already exists (safe default)
ALLOW_OVERWRITE_PNG <- FALSE

# ----------------------------
# 2) Helper functions
# ----------------------------

extract_yyyymm <- function(x) {
  m <- str_match(basename(x), "A(\\d{6})")[,2]
  if (is.na(m)) return(NA_character_)
  m
}

yyyymm_to_year  <- function(yyyymm) as.integer(substr(yyyymm, 1, 4))
yyyymm_to_month <- function(yyyymm) as.integer(substr(yyyymm, 5, 6))

safe_ggsave_versioned <- function(path, plot, width, height, dpi = 300) {
  dir <- dirname(path)
  base <- tools::file_path_sans_ext(basename(path))
  ext  <- tools::file_ext(path)
  
  # If file exists, append _vNN
  if (file.exists(path)) {
    i <- 2
    repeat {
      candidate <- file.path(dir, sprintf("%s_v%02d.%s", base, i, ext))
      if (!file.exists(candidate)) {
        path <- candidate
        break
      }
      i <- i + 1
      if (i > 99) stop("Too many versions exist for: ", base)
    }
  }
  
  ggsave(path, plot, width = width, height = height, dpi = dpi)
  if (!file.exists(path)) stop("ggsave reported success but file not found: ", path)
  path
}

# ----------------------------
# 3) Input checks
# ----------------------------

cat("\n=== TEMPORAL COVERAGE SANITY CHECK ===\n")
cat("DATA_DIR:", DATA_DIR, "\n")
cat("OUT_DIR: ", OUT_DIR,  "\n")

if (!dir.exists(DATA_DIR)) stop("DATA_DIR does not exist: ", DATA_DIR)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) stop("OUT_DIR could not be created: ", OUT_DIR)

files <- list.files(DATA_DIR, pattern = FILE_PATTERN, full.names = TRUE)
cat("Found files:", length(files), "\n")
if (length(files) == 0) stop("No files found. Check DATA_DIR and FILE_PATTERN.")

yyyymm <- vapply(files, extract_yyyymm, FUN.VALUE = character(1))
ok <- !is.na(yyyymm)
files <- files[ok]
yyyymm <- yyyymm[ok]

cat("Files with Ayyyymm parsed:", length(files), "\n")
if (length(files) == 0) stop("No files matched Ayyyymm pattern in filenames.")

years  <- vapply(yyyymm, yyyymm_to_year,  FUN.VALUE = integer(1))
months <- vapply(yyyymm, yyyymm_to_month, FUN.VALUE = integer(1))

# keep only within the check window
keep <- years >= START_YEAR & years <= END_YEAR
years <- years[keep]
months <- months[keep]

cat("Year range (kept):", min(years), "-", max(years), "\n")

# ----------------------------
# 4) Build expected vs observed table
# ----------------------------

observed <- tibble(year = years, month = months) %>%
  distinct() %>%
  mutate(present = 1L)

expected <- expand_grid(
  year = START_YEAR:END_YEAR,
  month = 1:12
) %>%
  left_join(observed, by = c("year", "month")) %>%
  mutate(
    present = ifelse(is.na(present), 0L, present),
    month_name = factor(month_names_short[month], levels = month_names_short)
  )

missing_tbl <- expected %>% filter(present == 0L) %>% arrange(year, month)

cat("\nExpected months:", nrow(expected), "\n")
cat("Observed distinct year-months:", nrow(observed), "\n")
cat("Missing months:", nrow(missing_tbl), "\n")

if (nrow(missing_tbl) > 0) {
  cat("\nMissing year-months:\n")
  print(missing_tbl %>% mutate(month = sprintf("%02d", month)) %>% select(year, month))
} else {
  cat("\n✅ No missing months detected in", START_YEAR, "-", END_YEAR, "\n")
}

# ----------------------------
# 5) Plot: Heatmap
# ----------------------------

p_heat <- ggplot(expected, aes(x = month_name, y = year, fill = factor(present))) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_y_reverse(breaks = seq(END_YEAR, START_YEAR, by = -2)) +
  scale_fill_manual(
    values = c("0" = "#c0392b", "1" = "#27ae60"),
    labels = c("Missing", "Present"),
    name = NULL
  ) +
  labs(
    title = "                      Temporal coverage check",
    subtitle = sprintf("Monthly completeness %d–%d (missing months highlighted)", START_YEAR, END_YEAR),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

out_heat <- file.path(OUT_DIR, sprintf("TemporalCoverage_%d_%d_heatmap.png", START_YEAR, END_YEAR))
saved_heat <- safe_ggsave_versioned(out_heat, p_heat, width = 8.5, height = 10.5, dpi = 300)
cat("Saved heatmap:", saved_heat, "\n")

print(p_heat)

cat("\n=== DONE ===\n")
