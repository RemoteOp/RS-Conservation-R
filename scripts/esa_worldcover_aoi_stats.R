# ============================================================
# ESA WorldCover land cover composition (AOI-only) — UTM raster
# Corrected to your paths + filenames AND avoids "Permission denied"
#
# What it does:
# - Reads your clipped AOI GeoTIFF (UTM 30N, meters)
# - Computes class pixel counts (freq)
# - Drops background value 0 (bounding-box NoData)
# - Joins ESA class names
# - Computes area (km²) + percent (AOI only)
# - Plots horizontal bar chart using ESA colours
#   * slimmer bars
#   * legend in ESA class order
# - Exports tidy CSV + PNG plot
# ============================================================

library(terra)
library(dplyr)
library(ggplot2)
library(tibble)

# ---- WORKING DIRECTORY (your folder) ----
setwd("C:/Users/...")

# ---- USER SETTINGS ----
in_raster <- "lc_roi_reproj_RasterExtraction.tif"
out_csv   <- "esa_lc_stats_AOI.csv"
out_plot  <- "esa_lc_composition_AOI.png"

# ---- SAFER OUTPUT PATHS (fixes Permission denied) ----
# Writes outputs to a subfolder you own (create if missing)
out_dir <- file.path(getwd(), "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_csv_path  <- file.path(out_dir, out_csv)
out_plot_path <- file.path(out_dir, out_plot)

# ---- READ RASTER ----
lc <- rast(in_raster)

# ---- OPTIONAL SANITY CHECKS ----
# res(lc)
# minmax(lc)
# freq(lc)

# ---- Safety: ensure projected CRS (meters). UTM 30N is projected ----
if (terra::is.lonlat(lc)) {
  stop("Raster appears to be lon/lat (degrees). Please use the UTM 30N GeoTIFF and retry.")
}

# ---- FREQUENCY TABLE (class pixel counts) ----
ff <- as.data.frame(freq(lc))

# ---- terra versions differ: sometimes (value,count), sometimes (layer,value,count) ----
if (ncol(ff) == 3) {
  ff <- ff %>% rename(class_code = value, pixels = count)
} else if (ncol(ff) == 2) {
  names(ff) <- c("class_code", "pixels")
} else {
  stop("Unexpected freq() output shape: ", ncol(ff), " columns")
}

# ---- Drop bounding-box background (NoData stored as 0) ----
ff <- ff %>% filter(class_code != 0)

# ---- Pixel area in m² (UTM units are meters) ----
pixel_area_m2 <- prod(res(lc))

# ---- ESA WORLD COVER LEGEND ----
class_names <- tribble(
  ~class_code, ~class_name,
  10,  "Tree cover",
  20,  "Shrubland",
  30,  "Grassland",
  40,  "Cropland",
  50,  "Built-up",
  60,  "Bare / sparse vegetation",
  70,  "Snow and ice",
  80,  "Water bodies",
  90,  "Herbaceous wetland",
  95,  "Mangroves",
  100, "Moss and lichen"
)

# ---- ESA WorldCover palette (hex) — keep in sync with class_names spelling ----
esa_cols <- c(
  "Tree cover"               = "#006400",
  "Shrubland"                = "#FFBB22",
  "Grassland"                = "#FFFF4C",
  "Cropland"                 = "#F096FF",
  "Built-up"                 = "#FA0000",
  "Bare / sparse vegetation" = "#B4B4B4",
  "Snow and ice"             = "#F0F0F0",
  "Water bodies"             = "#0064C8",
  "Herbaceous wetland"       = "#0096A0",
  "Mangroves"                = "#00CF75",
  "Moss and lichen"          = "#FAE6A0"
)

# ---- ESA order for legend (and optionally axis) ----
esa_order <- c(
  "Tree cover",
  "Shrubland",
  "Grassland",
  "Cropland",
  "Built-up",
  "Bare / sparse vegetation",
  "Water bodies",
  "Herbaceous wetland",
  "Snow and ice",
  "Mangroves",
  "Moss and lichen"
)

# ---- BUILD TIDY OUTPUT TABLE ----
out <- ff %>%
  left_join(class_names, by = "class_code") %>%
  mutate(
    class_name = ifelse(is.na(class_name), paste("Class", class_code), class_name),
    area_km2    = pixels * pixel_area_m2 / 1e6,
    percent     = pixels / sum(pixels) * 100
  )

# ---- Apply factor ordering (controls legend + axis order) ----
out$class_name <- factor(
  out$class_name,
  levels = rev(out$class_name[order(out$percent, decreasing = TRUE)])
)


# ---- Keep only classes that exist in the AOI (so legend doesn't show absent classes) ----
out <- out %>% filter(!is.na(class_name))

# ---- Sort table for readability in CSV (largest % first) ----
out <- out %>% arrange(desc(percent))

# ---- WRITE CSV (to output folder) ----
write.csv(out, out_csv_path, row.names = FALSE)
message("Saved CSV: ", out_csv_path)

# ---- PLOT (ESA colours, slimmer bars, clean grid, legend in ESA order) ----
p <- ggplot(out, aes(x = class_name, y = percent, fill = class_name)) +
  geom_col(width = 0.55) +
  coord_flip() +
  scale_fill_manual(values = esa_cols, breaks = esa_order, drop = TRUE) +
  labs(
    x = NULL,
    y = "Percent of AOI (%)",
    title = "ESA WorldCover composition (AOI only)",
    fill = "Land cover"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right"
  )

print(p)

# ---- SAVE PLOT (to output folder) ----
ggsave(out_plot_path, p, width = 7, height = 4.5, dpi = 300)
message("Saved plot: ", out_plot_path)

# ---- SANITY CHECKS ----
message("Sum of percent (should be ~100): ", round(sum(out$percent), 6))
message("Pixel area (m²): ", pixel_area_m2)
message("Raster resolution: ", paste(res(lc), collapse = " x "))
message("CRS: ", terra::crs(lc))
