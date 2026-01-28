library(terra)

# --- Load data ---
overlap <- rast("Lynx_Rabbit_Overlay_RasterExtraction.tif")
aoi     <- vect("AOI_reprojected.shp")

# Ensure same CRS (safety check)
aoi <- project(aoi, crs(overlap))

# --- Mask raster to AOI (if not already done in QGIS) ---
 overlap <- mask(overlap, aoi)

# --- Ensure proper binary logic ---
# Keep ONLY overlap pixels (1), everything else NA
overlap <- ifel(overlap == 1, 1, NA)

# --- Pixel area (m²) ---
pixel_area <- prod(res(overlap))  # e.g. 30m x 30m = 900

# --- Count overlap pixels ---
overlap_pixels <- global(overlap, "sum", na.rm = TRUE)[1,1]

overlap_area_m2 <- overlap_pixels * pixel_area

# --- Total AOI area ---
aoi_area_m2 <- expanse(aoi, unit = "m")

# --- Percentage overlap ---
overlap_percent <- 100 * overlap_area_m2 / aoi_area_m2

# --- Results ---
overlap_area_km2 <- overlap_area_m2 / 1e6

cat("Overlap area (km²):", round(overlap_area_km2, 2), "\n")
cat("Percentage of AOI:", round(overlap_percent, 2), "%\n")
