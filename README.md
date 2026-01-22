# Remote Sensing Scripts (R)

R scripts for geospatial analysis and remote sensing workflows in conservation.

## Scripts
- `esa_worldcover_aoi_stats.R` – ESA WorldCover class composition for an AOI (CSV + bar chart).
- `rzsm_noah_analysis_pipeline.R` – GLDAS NOAH RZSM climatologies, anomalies, and drought frequency products.
- `slope_class_area_stats.R` – Slope class areas and percent coverage.
- `temporal_coverage_check_rzsm.R` – Monthly coverage check from filenames (heatmap).
- `raster_overlap_area_stats.R` – Area/percent overlap from a binary raster within an AOI.

## Notes
- File paths are hard‑coded; update them before running.
- Outputs are written to local folders and may require write permissions.

## Requirements
- `terra`, `ggplot2`, `dplyr`, `tibble`, `stringr`, `tidyr` (varies by script)

## License
MIT (or your preferred license)
