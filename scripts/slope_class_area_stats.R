library(terra)
library(ggplot2)

# Load slope raster
slope <- rast("DEM_AOI_slope_Translated.tif")

# Define slope classes
classes <- c(0, 3, 6, 9, 11, 17, 22, 31, Inf)

# Build classification matrix for terra::classify()
# m <- cbind(classes[-length(classes)], classes[-1], seq_along(classes[-1]))

m <- matrix(c(
  0,   3,   1,
  3,   6,   2,
  6,   9,   3,
  9,  11,   4,
  11, 17,   5,
  17, 22,   6,
  22, 31,   7,
  31, Inf,  8
), ncol = 3, byrow = TRUE)


# Classify slope (correct raster method)
slope_class <- classify(slope, m, others = NA)

# Frequency table
freq_table <- as.data.frame(freq(slope_class))
freq_table <- freq_table[!is.na(freq_table$value), ]


# Area calculations
pixel_area_m2 <- abs(prod(res(slope)))
freq_table$area_m2 <- freq_table$count * pixel_area_m2

total_area <- sum(freq_table$area_m2)
freq_table$percent <- 100 * freq_table$area_m2 / total_area

print(freq_table)

# Labels
freq_table$class_label <- paste0(
  classes[freq_table$value], "–",
  classes[freq_table$value + 1], "°"
)

print(freq_table)

# Plot
ggplot(freq_table, aes(x = factor(value), y = percent)) +
  geom_bar(stat = "identity", fill = "#E69F00") +
  labs(x = "Slope class (degrees)",
       y = "Area (%)",
       title = "Percentage of AOI area by slope class") +
  scale_x_discrete(labels = labels) +
  theme_minimal()

write.csv(freq_table, "Slope_Class_Area_Statistics_DegreeClass.csv", row.names = FALSE)
