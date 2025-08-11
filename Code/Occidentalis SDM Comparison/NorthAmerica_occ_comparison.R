library(terra)
library(ggplot2)
library(dplyr)

#Load presence points and calculate thresholds based on suitability at presences

# Load presence points CSVs (update paths)
pre_occ <- read.csv("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Data/Occidentalis Comparison/Pre_presence_points.csv")   # historical presence points
post_occ <- read.csv("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Data/Occidentalis Comparison/Post_presence_points.csv") # current presence points

# Load predicted suitability rasters
historical_suit <- rast("/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_historical_vs_future/SDMHistoric_occ.tif")
current_suit <- rast("/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_historical_vs_future/SDMfuture_occ.tif")

# Convert presence points to SpatVector with correct CRS
pre_occ_v <- vect(pre_occ, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(historical_suit))
post_occ_v <- vect(post_occ, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(current_suit))

# Extract suitability values at presence points
pre_suit_vals <- terra::extract(historical_suit, pre_occ_v)[,2]
post_suit_vals <- terra::extract(current_suit, post_occ_v)[,2]

# Calculate threshold as 10th percentile of presence suitability values (captures 90% of presences)
thresh_historical <- quantile(pre_suit_vals, probs = 0.1, na.rm = TRUE)
thresh_current <- quantile(post_suit_vals, probs = 0.1, na.rm = TRUE)

cat("Historical threshold (10th percentile):", thresh_historical, "\n")
cat("Current threshold (10th percentile):", thresh_current, "\n")

# 1. Binarize suitability rasters using data-driven thresholds
historical_bin <- historical_suit >= thresh_historical
current_bin <- current_suit >= thresh_current

# 2. Calculate total suitable area (all raster extent)
calc_suitable_area <- function(bin_raster) {
  cell_area <- cellSize(bin_raster, unit = "km")   # area per cell in km²
  s <- global(bin_raster * cell_area, sum, na.rm = TRUE)
  return(as.numeric(s))
}

hist_area_total <- calc_suitable_area(historical_bin)
curr_area_total <- calc_suitable_area(current_bin)

cat("Total suitable area (historical):", hist_area_total, "km²\n")
cat("Total suitable area (current):", curr_area_total, "km²\n")

# 3. Find overlapping extent
ext1 <- ext(historical_suit)
ext2 <- ext(current_suit)

# Manually intersect their bounding boxes
xmin <- max(ext1[1], ext2[1])
xmax <- min(ext1[2], ext2[2])
ymin <- max(ext1[3], ext2[3])
ymax <- min(ext1[4], ext2[4])

# Create the SpatExtent manually
common_ext <- ext(xmin, xmax, ymin, ymax)

hist_crop <- crop(historical_bin, common_ext)
curr_crop <- crop(current_bin, common_ext)

# 4. Calculate suitable area inside overlap
hist_area_overlap <- calc_suitable_area(hist_crop)
curr_area_overlap <- calc_suitable_area(curr_crop)

cat("Suitable area inside overlap (historical):", hist_area_overlap, "km²\n")
cat("Suitable area inside overlap (current):", curr_area_overlap, "km²\n")

# 5. Calculate suitable area outside overlap
calc_outside_area <- function(full_raster, overlap_extent) {
  overlap_poly <- as.polygons(overlap_extent)
  crs(overlap_poly) <- crs(full_raster)
  
  outside <- mask(full_raster, overlap_poly, inverse = TRUE)
  
  calc_suitable_area(outside)
}

# 6. Optional: visualize
prepare_plot_df <- function(bin_raster, label) {
  df <- as.data.frame(bin_raster, xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("x", "y", "suitable")
  # Use logical levels, not numeric
  df$suitable <- factor(df$suitable, levels = c(FALSE, TRUE), labels = c("Unsuitable", "Suitable"))
  df$dataset <- label
  return(df)
}

hist_df <- prepare_plot_df(historical_bin, "Historical")
curr_df <- prepare_plot_df(current_bin, "Current")

plot_df <- bind_rows(hist_df, curr_df)

ggplot(plot_df) +
  geom_tile(aes(x = x, y = y, fill = suitable)) +
  facet_wrap(~ dataset, ncol = 1) +
  scale_fill_manual(values = c("grey90", "forestgreen")) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Changes in Suitable Habitat for Bombus occidentalis",
       fill = "Habitat Suitability")

# Export binarized suitability rasters with data-driven thresholds
writeRaster(historical_bin,
            filename = "/Users/whitneymaxfield/Downloads/Historical_Suitability_ThreshBased.tif",
            overwrite = TRUE)

writeRaster(current_bin,
            filename = "/Users/whitneymaxfield/Downloads/Current_Suitability_ThreshBased.tif",
            overwrite = TRUE)
