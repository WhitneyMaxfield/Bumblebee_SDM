library(terra)
library(ggplot2)
library(dplyr)

# ---- Load your predicted suitability rasters ----
# Replace these with your actual raster files or objects
historical_suit <- rast("/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_historical_vs_future/SDMHistoric_occ.tif")
current_suit <- rast("/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_historical_vs_future/SDMfuture_occ.tif")

# Thresholds for suitability (set these based on your model)
thresh_historical <- 0.5
thresh_current <- 0.5

# 1. Binarize suitability rasters
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
  # Create a polygon from overlap extent
  overlap_poly <- as.polygons(overlap_extent)
  crs(overlap_poly) <- crs(full_raster)  # ensure CRS matches
  
  # Mask full raster by overlap polygon, inverse = TRUE keeps outside
  outside <- mask(full_raster, overlap_poly, inverse = TRUE)
  
  calc_suitable_area(outside)
}

# 6. Optional: visualize

# Prepare data frames for plotting
prepare_plot_df <- function(bin_raster, label) {
  df <- as.data.frame(bin_raster, xy = TRUE, na.rm = TRUE)
  colnames(df) <- c("x", "y", "suitable")
  df$suitable <- factor(df$suitable, levels = c(0,1), labels = c("Unsuitable", "Suitable"))
  df$dataset <- label
  return(df)
}

hist_df <- prepare_plot_df(historical_bin, "Historical")
curr_df <- prepare_plot_df(current_bin, "Current")

plot_df <- bind_rows(hist_df, curr_df)

table(plot_df$suitable)
summary(plot_df)
sum(!is.na(values(historical_bin)))
sum(!is.na(values(current_bin)))

hist_df <- as.data.frame(historical_bin, xy = TRUE, na.rm = TRUE)
table(hist_df$layer)  # check if 'layer' column exists and what values it has
sum(!is.na(values(historical_bin)))

# Check summary of raster values
summary(values(historical_bin))

hist_df_all <- as.data.frame(historical_bin, xy = TRUE, na.rm = FALSE)
summary(hist_df_all)

# How many rows total?
nrow(hist_df_all)

# How many NAs in suitability column?
sum(is.na(hist_df_all[,3]))

# How many TRUE/FALSE in suitability column?
table(hist_df_all[,3], useNA = "ifany")

hist_df <- subset(hist_df_all, !is.na(hist_df_all[,3]))
colnames(hist_df) <- c("x", "y", "suitable")

# Convert logical to factor for plotting
hist_df$suitable <- factor(hist_df$suitable, levels = c(FALSE, TRUE),
                           labels = c("Unsuitable", "Suitable"))
table(hist_df$suitable)






library(ggplot2)

curr_df <- as.data.frame(current_bin, xy = TRUE, na.rm = FALSE) %>%
  subset(!is.na(focal_mean))

colnames(curr_df) <- c("x", "y", "suitable")
curr_df$suitable <- factor(curr_df$suitable,
                           levels = c(FALSE, TRUE),
                           labels = c("Unsuitable", "Suitable"))

plot_df <- rbind(
  transform(hist_df, dataset = "Historical"),
  transform(curr_df, dataset = "Current")
)

ggplot(plot_df) +
  geom_tile(aes(x = x, y = y, fill = suitable)) +
  facet_wrap(~ dataset, ncol = 1) +
  scale_fill_manual(values = c("grey90", "forestgreen")) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Changes in Suitable Habitat for Bombus occidentalis ",
       fill = "Habitat Suitability")

#export
# Export historical suitability raster
writeRaster(historical_bin,
            filename = "/Users/whitneymaxfield/Downloads/Historical_Suitability.tif",
            overwrite = TRUE)

# Export current suitability raster (use the aligned one!)
writeRaster(current_bin,
            filename = "/Users/whitneymaxfield/Downloads/Current_Suitability.tif",
            overwrite = TRUE)

