# Load required packages
library(terra)
library(tmap)
library(RColorBrewer)

# --- Step 1: Load historical and current suitability maps ---
historic <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Pre and Post 1998 Occ Figures/Occidentalis_<1998prediction_Oregon_geofixed.tif")  # Replace with your file
current  <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Pre and Post 1998 Occ Figures/Occidentalis_post1998prediction_Oregon.tif")   # Replace with your file

# --- Step 2: Resample historic to match current if needed ---
if (!compareGeom(historic, current, stopOnError = FALSE)) {
  historic <- resample(historic, current, method = "bilinear")
}

# --- Step 3: Calculate suitability change ---
suitability_change <- current - historic

# --- Step 4: Classify change into categories using app() ---
# Step 4: Reclassify into 1 = Decrease, 2 = No Change, 3 = Increase
change_class <- app(suitability_change, function(x) {
  ifelse(x < -0.1, 1,        # Decrease
         ifelse(x > 0.1, 3,  # Increase
                2))          # No change
})

# Set labels as factors using positive integers
change_class <- as.factor(change_class)
levels(change_class) <- data.frame(
  ID = c(1, 2, 3),
  Category = c("Decrease", "No Change", "Increase")
)

# Step 5: Plot with correct labels and colors
tm_shape(change_class) +
  tm_raster(title = "Change in Suitability",
            palette = c("red", "gray80", "blue"),
            labels = c("Decrease", "No Change", "Increase")) +
  tm_layout(title = "Change in Habitat Suitability for *Bombus occidentalis*",
            legend.outside = TRUE)

output_path <- "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Code/suitability_change_Occ_classified.tif"
writeRaster(change_class, filename = output_path, overwrite = TRUE)


