# Install required packages if not already installed
install.packages(c("sf", "FedData"))

# Load packages
library(terra)
library(sf)
library(tidyverse)
library(dismo)
library(geodata)
library(maxnet)
library(sp)
library(raster)

# ----------------------------------------
# 1. Define your study area (OR, CA, NV)
# ----------------------------------------
# Read in shapefiles
oregon <- st_read("/Users/whitneymaxfield/Downloads/Oregon/tl_2019_41_cousub.shp")
california <- st_read("/Users/whitneymaxfield/Downloads/Cali_shp/tl_2019_06_cousub.shp")
nevada <- st_read("/Users/whitneymaxfield/Downloads/Nevada_shape/tl_2019_32_cousub.shp")  # [NEW]

# Transform to EPSG:5070 (NLCD projection)
crs_target <- "EPSG:5070"
oregon <- st_transform(oregon, crs = crs_target)
california <- st_transform(california, crs = crs_target)
nevada <- st_transform(nevada, crs = crs_target)

# Combine shapefiles
or_ca_nv_combined <- rbind(oregon, california, nevada)

# Convert to SpatVector
aoi <- vect(or_ca_nv_combined)
crs(aoi) <- crs_target
or_ca_nv_vect <- vect(or_ca_nv_combined)

# ----------------------------------------
# 2. Load NLCD land cover raster
# ----------------------------------------
nlcd_1998_path <- "/Users/whitneymaxfield/Downloads/Annual_NLCD_LndCov_1998_CU_C1V1/Annual_NLCD_LndCov_1998_CU_C1V1.tif"
nlcd_1998 <- rast(nlcd_1998_path)

# Crop and mask to AOI
landcover_clipped <- crop(nlcd_1998, or_ca_nv_vect)
landcover_masked <- mask(landcover_clipped, or_ca_nv_vect)

# ----------------------------------------
# 3. Reclassify NLCD to floral habitat proxy
# ----------------------------------------
floral_classes <- c(52, 71, 81)
rcl <- matrix(c(
  0,  51, 0,
  52, 52, 1,
  53, 70, 0,
  71, 71, 1,
  72, 80, 0,
  81, 81, 1,
  82, 100, 0
), ncol = 3, byrow = TRUE)
floral_layer <- classify(landcover_masked, rcl)

# ----------------------------------------
# 4. Plot check
# ----------------------------------------
plot(floral_layer, main = "Floral Resource Proxy: Open Habitat (NLCD 1998)")
legend("topright", legend = c("Non-floral", "Open floral habitat"), fill = c("gray", "green"))

# ----------------------------------------
# 5. Resample to match climate raster
# ----------------------------------------
clim <- worldclim_global(var = "bio", res = 2.5, version = 2.1, path = "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old")
climList <- list.files(
  path = "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m",
  pattern = ".tif$",
  full.names = TRUE
)
clim_raster <- raster::stack(climList)
clim_raster_spat <- rast(clim_raster)

# Reproject and resample floral raster
floral_reprojected <- project(floral_layer, crs(clim_raster_spat), method = "near")
floral_resampled <- resample(floral_reprojected, clim_raster_spat, method = "near")

# Combine climate + floral
env_stack <- c(floral_resampled, clim_raster_spat)
or_ca_nv_vect <- project(or_ca_nv_vect, crs(env_stack))

# Crop and mask to AOI
env_stack_cropped <- crop(env_stack, or_ca_nv_vect)
env_stack_masked <- mask(env_stack_cropped, or_ca_nv_vect)

# Sample environmental space
sample_points <- spatSample(env_stack_masked, size = 2000, method = "random", na.rm = TRUE)
valid_mask <- !is.na(env_stack_masked[[1]])
for (i in 2:nlyr(env_stack_masked)) {
  valid_mask <- valid_mask & !is.na(env_stack_masked[[i]])
}
terra::global(valid_mask, "sum")
env_df <- as.data.frame(sample_points, na.rm = TRUE)

# ----------------------------------------
# 6. Multicollinearity filtering
# ----------------------------------------
library(corrplot)
cor_matrix <- cor(env_df, use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.4)

library(usdm)
vif_result <- vifstep(env_df, th = 10)
selected_vars <- vif_result@results$Variables
env_stack_filtered <- env_stack[[selected_vars]]

# ----------------------------------------
# 7. SDM with maxnet
# ----------------------------------------
shapefile_path <- "/Users/whitneymaxfield/Downloads/Franklins_historic_locations/Franklins_historic_locations.shp"
presence_points <- vect(shapefile_path)
presence_env <- extract(env_stack_filtered, presence_points)
presence_env$presence <- 1

# Generate background points
valid_cells <- which(!is.na(values(env_stack_filtered[[1]])))
set.seed(123)
sampled_cells <- sample(valid_cells, size = 5000)
background_points <- xyFromCell(env_stack_filtered[[1]], sampled_cells)
background_points_vect <- vect(background_points, type = "points", crs = crs(env_stack_filtered))
background_env <- extract(env_stack_filtered, background_points_vect)
background_env$presence <- 0

env_data <- rbind(presence_env, background_env)
env_data <- env_data[complete.cases(env_data), ]

response <- env_data$presence
predictors <- env_data[, setdiff(names(env_data), "presence")]
predictors_clean <- predictors[, !names(predictors) %in% c("ID", "presence")]
colnames(predictors_clean) <- make.names(colnames(predictors_clean))

maxnet_model <- maxnet(p = response, data = predictors_clean, f = maxnet.formula(response, predictors_clean))
suitability_raster <- terra::predict(env_stack_filtered, maxnet_model, type = "cloglog", na.rm = TRUE)

# Save and plot
writeRaster(suitability_raster, "/Users/whitneymaxfield/Desktop/Franklins_attempt_2_withLC.tif", overwrite = TRUE)
suit_df <- as.data.frame(suitability_raster, xy = TRUE, na.rm = TRUE)
colnames(suit_df)[3] <- "suitability"

# ----------------------------------------
# 8. Plot with counties and state outlines
# ----------------------------------------
library(tigris)
oregon_counties <- counties(state = "OR", cb = TRUE, class = "sf")
california_counties <- counties(state = "CA", cb = TRUE, class = "sf")
nevada_counties <- counties(state = "NV", cb = TRUE, class = "sf")
raster_crs <- crs(suitability_raster, proj = TRUE)

oregon_counties_proj <- st_transform(oregon_counties, crs = raster_crs)
california_counties_proj <- st_transform(california_counties, crs = raster_crs)
nevada_counties_proj <- st_transform(nevada_counties, crs = raster_crs)

or_ca_nv_counties <- rbind(oregon_counties_proj, california_counties_proj, nevada_counties_proj)

ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = or_ca_nv_counties, fill = NA, color = "white", size = 0.4) +
  scale_fill_viridis_c(option = "plasma") +
  coord_sf(crs = raster_crs, expand = FALSE) +
  theme_minimal() +
  labs(title = "Franklin's Bumblebee Habitat Suitability in OR, CA, and NV")

# Just state lines
states <- states(cb = TRUE, class = "sf")
target_states <- c("Oregon", "California", "Washington", "Idaho", "Nevada")
states_subset <- states %>% filter(NAME %in% target_states)
states_proj <- st_transform(states_subset, crs = raster_crs)

ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = states_proj, fill = NA, color = "white", size = 0.8) +
  scale_fill_viridis_c(option = "plasma") +
  coord_sf(crs = raster_crs, expand = FALSE) +
  theme_minimal() +
  labs(title = "Franklin's Bumblebee Habitat Suitability with State Boundaries")

# Save final output
writeRaster(suitability_raster, filename = "/Users/whitneymaxfield/Desktop/Franklins_withLC_final2.tif", overwrite = TRUE)

