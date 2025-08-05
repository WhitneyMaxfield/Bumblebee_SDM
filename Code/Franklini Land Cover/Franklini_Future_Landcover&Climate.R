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
# 1. Define your study area (e.g., Oregon)
# ----------------------------------------
# Read in shapefiles
library(sf)
library(terra)

# Step 1: Read shapefiles
oregon <- st_read("/Users/whitneymaxfield/Downloads/Oregon/tl_2019_41_cousub.shp")
california <- st_read("/Users/whitneymaxfield/Downloads/Cali_shp/tl_2019_06_cousub.shp")
nevada <- st_read("/Users/whitneymaxfield/Downloads/Nevada_shape/tl_2019_32_cousub.shp")

# Step 2: Transform both to EPSG:5070 (NLCD projection)
crs_target <- "EPSG:5070"
oregon <- st_transform(oregon, crs = crs_target)
california <- st_transform(california, crs = crs_target)
nevada <- st_transform(nevada, crs = crs_target)

# Step 3: Combine shapefiles
or_ca_combined <- rbind(oregon, california, nevada)

# Step 4: Convert to SpatVector for terra workflows
aoi <- vect(or_ca_combined)  # Converts directly from sf to terra format
crs(aoi) <- crs_target        # Ensure CRS is assigned (should already be)

# 2. Download NLCD 2024 land cover layer
# ----------------------------------------
#manually download LNCD land cover layers from https://www.mrlc.gov/data?f%5B0%5D=category%3ALand%20Cover&f%5B1%5D=project_tax_term_term_parents_tax_term_name%3AAnnual%20NLCD
# After downloading NLCD (e.g., 1992) from MRLC, unzip it and locate the .img or .tif file
nlcd_2024_path <- "/Users/whitneymaxfield/Downloads/Annual_NLCD_LndCov_2024_CU_C1V1/Annual_NLCD_LndCov_2024_CU_C1V1.tif"  # or .tif??
                 
# Read the raster using terra
nlcd_2024 <- rast(nlcd_2024_path)

# ----------------------------------------
# 3. Reclassify NLCD to open floral habitats
# ----------------------------------------
# NLCD legend: https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description

#making landcover clipped so its much smaller 
or_ca_vect <- vect(or_ca_combined)

# Crop and mask (this can take a bit)
landcover_clipped <- crop(nlcd_2024, or_ca_vect)
landcover_masked <- mask(landcover_clipped, or_ca_vect)

# Define classes that correspond to floral-rich, open habitats:
# (e.g., 52 = Shrubland, 71 = Grassland/Herbaceous, 81 = Pasture/Hay)
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
# 4. Plot to check
# ----------------------------------------
plot(floral_layer, main = "Floral Resource Proxy: Open Habitat (NLCD 2024)")
legend("topright", legend = c("Non-floral", "Open floral habitat"), fill = c("gray", "green"))

# ----------------------------------------
# 5.CLIMATE: This is very different from the original historical model 
# because we are predicting future model best practice from what i have read is to stack multiple models 
# from world clim. I decided to download MPI, MRI, and MIROC6 with: SSP245, 2021–2040, 2.5 minutes, "bioclim" (bc)
# ----------------------------------------

# Paths to each model folder
library(terra)
library(raster)

# Paths to manually downloaded models
library(terra)
library(raster)

# Paths to the stacked raster files (full filenames, not folders)
model_files <- c(
  "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Data/Franklini Data/Climate_for_future/wc2.1_2.5m_bioc_MIROC6_ssp245_2021-2040.tif",
  "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Data/Franklini Data/Climate_for_future/wc2.1_2.5m_bioc_MPI-ESM1-2-HR_ssp245_2021-2040.tif",
  "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Data/Franklini Data/Climate_for_future/wc2.1_2.5m_bioc_MRI-ESM2-0_ssp245_2021-2040.tif"
)

# Load all climate stacks as SpatRasters
model_stacks <- lapply(model_files, rast)

#aoi in wrong CRS so projecting
aoi_geo <- project(aoi, "EPSG:4326")

# Crop and mask each model stack to AOI
model_stacks_cropped <- lapply(model_stacks, function(m) {
  m_crop <- crop(m, aoi_geo)
  m_mask <- mask(m_crop, aoi_geo)
  return(m_mask)
})

# Combine all cropped model stacks into one big stack
all_climate <- rast(model_stacks_cropped)

# For each bioclim layer (1 to 19), extract layers from all models and calculate mean
bio_mean <- lapply(1:19, function(i) {
  layers_i <- all_climate[[seq(i, nlyr(all_climate), by = 19)]]
  app(layers_i, mean)
})

# Stack ensemble mean layers into one SpatRaster
clim_raster_spat <- rast(bio_mean)
names(clim_raster_spat) <- paste0("bio", 1:19)

# Convert to raster::stack for backward compatibility if needed
clim_raster <- stack()
for(i in 1:19){
  r <- raster(clim_raster_spat[[i]])
  clim_raster <- addLayer(clim_raster, r)
}
names(clim_raster) <- paste0("bio", 1:19)

# Check results by plotting one layer
plot(clim_raster_spat[[1]], main = "Ensemble Mean Bio1 (Clipped to AOI)")

# Now clim_raster and clim_raster_spat contain your ensemble mean layers
# Save if needed
# writeRaster(clim_raster_spat, "/path/to/save/clim_raster_spat.tif", overwrite=TRUE)

# Now resample floral_layer to clim_raster_spat
floral_layer_resampled <- resample(floral_layer, clim_raster_spat, method = "near")

# ----------------------------------------
# testing for multicolinearity 
# ----------------------------------------

# Resample floral raster to match resolution and extent of climate stack
floral_resampled <- resample(floral_layer, clim_raster_spat[[1]], method = "near")

#check crs 
crs(clim_raster_spat)
crs(floral_resampled)
#flora is wrong crs so convert: 
crs_clim <- crs(clim_raster_spat)
# Reproject floral raster to match climate raster CRS
floral_reprojected <- project(floral_layer, crs(clim_raster_spat), method = "near")
                 
#extents dont match (floral has to be clipped due to MRLC constraints)
#Have to resample: Resampling aligns all pixels to the climate raster grid, so each cell corresponds spatially.
#Cells outside the MRLC extent will be filled with NA in floral_resampled: 
floral_resampled <- resample(floral_reprojected, clim_raster_spat, method = "near")

# Combine into one stack
env_stack <- c(floral_resampled, clim_raster_spat)
# Reproject vector to match raster CRS
or_ca_vect <- project(or_ca_vect, crs(env_stack))

# Crop and mask the raster
env_stack_cropped <- crop(env_stack, or_ca_vect)
env_stack_masked <- mask(env_stack_cropped, or_ca_vect)

# Crop raster stack to Oregon extent
crs(env_stack)      # For SpatRaster
crs(or_ca_vect) 
env_stack_cropped <- crop(env_stack, or_ca_vect)

# Check size after crop
print(env_stack_cropped)
#NOW: 
sample_points <- spatSample(env_stack_cropped, size = 2000, method = "random", na.rm = TRUE)
#checking numbers: 
valid_mask <- !is.na(env_stack_cropped[[1]])
for (i in 2:nlyr(env_stack_cropped)) {
  valid_mask <- valid_mask & !is.na(env_stack_cropped[[i]])
}
terra::global(valid_mask, "sum")
                 
# Convert to a data frame
env_df <- as.data.frame(sample_points, na.rm = TRUE)

# Install if needed
library(corrplot)

# Correlation matrix
cor_matrix <- cor(env_df, use = "complete.obs")

# Plot
corrplot(cor_matrix, method = "color", tl.cex = 0.4)

# Install if needed
library(usdm)

# Run VIF analysis (threshold 10 by default, or try 5 for stricter filtering)
vif_result <- vifstep(env_df, th = 10)

# View retained and excluded variables
vif_result@results     # Retained
vif_result@excluded    # Dropped due to high collinearity
                 
# Keep only retained variables
selected_vars <- vif_result@results$Variables
env_stack_filtered <- env_stack[[selected_vars]]
                 
# ----------------------------------------
# SDM !!!
# ----------------------------------------

# Extract environmental values at presence points
shapefile_path <- "/Users/whitneymaxfield/Downloads/Franklins_historic_locations/Franklins_historic_locations.shp"           # Replace with actual path
presence_points <- vect(shapefile_path)
presence_env <- extract(env_stack_filtered, presence_points)
presence_env$presence <- 1

#trying this instead (previous background sampling not working)
env_stack_filtered <- crop(env_stack_filtered, or_ca_vect)
# 1. Identify all valid cell indices (not NA)
valid_cells <- which(!is.na(values(env_stack_filtered[[1]])))

# 2. Randomly sample 5000 of those cell indices
set.seed(123)
sampled_cells <- sample(valid_cells, size = 5000)

# 3. Convert those cells to spatial points
background_points <- xyFromCell(env_stack_filtered[[1]], sampled_cells)

# 4. Convert to SpatVector so you can extract full env data
background_points_vect <- vect(background_points, type = "points", crs = crs(env_stack_filtered))

# 5. Extract environmental values from full stack
background_env <- extract(env_stack_filtered, background_points_vect)
background_env$presence <- 0

# Combine presence + background
env_data <- rbind(presence_env, background_env)
env_data <- env_data[complete.cases(env_data), ]  # Remove any rows with NA

# Response vector: presence (1) and background (0)
response <- env_data$presence

# Predictor variables (drop the 'presence' column)
predictors <- env_data[, setdiff(names(env_data), "presence")]

# Fit maxnet model
# Remove non-predictor columns like 'ID'
predictors_clean <- predictors[, !names(predictors) %in% c("ID", "presence")]
                 
# Clean column names
colnames(predictors_clean) <- make.names(colnames(predictors_clean))
                 
# Re-run maxnet
maxnet_model <- maxnet(p = response, data = predictors_clean, f = maxnet.formula(response, predictors_clean))
# Predict suitability raster using cloglog output (interpretable as habitat suitability)
suitability_raster <- terra::predict(env_stack_filtered, maxnet_model, type = "cloglog", na.rm = TRUE)
                 
# Plot suitability
plot(suitability_raster, main = "MaxNet Habitat Suitability for Franklin's Bumblebee")
writeRaster(suitability_raster, 
            filename = "/Users/whitneymaxfield/Desktop/Future_Franklins_LC.tif", 
            overwrite = TRUE)
suit_df <- as.data.frame(suitability_raster, xy = TRUE, na.rm = TRUE)
colnames(suit_df)[3] <- "suitability"
                 
library(maps)
library(tigris)
# Get US states polygons as a dataframe
states_df <- map_data("state")
oregon_counties <- counties(state = "OR", cb = TRUE, class = "sf")
 
raster_crs <- crs(suitability_raster, proj = TRUE)  # should be EPSG:5070

# Transform the counties to the same CRS
oregon_counties_proj <- st_transform(oregon_counties, crs = raster_crs)

library(tigris)
library(sf)

# 1. Load California counties (using tigris)
california_counties <- counties(state = "CA", cb = TRUE, class = "sf")

# 2. Transform to raster CRS (assuming EPSG:5070, same as suit_df)
raster_crs <- crs(suitability_raster, proj = TRUE)
oregon_counties_proj <- st_transform(oregon_counties, crs = raster_crs)
california_counties_proj <- st_transform(california_counties, crs = raster_crs)

# 3. Combine Oregon + California counties
or_ca_counties <- rbind(oregon_counties_proj, california_counties_proj)

# 4. Plot combined counties with raster
ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = or_ca_counties, fill = NA, color = "white", size = 0.4) +
  scale_fill_viridis_c(option = "plasma") +
  coord_sf(crs = raster_crs, expand = FALSE) +
  theme_minimal() +
  labs(title = "Franklin's Bumblebee Habitat Suitability in OR + CA")

#just state lines: 
library(tigris)
library(sf)

# Download all US states polygons
states <- states(cb = TRUE, class = "sf")
target_states <- c("Oregon", "California", "Washington", "Idaho", "Nevada")  # add/remove as needed

states_subset <- states %>%
  filter(NAME %in% target_states)

raster_crs <- crs(suitability_raster, proj = TRUE)
states_proj <- st_transform(states_subset, crs = raster_crs)

ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = states_proj, fill = NA, color = "white", size = 0.8) +  # state borders only
  scale_fill_viridis_c(option = "plasma") +
  coord_sf(crs = raster_crs, expand = FALSE) +
  theme_minimal() +
  labs(title = "Franklin's Bumblebee Habitat Suitability with State Boundaries")
 
# Save as GeoTIFF
output_path <- "/Users/whitneymaxfield/Desktop/Future_Franklins_LC.tif"
writeRaster(suitability_raster, filename = output_path, overwrite = TRUE)
