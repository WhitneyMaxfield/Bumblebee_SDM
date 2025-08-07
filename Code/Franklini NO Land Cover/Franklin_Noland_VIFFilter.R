# Install required packages if not already installed
install.packages(c("sf", "terra", "tidyverse", "dismo", "geodata", "maxnet", "sp", "raster", "corrplot", "usdm", "tigris"))

# Load libraries
library(terra)
library(sf)
library(tidyverse)
library(dismo)
library(geodata)
library(maxnet)
library(sp)
library(raster)
library(corrplot)
library(usdm)
library(tigris)

# ----------------------------------------
# 1. Define your study area (OR, CA, NV)
# ----------------------------------------
oregon <- st_read("/Users/whitneymaxfield/Downloads/Oregon/tl_2019_41_cousub.shp")
california <- st_read("/Users/whitneymaxfield/Downloads/Cali_shp/tl_2019_06_cousub.shp")
nevada <- st_read("/Users/whitneymaxfield/Downloads/Nevada_shape/tl_2019_32_cousub.shp")

crs_target <- "EPSG:5070"
oregon <- st_transform(oregon, crs = crs_target)
california <- st_transform(california, crs = crs_target)
nevada <- st_transform(nevada, crs = crs_target)

or_ca_nv_combined <- rbind(oregon, california, nevada)
or_ca_nv_vect <- vect(or_ca_nv_combined)

# ----------------------------------------
# 2. Load Climate Variables
# ----------------------------------------
clim <- worldclim_global(var = "bio", res = 2.5, version = 2.1, path = "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old")
climList <- list.files(
  path = "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m",
  pattern = ".tif$",
  full.names = TRUE
)
clim_raster <- raster::stack(climList)
clim_raster_spat <- rast(clim_raster)

# Reproject AOI to match climate CRS
or_ca_nv_vect <- project(or_ca_nv_vect, crs(clim_raster_spat))

# Crop and mask climate data to AOI
env_stack_cropped <- crop(clim_raster_spat, or_ca_nv_vect)
env_stack_masked <- mask(env_stack_cropped, or_ca_nv_vect)

# Sample environmental space
sample_points <- spatSample(env_stack_masked, size = 2000, method = "random", na.rm = TRUE)
env_df <- as.data.frame(sample_points, na.rm = TRUE)

# ----------------------------------------
# 3. Multicollinearity filtering
# ----------------------------------------
cor_matrix <- cor(env_df, use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.4)

vif_result <- vifstep(env_df, th = 10)
selected_vars <- vif_result@results$Variables
env_stack_filtered <- env_stack_masked[[selected_vars]]

#variables retained: 
print(vif_result)

#variables removed:
print(vif_result@excluded)

# ----------------------------------------
# 4. SDM with maxnet (climate only)
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

# Combine data
env_data <- rbind(presence_env, background_env)
env_data <- env_data[complete.cases(env_data), ]

response <- env_data$presence
predictors <- env_data[, setdiff(names(env_data), "presence")]
predictors_clean <- predictors[, !names(predictors) %in% c("ID", "presence")]
colnames(predictors_clean) <- make.names(colnames(predictors_clean))

maxnet_model <- maxnet(p = response, data = predictors_clean, f = maxnet.formula(response, predictors_clean))
suitability_raster <- terra::predict(env_stack_filtered, maxnet_model, type = "cloglog", na.rm = TRUE)

# Save output
writeRaster(suitability_raster, "/Users/whitneymaxfield/Desktop/Franklins_climateOnly.tif", overwrite = TRUE)

# Convert to dataframe for plotting
suit_df <- as.data.frame(suitability_raster, xy = TRUE, na.rm = TRUE)
colnames(suit_df)[3] <- "suitability"

# ----------------------------------------
# 5. Plot with Counties and States
# ----------------------------------------
oregon_counties <- counties(state = "OR", cb = TRUE, class = "sf")
california_counties <- counties(state = "CA", cb = TRUE, class = "sf")
nevada_counties <- counties(state = "NV", cb = TRUE, class = "sf")
raster_crs <- crs(suitability_raster, proj = TRUE)

or_ca_nv_counties <- rbind(
  st_transform(oregon_counties, crs = raster_crs),
  st_transform(california_counties, crs = raster_crs),
  st_transform(nevada_counties, crs = raster_crs)
)

ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = or_ca_nv_counties, fill = NA, color = "white", size = 0.4) +
  scale_fill_viridis_c(option = "plasma") +
  coord_sf(crs = raster_crs, expand = FALSE) +
  theme_minimal() +
  labs(title = "Franklin's Bumblebee Habitat Suitability (Climate Only)")

# State boundaries
states <- states(cb = TRUE, class = "sf")
target_states <- c("Oregon", "California", "Washington", "Idaho", "Nevada")
states_proj <- st_transform(filter(states, NAME %in% target_states), crs = raster_crs)

ggplot() +
  geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
  geom_sf(data = states_proj, fill = NA, color = "white", size = 0.8) +
  scale_fill_viridis_c(option = "plasma") +
  coord_sf(crs = raster_crs, expand = FALSE) +
  theme_minimal() +
  labs(title = "Franklin's Bumblebee Suitability (Climate Only, State Borders)")

