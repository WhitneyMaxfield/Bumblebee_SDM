#==============================================================================
# Species Distribution Modeling with Landcover and Multicollinearity Filtering
#==============================================================================

# Dependencies
library(tidyverse)
library(dismo)
library(geodata)
library(maxnet)
library(sp)
library(raster)
library(terra)
library(usdm)  # For VIF

#==============================================================================
# Step 1: Load Occurrence Data 
#==============================================================================
shapefile_path <- "/Users/whitneymaxfield/Downloads/Franklins_historic_locations/Franklins_historic_locations.shp"           # Replace with actual path
franklini_occur <- vect(shapefile_path)

#==============================================================================
# Step 2: Load Climate and Landcover Data
#==============================================================================
# List .tif files and load them as SpatRaster (not raster::stack)
# Load climate layers as SpatRaster
# don't forget to CHANGE PATH (and make sure you have climate data downloaded, check data folder!)
climList <- list.files(
  path = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2/climate/wc2.1_2.5m",
  pattern = ".tif$",
  full.names = TRUE
)
clim <- rast(climList)  # Multi-layer SpatRaster

# Load and resample landcover
# Download NLCD landcover from git Data folder 
# CHANGE PATH 
landcover <- rast("/Users/whitneymaxfield/Downloads/Annual_NLCD_LndCov_1998_CU_C1V1/Annual_NLCD_LndCov_1998_CU_C1V1.tif")
landcover <- resample(landcover, clim[[1]], method = "near")

# Combine all into one multilayer SpatRaster
landcover <- project(landcover, clim, method = "near")
full_stack <- c(clim, landcover)

#==============================================================================
# Step 3: Multicollinearity Check with VIF
#==============================================================================
numeric_vars <- subset(full_stack, names(full_stack) != "landcover")
sample_df <- terra::spatSample(numeric_vars, size = 10000, method = "random", na.rm = TRUE, as.df = TRUE)
vif_result <- vifstep(sample_df, th = 10)
selected_vars <- vif_result@results$Variables
env_stack_filtered <- subset(full_stack, c(selected_vars, "landcover"))

#==============================================================================
# Step 4: Generate Pseudo-Absence Points
#==============================================================================
mask <- raster(clim[[1]])
geographicExtent <- extent(OccDataSpatialPts)
set.seed(7536)
backgroundPoints <- randomPoints(mask = mask, n = nrow(OccDataNotCoords), ext = geographicExtent, extf = 1.25, warn = 0)
colnames(backgroundPoints) <- c("longitude", "latitude")

#==============================================================================
# Step 5: Extract Environmental Data at Points
#==============================================================================
occEnv <- na.omit(raster::extract(env_stack_filtered, OccDataNotCoords))
absenceEnv <- na.omit(raster::extract(env_stack_filtered, backgroundPoints))
presenceAbsenceV <- c(rep(1, nrow(occEnv)), rep(0, nrow(absenceEnv)))
presenceAbsenceEnvDf <- as.data.frame(rbind(occEnv, absenceEnv))
presenceAbsenceEnvDf$landcover <- as.factor(presenceAbsenceEnvDf$landcover)

#==============================================================================
# Step 6: Build Maxnet Model
#==============================================================================
occSDM <- maxnet(
  p = presenceAbsenceV,
  data = presenceAbsenceEnvDf,
  f = maxnet.formula(presenceAbsenceV, presenceAbsenceEnvDf, classes = "default")
)

#==============================================================================
# Step 7: Predict Across Oregon
#==============================================================================
clim <- terra::rast(climList)
oregon <- geodata::gadm(country = "USA", level = 1, path = tempdir())
oregon_boundary <- oregon[oregon$NAME_1 == "Oregon", ]
oregon_boundary <- terra::project(oregon_boundary, clim)
geographicArea <- terra::crop(clim, oregon_boundary)
geographicArea <- terra::mask(geographicArea, oregon_boundary)
landcover <- rast("/path/to/your/landcover.tif")
landcover <- terra::resample(landcover, geographicArea[[1]], method = "near")
landcover <- terra::crop(landcover, oregon_boundary)
landcover <- terra::mask(landcover, oregon_boundary)
landcover <- as.factor(landcover)
final_stack <- c(geographicArea[[selected_vars]], landcover)
names(final_stack) <- names(presenceAbsenceEnvDf)
occpredictplot <- terra::predict(final_stack, model = occSDM, type = "cloglog", na.rm = TRUE)
writeRaster(occpredictplot, filename = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2/Occidentalis_pred_wVIF_LC.tif", overwrite = TRUE)

#==============================================================================
# Step 8: Plot Results with ggplot
#==============================================================================
occPredictDf <- as.data.frame(occpredictplot, xy = TRUE, na.rm = TRUE)
colnames(occPredictDf)[3] <- "layer"
wrld <- map_data("world")
buffer <- 0.2
xmax <- max(occPredictDf$x)
xmin <- min(occPredictDf$x)
ymax <- max(occPredictDf$y)
ymin <- min(occPredictDf$y)

# Plot
ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = occPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(xlim = c(xmin - buffer, xmax + buffer), ylim = c(ymin - buffer, ymax + buffer), expand = FALSE) +
  borders("state") +
  geom_point(data = OccDataNotCoords, aes(x = decimalLongitude, y = decimalLatitude), alpha = 0.5)
