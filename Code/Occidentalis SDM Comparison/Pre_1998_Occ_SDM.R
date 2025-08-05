#==============================================================================
#
#                    Welcome to species distrubution modeling! 
#           To run this code you need to use the RStudio desktop version 
#       While R studio Cloud (online version) is great, it won't work for this! (issues with Terra rn) 
#   If you need to download R studio visit https://posit.co/download/rstudio-desktop/ 
#
#==============================================================================

# dependencies
library(tidyverse)
library(dismo)
library(geodata)
library(maxnet)
library(sp)
library(raster)
library(terra)
#sometimes R gets angry with terra and raster. If it does ask whitney or try chatgpt to problem solve it!

#==============================================================================
# filtering using IUCN shp file for range of western bumblebee 
#==============================================================================
# Step 1: Load the shapefile (IUCN range map)
range_shp <- vect("/Users/whitneymaxfield/Downloads/redlist_species_data_26b44e70-1fb1-4d20-91ab-3d384a808d87")

# Step 2: Convert your tabular data to spatial points
# Use the pre-1998 data
occ_data_pre_1998 <- pre1998_occ_occur

# Convert to SpatVector
occ_spat <- vect(occ_data_pre_1998, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

# Step 3: Filter points INSIDE the IUCN range shapefile
occ_filtered <- occ_spat[range_shp, ]

# Step 5: Convert to data frame with coordinates
filtered_df <- as.data.frame(occ_filtered, geom = "xy")

# Step 6: Extract just coordinates and convert to SpatialPoints
OccDataNotCoords <- dplyr::select(filtered_df, decimalLongitude = x, decimalLatitude = y)

OccDataSpatialPts <- SpatialPoints(
  OccDataNotCoords,
  proj4string = CRS("+proj=longlat +datum=WGS84")
)

#==============================================================================

#climate data from worldclim v.2 for 1970-2000
clim <- worldclim_global(var = "bio", res = 2.5, version = 2.1, path = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2")
climList <- list.files(path = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2/climate/wc2.1_2.5m", pattern = ".tif", 
                       full.names = T)

clim <- raster::stack(climList)

plot(clim[[12]])
plot(OccDataSpatialPts, add=T)

#==============================================================================

#             Section 2: Adding Pseudo-Absence Points 
#==============================================================================

mask <- raster(clim[[1]]) 
geographicExtent <- extent(x = OccDataSpatialPts)

set.seed(7536) 
backgroundPoints <- randomPoints(mask = mask, 
                                 n = nrow(OccDataNotCoords), 
                                 ext = geographicExtent, 
                                 extf = 1.25,
                                 warn = 0) 

colnames(backgroundPoints) <- c("longitude", "latitude")

#==============================================================================

#    Section 3: Collate Env Data and Point Data into Proper Model Formats
#==============================================================================

occEnv <- na.omit(raster::extract(x = clim, y = OccDataNotCoords)) 
absenceEnv<- na.omit(raster::extract(x = clim, y = backgroundPoints)) 

presenceAbsenceV <- c(rep(1, nrow(occEnv)), rep(0, nrow(absenceEnv)))
presenceAbsenceEnvDf <- as.data.frame(rbind(occEnv, absenceEnv)) 

#==============================================================================

#                Section 4: Create SDM with Maxnet 
#==============================================================================

occSDM <- maxnet(
  p = presenceAbsenceV,
  data = presenceAbsenceEnvDf,
  f = maxnet.formula(presenceAbsenceV, presenceAbsenceEnvDf)
)

plot_maxnet_responses_to_png <- function(model, data, output_dir = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2/maxent_curves_just_OR_pre1998", type = "cloglog") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  vars <- names(data)
  for (var in vars) {
    xseq <- seq(min(data[[var]], na.rm=TRUE), max(data[[var]], na.rm=TRUE), length.out = 100)
    newdata <- data.frame(matrix(nrow=100, ncol=length(vars)))
    names(newdata) <- vars
    for (v in vars) {
      newdata[[v]] <- if (v == var) xseq else mean(data[[v]], na.rm=TRUE)
    }
    preds <- predict(model, newdata, type = type)
    png_filename <- file.path(output_dir, paste0("response_", var, ".png"))
    png(png_filename, width = 800, height = 600)
    plot(xseq, preds, type = "l", xlab = var, ylab = "Suitability", main = paste("Response:", var))
    dev.off()
  }
  message("Saved response plots to: ", normalizePath(output_dir))
}

plot_maxnet_responses_to_png(occSDM, presenceAbsenceEnvDf)

#==============================================================================

#                         Section 5: Plot the Model (Oregon only!)
#==============================================================================

# Load Oregon boundary as a SpatVector
clim <- terra::rast(climList)
oregon <- geodata::gadm(country = "USA", level = 1, path = tempdir())
oregon_boundary <- oregon[oregon$NAME_1 == "Oregon", ]

# Reproject Oregon to match the climate raster CRS
oregon_boundary <- terra::project(oregon_boundary, clim)

# Crop and mask the climate data to Oregon
geographicArea <- terra::crop(clim, oregon_boundary)
geographicArea <- terra::mask(geographicArea, oregon_boundary)

#rename so maxent doesn't get mad 
model_vars <- colnames(presenceAbsenceEnvDf)

# Check the current names of your prediction rasters
names(geographicArea)

# Rename the prediction rasters to match the model's expected names
names(geographicArea) <- model_vars

# Predict using the MaxNet model
occpredictplot <- terra::predict(
  geographicArea,
  model = occSDM,
  type = "cloglog",
  na.rm = TRUE
)

# Save the prediction raster as PNG
png("/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2/maxent_map_Oregon<1998_geofixed.png", width = 800, height = 600)
plot(occpredictplot, main = "Maxnet Predicted Suitability for Bombus occidentalis in Oregon")
dev.off()

# Convert prediction to a data frame for ggplot
occPredictDf <- as.data.frame(occpredictplot, xy = TRUE, na.rm = TRUE)
colnames(occPredictDf)
colnames(occPredictDf)[3] <- "layer"
# Plot with ggplot
wrld <- ggplot2::map_data("world")

xmax <- max(occPredictDf$x)
xmin <- min(occPredictDf$x)
ymax <- max(occPredictDf$y)
ymin <- min(occPredictDf$y)

buffer <- 0.2  # degrees; adjust as needed


ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = occPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(
    xlim = c(xmin - buffer, xmax + buffer), 
    ylim = c(ymin - buffer, ymax + buffer), 
    expand = FALSE
  ) + borders("state") + 
  geom_point(data = OccDataNotCoords, aes(x = decimalLongitude, y=decimalLatitude, alpha = 0.5)) 

# Version without points for export
ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = occPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(
    xlim = c(xmin - buffer, xmax + buffer), 
    ylim = c(ymin - buffer, ymax + buffer), 
    expand = FALSE
  ) + borders("state")

# Save as GeoTIFF (Oregon only)
output_path <- "/Users/whitneymaxfield/Desktop/Bee_SDMs/Westerns_past&future_oregon2/Occidentalis_<1998prediction_Oregon_geofixed.tif"
writeRaster(occpredictplot, filename = output_path, overwrite = TRUE)
