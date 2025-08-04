#==============================================================================
#==============================================================================


# dependencies
install.packages("tidyverse")
library(tidyverse)

install.packages("dismo")
library(dismo)

install.packages("geodata")
library(geodata)

install.packages("maxnet")
library(maxnet)

install.packages("sp")
library(sp)

install.packages("raster")
library(raster)

install.packages("terra")
library(terra)
#sometimes R gets angry with terra and raster. If it does ask whitney or try chatgpt to problem solve it!

#==============================================================================
#==============================================================================
library(sf)
points_sf <- st_read("/Users/whitneymaxfield/Downloads/Franklins_historic_locations/Franklins_historic_locations.shp")
#turning into csv for future 
write.csv(points_sf, "/Users/whitneymaxfield/Downloads/Franklini_csv.csv", row.names = FALSE)
Frank_points <- points_sf
FrankDataNotCoords <- dplyr::select(Frank_points, LONG, LAT)

# convert to spatial points, necessary for modelling and mapping
FrankDataSpatialPts <- as(FrankDataNotCoords, "Spatial")

#==============================================================================
#==============================================================================
#wc2 <- worldclim_global(var = "bio", res = 2.5, path = "/Users/whitneymaxfield/Desktop/Bee_SDMS/Occidentalis")
# For franklins I want to use SAME climate data from the original landcover climate SDM (not loading new climate data) 
climList <- list.files(
  path = "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m",
  pattern = ".tif$",
  full.names = TRUE
)

#==============================================================================
#==============================================================================

# stacking the bioclim variables to process them at one go
clim <- raster::stack(climList)

plot(clim[[12]])
plot(FrankDataSpatialPts, add=T)
#should see a map on the right with a big mass of black points lol 

#==============================================================================
#
#             Section 2: Adding Pseudo-Absence Points 
# Create pseudo-absence points (making them up, using 'background' approach)
# first we need a raster layer to make the points up on, just picking 1
#
#==============================================================================

mask <- raster(clim[[1]]) 
# mask is the raster object that determines the area where we are generating pts

# determine geographic extent of our data (so we generate random points reasonably nearby)
geographicExtent <- extent(x = FrankDataSpatialPts)

# Random points for background (same number as our observed points we will use )
set.seed(7536) # seed set so we get the same background points each time we run this code! 

crs(mask)
crs(FrankDataSpatialPts)
#they have different CRS so need to reproject 
# Reproject points to match raster CRS
FrankDataSpatialPts_longlat <- spTransform(FrankDataSpatialPts, crs(mask))

# Now get extent in lon/lat
geographicExtent <- extent(FrankDataSpatialPts_longlat)

# Generate background points
backgroundPoints <- randomPoints(mask = mask, 
                                 n = length(FrankDataSpatialPts_longlat), 
                                 ext = geographicExtent, 
                                 extf = 1.25, 
                                 warn = 0)
# add col names (can click and see right now they are x and y)
colnames(backgroundPoints) <- c("longitude", "latitude")

#==============================================================================
#
#    Section 3: Collect Env Data and Point Data into Proper Model Formats
#   Data for observation sites (presence and background), with climate data
#
#==============================================================================

FrankEnv <- na.omit(raster::extract(x = clim, y = FrankDataNotCoords)) 
absenceEnv<- na.omit(raster::extract(x = clim, y = backgroundPoints)) # again, many NA values

# Create data frame with presence training data and backround points (0 = abs, 1 = pres)
presenceAbsenceV <- c(rep(1, nrow(FrankEnv)), rep(0, nrow(absenceEnv)))
presenceAbsenceEnvDf <- as.data.frame(rbind(FrankEnv, absenceEnv)) 

#==============================================================================
#
#                Section 4: Create SDM with Maxnet ### 
#
#==============================================================================


FrankSDM <- maxnet(
  p = presenceAbsenceV,
  data = presenceAbsenceEnvDf,
  f = maxnet.formula(presenceAbsenceV, presenceAbsenceEnvDf)
)

plot_maxnet_responses_to_png <- function(model, data, output_dir = "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Franklini no Landcover/Maxent_curves", type = "cloglog") {
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


plot_maxnet_responses_to_png(FrankSDM, presenceAbsenceEnvDf)

#==============================================================================
#
#                         Section 5: Plot the Model 
#   clim is huge and it isn't reasonable to predict over whole world!
#   first we will make it smaller
#
#==============================================================================

# Approximate bounding box (lon/lat)
ca_nv_or_extent <- raster::extent(-125, -114, 32, 48)
clim_ca_nv_or <- crop(clim, ca_nv_or_extent)

Frankpredictplot <- raster::predict(
  clim_ca_nv_or,
  model = FrankSDM,
  type = "cloglog"
)

plot(Frankpredictplot, main = "Maxnet Predicted Suitability for Bombus franklini")
png("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Franklini no Landcover/maxent_map_franklini.png", width = 800, height = 600)
plot(Frankpredictplot, main = "Maxnet Predicted Suitability")
dev.off()

# for ggplot, we need the prediction to be a data frame 
raster.spdf <- as(Frankpredictplot, "SpatialPixelsDataFrame")
FrankPredictDf <- as.data.frame(raster.spdf)

# plot in ggplot
wrld <- ggplot2::map_data("world")

xmin <- -125   # west coast (Pacific Ocean)
xmax <- -114   # east edge of Nevada
ymin <- 34     # southern California border
ymax <- 51     # northern Washington / southern BC border

ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = FrankPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) + 
  borders("state") + 
  geom_point(data = FrankDataNotCoords, aes(x = LONG, y = LAT), alpha = 0.5)

#look at that overlap! this means our prediction was good!
#can re-run to remove points and then save as .tif for QGIS: 

ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = FrankPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) + 
  borders("state") 

# Set output path — change to your desired folder
output_path <- "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Franklini_prediction_US.tif"

# Save the raster as GeoTIFF
writeRaster(Frankpredictplot, filename = output_path, format = "GTiff", overwrite = FALSE)



