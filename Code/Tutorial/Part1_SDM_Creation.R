#==============================================================================
#
#                    Welcome to species distrubution modeling! 
#           To run this code you need to use the RStudio desktop version 
#       While R studio Cloud (online version) is great, it won't work for this! (issues with Terra rn) 
#   If you need to download R studio visit https://posit.co/download/rstudio-desktop/ 
#
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

#   Section 1: Obtaining and Formatting Occurrence / Climate Data 
#
# Download the occidentalis_data file from the "data" folder on the github 
# Add this file to the environment by selecting "import dataset" in the top right Environment box
# Select "from text (base)" then find the .csv file in your downloads 
# Make sure "heading" is selected to YES when adding! 

#==============================================================================

occ_data <- Occidentalis_data
OccDataNotCoords <- dplyr::select(occ_data, decimalLongitude, decimalLatitude)

# convert to spatial points, necessary for modelling and mapping
OccDataSpatialPts <- SpatialPoints(OccDataNotCoords , proj4string = CRS("+proj=longlat"))

#==============================================================================
#
# before running the code below we need to sub in our own path info! 
#
# set your own path (a folder you want this to save to!)
# To set the path: 
# 1. create a folder in your desktop and name it
# 2. Right click the folder and select "get info" 
# 3. in the pop up window copy everything after "where" 
# (it should be something like Macintosh... > Users > ...)
# 4. paste everything in the "where" line into the path = "" below! 
# 5. dont forget the "" around the text and the closing parentheses! 
#
#==============================================================================
wc2 <- worldclim_global(var = "bio", res = 2.5, path = "/Users/whitneymaxfield/Desktop/Bee_SDMS/Occidentalis")
#res can be changed depending on what you are doing (can google this or ask whit for more info)

#==============================================================================
#                           *Note on Rasters* 
#
# raster consists of a matrix of cells (or pixels) organized into rows and columns
# (or a grid) where each cell contains a value representing information, such as 
# temperature. Rasters are digital aerial photographs, imagery from satellites, 
# digital pictures, or even scanned maps.
#
#==============================================================================

#use the same path we just set as above!!!! 
### for whitney: old worldclim used .bif but new uses .tif!

climList <- list.files(path = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Occidentalis/climate/wc2.1_2.5m", pattern = ".tif", 
                       full.names = T)  # '..' leads to the path above the folder where the .rmd file is located

# stacking the bioclim variables to process them at one go
clim <- raster::stack(climList)

plot(clim[[12]])
plot(OccDataSpatialPts, add=T)
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
geographicExtent <- extent(x = OccDataSpatialPts)

# Random points for background (same number as our observed points we will use )
set.seed(7536) # seed set so we get the same background points each time we run this code! 
backgroundPoints <- randomPoints(mask = mask, 
                                 n = nrow(OccDataNotCoords), # n should be same n as in the pts to be used to test
                                 ext = geographicExtent, 
                                 extf = 1.25, # draw a slightly larger area than where our sp was found (ask katy what is appropriate here)
                                 warn = 0) # don't complain about not having a coordinate reference system

# add col names (can click and see right now they are x and y)
colnames(backgroundPoints) <- c("longitude", "latitude")

#==============================================================================
#
#    Section 3: Collate Env Data and Point Data into Proper Model Formats
#   Data for observation sites (presence and background), with climate data
#
#==============================================================================

occEnv <- na.omit(raster::extract(x = clim, y = OccDataNotCoords)) 
absenceEnv<- na.omit(raster::extract(x = clim, y = backgroundPoints)) # again, many NA values

# Create data frame with presence training data and backround points (0 = abs, 1 = pres)
presenceAbsenceV <- c(rep(1, nrow(occEnv)), rep(0, nrow(absenceEnv)))
presenceAbsenceEnvDf <- as.data.frame(rbind(occEnv, absenceEnv)) 

#==============================================================================
#
#                Section 4: Create SDM with Maxnet ### 
#
#==============================================================================


occSDM <- maxnet(
  p = presenceAbsenceV,
  data = presenceAbsenceEnvDf,
  f = maxnet.formula(presenceAbsenceV, presenceAbsenceEnvDf)
)
#THIS CAN TAKE A LONG TIME (like a few minutes!!) depending on how many occurences/layers 
#response(occSDM) but tailored for maxnet:
# !!!!you do not need to understand what this is doing lol it is very complicated!!!!
#do need to change output/dir to where you want your curves saved :P

plot_maxnet_responses_to_png <- function(model, data, output_dir = "/Users/whitneymaxfield/Desktop/Bee_SDMs/Occidentalis/maxent_curves", type = "cloglog") {
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

# Run it; 
#note we had to save these plots to a folder on your computer because they are too large for R to show us!
#to view the plot travel to the location stated after the "saved response plots to:..." text
plot_maxnet_responses_to_png(occSDM, presenceAbsenceEnvDf)

#==============================================================================
#
#                         Section 5: Plot the Model 
#   clim is huge and it isn't reasonable to predict over whole world!
#   first we will make it smaller
#
#==============================================================================

predictExtent <- 1.25 * geographicExtent
geographicArea <- crop(clim,predictExtent) #crop clim data to predict extent 

# Make sure geographicArea is a RasterStack or RasterBrick
class(geographicArea)
# Should return RasterStack or RasterBrick

# Now predict across the raster area using the maxnet model
occpredictplot <- raster::predict(
  geographicArea,      # RasterStack or RasterBrick
  model = occSDM,       # maxnet model
  type = "cloglog"      # Match training type
)

# Plot the result!
#the plots are too big so we have to save them in a png, we can view them in the files area!

#change where you want your png to save! can also manually export in the plots area to the right. 
plot(occpredictplot, main = "Maxnet Predicted Suitability for Bombus occidentalis")
png("/Users/whitneymaxfield/Desktop/Bee_SDMs/Occidentalis/maxent_map.png", width = 800, height = 600)
plot(occpredictplot, main = "Maxnet Predicted Suitability")
dev.off()

# for ggplot, we need the prediction to be a data frame 
raster.spdf <- as(occpredictplot, "SpatialPixelsDataFrame")
occPredictDf <- as.data.frame(raster.spdf)

# plot in ggplot
wrld <- ggplot2::map_data("world")

xmax <- max(occPredictDf$x)
xmin <- min(occPredictDf$x)
ymax <- max(occPredictDf$y)
ymin <- min(occPredictDf$y)

ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = occPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) + borders("state") + 
  geom_point(data = OccDataNotCoords, aes(x = decimalLongitude, y=decimalLatitude, alpha = 0.5)) 

#look at that overlap! this means our prediction was good!
#can re-run to remove points and then save as .tif for QGIS: 

ggplot() + 
  geom_polygon(data = wrld, aes(x = long, y = lat, group = group), fill = "grey60") +
  geom_raster(data = occPredictDf, aes(x = x, y = y, fill = layer)) +
  scale_fill_gradientn(colors = terrain.colors(10, rev = TRUE)) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) + borders("state")

# Set output path — change to your desired folder
output_path <- "/Users/whitneymaxfield/Desktop/Bee_SDMs/Occidentalis/Occidentalis_prediction_US.tif"

# Save the raster as GeoTIFF
writeRaster(occpredictplot, filename = output_path, format = "GTiff", overwrite = FALSE)



