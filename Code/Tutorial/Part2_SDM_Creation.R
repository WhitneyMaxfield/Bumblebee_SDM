# future SDM 
# using CMIP6

#==============================================================================
#                        Welcom to Future SDM! 
#                This code should be run after SDM.R! 
#   Everything from SDM.R should still be in your upper right environemnt 
#     if it is NOT, you will need to re-run SDM.R and then do this code. 
#      you can open both R files in the same window by simply clicking 
#     on the other .R file you want to open while one file is loaded into R

#==============================================================================

# Install packages & load libraries IF they aren't yet part of your project environment
# They may have already been loaded in sdm_maxent_lesson.R. 
# it is good to check by running JUST the "library" function 

#install.packages("dismo")
library(dismo)

#install.packages("terra")
library(terra)

#install.packages("tidyverse")
library(tidyverse)

#might not need this one!
#install.packages("remotes")
#library(remotes)
#remotes::install_github("rspatial/geodata") #force download new geodata package


#==============================================================================
#
#                      Section 1: Obtaining future climate data 
#   This data is predicted based on current emissions and env. data (big maths lol)
#
#==============================================================================

# first we are renaming our climate data from the last sdm because "clim" isn't very specific 
currentEnv <- clim 

#manual download attempt if this doesn't work see below 
futureEnv <- cmip6_world(
  var = "bio",
  res = 2.5,
  model = "IPSL-CM5A-LR",
  ssp = 45,
  time = 70,
  path = "data"
)

#==============================================================================
#
#with manual download: if above code worked SKIP THIS: 
#if the code above DID NOT WORK you will need to manually download the climate data. 
#to do this visit: https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html 
#at this link scroll down to "2081-2100" find IPSL-CM6A-LR and click on the bc option 
#like we did in the SDM.R file, you will need to change the path below to wherever this new data ends up! 
#
#==============================================================================
data_path <- "/Users/whitneymaxfield/Downloads/wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2081-2100.tif"  # or wherever you saved the files
clim_files <- list.files(data_path, pattern = "\\.tif$", full.names = TRUE)

# Check how many layers it contains
futureEnv <- rast("/Users/whitneymaxfield/Downloads/wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2081-2100.tif")

#==============================================================================
#                         ***NOTE***
#   when the getdata function starts working again
#   you can specific with rcp, which is the greenhouse gas emission prediction,
#   can also specific which model https://rdrr.io/cran/raster/man/getData.html
#   and which year: 50 or 70 years from now
#==============================================================================

#making sure the names match up
names(futureEnv) = names(currentEnv)
# look at current vs future climate vars
plot(currentEnv[[1]])
plot(futureEnv[[1]])

predictExtent_terra <- ext(predictExtent) 

## note for Whitney: predictExtent is defined in SDM.R
# crop clim to the extent of the map you want
geographicAreaFutureC5 <- crop(futureEnv, predictExtent_terra)

geographicAreaFutureC5_raster <- stack(geographicAreaFutureC5)  # convert to RasterStack
occPredictPlotFutureC5 <- raster::predict(geographicAreaFutureC5_raster, occSDM)

#==============================================================================
#
#     Step 6: Convert for ggplot (convert raster prediction to data frame)
#
#==============================================================================

raster.spdfFutureC5 <- as(occPredictPlotFutureC5, "SpatialPixelsDataFrame")
occPredictDfFutureC5 <- as.data.frame(raster.spdfFutureC5)
colnames(occPredictDfFutureC5) <- c("value", "x", "y")

# plot with ggplot
ggplot() +
  geom_tile(data = occPredictDfFutureC5, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal()

# get our lat/lon boundaries
xmax <- max(occPredictDfFutureC5$x)+10
xmin <- min(occPredictDfFutureC5$x)-10
ymax <- max(occPredictDfFutureC5$y)+2
ymin <- min(occPredictDfFutureC5$y)-2

head(occPredictDfFutureC5)
ggplot() +
  geom_polygon(data = wrld, mapping = aes(x = long, y = lat, group = group),
               fill = "grey75") +
  geom_raster(data = occPredictDfFutureC5, aes(x = x, y = y, fill = value)) + 
  scale_fill_gradientn(colors = terrain.colors(10, rev = T)) +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = F) +
  scale_size_area() +
  borders("world") +
  borders("state") +
  labs(title = "SDM of ___ Under CMIP 5 Climate Conditions",
       x = "longitude",
       y = "latitude",
       fill = "Env Suitability") +
  theme(legend.box.background=element_rect(),legend.box.margin=margin(5,5,5,5)) 


#if you wnat to save this map we can use the ggsave function, dont forget to change the path!
ggsave(filename = "futureSDM.jpg", plot=last_plot(),path="/Users/whitneymaxfield/Desktop/ENM tutorials", 
       width=2000, height=900, units="px")

