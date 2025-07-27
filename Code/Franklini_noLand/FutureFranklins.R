# future SDM 
# using CMIP6

#==============================================================================
#==============================================================================

#install.packages("dismo")
library(dismo)

#install.packages("terra")
library(terra)

#install.packages("tidyverse")
library(tidyverse)

#==============================================================================
#==============================================================================

# first we are renaming our climate data from the last sdm because "clim" isn't very specific 
currentEnv <- clim 

#manual download attempt if this doesn't work see below 
futureEnv <- ("/Users/whitneymaxfield/Downloads/wc2.1_2.5m_bioc_IPSL-CM6A-LR_ssp245_2021-2040.tif")

#==============================================================================
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

predictExtent_ca_nv_or <- ext(-125, -114, 32, 48)
geographicAreaFutureC5 <- crop(futureEnv, predictExtent_ca_nv_or)
geographicAreaFutureC5_raster <- raster::stack(geographicAreaFutureC5)

FrankPredictPlotFutureC5 <- raster::predict(
  geographicAreaFutureC5_raster,
  model = FrankSDM,
  type = "cloglog"
)

#==============================================================================
#==============================================================================

raster.spdfFutureC5 <- as(FrankPredictPlotFutureC5, "SpatialPixelsDataFrame")
FrankPredictDfFutureC5 <- as.data.frame(raster.spdfFutureC5)
colnames(FrankPredictDfFutureC5) <- c("value", "x", "y")

# plot with ggplot
ggplot() +
  geom_tile(data = FrankPredictDfFutureC5, aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_c() +
  coord_fixed(xlim = c(-125, -114), ylim = c(32, 45)) +  # match crop extent
  theme_minimal()

###??
head(FrankPredictDfFutureC5)
ggplot() +
  geom_polygon(data = wrld, mapping = aes(x = long, y = lat, group = group),
               fill = "grey75") +
  geom_raster(data = FrankPredictDfFutureC5, aes(x = x, y = y, fill = value)) + 
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

# Set output path — change to your desired folder
output_path <- "/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Franklini_Futureprediction_US.tif"

# Save the raster as GeoTIFF
writeRaster(FrankPredictPlotFutureC5, filename = output_path, format = "GTiff", overwrite = FALSE)




