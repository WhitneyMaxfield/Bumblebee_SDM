#===============================
# Setup
#===============================
install.packages(c("terra", "raster", "tidyverse"), dependencies = TRUE)

library(terra)
library(raster)
library(tidyverse)

#===============================
# Load your rasters
#===============================
r_LC <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Franklini with Landcover/Franklins_historical.tif")     # SDM with Landcover
r_noLC <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Franklini no Landcover/Franklins_climateOnlyVIF.tif")         # SDM without Landcover

# Align extents, resolution, and NA handling (if needed)
r_LC <- resample(r_LC, r_noLC, method = "bilinear")  # Makes sure rasters align
r_LC[is.na(r_noLC)] <- NA                             # Mask areas not shared

#===============================
# 1. Plot habitat suitability
#===============================
plot(r_LC, main = "SDM with land cover")
plot(r_noLC, main = "SDM without land cover")

#===============================
# 2. Difference map
#===============================
r_diff <- r_LC - r_noLC
plot(r_diff, main = "Difference (land cover - No land cover)")

#===============================
# 3. Correlation between rasters
#===============================
cor_val <- cor(values(r_LC), values(r_noLC), use = "complete.obs")
print(paste("Pearson correlation:", round(cor_val, 3)))
#"Pearson correlation: 0.982"

#===============================
# 4. Schoener’s D (niche overlap)
#===============================
# Extract values
v1 <- values(r_LC)
v2 <- values(r_noLC)

# Remove NA values
valid <- complete.cases(v1, v2)
v1 <- v1[valid]
v2 <- v2[valid]

# Normalize to probability distributions
v1_norm <- v1 / sum(v1)
v2_norm <- v2 / sum(v2)

# Schoener’s D
D <- 1 - 0.5 * sum(abs(v1_norm - v2_norm))
print(paste("Schoener’s D:", round(D, 3)))
#"Schoener’s D: 0.9"



