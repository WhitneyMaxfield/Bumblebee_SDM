library(terra)

# Choose a template raster (one of the WorldClim BIO variables)
template <- rast("/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif")

# List your predictor files in desired order
files <- c(
  "/Users/whitneymaxfield/Downloads/NLCD_JtdPwMwal77VH8CTmKxZ/Annual_NLCD_LndCov_1992_CU_C1V1_JtdPwMwal77VH8CTmKxZ.tiff",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_3.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_2.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_9.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_8.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_13.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_14.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif",
  "/Users/whitneymaxfield/Desktop/Franklini_climate_data_old/climate/wc2.1_2.5m/wc2.1_2.5m_bio_15.tif"
)

# Align each raster to the template
aligned_files <- lapply(files, function(f) {
  r <- rast(f)
  r <- project(r, template)      # ensure projection matches
  r <- resample(r, template)     # ensure resolution & extent match
  return(r)
})

# Stack them
pred_stack <- rast(aligned_files)

# Save
writeRaster(pred_stack, "/Users/whitneymaxfield/Desktop/predictor_stack.tif", overwrite = TRUE)


