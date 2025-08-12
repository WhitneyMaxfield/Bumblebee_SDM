# --- Load packages ---
library(terra)
library(ecospat)
library(ade4)
library(dplyr)

# -------------------------
# 1. Load data
# -------------------------
# Predictor stack (same variables used for SDMs)
pred_stack <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/predictor_stack.tif") 

# Occurrence CSVs
vos_pts <- Vosnesenskii_occurence
west_pts <- Western_to2006

# -------------------------
# 2. PRESENCE-ONLY NICHE OVERLAP
# -------------------------
# Create SpatVectors of occurrences (make sure CRS matches pred_stack)
vos_vect <- vect(vos_pts, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(pred_stack))
west_vect <- vect(west_pts, geom = c("decimalLongitude", "decimalLatitude"), crs = crs(pred_stack))

# Extract environmental values at occurrence points
vos_env <- extract(pred_stack, vos_vect)
west_env <- extract(pred_stack, west_vect)

# Combine for PCA 
combined_env <- rbind(vos_env, west_env)

# Remove ID and categorical columns
combined_env_clean <- combined_env %>%
  dplyr::select(-ID, -`NLCD Land Cover Class`) %>%  # remove non-numeric columns
  dplyr::select(where(is.numeric)) %>%             # keep only numeric columns
  na.omit()

# Run PCA on combined environmental variables
pca_env <- dudi.pca(combined_env_clean, scannf = FALSE, nf = 2)

# Clean supplementary data to match PCA input columns
vos_env_clean <- vos_env %>%
  dplyr::select(-ID, -`NLCD Land Cover Class`) %>%
  dplyr::select(where(is.numeric)) %>%
  na.omit()

west_env_clean <- west_env %>%
  dplyr::select(-ID, -`NLCD Land Cover Class`) %>%
  dplyr::select(where(is.numeric)) %>%
  na.omit()

# Project supplementary points on PCA space
vos_scores <- suprow(pca_env, vos_env_clean)$li
west_scores <- suprow(pca_env, west_env_clean)$li

# Create PCA grids for niche overlap
vos_grid <- ecospat.grid.clim.dyn(glob = pca_env$li,
                                  glob1 = pca_env$li,
                                  sp = vos_scores,
                                  R = 100)

west_grid <- ecospat.grid.clim.dyn(glob = pca_env$li,
                                   glob1 = pca_env$li,
                                   sp = west_scores,
                                   R = 100)

# Calculate overlap metrics
presence_D <- ecospat.niche.overlap(vos_grid, west_grid, cor = TRUE)
presence_equiv <- ecospat.niche.equivalency.test(vos_grid, west_grid, rep = 100)
presence_sim <- ecospat.niche.similarity.test(vos_grid, west_grid, rep = 100)

# Plot presence-only niche overlap
ecospat.plot.niche(vos_grid, west_grid, name.axis1 = "PCA1", name.axis2 = "PCA2")

# -------------------------
# 3. MODEL-BASED NICHE OVERLAP
# -------------------------
# Load modeled suitability rasters
vos_rast <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Vosnesenskii with Landcover/Vos_2006.tif")
west_rast <- rast("/Users/whitneymaxfield/Desktop/Git_attempt_SDM/Bumblebee_SDM/Figures/Western with Land Cover/Western_2006.tif")

# Convert rasters to data.frames
vos_df <- as.data.frame(vos_rast, xy = TRUE)
west_df <- as.data.frame(west_rast, xy = TRUE)

# Merge on coordinates to keep common cells
merge_df <- merge(vos_df, west_df, by = c("x", "y"), suffixes = c("_vos", "_west"))

# Remove NAs
merge_df <- na.omit(merge_df)

# PCA on modeled suitability values
model_env <- merge_df[, c("lyr1_vos", "lyr1_west")]
pca_model <- dudi.pca(model_env, scannf = FALSE, nf = 2)

# Scores weighted by suitability
vos_scores_model <- suprow(pca_model, model_env)$li
west_scores_model <- vos_scores_model  # same cells, can adjust if needed

# Create PCA grids for model niche overlap
vos_grid_model <- ecospat.grid.clim.dyn(glob = pca_model$li,
                                        glob1 = pca_model$li,
                                        sp = vos_scores_model,
                                        R = 100)

west_grid_model <- ecospat.grid.clim.dyn(glob = pca_model$li,
                                         glob1 = pca_model$li,
                                         sp = west_scores_model,
                                         R = 100)

# Overlap metrics
model_D <- ecospat.niche.overlap(vos_grid_model, west_grid_model, cor = TRUE)
model_equiv <- ecospat.niche.equivalency.test(vos_grid_model, west_grid_model, rep = 100)
model_sim <- ecospat.niche.similarity.test(vos_grid_model, west_grid_model, rep = 100)

# Plot model-based niche overlap
ecospat.plot.niche(vos_grid_model, west_grid_model,
                   name.axis1 = "PCA1", name.axis2 = "PCA2")

# -------------------------
# 4. Output results
# -------------------------
cat("Presence-only Schoener's D:", presence_D$D, "\n")
print(presence_equiv)
print(presence_sim)

cat("Model-based Schoener's D:", model_D$D, "\n")
print(model_equiv)
print(model_sim)
