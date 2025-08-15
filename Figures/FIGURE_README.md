## Oregon: Western Bumblebee Habitat Before and After Population Decline

<img src="https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Pre%20and%20Post%201998%20Occidentalis/Portrait_bothMapsJPEG.jpeg" width="500">

### Methods Summary

Species distribution models for *Bombus occidentalis* in Oregon were developed using presence-background data and the **maxnet** package implementation of Maximum Entropy modeling (Phillips et al., 2006). Occurrence records were filtered using IUCN range boundaries and divided into pre-1998 and post-1998 datasets based on documented population decline timing. Environmental predictors consisted of WorldClim v2.1 bioclimatic variables at 2.5-minute resolution (Fick & Hijmans, 2017).

For each time period, pseudo-absence points were generated within the species' geographic extent using **randomPoints** from the **dismo** package. MaxEnt models were trained using presence-background datasets, with model predictions generated across Oregon using **terra** spatial operations (Hijmans, 2023). Oregon boundaries were obtained via **geodata** and used to crop and mask environmental layers to the study region.

Final habitat suitability maps display continuous probability values using a terrain color palette, with visualization conducted in **ggplot2** (Wickham, 2020). Models were saved as GeoTIFF format for further spatial analysis and comparison in QGIS. 

### Key Findings

- Historic model predicts widespread suitable habitat across Oregon, particularly in eastern Cascades and Blue Mountains
- Current model shows significant habitat contraction with reduced suitability in formerly optimal areas
- Spatial comparison reveals potential habitat loss corresponding to documented population declines since 1998

### References

- **Fick, S. E., & Hijmans, R. J.** (2017). WorldClim 2: New 1km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37(12), 4302–4315. https://doi.org/10.1002/joc.5086
- **Hijmans, R. J.** (2023). *terra: Spatial Data Analysis* (R package version 1.8-61) [Computer software]. https://github.com/rspatial/terra
- **Phillips, S. J., Anderson, R. P., & Schapire, R. E.** (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190(3-4), 231–259. https://doi.org/10.1016/j.ecolmodel.2005.03.026
- **Wickham, H.** (2020). *ggplot2: Elegant graphics for data analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org

---

## North America: Western Bumblebee Habitat Before and After Population Decline
<img src="https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Pre%20and%20Post%201998%20Occidentalis/North_America_suitability.jpeg" width="500">

### Methods Summary

Historical and current habitat suitability for *Bombus occidentalis* were analyzed by integrating presence point data with predicted suitability rasters using the **terra** package (Hijmans, 2023). Presence records were loaded as spatial vectors and aligned with their corresponding raster CRS to accurately extract predicted suitability values at known occurrence locations.

Suitability thresholds were calculated as the 10th percentile of suitability values at presence points for both historical and current periods, following recommended practices to balance omission and commission errors in species distribution modeling (Liu et al., 2013; Pearson et al., 2006). These thresholds were then used to binarize the continuous suitability rasters, distinguishing suitable from unsuitable habitat areas.

To enable fair spatial comparisons, the raster extents were intersected to define a common geographic extent, over which suitable habitat area was quantified using cell-specific area calculations (km²), accounting for spatial resolution and projection. Suitable habitat outside this overlapping extent was also evaluated to capture potential range shifts. Visualization of binarized suitability was performed using QGIS. 

### Key Findings
- Illustrates modeled habitat suitability for Bombus occidentalis across western North America during two time periods: before (top) and after (bottom) the species' widespread decline around 1998 (Xerces, 2015).
- Dark blue areas indicate regions classified as suitable habitat, while light gray represents unsuitable habitat.
- The historic model (pre-1998) shows a broader range of suitable habitat, which contracts significantly in the post-1998 model highlighting the species' range reduction over time.

### References
- **Hijmans, R. J.** (2023). *terra: Spatial Data Analysis* (R package version 1.8-61) [Computer software]. https://github.com/rspatial/terra
- **Liu, C., White, M., & Newell, G.** (2013). Selecting thresholds for the prediction of species occurrence with presence-only data. *Journal of Biogeography, 40*(4), 778–789. https://doi.org/10.1111/jbi.12058
- **Pearson, R. G., Raxworthy, C. J., Nakamura, M., & Townsend Peterson, A.** (2006). Predicting species distributions from small numbers of occurrence records: A test case using cryptic geckos in Madagascar. *Journal of Biogeography, 34*(1), 102–117. https://doi.org/10.1111/j.1365-2699.2006.01594.x

---

## Changes in Habitat Suitability for Western Bumblebee
<img src="https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Pre%20and%20Post%201998%20Occidentalis/Continous_Bocc.jpeg" width="500">

### Methods Summary

Habitat suitability for Bombus occidentalis was compared between historic and current periods by processing raster layers in R with the **terra** package (Hijmans, 2023). The historic raster was resampled to match the current raster’s resolution and extent, and suitability change was calculated by subtracting historic from current values, producing a continuous raster of change scores from –1 to +1 (Guisan & Zimmermann, 2000). For visualization, continuous change maps highlighted areas of gain and loss, with shifts smaller than ±0.1 masked as “no change.” Values above +0.1 indicated increases in suitability, while values below –0.1 indicated decreases—thresholds chosen to capture meaningful changes while remaining sensitive to subtle shifts.

To quantify total suitable habitat, we applied a 10th percentile training presence threshold to binarize suitability maps (Pearson et al., 2006; Liu et al., 2013) and calculated area by multiplying suitable cells by their per-cell area (Elith et al., 2011; Guisan & Zimmermann, 2000). Thematic maps showing both continuous and classified change were created with tmap (Tennekes, 2018), overlaid with county and state boundaries from **tigris** (Walker, 2023), and finalized in QGIS.

### Key Findings

- Extensive white regions stretching through the central part of the state indicate areas where habitat suitability remained largely stable between the pre-1998 and post-1998 models.
- Southwestern Oregon shows notable decreases in habitat suitability

### References

- **Elith, J., Phillips, S. J., Hastie, T., Dudík, M., Chee, Y. E., & Yates, C. J.** (2011). A statistical explanation of MaxEnt for ecologists. *Diversity and Distributions*, 17(1), 43–57. https://doi.org/10.1111/j.1472-4642.2010.00725.x
- **Guisan, A., & Zimmermann, N. E.** (2000). Predictive habitat distribution models in ecology. *Ecological Modelling*, 135(2-3), 147–186. https://doi.org/10.1016/s0304-3800(00)00354-9
- **Hijmans, R. J.** (2023). *terra: Spatial Data Analysis* (R package version 1.8-61) [Computer software]. https://github.com/rspatial/terra
- **Liu, C., White, M., & Newell, G.** (2013). Selecting thresholds for the prediction of species occurrence with presence-only data. *Journal of Biogeography*, 40(4), 778–789. https://doi.org/10.1111/jbi.12058
- **Pearson, R. G., Raxworthy, C. J., Nakamura, M., & Townsend Peterson, A.** (2006). Predicting species distributions from small numbers of occurrence records: A test case using cryptic geckos in Madagascar. *Journal of Biogeography*, 34(1), 102–117. https://doi.org/10.1111/j.1365-2699.2006.01594.x
- **Tennekes, M.** (2018). tmap: Thematic maps in R. *Journal of Statistical Software*, 84(6). https://doi.org/10.18637/jss.v084.i06
- **Tronstad, L. M., Bell, C., Cook, K., & Dillon, M. E.** (2024). Using species distribution models to assess the status of the declining Western bumble bee (Hymenoptera: Apidae: Bombus occidentalis) in Wyoming, USA. *Environments*, 12(1), 2. https://doi.org/10.3390/environments12010002

--- 

## Habitat Suitability Models for Franklin's Bumblebee

<img src="https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Franklini%20with%20Land%20Cover/Franklini_Historic_withLand.jpeg" width="500">

### Methods Summary

Environmental and land cover data for Oregon, California, and Nevada were processed using the **sf** and **terra** packages in R to define the study area, reproject, crop, and mask raster layers (Hijmans, 2023; Pebesma, 2023). State boundaries were transformed to EPSG:5070 projection to match the National Land Cover Database.

The 1998 NLCD land cover layer was reclassified to a floral habitat proxy by identifying open habitat classes (Shrub/Scrub, Grassland/Herbaceous, Pasture/Hay) representing potential bumblebee foraging areas. Climate variables from WorldClim v2.1 (Fick & Hijmans, 2017) were combined with landcover layers, reprojected and resampled, then masked to the three-state area of interest.

Multicollinearity among predictors was assessed and reduced using variance inflation factor (VIF) procedures with the **usdm** package (Naimi, 2017). Species distribution modeling was conducted with the **maxnet** package, implementing a Maxent approach with 5000 background points randomly sampled across environmentally suitable areas (Phillips et al., 2006; Breiner et al., 2018). Resulting habitat suitability maps were visualized in **ggplot2**, overlaid with county and state boundaries from **tigris** (Walker, 2023), and then exported for final map production in QGIS. 

Code for generating this map was also adapted to produce comparison maps for two additional bumblebee species—B. vosnesenskii and B. occidentalis—by substituting the occurrence data while keeping the climate and land cover variables constant. The same code was further modified to create maps assessing the significance of including land cover data and the effects of VIF filtering. These outputs, available in the Figures folder (with corresponding scripts in the Code folder), were not included here due to their redundancy.

### Key Findings

- Model predicts the highest suitability in the Klamath-Siskiyou region, particularly around Mt. Ashland in southwestern Oregon
- When comparing various years, overall spatial distribution remains largely unchanged suggesting stable environmental conditions in the species' preferred  habitat.

### References

- **Breiner, F. T., Nobis, M. P., Bergamini, A., & Guisan, A.** (2018). Optimizing ensembles of small models for predicting the distribution of species with few occurrences. *Methods in Ecology and Evolution*, 9(4), 802–808. https://doi.org/10.1111/2041-210x.12957
- **Fick, S. E., & Hijmans, R. J.** (2017). WorldClim 2: New 1 km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37(12), 4302–4315. https://doi.org/10.1002/joc.5086
- **Hijmans, R. J.** (2023). *terra: Spatial Data Analysis* (R package version 1.8-61) [Computer software]. https://github.com/rspatial/terra
- **Naimi, B., Hamm, N. A., Groen, T. A., Skidmore, A. K., & Toxopeus, A. G.** (2014). Where is positional uncertainty a problem for species distribution modelling? *Ecography*, 37, 191–203. https://doi.org/10.1111/j.1600-0587.2013.00205.x
- **Pebesma, E., & Bivand, R.** (2023). *Spatial data science: With applications in R*. Chapman and Hall/CRC. https://r-spatial.org/book/
- **Phillips, S. J., Anderson, R. P., & Schapire, R. E.** (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190(3-4), 231–259. https://doi.org/10.1016/j.ecolmodel.2005.03.026
- **Walker, K.** (2023). *tigris: Load Census TIGER/Line shapefiles* (R package version 2.2) [Computer software]. https://CRAN.R-project.org/package=tigris
- **Hijmans, R. J., Phillips, S., Leathwick, J., & Elith, J.** (2024). dismo: Species Distribution Modeling (R package, version 1.3-15). https://github.com/rspatial/dismo
- **Hijmans, R. J.** (2025). raster: Geographic Data Analysis and Modeling (R package, version 3.6-32). https://rspatial.org/raster

