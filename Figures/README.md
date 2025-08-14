## Oregon: Western Bumblebee Habitat Before and After Population Decline

<img src="https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Pre%20and%20Post%201998%20Occidentalis/Portrait_bothMapsJPEG.jpeg" width="500">

### Methods Summary

Species distribution models for *Bombus occidentalis* in Oregon were developed using presence-background data and the **maxnet** package implementation of Maximum Entropy modeling (Phillips et al., 2006). Occurrence records were filtered using IUCN range boundaries and divided into pre-1998 and post-1998 datasets based on documented population decline timing. Environmental predictors consisted of WorldClim v2.1 bioclimatic variables at 2.5-minute resolution (Fick & Hijmans, 2017).

For each time period, pseudo-absence points were generated within the species' geographic extent using **randomPoints** from the **dismo** package. MaxEnt models were trained using presence-background datasets, with model predictions generated across Oregon using **terra** spatial operations (Hijmans, 2023). Oregon boundaries were obtained via **geodata** and used to crop and mask environmental layers to the study region.

Final habitat suitability maps display continuous probability values using a terrain color palette, with visualization conducted in **ggplot2** (Wickham, 2020). Models were saved as GeoTIFF format for further spatial analysis and comparison.

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

To enable fair spatial comparisons, the raster extents were intersected to define a common geographic extent, over which suitable habitat area was quantified using cell-specific area calculations (km²), accounting for spatial resolution and projection. Suitable habitat outside this overlapping extent was also evaluated to capture potential range shifts. Visualization of binarized suitability was performed using **ggplot2** to clearly illustrate spatial changes in habitat suitability between the two time periods (Wickham, 2020).

### Key Findings
- ADD HERE 

### References
- **Hijmans, R. J.** (2023). *terra: Spatial Data Analysis* (R package version 1.8-61) [Computer software]. https://github.com/rspatial/terra
- **Liu, C., White, M., & Newell, G.** (2013). Selecting thresholds for the prediction of species occurrence with presence-only data. *Journal of Biogeography, 40*(4), 778–789. https://doi.org/10.1111/jbi.12058
- **Pearson, R. G., Raxworthy, C. J., Nakamura, M., & Townsend Peterson, A.** (2006). Predicting species distributions from small numbers of occurrence records: A test case using cryptic geckos in Madagascar. *Journal of Biogeography, 34*(1), 102–117. https://doi.org/10.1111/j.1365-2699.2006.01594.x
- **Wickham, H.** (2020). *ggplot2: Elegant graphics for data analysis*. Springer-Verlag New York. https://ggplot2.tidyverse.org

