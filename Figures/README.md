## North American Historical vs. Current Habitat Suitability Analysis
![Image](https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Pre%20and%20Post%201998%20Occidentalis/North_America_suitability.jpeg|width=100)
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

