## Oregon: Western Bumblebee Habitat Before and After Population Decline

<img src="https://github.com/WhitneyMaxfield/Bumblebee_SDM/blob/main/Figures/Pre%20and%20Post%201998%20Occidentalis/Portrait_bothMapsJPEG.jpeg" width="500">

### What This Figure Shows

This figure compares modeled habitat suitability for *Bombus occidentalis* (Western bumblebee) in Oregon across two critical time periods: before and after the species experienced a sharp population decline around 1998. The color gradient represents relative habitat suitability, with dark purple indicating high-quality habitat and light yellow showing areas of low suitability. County boundaries provide spatial reference points across Oregon.

### Our Modeling Approach

We created separate species distribution models for each time period using **MaxEnt** (Maximum Entropy) modeling implemented through the **maxnet** R package. Here's how we built these models:

**Data Preparation:** We filtered occurrence records using IUCN range maps to focus on the species' known geographic range, then split the data into pre-1998 and post-1998 periods. For each time period, we used WorldClim v2.1 bioclimatic variables as environmental predictors.

**Model Training:** We generated pseudo-absence points (background points) within the species' range to create presence-background datasets. The MaxEnt algorithm compared environmental conditions at known occurrence sites with randomly selected background locations to identify suitable habitat conditions.

**Spatial Predictions:** Both models were trained using the same environmental variables and then projected across Oregon's landscape to create continuous habitat suitability maps. We focused specifically on Oregon to examine regional patterns in this ecologically important state.

### What the Results Tell Us

The comparison reveals a dramatic story of habitat change:

**Historic Period (top map):** The model predicts widespread suitable habitat across Oregon, with particularly high suitability in the eastern Cascades, Blue Mountains, and portions of southwestern Oregon. This suggests the species historically had access to extensive, well-connected habitat.

**Current Period (bottom map):** The model shows a significant contraction in suitable habitat quality. Much of the previously suitable range now shows reduced suitability or has become unsuitable altogether, potentially explaining the species' documented population decline.

This pattern aligns with known conservation concerns for Western bumblebees, which have experienced severe population declines across their range since the late 1990s.

### Methods Details

**Occurrence Data:** We used spatially filtered occurrence records within IUCN range boundaries, split into pre-1998 and post-1998 datasets based on the documented population decline timing.

**Environmental Variables:** WorldClim v2.1 bioclimatic variables at 2.5-minute resolution provided climate data for model training.

**Modeling Algorithm:** MaxEnt modeling via the **maxnet** package, using presence-background data with pseudo-absence points generated within the species' geographic extent.

**Spatial Analysis:** All spatial operations conducted using **terra** and **raster** packages in R, with **geodata** for administrative boundaries.

### Key References

- **Xerces Society.** (2015). Western bumble bee (*Bombus occidentalis*): Species profile and conservation status review.
- **Phillips, S. J., et al.** (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190(3-4), 231-259.
- **Fick, S. E., & Hijmans, R. J.** (2017). WorldClim 2: New 1km spatial resolution climate surfaces for global land areas. *International Journal of Climatology*, 37(12), 4302-4315.

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

