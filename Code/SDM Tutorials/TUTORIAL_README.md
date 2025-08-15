# Bombus occidentalis Habitat Modeling Tutorial

*Edited from a class taught by Dr. Greta Binford, Jeremy McWilliams, & Katy Prudic (2023).*

This tutorial demonstrates species distribution modeling (SDM) for *Bombus occidentalis* using presence-background data and environmental predictors.

## Workflow

1. **Data Preparation**  
   Occurrence data for *B. occidentalis* were combined with background points generated within the species’ geographic extent.

2. **Environmental Variables**  
   Bioclimatic variables from **WorldClim v2.1** (Fick & Hijmans, 2017) were processed using the `raster` and **terra** packages (Hijmans, 2023) to crop to the study region and extract environmental values at occurrence and background locations.

3. **Modeling**  
   Habitat suitability was estimated using the Maxent algorithm via the **maxnet** package (Phillips et al., 2010; Breiner et al., 2018) with presence-background data.

4. **Visualization**  
   Predictions were visualized with base R plotting and **ggplot2** (Wickham, 2020) including geographic basemaps for context.

5. **Export**  
   Final suitability maps were saved as GeoTIFFs for further use in GIS software.

## References

- Fick, S. E., & Hijmans, R. J. (2017). WorldClim 2: New 1 km spatial resolution climate surfaces for global land areas. *International Journal of Climatology, 37*(12), 4302–4315. https://doi.org/10.1002/joc.5086  
- Hijmans, R. J. (2023). *terra: Spatial Data Analysis* (R package version 1.8-61) [Computer software]. https://github.com/rspatial/terra  
- Elith, J., Phillips, S. J., Hastie, T., Dudík, M., Chee, Y. E., & Yates, C. J. (2010). A statistical explanation of MaxEnt for ecologists. *Diversity and Distributions, 17*(1), 43–57. https://doi.org/10.1111/j.1472-4642.2010.00725.x  
- Breiner, F. T., Nobis, M. P., Bergamini, A., & Guisan, A. (2018). Optimizing ensembles of small models for predicting the distribution of species with few occurrences. *Methods in Ecology and Evolution, 9*(4), 802–808. https://doi.org/10.1111/2041-210x.12957  
- Wickham, H. (2020). *ggplot2: Elegant graphics for data analysis.* Springer-Verlag New York. https://ggplot2.tidyverse.org
