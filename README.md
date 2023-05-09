# The Spatial Distribution of Zooplankton Biomass

This repository contains the statistical models used for the analysis of global zooplankton biomass patterns in Heneghan et al. (2020) and Petrik et al. (2022). 

The COPEPOD dataset (O'Brien 2014) was filtered for duplicate or incorrect values before being used in this analysis. The filtered dataset can be found in the `Data/` folder. This data was generated as part of the Heneghan et al. (2020) work, and further information about the process is available from the lead authors of that paper. Our goal is to publish this process as a separate analysis in the near future.
  
The statistical model used in this research was generated in an iterative process whereby the fit, error structure and biological meaning were evaluated at each step. Model 7 (m7) was selected as the best model and used in this analysis. This selection process can be recreated following the code in `GlobalZooBio_01_Model.R`. The output from each model is stored in `Output/` and the corresponding plots are in `Figures/`, but please note, Model 7 (m7) was the chosen model for further analysis.

The gridded global zooplankton biomass used in Petrik et al. (2022) was generated using the code in `GlobalZooBio_02_Predict.R`. The one-degree gridded satellite SST, satellite Chlorophyll *a* and bathymetry data used in the mapping can be found in `one/`. The output of the prediction is stored in `Output/`

## References

Heneghan, R.F., Everett, J.D., Sykes, P., Batten, S.D., Edwards, M., Takahashi, K., Suthers, I.M., Blanchard, J.L. and Richardson, A.J., 2020. A functional size-spectrum model of the global marine ecosystem that resolves zooplankton composition. Ecological Modelling, 435, p.109265. https://doi.org/10.1016/j.ecolmodel.2020.109265

O'Brien, T.D. 2014.   COPEPOD:   The Global Plankton Database. An overview of the 2014 database contents, processing methods, and access interface. U.S. Dep. Commerce, NOAA Tech. Memo. NMFS-F/ST-37, 29p.

Petrik, C. M., Luo, J. Y., Heneghan, R. F., Everett, J. D., Harrison, C. S., & Richardson, A. J. (2022). Assessment and constraint of mesozooplankton in CMIP6 Earth system models. Global Biogeochemical Cycles, 36, e2022GB007367. https://doi.org/10.1029/2022GB007367
