# Australian marine heatwave data (mhw_australia)

Code for generating marine heatwave (MHW) metrics from observational data for Kajtar et al. (2021). Marine heatwaves are calculated following the Hobday et al. (2016) definition, using the marineHeatWaves python module (https://github.com/ecjoliver/marineHeatWaves). Marine heatwave categories follow the definition given by Hobday et al. (2018).

This code has been written to generate MHW metrics for the marine region around Australia (100°E-170°E, 50°S-0°). The metrics are computed from daily sea surface temperature (SST) data, from both observational data. The observed marine heatwave metrics are calculated from NOAA 0.25° daily Optimum Interpolation Sea Surface Temperature (OISST) over the period 1982-2020. Marine heatwaves are computed with respect to the 1983-2012 climatology. The marine heatwave data are provided on a grid point basis across the domain. Area-averaged SST timeseries and marine heatwave metrics are also provided for five case study regions: the Western Australia coast, Torres Strait, Great Barrier Reef, Tasman Sea, and South Australian Basin.

# Contents

|Directory         |Description|
|------------------|-----------|
|code              |Code for computing MHW metrics|
|marineHeatWaves   |Customised MHW detection code|
|post-processing   |Code for post-processing data|
|data              |Output and source data|

# code

Python code for computing MHW metrics. Data source paths must first be specified in each of the scripts.

|File              |Description|
|------------------|-----------|
|compute_australian_mhw_stats.py  |Compute MHW metrics across the Australian region.|
|store_aus_case_study_regions.py  |Store area-averaged SST timeseries, for further processing.|
|plot_fig_case_studies.py         |Processing and plotting of case study region MHW metrics. Code for generating Figs. 5-9.|

# marineHeatWaves

Marine heatwaves detection code, available at: https://github.com/ecjoliver/marineHeatWaves. The code provided here has been slightly adapted, but the original module should also work.


# References

Hobday, A.J., L. V. Alexander, S.E. Perkins-Kirkpatrick, D.A. Smale, S.C. Straub, E.C.J. Oliver, J.A. Benthuysen, M.T. Burrows, M.G. Donat, M. Feng, N.J. Holbrook, P.J. Moore, H.A. Scannell, A. Sen Gupta, T. Wernberg, A hierarchical approach to defining marine heatwaves, Prog. Oceanogr. 141 (2016) 227–238. https://doi.org/10.1016/j.pocean.2015.12.014.

Hobday, A.J., E.C.J. Oliver, A. Sen Gupta, J.A. Benthuysen, M.T. Burrows, M.G. Donat, N.J. Holbrook, P.J. Moore, M.S. Thomsen, T. Wernberg, D.A. Smale, Categorizing and naming marine heatwaves, Oceanography. 31 (2018) 162–173. https://doi.org/10.5670/oceanog.2018.205.

Kajtar, J.B., N.J. Holbrook, V. Hernaman, A catalogue of marine heatwave characteristics and trends for the Australian region. Under review at the Journal of Southern Hemisphere Earth Systems Science. (2021)

# Acknowledgements

We acknowledge the World Climate Research Programme's Working Group on Coupled Modelling, which is responsible for CMIP, and we thank the climate modelling groups for producing and making available their model output. CMIP6 model outputs were made available with the assistance of resources from the National Computational Infrastructure (NCI), which is supported by the Australian Government. NOAA High Resolution SST data were provided by the NOAA National Centers for Environmental Information. We acknowledge Eric Oliver for the use of his marineHeatWaves python module, freely available at https://github.com/ecjoliver/marineHeatWaves. 
