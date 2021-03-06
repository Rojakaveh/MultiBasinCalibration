Multi-Basin Calibration (MBC) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6609237.svg)](https://doi.org/10.5281/zenodo.6609237)
=================

Multi-Basin Calibration method (MBC) to calibrate several neighboring and similar watersheds simultaneously in Low Data Environment (LDE). Running the codes in this repository requirs `R` and `RStudio`.

## Links
See the following links for more information on  `R` and `RStudio` download and installation:

- An introduction to `R`: <https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf>
- `R` download: <https://www.r-project.org/>
- `RStudio` download: <https://www.rstudio.com/>

There is also a cloud-based `RStudio` sever at the following location:

- Cloud-based `RStudio` server: <https://rstudio.cloud/>
## Description
This repository contains R codes for a project called **Watershed Model Parameter Estimation in Low Data Environments**. The study area is located within the Lake Champlain Basin (LCB) of Vermont (VT), USA. "USGS 04282650 - LOC-Ferrisburg", "USGS 04282629 - Little Otter Creek-Monkton" (LOC-Monkton), "USGS 04282586 - West Branch Dead Creek" (WBDC), and "USGS 04282581 - East Branch Dead Creek" (EBDC) are four watersheds located in the southern portion of LCB that were used to develop the model. LOC-Monkton, WBDC, and EBDC are three newly instrumented watersheds with short recorded streamflow data (26 months), and LOC-Ferrisburg has a longer recorded data (31 years of daily streamflow data).

## Quick start

### R packages that need to be installed:
•   httr
•   EcoHydRology
•   curl
•   tidyverse
•   lubridate
•   data.table
•   forecast
•   rnoaa
•   dplyr
•   SWATmodel
•   ggplot2

        if (!require("pacman")) install.packages("pacman")
        pacman::p_load(httr,EcoHydRology,curl,tidyverse,lubridate,data.table,forecast,rnoaa,dplyr,SWATmodel,ggplot2)

##### SWATInitSingleCalib.R
SWAT initialization for all 4 watersheds and single single basin calibration code.
For EBDC watershed we used rating curve to generate measured daily streamflows that can be found in **EBDCflowdata.RData** file in this repository.
In order to run this script in `RStudio`:

        download.file("https://raw.githubusercontent.com/Rojakaveh/MultiBasinCalibration/main/SWATInitSingleCalib.R","SWATInitSingleCalib.R")
        file.edit("SWATInitSingleCalib.R")

##### SensitivityAnalysis.R
New approach for sesnitivity analysis. We used **DEoptim** function within R DEoptim package to generate five datasets of 300 evolutionary optimization iteration parameter sets aligned with the progressively increasing Nash-Sutcliffe model efficiency coefficients (NSEs) for each period of 2015-2021 and 2010-2021 (10 datasets in total). Datasets can be found in **optimizitiondatasets** file in this repository. This R code plots all parameters and NSE values, parameters and NSE percent changes between each two itterations, remove non-influential parameters, and rank the sensitive parameters using mean of reletive sensitivity and visualize rank of the parameters using a violin-boxplot plot.
In order to run this script in `RStudio`:

        download.file("https://raw.githubusercontent.com/Rojakaveh/MultiBasinCalibration/main/SensitivityAnalysis.R","SensitivityAnalysis.R")
        file.edit("SensitivityAnalysis.R")

##### MBCfunc.R
This code contains an R function that consists of two steps, first streamflow data from each of three LDE watersheds (LOC-Monkton, EBDC, and WBDC) are aggregated and second these data are then calibrated as a single dataset.In order to run this script in `RStudio`:

        download.file("https://raw.githubusercontent.com/Rojakaveh/MultiBasinCalibration/main/MBCfunc.R","MBCfunc.R")
        file.edit("MBCfunc.R")
        
##### LOOCV.R
The code for implementation Leave One Out Cross-Validation.In order to run this script in `RStudio`:

        download.file("https://raw.githubusercontent.com/Rojakaveh/MultiBasinCalibration/main/LOOCV.R","LOOCV.R")
        file.edit("LOOCV.R")

# License
Please see the LICENSE.md file for license information.
