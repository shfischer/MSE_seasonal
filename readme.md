Exploration of seasonal Management Strategy evaluation (MSE) modelling
with FLR
================

## Introduction

This repository contains the code for exploratory seasonal Management
Strategy Evaluation (MSE) modelling. The simulation is based on the
Fisheries Library in R ([FLR](http://www.flr-project.org/)) and use the
([`FLasher`](github.com/FLR/FLasher)) R package for the seasonal
projections.

## Repository structure

The root folder contains the following R scripts:

-   `OM_san.R`: This script creates the operating models (OMs),
-   `funs.R` contains functions and methods used for the creation of the
    operating models and for running the MSE,
-   `MSY_HR.R`: contains the code for finding deterministic long-term
    optimum values for the harvest rate and escapement strategy,
-   `MP_run.R`: is script for running stochastic simulations,
-   `MP_run.pbs`: is job submission script for an HPC system to run
    `MP_run.R`,
-   `MP.R`: contains the code for running stochastic simulations
    locally.

The following input files are provided:

-   `input/stocks.csv` contains the stock definitions and life-history
    parameters

## R, R packages and version info

The simulations were run with R 4.1:

``` r
> sessionInfo()
R version 4.1.0 (2021-05-18)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19041)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252   
[3] LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] cowplot_1.1.1        forcats_0.5.1        stringr_1.4.0        purrr_0.3.4         
 [5] readr_2.0.1          tibble_3.1.6         tidyverse_1.3.1      doParallel_1.0.16   
 [9] foreach_1.5.1        FLasher_0.6.8        FLFishery_0.3.7.9003 tidyr_1.1.3         
[13] dplyr_1.0.7          FLife_3.4.0          FLBRP_2.5.8          ggplotFL_2.6.10.9001
[17] ggplot2_3.3.5        FLCore_2.6.18.9005   iterators_1.0.13     lattice_0.20-44     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.7        lubridate_1.7.10  assertthat_0.2.1  utf8_1.2.2        R6_2.5.1         
 [6] cellranger_1.1.0  backports_1.2.1   reprex_2.0.1      stats4_4.1.0      httr_1.4.2       
[11] pillar_1.6.4      rlang_0.4.12      readxl_1.3.1      rstudioapi_0.13   data.table_1.14.2
[16] Matrix_1.3-3      munsell_0.5.0     broom_0.7.9       compiler_4.1.0    modelr_0.1.8     
[21] pkgconfig_2.0.3   tidyselect_1.1.1  gridExtra_2.3     codetools_0.2-18  fansi_0.5.0      
[26] crayon_1.4.2      tzdb_0.1.2        dbplyr_2.1.1      withr_2.4.2       MASS_7.3-54      
[31] grid_4.1.0        jsonlite_1.7.2    gtable_0.3.0      lifecycle_1.0.1   DBI_1.1.1        
[36] magrittr_2.0.1    scales_1.1.1      cli_3.1.0         stringi_1.7.4     renv_0.14.0      
[41] fs_1.5.0          xml2_1.3.2        ellipsis_0.3.2    generics_0.1.0    vctrs_0.3.8      
[46] tools_4.1.0       glue_1.5.0        hms_1.1.0         colorspace_2.0-2  rvest_1.0.1      
[51] haven_2.4.3
```

The versions of the R packages are recorded with `renv` in the file
[`renv.lock`](renv.lock). In order to reproduce the analyses of this
repository, exactly these R package versions (particularly for the FLR
packages) are required.
