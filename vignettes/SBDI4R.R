## ----setup, include = FALSE-------------------------------------------------------------------------------------------
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(width = 120)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
library(SBDI4R)
sbdi_config(caching="off")

## ---- message=FALSE---------------------------------------------------------------------------------------------------
to_install <- c("ape", "dplyr", "ggplot2", "jpeg", "leaflet","maps", "mapdata",
                "maptools", "phytools", "sp", "rgeos", "tidyr", "vegan")
to_install <- to_install[!sapply(to_install, requireNamespace, quietly=TRUE)]
if(length(to_install)>0)
    install.packages(to_install, repos="http://cran.us.r-project.org")

