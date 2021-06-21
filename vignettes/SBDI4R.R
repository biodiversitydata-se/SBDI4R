## ----setup, include = FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# options(width = 120)

## ---- message=FALSE-----------------------------------------------------------
library(SBDI4R)
sbdi_config(caching="off")
# sbdi_config(cache_directory="Z:/mydir/sbdi-cache")

## ---- message=FALSE-----------------------------------------------------------
to_install <- c("ape", "dplyr", "ggplot2", "jpeg", "leaflet","maps", "mapdata",
                "maptools", "phytools", "sf", "rgeos", "tidyr", "vegan", 
                "phytools", "BIRDS", "leaflet", "rgdal")
to_install <- to_install[!sapply(to_install, requireNamespace, quietly=TRUE)]
if(length(to_install)>0)
    install.packages(to_install, repos="http://cran.us.r-project.org")

## ---- warning=FALSE, message=FALSE--------------------------------------------
sx <- search_fulltext("parus")
sx$data[,c( "name","species", "speciesGuid", "rank")]

## ---- message=FALSE-----------------------------------------------------------
sx <- search_fulltext("parus", fq="family_s:Paridae")
sx$data[,c( "name","species", "speciesGuid", "rank")]

## ---- message=FALSE-----------------------------------------------------------
sx <- search_fulltext("parus", fq="class_s:Aves", page_size=100)
head(sx$data[,c( "name","species", "speciesGuid", "rank")])

## ---- message=FALSE-----------------------------------------------------------
tx <- taxinfo_download("family_s:Paridae", 
                       fields = c("guid", "genus_s", "specificEpithet_s", "scientificName",  "canonicalName_s", "rank"), 
                       verbose = FALSE)
tx <- tx[tx$rank == "species" & tx$genusS != "",] ## restrict to species and not hybrids

## ---- message=FALSE, fig.width=8, fig.height=6--------------------------------
library(phytools)
## as.phylo requires the taxonomic columns to be factors
tx$genusS <- as.factor(tx$genusS)
tx$scientificName <- as.factor(tx$scientificName)
tx$canonicalNameS <- as.factor(tx$canonicalNameS)
## create phylo object of canonical name nested within Genus
ax <- as.phylo(~genusS/canonicalNameS, data=tx[1:50,])
plotTree(ax, fsize=0.7, ftype="i") ## plot it

## ---- eval=FALSE--------------------------------------------------------------
#  x <- occurrences(taxon="Callitriche cophocarpa",
#                   email="sbdi4r-test@biodiversitydata.se",
#                   download_reason_id=10)
#  head(x$data)
#  table(x$data$dataResourceName)
#  table(x$data$dataResourceID)

## ---- eval=FALSE--------------------------------------------------------------
#  x <- occurrences(taxon="sommarlånke",
#                   email="sbdi4r-test@biodiversitydata.se",
#                   download_reason_id=10,
#                   verbose = FALSE)
#  head(x$data)
#  table(x$data$dataResourceName)
#  table(x$data$dataResourceID)

