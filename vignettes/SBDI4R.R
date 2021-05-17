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
# sbdi_config(cache_directory="Z:/mydir/sbdi-cache")

## ---- message=FALSE---------------------------------------------------------------------------------------------------
to_install <- c("ape", "dplyr", "ggplot2", "jpeg", "leaflet","maps", "mapdata",
                "maptools", "phytools", "sp", "rgeos", "tidyr", "vegan", 
                "phytools", "BIRDS", "leaflet", "rgdal")
to_install <- to_install[!sapply(to_install, requireNamespace, quietly=TRUE)]
if(length(to_install)>0)
    install.packages(to_install, repos="http://cran.us.r-project.org")

## ---- warning=FALSE, message=FALSE------------------------------------------------------------------------------------
sx <- search_fulltext("parus")
sx$data[,c( "name","species", "speciesGuid", "rank")]

## ---- message=FALSE---------------------------------------------------------------------------------------------------
sx <- search_fulltext("parus", fq="family_s:Paridae")
sx$data[,c( "name","species", "speciesGuid", "rank")]

## ---- message=FALSE---------------------------------------------------------------------------------------------------
sx <- search_fulltext("parus", fq="class_s:Aves", page_size=100)
head(sx$data[,c( "name","species", "speciesGuid", "rank")])

## ---- message=FALSE, results=FALSE------------------------------------------------------------------------------------
tx <- taxinfo_download("family_s:Paridae", 
                       fields = c("guid", "genus_s", "specificEpithet_s", "scientificName",  "canonicalName_s", "rank"), 
                       verbose = FALSE)
tx <- tx[tx$rank == "species" & tx$genusS != "",] ## restrict to species and not hybrids

## ---- message=FALSE---------------------------------------------------------------------------------------------------
library(phytools)
## as.phylo requires the taxonomic columns to be factors
tx$genusS <- as.factor(tx$genusS)
tx$scientificName <- as.factor(tx$scientificName)
tx$canonicalNameS <- as.factor(tx$canonicalNameS)
## create phylo object of canonical name nested within Genus
ax <- as.phylo(~genusS/canonicalNameS, data=tx[1:50,])
plotTree(ax, fsize=0.7, ftype="i") ## plot it

## ---- message=FALSE, eval=FALSE---------------------------------------------------------------------------------------
#  x <- occurrences(taxon="Callitriche cophocarpa",
#                   email="test@test.org",
#                   download_reason_id=10)
#  head(x$data)
#  table(x$data$dataResourceName)

## ---- message=FALSE, eval=FALSE---------------------------------------------------------------------------------------
#  x <- occurrences(taxon="sommarlånke",
#                   email="test@test.org",
#                   download_reason_id=10,
#                   verbose = FALSE)
#  head(x$data)
#  table(x$data$dataResourceName)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
taxa <- c("Callitriche", "Anarrhinum")
fq_str <- paste0("raw_name:&quot;", taxa, "&quot;")
fq_str <- paste0(fq_str, collapse = " OR ")
xbatch <- occurrences(fq=fq_str, 
                 email="test@test.org", 
                 download_reason_id=10, 
                 verbose = FALSE)
head(xbatch$data)
table(xbatch$data$dataResourceName)
table(xbatch$data$basisOfRecord)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
xf <- occurrences(taxon="Callitriche cophocarpa", 
                 fq = "data_resource_uid:dr5",
                 email="test@test.org", 
                 download_reason_id=10, verbose = FALSE)
table(xf$data$dataResourceName)

## ---- message=FALSE, eval= FALSE--------------------------------------------------------------------------------------
#  fq_str <- pick_filter("resource")
#  ## follow the instructions
#  
#  xf <- occurrences(taxon="Callitriche cophocarpa",
#                   fq = fq_str,
#                   email="test@test.org",
#                   download_reason_id=10)

## ---- message=FALSE, eval= FALSE--------------------------------------------------------------------------------------
#  # fq_str <- pick_filter("layer")
#  # Follow the instructions, but here we just use the county Uppsala
#  fq_str <- "cl10097:Uppsala"
#  
#  xf <- occurrences(taxon="Callitriche cophocarpa",
#                   fq = fq_str,
#                   email="test@test.org",
#                   download_reason_id=10)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
xf <- occurrences(taxon="Callitriche cophocarpa", 
                 fq="coordinate_uncertainty:[0 TO 100]",
                 email="test@test.org", 
                 download_reason_id=10)

range(xf$data$coordinateUncertaintyInMetres)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
# year = 2019
x2019 <- occurrences(taxon="Reynoutria japonica", 
                     fq="year:2019",
                     email="test@test.org", 
                     download_reason_id=10)
nrow(x2019$data)

x2yr <- occurrences(taxon="Reynoutria japonica", 
                    fq=c("year:2018 OR year:2019"),
                    email="test@test.org", 
                    download_reason_id=10)
nrow(x2yr$data)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
xf <- occurrences(taxon="Callitriche cophocarpa", 
                 fq="year:[2010 TO 2020]",
                 email="test@test.org", 
                 download_reason_id=10)

hist(xf$data$year, xlab = "Year", main = "")

## ---- message=FALSE---------------------------------------------------------------------------------------------------
xf <- occurrences(taxon="Callitriche cophocarpa", 
                 fq=c("year:[2010 TO 2020]", "month:[06 TO 08]"),
                 email="test@test.org", 
                 download_reason_id=10)

hist(xf$data$year, xlab = "Year", main = "")

## ---- message=FALSE---------------------------------------------------------------------------------------------------
xf <- occurrences(taxon="Callitriche cophocarpa", 
                 fq="basis_of_record:HumanObservation",
                 email="test@test.org", 
                 download_reason_id=10)

unique(xf$data$basisOfRecord)

## ---- message=FALSE---------------------------------------------------------------------------------------------------
x <- occurrences(taxon="Callitriche cophocarpa",
                 fq = "data_resource_uid:dr5",
                 email="test@test.org",
                 download_reason_id=10, 
                 verbose = FALSE)
summary(x)

## ---------------------------------------------------------------------------------------------------------------------
assert <- sbdi_fields("assertions")
assertFatal <- assert[assert$isFatal==TRUE,"name"]
wAssertInX <- assertFatal %in% colnames(x$data)
colSums(x$data[,assertFatal[wAssertInX]])

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  occurrences_plot(x, "obsPlot.pdf", qa="error",
#                    grouped=FALSE, taxon_level="species",
#                    pch='+')

## ---- message=FALSE---------------------------------------------------------------------------------------------------
library(leaflet)
## drop any records with missing lat/lon values
x$data <- x$data[!is.na(x$data$longitude) & !is.na(x$data$latitude),] 
xa <- check_assertions(x)
## columns of x corresponding to a fatal assertion
x_afcols <- which(names(x$data) %in% xa$occurColnames)
## rows of x that have a fatal assertion
x_afrows <- apply(x$data[,x_afcols], 1, any)
## which taxonIdentificationIssue assertions are present in this data?
these_assertions <- names(x$data)[x_afcols]
## make a link to the web page for each occurrence
popup_link <- paste0("<a href=\"https://records.biodiversitydata.se/occurrences/",
                      x$data$id,"\">Link to occurrence record</a>")
## colour palette
pal <- c(sub("FF$","", heat.colors(length(these_assertions))))
## map each data row to colour, depending on its assertions
marker_colour <- rep("#00FF00", nrow(x$data))
if(length(these_assertions)>0){
  for (k in 1:length(these_assertions)){
    marker_colour[x$data[,x_afcols[k]]] <- pal[k]
  } 
}

## blank map, with imagery background
m <- addProviderTiles(leaflet(),"Esri.WorldImagery") %>% 
  ## add markers
  addCircleMarkers(x$data$longitude, x$data$latitude,  
                   radius = 2, fillOpacity =.5, opacity = 1,
                   col=marker_colour, popup=popup_link) %>% 
  addLegend(colors = pal, opacity = 1, labels = these_assertions)
m

## ---- eval=FALSE------------------------------------------------------------------------------------------------------
#  # save as data.frame
#  Callitriche <- as.data.frame(x$data)
#  
#  # simplyfy data frame
#  calli <- data.frame(Callitriche$scientificName,
#                     Callitriche$latitude,
#                     Callitriche$longitude)
#  # simplify column names
#  colnames(calli) <- c("species","latitude","longitude")
#  # remove rows with missing values (NAs)
#  calli <- na.omit(calli)
#  
#  # save new dataframe
#  write.csv(calli,"Callitriche.csv")

## ---- message=FALSE, warning=FALSE------------------------------------------------------------------------------------

x <- occurrences(taxon="Callitriche cophocarpa",
                 fq = "data_resource_uid:dr5",
                 email="test@test.org",
                 download_reason_id=10, 
                 verbose = FALSE)

library(sp) # the function coordinates() and proj4string() are in sp
library(rgeos) #  the function over() is in package rgeos
# load some shapes over Sweden
# Political borders
data("swe_wgs84", package="SBDI4R", envir=environment()) 
# A standard 50km grid
data("Sweden_Grid_50km_Wgs84", package="SBDI4R", envir=environment()) 

grid <- Sweden_Grid_50km_Wgs84
# grid <- spTransform(grid, CRS("+init=epsg:4326")) ## it has the same CRS 
# but changes are undergoing in the sp package and this step is needed

# make the observations spatial
# NOTE: make sure there are no NAs on either column defining the coordinates 
# see example 2 for cleaning your dataset.

obs <- as.data.frame(x$data)
coordinates(obs) <- obs[,c("longitude","latitude")]
wkt <- sf::st_crs(4326)[[2]]
proj4string(obs) <- sp::CRS(wkt)

nObs <- nrow(obs)

## overlay the data with the grid
ObsInGridList <- over(grid, obs, returnList=TRUE)
wNonEmpty <- unname( which( unlist(lapply(ObsInGridList, nrow)) != 0) )
if(length(wNonEmpty)==0) message("Observations don't overlap any grid cell.")

## check nObs
nObsInGrid <- sum(unlist(lapply(ObsInGridList, nrow)))

## ---------------------------------------------------------------------------------------------------------------------
## apply a summary over the grid
nCells <- length(ObsInGridList)

res <- data.frame("nObs"=as.numeric(rep(NA,nCells)),
                  "nYears"=as.numeric(rep(NA,nCells)),
                  row.names = row.names(grid),
                  stringsAsFactors = FALSE)

cols2use <- c("scientificName", "year")

dataRes <- lapply(ObsInGridList[wNonEmpty], function(x){
  x <- x[,cols2use]
  colnames(x) <- c("scientificName", "year")
  
  return(c("nObs" = as.numeric(length(x[,"scientificName"])),
           "nYears" = length(unique(x[,"year"]))
  ))
})

dataRes <- as.data.frame(dplyr::bind_rows(dataRes, .id = "gridID"))

res[wNonEmpty,] <- dataRes[,-1]
res$nObs <- as.numeric(res$nObs)
resSp <- sp::SpatialPolygonsDataFrame(grid, res)

## ----grid, warning=FALSE, fig.width=6, fig.height=6, message=FALSE----------------------------------------------------
palBW <- leaflet::colorNumeric(palette = c("white", "navyblue"),
                               domain = c(0, max(resSp@data$nObs, na.rm = TRUE)), 
                               na.color = "transparent")
oldpar <- par()
par(mar = c(1,1,0,0))
plot(resSp, col=palBW(resSp@data$nObs), border = NA)
plot(swe_wgs84$Border, border=1, lwd=1, add=T)
legend("bottomleft", 
       legend = round(seq(0, max(resSp@data$nObs, na.rm = TRUE), length.out = 5)),
       col = palBW(seq(0, max(resSp@data$nObs, na.rm = TRUE), length.out = 5)),
       title = "Number of \nobservations", pch = 15, bty="n")
suppressWarnings(par(oldpar))

