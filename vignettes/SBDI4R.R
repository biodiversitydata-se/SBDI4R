## ----setup, include = FALSE---------------------------------------------------
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
# options(width = 120)

## ---- message=FALSE-----------------------------------------------------------
library(SBDI4R)
sbdi_config(caching = "off")
sbdi_config(email = "sbdi4r-test@biodiversitydata.se", 
            download_reason_id = 10)
# sbdi_config(cache_directory="Z:/mydir/sbdi-cache")

## ---- message=FALSE-----------------------------------------------------------
to_install <- c("ape", "dplyr", "ggplot2", "jpeg", "leaflet","maps", "mapdata",
                "maptools", "phytools", "sf",  "tidyr", "vegan") #, "rgdal", "rgeos",
to_install <- to_install[!sapply(to_install, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0)
    install.packages(to_install, repos = "http://cran.us.r-project.org")

## ---- warning=FALSE, message=FALSE--------------------------------------------
sx <- search_fulltext("parus")
sx$data[,c( "name","species", "speciesGuid", "rank")]

## ---- message=FALSE-----------------------------------------------------------
sx <- search_fulltext("parus", fq = "family_s:Paridae")
sx$data[,c( "name","species", "speciesGuid", "rank")]

## ---- message=FALSE-----------------------------------------------------------
sx <- search_fulltext("parus", fq = "class_s:Aves", page_size=100)
head(sx$data[,c( "name","species", "speciesGuid", "rank")])

## ---- message=FALSE-----------------------------------------------------------
tx <- taxinfo_download("family_s:Paridae", 
                       fields = c("guid", "genus_s", "specificEpithet_s", 
                                  "scientificName",  "canonicalName_s", "rank"), 
                       verbose = FALSE)
tx <- tx[tx$rank == "species" & tx$genusS != "",] ## restrict to species and not hybrids

## ---- message=FALSE, fig.width=8, fig.height=6--------------------------------
library(phytools)
## as.phylo requires the taxonomic columns to be factors
tx$genusS <- as.factor(tx$genusS)
tx$scientificName <- as.factor(tx$scientificName)
tx$canonicalNameS <- as.factor(tx$canonicalNameS)
## create phylo object of canonical name nested within Genus
ax <- as.phylo(~genusS/canonicalNameS, data = tx[1:50,])
plotTree(ax, fsize = 0.7, ftype="i") ## plot it

## ---- eval=FALSE--------------------------------------------------------------
#  x <- occurrences(taxon = "Callitriche cophocarpa")
#  head(x$data)
#  table(x$data$dataResourceName)
#  table(x$data$dataResourceID)

## ---- eval=FALSE--------------------------------------------------------------
#  x <- occurrences(taxon = "sommarlånke",
#                   verbose = FALSE)
#  head(x$data)
#  table(x$data$dataResourceName)
#  table(x$data$dataResourceID)

## ---- message=FALSE-----------------------------------------------------------
taxa <- c("Callitriche", "Anarrhinum")
fq_str <- paste0("raw_name:&quot;", taxa, "&quot;")
# fq_str <- paste0(fq_str, collapse = " OR ")
xbatch <- occurrences(fq = fq_str,
                      verbose = FALSE)
table(xbatch$data$dataResourceName)
table(xbatch$data$dataResourceID)
table(xbatch$data$basisOfRecord)

## ---- message=FALSE-----------------------------------------------------------
xf <- occurrences(taxon = "Callitriche cophocarpa", 
                  fq = "data_resource_uid:dr2",
                  verbose = FALSE)
table(xf$data$dataResourceName)
table(xf$data$dataResourceID)

## ---- eval= FALSE-------------------------------------------------------------
#  fq_str <- pick_filter("resource")
#  ## follow the instructions
#  
#  xf <- occurrences(taxon = "Callitriche cophocarpa",
#                    fq = fq_str)

## ---- eval= FALSE-------------------------------------------------------------
#  # fq_str <- pick_filter("layer")
#  # Follow the instructions, but here we just use the county Uppsala
#  fq_str <- "cl10097:Uppsala"
#  
#  xf <- occurrences(taxon = "Callitriche cophocarpa",
#                    fq = fq_str)

## ---- message=FALSE-----------------------------------------------------------
xf <- occurrences(taxon = "Callitriche cophocarpa", 
                  fq = "coordinate_uncertainty:[0 TO 100]")

range(xf$data$coordinateUncertaintyInMetres)

## ---- message=FALSE-----------------------------------------------------------
# year = 2019
x2011 <- occurrences(taxon = "Callitriche cophocarpa", 
                     fq = "year:2011")
nrow(x2011$data)

x2yr <- occurrences(taxon = "Callitriche cophocarpa", 
                    fq = c("year:2011 OR year:2014"))
nrow(x2yr$data)

## ---- message=FALSE, fig.width=8, fig.height=6--------------------------------
xf <- occurrences(taxon = "Callitriche cophocarpa", 
                  fq = "year:[2010 TO 2020]")

hist(xf$data$year, xlab = "Year", main = "")

## ---- message=FALSE, fig.width=8, fig.height=6--------------------------------
xf <- occurrences(taxon = "Callitriche cophocarpa", 
                 fq = c("year:[2010 TO 2020]", "month:[06 TO 08]"))

hist(xf$data$year, xlab = "Year", main = "")

## ---- message=FALSE-----------------------------------------------------------
xf <- occurrences(taxon = "Callitriche cophocarpa", 
                  fq = "basis_of_record:PreservedSpecimen")

unique(xf$data$basisOfRecord)

## ---- message=FALSE-----------------------------------------------------------
x <- occurrences(taxon = "Callitriche cophocarpa",
                 fq = "data_resource_uid:dr2",
                 verbose = FALSE)
summary(x)

## ----assertions---------------------------------------------------------------
assert <- sbdi_fields("assertions")
assertFatal <- assert[assert$isFatal == TRUE,"name"]
wAssertInX <- assertFatal %in% colnames(x$data)
colSums(x$data[,assertFatal[wAssertInX]])

## ---- eval=FALSE--------------------------------------------------------------
#  occurrences_plot(x, "obsPlot.pdf", qa = "error",
#                   grouped = FALSE, taxon_level = "species",
#                   pch = '+')

## ---- message=FALSE, fig.width=9, fig.height=9--------------------------------
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
m <- leaflet() |>  
  addProviderTiles("Esri.WorldImagery") |>  
  ## add markers
  addCircleMarkers(x$data$longitude, x$data$latitude,  
                   radius = 2, fillOpacity = .5, opacity = 1,
                   col = marker_colour, popup = popup_link) |> 
  addLegend(colors = pal, opacity = 1, labels = these_assertions)
m

## ---- eval=FALSE--------------------------------------------------------------
#  # save as data.frame
#  Callitriche <- as.data.frame(x$data)
#  
#  # simplyfy data frame
#  calli <- data.frame(Callitriche$scientificName,
#                      Callitriche$latitude,
#                      Callitriche$longitude)
#  # simplify column names
#  colnames(calli) <- c("species","latitude","longitude")
#  # remove rows with missing values (NAs)
#  calli <- na.omit(calli)
#  
#  # save new dataframe
#  write.csv(calli,"Callitriche.csv")

## ---- message=FALSE, warning=FALSE--------------------------------------------

x <- occurrences(taxon = "Callitriche cophocarpa",
                 fq = "data_resource_uid:dr2", 
                 verbose = FALSE)

library(sf)
# load some shapes over Sweden
# Political borders
data("swe_wgs84", package="SBDI4R", envir=environment()) 
# A standard 50km grid
data("Sweden_Grid_50km_Wgs84", package="SBDI4R", envir=environment()) 

grid <- Sweden_Grid_50km_Wgs84
grid <- st_transform(grid, crs = st_crs(3006))

# make the observations spatial
# NOTE: make sure there are no NAs on either column defining the coordinates 
# see example 2 for cleaning your dataset.
wNAcoor <- which(is.na(x$data$latitude))

obs <- st_as_sf(x = x$data[-wNAcoor,], 
                coords = c("longitude","latitude"))
st_crs(obs) <- st_crs(4326)
obs <- st_transform(obs, crs = st_crs(grid))

nObs <- nrow(obs)

## overlay the data with the grid
listGrid <- st_intersects(grid, obs)

ObsInGridList <- list()
for (i in seq(length(listGrid))) {
  if (length(listGrid[[i]]) == 0) {
    ObsInGridList[[i]] <- NA
  } else {
    ObsInGridList[[i]] <- st_drop_geometry(obs[listGrid[[i]],])
  }
}
wNonEmpty <- which( unlist(lapply(ObsInGridList, function(x) !all(is.na(x)))) )
if (length(wNonEmpty) == 0) message("Observations don't overlap any grid cell.")

## check nObs
sum(unlist(lapply(ObsInGridList, nrow))) == nObs
length(ObsInGridList) == nrow(grid)


## -----------------------------------------------------------------------------
## apply a summary over the grid
nCells <- length(ObsInGridList)

res <- data.frame("nObs" = as.numeric(rep(NA, nCells)),
                  "nYears" = as.numeric(rep(NA, nCells)),
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

dataRes <- as.data.frame(dplyr::bind_rows(dataRes, .id = "id"))

res[wNonEmpty,] <- dataRes[,-1]
res$nObs <- as.numeric(res$nObs)
resSf <- st_as_sf(cbind(res, st_geometry(grid)) )
rownames(resSf) <- grid$id

## ----grid, message=FALSE, warning=FALSE, fig.width=6, fig.height=6------------
palBW <- leaflet::colorNumeric(palette = c("white", "navyblue"),
                               domain = c(0, max(resSf$nObs, na.rm = TRUE)), 
                               na.color = "transparent")
oldpar <- par()
par(mar = c(1,1,0,0))
plot(resSf$geometry, col = palBW(resSf$nObs), border = NA)
plot(swe$Border, border = 1, lwd = 1, add = T)
legend("bottomleft", 
       legend = round(seq(0, max(resSf$nObs, na.rm = TRUE), length.out = 5)),
       col = palBW(seq(0, max(resSf$nObs, na.rm = TRUE), length.out = 5)),
       title = "Number of \nobservations", pch = 15, bty = "n")
suppressWarnings(par(oldpar))

## ---- message=FALSE, warning=FALSE--------------------------------------------
counties <- swe_wgs84$Counties
obs <- st_transform(obs, crs = st_crs(counties))

## overlay the data with the counties
listGrid <- st_intersects(counties, obs)

ObsInCountyList <- list()
for (i in seq(length(listGrid))) {
  if (length(listGrid[[i]]) == 0) {
    ObsInCountyList[[i]] <- NA
  } else {
    ObsInCountyList[[i]] <- st_drop_geometry(obs[listGrid[[i]],])
  }
}
wNonEmpty <- which( unlist(lapply(ObsInCountyList, function(x) !all(is.na(x)))) )
if (length(wNonEmpty) == 0) message("Observations don't overlap any grid cell.")

## check nObs
sum(unlist(lapply(ObsInCountyList, nrow))) == nObs # some observations are not in the counties territory
length(ObsInCountyList) == nrow(counties)

## apply a summary over the grid
nCells <- length(ObsInCountyList)

res <- data.frame("nObs" = as.numeric(rep(NA, nCells)),
                  "nYears" = as.numeric(rep(NA, nCells)),
                  stringsAsFactors = FALSE)

cols2use <- c("scientificName", "year")

dataRes <- lapply(ObsInCountyList[wNonEmpty], function(x){
  x <- x[,cols2use]
  colnames(x) <- c("scientificName", "year")
  
  return(c("nObs" = length(x[,"scientificName"]),
           "nYears" = length(unique(x[,"year"]))
  ))
})


dataRes <- as.data.frame(dplyr::bind_rows(dataRes, .id = "id"))
res[wNonEmpty,] <- dataRes[,-1]
res$nObs <- as.numeric(res$nObs)

resSf <- st_as_sf(cbind(res, st_geometry(counties)))
rownames(resSf) <- counties$LnNamn

## ----counties, warning=FALSE, fig.width=6, fig.height=6-----------------------
palBW <- leaflet::colorNumeric(c("white", "navyblue"), 
                               c(0, max(resSf$nObs, na.rm = TRUE)), 
                               na.color = "transparent")
oldpar <- par()
par(mar = c(1,1,0,0))
plot(resSf$geometry, col = palBW(resSf$nObs), border = NA)
plot(swe_wgs84$Border, border = 1, lwd = 1, add = T)
text(st_coordinates(st_centroid(counties)), 
    labels = as.character(counties$LnNamn), font = 2, cex = .5 )
legend("bottomleft", 
       legend = round(seq(0, max(resSf$nObs, na.rm = TRUE), length.out = 5)),
       col = palBW(seq(0, max(resSf$nObs, na.rm = TRUE), length.out = 5)),
       title = "Number of \nobservations", pch = 15, bty = "n")
suppressWarnings(par(oldpar))

## ---- warning=FALSE, fig.width=8, fig.height=6--------------------------------
countiesLab <- as.character(counties$LnNamn)

## Add a column to the obs data.frame to hold the id of the overlapped polygon, 
## in this case, Län (county) and plot which observation didn't fall with any 
## territory.
obs$overId <- NA 

for (c in 1:length(ObsInCountyList)) {
  if (nrow(ObsInCountyList[[c]]) > 0) {
    idsC <- ObsInCountyList[[c]]$id
    wObs <- match(idsC, obs$id)
    obs$overId[wObs] <- rep(countiesLab[c], length(wObs))
  }
}

oldpar <- par()
par(mar = c(1,1,0,0))
plot(counties$geometry, border = 1, lwd = 1)
plot(obs$geometry[which(is.na(obs$overId))], 
     pch = 19, cex = .5, col = "red", add = T)
suppressWarnings(par(oldpar))

## ---- eval=FALSE--------------------------------------------------------------
#  shape <- st_read(dsn=file.path("your/path/to/file", "Kommun_Sweref99TM_region.shp"))

## -----------------------------------------------------------------------------
shape <- swe$Municipalities
## extract just the Municipality of Örebro
shape <- shape[shape$KnNamn == "Örebro", ]

## ---- eval=FALSE--------------------------------------------------------------
#  library(rgeos)
#  wkt <- writeWKT(as_Spatial(st_transform(shape, crs = st_crs(4326))))

## ---- warning=FALSE-----------------------------------------------------------
shape <- st_transform(shape,
                      crs = st_crs(4326))
shape <- st_union(shape)

## extract the convex hull of the polygon to reduce the length of the WKT string
shape_ch <- st_convex_hull(shape)

# cast it as MULTIPOLYGON as this is what SBDIs API need
# NOTE: as of today, the SBDI APIs will only work properly if the polygon is 
# submitted as a MULTIPOLYGON
shape_ch <- st_cast(shape_ch, to = "MULTIPOLYGON")

# create WKT string
wkt <- st_as_text(shape_ch)

## ----eval=FALSE---------------------------------------------------------------
#  species_list(wkt=wkt, fq="rank:species") %>%
#      dplyr::arrange(desc(occurrenceCount)) %>%
#      dplyr::select(species, family, occurrenceCount) %>%
#      head(10)

## ----message=FALSE, echo=FALSE------------------------------------------------
tryCatch({
  species_list(wkt=wkt, fq="rank:species") %>%
      dplyr::arrange(desc(occurrenceCount)) %>%
      dplyr::select(species, family, occurrenceCount) %>%
      head(10)
}, error = function(e) { print(e$message)})

## ----message=FALSE, warning=FALSE---------------------------------------------
library(vegan)

## ----eval=FALSE---------------------------------------------------------------
#  ## A rough polygon around the Mällardalen
#  wkt <- "MULTIPOLYGON(((14.94 58.88, 14.94 59.69, 18.92 59.69, 18.92 58.88, 14.94 58.88)))"
#  
#  ## define some environmental layers of interest [see sbdi_fields(fields_type = "occurrence")]
#  # el10011 https://spatial.biodiversitydata.se/ws/layers/view/more/worldclim_bio_12
#  # el10009 https://spatial.biodiversitydata.se/ws/layers/view/more/worldclim_bio_10
#  env_layers <- c("el10009","el10011")
#  
#  ## Download the data.  We use the `occurrences()` function, adding environmental
#  ## data via the 'extra' parameter.
#  x <- occurrences(fq="family:Fabaceae",
#                   wkt=wkt, qa="none",
#                   email="sbdi4r-test@biodiversitydata.se",
#                   download_reason_id="testing",
#                   extra=env_layers)

## ---- eval=FALSE--------------------------------------------------------------
#  library(dplyr)
#  library(tidyr)
#  xgridded <- x$data %>%
#      ## discard genus- and higher-level records
#      filter(rank %in%
#                    c("species", "subspecies", "variety", "form", "cultivar")) %>%
#      mutate(longitude=round(longitude*4)/4,
#             latitude=round(latitude*4)/4,
#             worldClimMeanTemperatureOfWarmestQuarter = worldClimMeanTemperatureOfWarmestQuarter /10) %>%
#      ## average environmental vars within each bin
#      group_by(longitude,latitude) %>%
#      mutate(worldClimAnnualPrecipitation = mean(worldClimAnnualPrecipitation, na.rm=TRUE),
#             worldClimMeanTemperatureOfWarmestQuarter = mean(worldClimMeanTemperatureOfWarmestQuarter, na.rm=TRUE)) %>%
#      ## subset to vars of interest
#      select(longitude, latitude, species,
#              worldClimAnnualPrecipitation,
#              worldClimMeanTemperatureOfWarmestQuarter) %>%
#      ## take one row per cell per species (presence)
#      distinct() %>%
#      ## calculate species richness
#      mutate(richness=n()) %>%
#      ## convert to wide format (sites by species)
#      mutate(present=1) %>%
#      do(tidyr::spread(data=., key=species, value=present, fill=0)) %>%
#      ungroup()
#  ## where a species was not present, it will have NA: convert these to 0
#  sppcols <- setdiff(names(xgridded),
#                     c("longitude", "latitude",
#                       "worldClimAnnualPrecipitation",
#                       "worldClimMeanTemperatureOfWarmestQuarter",
#                       "richness"))
#  xgridded <- xgridded %>%
#    mutate_at(sppcols, function(z) ifelse(is.na(z), 0, z))

## ---- include=FALSE-----------------------------------------------------------
## load data from a local copy so that vignette building doesn't require downloading a big chunk of data and slow sites-by-species processing
## this file generated by running the above unevaluated code blocks, then
## saveRDS(xgridded, file="vignette_fabaceae.rds")
xgridded <- readRDS("vignette_fabaceae.rds")
sppcols <- setdiff(names(xgridded), c("longitude", "latitude", 
                                      "worldClimAnnualPrecipitation", 
                                      "worldClimMeanTemperatureOfWarmestQuarter", 
                                      "richness"))

## ---- message=FALSE, warning=FALSE--------------------------------------------
xgridded[, 1:10]

## ---- warning=FALSE, fig.width=8, fig.height=6--------------------------------
library(ggplot2)
ggplot(xgridded, aes(longitude, richness)) + 
  labs(x = "Longitud (º)", 
       y = "Species richness") +
  lims(y = c(0,100)) +
  geom_point() + 
  theme_bw()

## ---- warning=FALSE, fig.width=8, fig.height=6--------------------------------
ggplot(xgridded, aes(worldClimMeanTemperatureOfWarmestQuarter , 
                     worldClimAnnualPrecipitation, 
                     colour=richness)) +
  labs(x = "Mean temperature of warmest quarter (ºC)" , 
       y = "Annual precipitation (mm)",
       colour = "Species \nrichness") + 
  scale_colour_distiller(palette="Spectral") +
  geom_point(size=3) + 
  theme_bw()

## ---- fig.width=6, fig.height=6, message=FALSE, warning=FALSE-----------------
library(vegan)
## Bray-Curtis dissimilarity
D <- vegdist(xgridded[, sppcols], "bray")
## UPGMA clustering
cl <- hclust(D, method="ave")
## plot the dendrogram
plot(cl)
## extract group labels at the 5-group level
grp <- cutree(cl, 5)
## coalesce small (outlier) groups into a single catch-all group
sing <- which(table(grp)<5)
# grp[grp %in% sing] <- 6 ## put these in a new combined group
grp <- sapply(grp, function(z)which(unique(grp)==z)) ## renumber groups
xgridded$grp <- as.factor(grp)
## plot
## colours for clusters
thiscol <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
             "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
ggplot(xgridded, aes(longitude, latitude, colour=grp)) + 
  labs(x="Longitude", y="Latitude", colour="Group") + 
  geom_point(size=3) +
  scale_colour_manual(values=thiscol) + 
  theme_bw()

## or a slightly nicer map plot
library(maps)
library(mapdata)
oldpar <- par()
par(mar = c(1,1,0,0))
map("worldHires", "Sweden", 
    xlim=c(14.5, 20), ylim=c(58.8, 59.95), 
    col="gray90", fill=TRUE)
with(xgridded, points(longitude, latitude, 
                      pch=21, col=thiscol[grp], 
                      bg=thiscol[grp], cex=0.75))
suppressWarnings(par(oldpar))

## ----Datahandling, message=FALSE, fig.width=8, fig.height=6-------------------
x <- occurrences(taxon="Callitriche cophocarpa",
                  fq = "data_resource_uid:dr2",
                  verbose = FALSE)
#keep spatially unique data at 0.01 degrees (latitude and longitude)
ll001 <- unique(x, spatial = 0.01)

# #keep only information for which fatal or "error" assertions do not exist
nofat <- subset(x, remove.fatal = TRUE)

#keep only observations with a maximum spatial uncertainty of 50m
SpatCert <- subset(x, max.spatial.uncertainty = 50)

# quickly get some more info about the data:
# no. observations (records)
nrow(x$data)
 
# no. obs/year
freq_year <- table(x$data$year)

# no. obs across years:
plot(freq_year, bty = "n", ylab="Frequency")
# or
hist(x$data$year, 20, main = "", xlab = "Year")

# Subsetting is done using '[ ]'
x10yr <- x$data[(x$data$year>=2010 & x$data$year<=2019),] 
table(x10yr$year)

## ----BIRDSspatial, message=FALSE, warning=FALSE-------------------------------
library(BIRDS)
x <- occurrences(taxon="Callitriche cophocarpa",
                  fq = "data_resource_uid:dr2",
                  verbose = FALSE)
wNAcoor <- which(is.na(x$data$latitude))
x$data <- x$data[-wNAcoor,]

### fill in NAs
x$data <- x$data[-which(is.na(x$data$year)),]
x$data$month <- ifelse(is.na(x$data$month), 1, x$data$month)
x$data$day <- ifelse(is.na(x$data$day), 1, x$data$day)

## Define the visit
OB <- organiseBirds(x$data[-wNAcoor,], 
                    sppCol = "scientificName", 
                    idCols = c("locality"), 
                    timeCols = c("year", "month","day"),
                    xyCols = c("longitude", "latitude"))

SB <- summariseBirds(OB, grid = Sweden_Grid_25km_Wgs84)

## ----plotBIRDSspatial, message=FALSE, warning=FALSE, fig.width=8, fig.height=6----
SBspatial <- SB$spatial
maxC <- max(SBspatial$nObs, na.rm = TRUE)
palBW <- leaflet::colorNumeric(c("white", "navyblue"), 
                               c(0, maxC), 
                               na.color = "transparent")
oldpar <- par()
par(mar = c(4,0,4,0), mfrow=c(1,3))
plot(SBspatial$geometry, col=palBW(SBspatial$nObs),
     border = "grey", main="All years") ## with palette
legend("topleft", inset = c(0,0.05),
       legend = round(seq(0, maxC, length.out = 5)),
       col = palBW(seq(0, maxC, length.out = 5)),
       title = "Number of \nobservations", pch = 15, bty="n")

## or export other combinations, e.g. one map per observed year
yearlySf <- exportBirds(SB, 
                        dimension = "spatial", 
                        timeRes = "yearly", 
                        variable = "nObs", 
                        method = "sum")

maxC <- max(st_drop_geometry(yearlySf), na.rm = TRUE)
palBW <- leaflet::colorNumeric(c("white", "navyblue"), 
                               c(0, maxC), 
                               na.color = "transparent")

plot(yearlySf$geometry, col = palBW(yearlySf$'1987'), 
     border = "grey", main = "1987")
legend("topleft", inset = c(0,0.05),
       legend = round(seq(0, maxC, length.out = 5)),
       col = palBW(seq(0, maxC, length.out = 5)),
       title = "Number of \nobservations", pch = 15, bty = "n")

plot(yearlySf$geometry, col = palBW(yearlySf$'2017'), 
     border = "grey", main = "2017")
suppressWarnings(par(oldpar))

## ----BIRDSsave, message=FALSE, warning=FALSE, eval=FALSE----------------------
#  gridSummary <- st_drop_geometry(SB$spatial)
#  write.csv(gridSummary, "Callitriche_grid_summary.csv")
#  

## ----BIRDStemporal, message=FALSE, warning=FALSE, fig.width=8, fig.height=6----
# But the function exportBIrds() offers planty of combinations
yearlyXTS <- exportBirds(SB, 
                         dimension = "temporal", 
                         timeRes = "yearly", 
                         variable = "nObs",
                         method = "sum")

plot(yearlyXTS, col = "darkblue", 
     grid.ticks.on = "year",  
     grid.col = "lightgrey",  
     main = "Number of observations")


