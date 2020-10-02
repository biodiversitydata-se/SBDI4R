# SBDI4R

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![Build Status](https://travis-ci.com/biodiversitydata-se/SBDI4R.svg?branch=master)](https://travis-ci.com/biodiversitydata-se/SBDI4R)

R functionality for SBDI data portal

The Swedish Biodiversity Data Infrastructure (SBDI) provides tools to enable users of biodiversity information to find, access, combine and visualise data on Swedish plants and animals; available through [Biodiversity Atlas Sweden](https://bioatlas.se/). The R package SBDI4R provides a subset of the tools, and some extension tools (found previously in Analysportalen.se), to be directly used within R.

SBDI4R enables the R community to directly access data and resources hosted by SBDI. Our goal is to enable observations of species to be queried and output in a range of standard formats. This tool is built on the Atlas of Living Australia [ALA4R](https://github.com/AtlasOfLivingAustralia/ALA4R) package which provides similar services for the ALA. Similar to the [NBN4R](https://github.com/fozy81/NBN4R) package SBDI4R wraps ALA4R functions but redirects requests to local web servers. All SBDI, NBN and ALA share similar Application Protocol Interface (API) web services. 

The use-examples based on ALA4R are presented at the [2014 ALA Science Symposium](http://www.ala.org.au/blogs-news/2014-atlas-of-living-australia-science-symposium/) are available in the package vignette, via (in R): `vignette("SBDI4R")`, and a draft modifed version using NBN data is found below.

## Installing SBDI4R

### Windows

In R:

Or the development version from GitHub:

```R
install.packages("remotes")
library(remotes)
install_github("biodiversitydata-se/SBDI4R")
```

If you see an error about a missing package, you will need to install it manually, e.g.:
```R
install.packages(c("stringr","sp"))
```
and then `install_github("biodiversitydata-se/SBDI4R")` again.

If you see an error about "ERROR: lazy loading failed for package 'SBDI4R'", this may be due to you trying to install on a network location. 
Try instead to install on a local location: first create the local location you want to use, and then specify this location for isntalling, and later loading the package:
```R
install_github("biodiversitydata-se/SBDI4R", lib = "C:/pathname/MyLibrary")
library(SBDI4R, lib.loc = "C:/pathname/MyLibrary")
```

If you wish to use the `data.table` package for potentially faster loading of data matrices (optional), also do:
```R
install.packages("data.table")
```

### Linux

First, ensure that `libcurl` is installed on your system --- e.g. on Ubuntu, open a terminal and do:
```BASH
sudo apt-get install libcurl4-openssl-dev
```
or install `libcurl4-openssl-dev` via the Software Centre.

Then, in R:

Or the development version from GitHub:

```{r eval=FALSE}
install.packages("remotes")
library(remotes)
install_github("biodiversitydata-se/SBDI4R")
```

If you see an error about a missing package, you will need to install it manually, e.g.:
```R
install.packages(c("stringr","sp"))
```
and then try installing SBDI4R again.


If you wish to use the `data.table` package for potentially faster loading of data matrices (optional), also do:
```R
install.packages("data.table")
```

## Using SBDI4R  
The SBDI4R package must be loaded for each new R session with `library(SBDI4R)`,
or specifying your local location with `library(SBDI4R, lib.loc = "C:/pathname/MyLibrary")`.

## Customizing SBDI4R  
Various aspects of the SBDI4R package can be customized.

### Caching  
SBDI4R can cache most results to local files. This means that if the same code is run multiple times, the second and subsequent iterations will be faster. This will also reduce load on the web servers.
By default, this caching is session-based, meaning that the local files are stored in a temporary directory that is automatically deleted when the R session is ended. This behaviour can be altered so that caching is permanent, by setting the caching directory to a non-temporary location. For example, under Windows, use something like:
```R
sbdi_config(cache_directory = file.path("c:","mydata","sbdi_cache")) ## Windows
```
or for Linux:
```R
sbdi_config(cache_directory = "~/mydata/sbdi_cache") ## Linux
```
Note that this directory must exist (you need to create it yourself).

All results will be stored in that cache directory and will be used from one session to the next. They won’t be re-downloaded from the server unless the user specifically deletes those files or changes the caching setting to “refresh”.

If you change the cache_directory to a permanent location, you may wish to add something like this to your .Rprofile file, so that it happens automatically each time the SBDI4R package is loaded:
```R
setHook(packageEvent("SBDI4R", "onLoad"), 
        function(...) sbdi_config(cache_directory=file.path("~","mydata","sbdi_cache")))
```
Caching can also be turned off entirely by:
```R
sbdi_config(caching="off")
```
or set to “refresh”, meaning that the cached results will re-downloaded from the SBDI servers and the cache updated. (This will happen for as long as caching is set to “refresh” — so you may wish to switch back to normal “on” caching behaviour once you have updated your cache with the data you are working on).

### User-agent string  
Each request to SBDI servers is accompanied by a “user-agent” string that identifies the software making the request. This is a standard behaviour used by web browsers as well. The user-agent identifies the user requests to SBDI, helping SBDI to adapt and enhance the services that it provides. By default, the SBDI4R user-agent string is set to “SBDI4R” plus the SBDI4R version number (e.g. “SBDI4R 1.0”).

### E-mail address  
Each request to SBDI servers is also accompanied by an “e-mail address” string that identifies the user making the request. This is a standard behaviour used by web browsers as well. There is no default for this field. You can provide your e-mail address as a parameter directly to each call of the function occurrences(), or you can set it once per session specifying it in the package configuration:
```R
sbdi_config(email="your.valid@emailaddress.com")
```

*NO* other personal identification information is sent. You can see all configuration settings, including the the user-agent string that is being used, with the command:
```R
sbdi_config()
```

### Debugging  
If things aren’t working as expected, more detail (particularly about web requests and caching behaviour) can be obtained by setting the verbose configuration option:
```R
sbdi_config(verbose=TRUE)
```

### Setting the download reason  
SBDI requires that you provide a reason when downloading occurrence data (via the SBDI4R `occurrences()` function). You can provide this as a parameter directly to each call of `occurrences()`, or you can set it once per session using:
```R
sbdi_config(download_reason_id=your_reason_id)
```

(See `sbdi_reasons()` for valid download reasons, 
e.g. download_reason_id=10 for "testing", or 7 for "ecological research", 8 for "systematic research/taxonomy", 3 for "education")


### Other options
If you make a request that returns an empty result set (e.g. an un-matched name), by default you will simply get an empty data structure returned to you without any special notification. If you would like to be warned about empty result sets, you can use:

```R
sbdi_config(warn_on_empty=TRUE)
```

## Examples  
First, check that we have some additional packages that we'll use in the examples, and install them if necessary.
```R
to_install <- c("plyr","phytools","ape","leaflet")
to_install <- to_install[!sapply(to_install, requireNamespace, quietly=TRUE)]
if (length(to_install) > 0) install.packages(to_install)
```

```R
library(plyr) 
library(ape)
library(phytools)
```

### Example 1: Name searching and taxonomic trees
We want to look at the taxonomy of titmice, but we don’t know what the correct scientific name is, so let’s search for it:
```R
sx <- search_fulltext("parus")
(sx$data[,c( "name","species", "speciesGuid", "rank")])
```
```
#                    name             species speciesGuid    rank
# 1  Parus Linnaeus, 1758                <NA>        <NA>   genus
# 2        Parus dichrous      Parus dichrous     5788759 species
# 3     Neuroctenus parus   Neuroctenus parus    10213928 species
# 4         Parus elegans       Parus elegans     5788758 species
# 5   Parus rubidiventris Parus rubidiventris     5788753 species
# 6        Parus borealis      Parus borealis     8786895 species
# 7      Parus salicarius    Parus salicarius     9104351 species
# 8   Parus superciliosus Parus superciliosus     5788752 species
# 9        Parus amabilis      Parus amabilis     5788755 species
# 10   Parus rufonuchalis  Parus rufonuchalis     5788756 species
```

But we see some e.g. insects (**Neuroctenus parus**) are also returned. We want to restrict the search to Paridae.
```R
sx <- search_fulltext("parus", fq="family_s:Paridae")
(sx$data[,c( "name","species", "speciesGuid", "rank")])
```
To restrict the query specifically to birds we can also use the fq argument to filter the query (see sbdi_fields("general",as_is=TRUE) for all the fields that are queryable), and increase page_size to include more records (default=10):
```R
sx <- search_fulltext("parus", fq="class_s:Aves", page_size=100)
(sx$data[,c( "name","species", "speciesGuid", "rank")])
```

Now we can download the taxonomic data (note that the search is case-sensitive):
```R
tx <- taxinfo_download("family_s:Paridae", fields=c("guid", "genus_s", "scientificName", "rank"))
tx <- tx[tx$rank == "species",] ## restrict to species
```
We can make a taxonomic tree plot using the `phytools` package:
```R
## as.phylo requires the taxonomic columns to be factors
tx$genusS<-as.factor(tx$genusS)
tx$scientificName<-as.factor(tx$scientificName)
## create phylo object of Scientific.Name nested within Genus
ax <- as.phylo(~genusS/scientificName, data=tx)
plotTree(ax, fsize=0.7, ftype="i") ## plot it
```

### Example 2: Get some data, quality assertions, plotting data on a map and save data  
Download occurrence data for the Blunt-fruited Water-starwort and view top of the data table:
```R
x <- occurrences(taxon="Callitriche cophocarpa", download_reason_id=10)
head(x$data)
```

#### Quality assertions
Data quality assertions are a suite of fields that are the result of a set of tests performed on data. We continue using the data for the Blunt-fruited Water-starwort and get a summary of the data quality assertions:
```R
summary(x)

# number of original names: 3 
# number of taxonomically corrected names: 1 
# number of observation records: 5634 
# number of assertions listed: 16  -- ones with flagged issues are listed below
# 	invalidCollectionDate: 45 records 
# 	zeroCoordinates: 2 records -- considered fatal
# 	incompleteCollectionDate: 271 records 
# 	unrecognisedInstitutionCode: 1732 records 
# 	uncertaintyRangeMismatch: 1 records 
# 	firstOfYear: 24 records 
# 	zeroLatitude: 2 records -- considered fatal
# 	geodeticDatumAssumedWgs84: 5034 records 
# 	unrecognisedCollectionCode: 1732 records 
# 	uncertaintyInPrecision: 81 records 
# 	countryCoordinateMismatch: 2 records 
# 	assumedPresentOccurrenceStatus: 1180 records 
# 	zeroLongitude: 2 records -- considered fatal
# 	idPreOccurrence: 4 records 
# 	recordedByUnparsable: 3 records 
# 	firstOfMonth: 131 records
```  

You can see a list of all record issues using `sbdi_fields("assertions",as_is=TRUE)` and see what is considered as fatal quality issues.

#### Plotting data on a map  
You can quickly plot all the observations with the function `ocurrence_plot()`, here we specify to map all fatal issues:
```R
occurrences_plot(x, "obsPlot.pdf", qa="fatal", 
                  grouped=FALSE, taxon_level="species", 
                  pch='+')
```
Note that the plot is saved to a pdf file in the current working directory. You can find that by `getwd()`.  

<img src=https://github.com/biodiversitydata-se/SBDI4R/tree/master/man/figures/obsPlot_CallitricheCophocarpa.png />

There are many other ways of producing spatial plots in R. The `leaflet` package provides a simple method of producing browser-based maps with panning, zooming, and background layers:
```R
library(leaflet)
## drop any records with missing lat/lon values
x$data <- x$data[!is.na(x$data$longitude) & !is.na(x$data$latitude),] 
xa <- check_assertions(x)
## columns of x corresponding to a fatal assertion
x_afcols <- which(names(x$data) %in% xa$occurColnames[xa$isFatal])
## rows of x that have a fatal assertion
x_afrows <- apply(x$data[,x_afcols], 1, any)
## which taxonIdentificationIssue assertions are present in this data?
these_assertions <- names(x$data)[x_afcols]
## make a link to the web page for each occurrence
popup_link <- paste0("<a href=\"https://records.bioatlas.se/occurrences/",
                      x$data$id,"\">Link to occurrence record</a>")
## colour palette
pal <- c(sub("FF$","", heat.colors(length(these_assertions))))
## map each data row to colour, depending on its assertions
marker_colour <- rep("#00FF00", nrow(x$data))
for (k in 1:length(these_assertions)) marker_colour[x$data[,x_afcols[k]]] <- pal[k]
## blank map, with imagery background
m <- addProviderTiles(leaflet(),"Esri.WorldImagery")
## add markers
m <- addCircleMarkers(m, x$data$longitude, x$data$latitude,  radius = 5,
                      col=marker_colour, popup=popup_link)
m
```  

#### Save data
```R
# save as data.frame
Callitriche <- as.data.frame(x$data)

# simplyfy data frame
calli <- data.frame(Callitriche$scientificName,
                   Callitriche$latitude,
                   Callitriche$longitude)
# simplify column names
colnames(calli) <- c("species","latitude","longitude")
# remove rows with missing values (NAs)
calli <- na.omit(calli)

# save new dataframe
write.csv(calli,"Callitriche.csv")
```


### Example 3: Summarise occurrences over a defined grid
Now we want to summarise occurrences over a defined grid instead of plotting every 
observation point. For this example you will need the `sp` and `rgeso` packages. 
Install them if you haven´t with `install.packages("sp", "rgeos")`.
First we need to overlay the observations with the grid:  

```R
library(sp) # the function coordinates() and proj4string() are in sp
library(rgeos) #  the function over() is in package rgeos
## load some shapes over Sweden
data("swe_wgs84", package="SBDI4R", envir=environment()) # Political borders
data("swe100kmGrid", package="SBDI4R", envir=environment()) # A grid
grid <- swe100kmGrid

## make the observations spatial
obs <- as.data.frame(x$data)
coordinates(obs) <- obs[,c("longitude","latitude")]
proj4string(obs) <- "+init=epsg:4326"

nObs <- nrow(obs)

## overlay the data with the grid
ObsInGridList <- over(grid, obs, returnList=TRUE)
wNonEmpty <- unname( which( unlist(lapply(ObsInGridList, nrow)) != 0) )
if(length(wNonEmpty)==0) message("Observations don't overlap any grid cell.")

## check nObs
nObsInGrid <- sum(unlist(lapply(ObsInGridList, nrow)))
```

Now summarise occurrences within grid cells:
```R
## apply a summary over the grid
nCells <- length(ObsInGridList)

res <- data.frame("nObs"=as.numeric(rep(NA,nCells)),
                  "nYears"=as.numeric(rep(NA,nCells)),
                  stringsAsFactors = FALSE)

cols2use<-c("scientificName", "year")

dataRes<-lapply(ObsInGridList[wNonEmpty], function(x){
  x<-x[,cols2use]
  colnames(x) <- c("scientificName", "year")
  
  return(c("nObs" = length(x[,"scientificName"]),
           "nYears" = length(unique(x[,"year"]))
  ))
})

dataRes <- data.frame(matrix(unlist(dataRes),
                             nrow=length(dataRes), 
                             byrow=TRUE),
                    stringsAsFactors = FALSE)

dataRes$X1 <- as.numeric(dataRes$X1)
dataRes$X2 <- as.numeric(dataRes$X2)

res[wNonEmpty,] <- dataRes
rownames(res) <- row.names(grid)

resSp <- sp::SpatialPolygonsDataFrame(grid, res)
```

Finally plot the grid summary as a map:
```R
palBW <- leaflet::colorNumeric(c("white", "navyblue"), 
                               c(0, max(resSp@data$nObs, na.rm = TRUE)), 
                               na.color = "transparent")
plot(resSp, col=palBW(resSp@data$nObs), border = NA)
plot(swe_wgs84$Border, border=1, lwd=1, add=T)
legend("bottomleft", 
       legend = round(seq(0, max(resSp@data$nObs, na.rm = TRUE), length.out = 5)),
       col = palBW(seq(0, max(resSp@data$nObs, na.rm = TRUE), length.out = 5)),
       title = "Number of observations", pch = 15, bty="n")
```

## More examples - case studies

### Case 1: 
coming soon...

