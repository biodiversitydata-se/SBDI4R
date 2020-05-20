# SBDI4R

R functionality for SBDI data portal

The Swedish Biodiversity Data Infrastructure (SBDI) provides tools to enable users of biodiversity information to find, access, combine and visualise data on Swedish plants and animals; these have been made available from [Biodiversity Atlas Sweden](https://bioatlas.se/). Here we provide a subset of the tools, and some extension tools (found previously in Analysportalen.se) to be directly used within R.

SBDI4R enables the R community to directly access data and resources hosted by the Biodiversity Atlas Sweden. Our goal is to enable outputs (e.g. observations of species) to be queried and output in a range of standard formats. This tool is built on the Atlas of Living Australia [ALA4R](https://github.com/AtlasOfLivingAustralia/ALA4R) package which provides similar services for the ALA and. This package also builds upon a similar package [NBN4R](https://github.com/fozy81/NBN4R) that wraps ALA4R functions but redirects requests to different web servers, in this case, SBDI servers. Both SBDI, NBN and ALA share similar Application Protocol Interface (API) web services. 

The use-examples based on ALA4R are presented at the [2014 ALA Science Symposium](http://www.ala.org.au/blogs-news/2014-atlas-of-living-australia-science-symposium/) are available in the package vignette, via (in R): `vignette("SBDI4R")`, and a draft modifed version using NBN data is below.

## Installing SBDI4R

### Windows

In R:

Or the development version from GitHub:

```R
install.packages("remotes")
library(remotes)
install_github("bioatlas/SBDI4R")
```

If you see an error about a missing package, you will need to install it manually, e.g.:
```R
install.packages(c("stringr","sp"))
```
and then `install_github("bioatlas/SBDI4R")` again.


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
install_github("bioatlas/SBDI4R")
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
The SBDI4R package must be loaded for each new R session:
```R
# library(SBDI4R)
```

##Customizing SBDI4R
Various aspects of the SBDI4R package can be customized.

###Caching
SBDI4R can cache most results to local files. This means that if the same code is run multiple times, the second and subsequent iterations will be faster. This will also reduce load on the NBN servers.
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
setHook(packageEvent("SBDI4R", "attach"), 
function(...) sbdi_config(cache_directory=file.path("~","mydata","sbdi_cache")))
```
Caching can also be turned off entirely by:
```R
sbdi_config(caching="off")
```
or set to “refresh”, meaning that the cached results will re-downloaded from the SBDI servers and the cache updated. (This will happen for as long as caching is set to “refresh” — so you may wish to switch back to normal “on” caching behaviour once you have updated your cache with the data you are working on).

###User-agent string
Each request to SBDI servers is accompanied by a “user-agent” string that identifies the software making the request. This is a standard behaviour used by web browsers as well. The user-agent identifies the user requests to SBDI, helping SBDI to adapt and enhance the services that it provides. By default, the SBDI4R user-agent string is set to “SBDI4R” plus the SBDI4R version number (e.g. “SBDI4R 1.0”).

###E-mail address
Each request to SBDI servers is also accompanied by an “e-mail address” string that identifies the user making the request. This is a standard behaviour used by web browsers as well. There is no default for this field, and unless specified in the package configuration, needs to be specified for each request using the function 'occurrences()'.
You can set this up by:
```R
sbdi_config(email="your.valid@emailaddress.com")
```

*NO* other personal identification information is sent. You can see all configuration settings, including the the user-agent string that is being used, with the command:
```R
sbdi_config()
```

###Debugging
If things aren’t working as expected, more detail (particularly about web requests and caching behaviour) can be obtained by setting the verbose configuration option:
```R
sbdi_config(verbose=TRUE)
```

### Setting the download reason
SBDI requires that you provide a reason when downloading occurrence data (via the SBDI4R `occurrences()` function). You can provide this as a parameter directly to each call of `occurrences()`, or you can set it once per session using:

```R
sbdi_config(download_reason_id=your_reason_id)
```

(See `sbdi_reasons()` for valid download reasons)


### Other options
If you make a request that returns an empty result set (e.g. an un-matched name), by default you will simply get an empty data structure returned to you without any special notification. If you would like to be warned about empty result sets, you can use:

```R
sbdi_config(warn_on_empty=TRUE)
```

##Examples
First, check that we have some additional packages that we'll use in the examples, and install them if necessary.
```R
to_install <- c("plyr","jpeg","phytools","ape","leaflet","vegan","mgcv","geosphere","maps","mapdata","maptools")
to_install <- to_install[!sapply(to_install, requireNamespace, quietly=TRUE)]
if (length(to_install) > 0) install.packages(to_install)
```

We’ll use the `plyr` package throughout these examples, so load that now:
```R
library(plyr) 
```

###Example 1: Name searching and taxonomic trees

```R
library(ape)
library(phytools)
```

We want to look at the taxonomic finches, but we don’t know what the correct scientific name is, so let’s search for it:
```R
sx <- search_fulltext("parus", )
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
##TO UPDATE
And we can see that violets correspond to a number of families. We want to use Estrildidae. Now we can download the taxonomic data (note that the search is case-sensitive):
```R
tx <- taxinfo_download("family_s:Paridae", fields=c("guid", "genus_s", "scientificName", "rank"))
tx <- tx[tx$rank %in% c("species","subspecies"),] ## restrict to species and subspecies
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

###Example 2: Quality assertions
Data quality assertions are a suite of fields that are the result of a set of tests peformed on NBN data. Download occurrence data for the Blunt-fruited Water-starwort:
```R
x <- occurrences(taxon="Callitriche obtusangula", download_reason_id=10)
summary(x)
```
```
# number of original names: 2 
# number of taxonomically corrected names: 1 
# number of observation records: 1319 
# number of assertions listed: 9  -- ones with flagged issues are listed below
# 	invalidCollectionDate: 1 records 
# 	incompleteCollectionDate: 1 records 
# 	precisionRangeMismatch: 1040 records 
# 	firstOfYear: 338 records 
# 	unknownCountry: 474 records 
# 	countryInferredByCoordinates: 806 records 
# 	decimalLatLongCalculatedFromGridReference: 1298 records 
# 	assumedPresentOccurrenceStatus: 992 records 
# 	firstOfMonth: 479 records 
```

There are many other ways of producing spatial plots in R. The `leaflet` package provides a simple method of producing browser-based maps iwth panning, zooming, and background layers:

```R
library(leaflet)
## drop any records with missing lat/lon values
x$data <- x$data[!is.na(x$data$longitude) & !is.na(x$data$latitude),] 
xa <- check_assertions(x)
## columns of x corresponding to a fatal assertion
x_afcols <- which(names(x$data) %in% xa$occurColnames[xa$taxonIdentificationIssue])
## rows of x that have a fatal assertion
x_afrows <- apply(x$data[,x_afcols],1,any)
## which taxonIdentificationIssue assertions are present in this data?
these_assertions <- names(x$data)[x_afcols]
## make a link to the web page for each occurrence
popup_link <- paste0("<a href=\"https://records.nbnatlas.org/occurrences/",x$data$id,"\">Link to occurrence record</a>")
## colour palette
pal <- c(sub("FF$","",heat.colors(length(these_assertions))))
## map each data row to colour, depending on its assertions
marker_colour <- rep("#00FF00",nrow(x$data))
for (k in 1:length(these_assertions)) marker_colour[x$data[,x_afcols[k]]] <- pal[k]
## blank map, with imagery background
m <- addProviderTiles(leaflet(),"Esri.WorldImagery")
## add markers
m <- addCircleMarkers(m,x$data$longitude,x$data$latitude,col=marker_colour,popup=popup_link)
print(m)
```



