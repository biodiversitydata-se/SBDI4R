# SBDI4R <img src="https://github.com/biodiversitydata-se/SBDI4R/raw/master/man/figures/SBDI.png" align="right" width="120"/>

[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R-CMD-check](https://github.com/biodiversitydata-se/SBDI4R/workflows/R-CMD-check/badge.svg)](https://github.com/biodiversitydata-se/SBDI4R/actions)


### R functionality for the SBDI data portal

The Swedish Biodiversity Data Infrastructure (SBDI) provides tools to
enable users of biodiversity information to find, access, combine and
visualize data on Swedish plants and animals; available through [Swedish
Biodiversity Data Infrastructure](https://biodiversitydata.se/). The R
package SBDI4R provides a subset of the tools, and some extension tools
(found previously in Analysportalen.se), to be directly used within R.

SBDI4R enables the R community to directly access data and resources
hosted by SBDI. Our goal is to enable observations of species to be
queried and output in a range of standard formats. This tool is built on
the Atlas of Living Australia
[ALA4R](https://github.com/AtlasOfLivingAustralia/ALA4R) package which
provides similar services for the ALA. Similar to the
[NBN4R](https://github.com/fozy81/NBN4R) package SBDI4R wraps ALA4R
functions but redirects requests to local web servers. All SBDI, NBN and
ALA share similar Application Protocol Interface (API) web services.

Use-examples are available in the package [vignette
here](https://biodiversitydata-se.github.io/SBDI4R/articles/SBDI4R.html),
or via (in R): `vignette("SBDI4R")`. If you have any questions please
contact the author and maintainer Debora Arlt via
[email](mailto:debora.arlt@slu.se?subject=%5BGitHub%5D%20SBDI4R).

## Installing SBDI4R

### Windows

In R:

Or the development version from GitHub:

```{r}
install.packages("remotes")
library(remotes)
install_github("biodiversitydata-se/SBDI4R")
```

If you see an error about a missing package, you will need to install it
manually, e.g.:

```{r}
install.packages(c("stringr","sp"))
```

and then `install_github("biodiversitydata-se/SBDI4R")` again.

If you see an error about "ERROR: lazy loading failed for package
'SBDI4R'", this may be due to you trying to install on a network
location. Try instead to install on a local location: first create the
local location you want to use, and then specify this location for
installing, and later loading the package:

```{r}
install_github("biodiversitydata-se/SBDI4R", lib = "C:/pathname/MyLibrary")
library(SBDI4R, lib.loc = "C:/pathname/MyLibrary")
```

If you wish to use the `data.table` package for potentially faster
loading of data matrices (optional), also do:

```{r}
install.packages("data.table")
```

### Mac

Follow the instructions for Windows.

If you see an error about a failure to set default locale, you will need
to manually set this:

```{r}
system('defaults write org.R-project.R force.LANG en_US.UTF-8')
```

and restart R.

More information can be found on the [CRAN R for Mac
page](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Internationalization-of-the-R_002eapp).

### Linux

First, ensure that `libcurl` is installed on your system --- e.g. on
Ubuntu, open a terminal and do:

``` bash
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

If you see an error about a missing package, you will need to install it
manually, e.g.:

```{r}
install.packages(c("stringr","sp"))
```

and then try installing SBDI4R again.

If you wish to use the `data.table` package for potentially faster
loading of data matrices (optional), also do:

```{r}
install.packages("data.table")
```

## Using SBDI4R

The SBDI4R package must be loaded for each new R session with
`library(SBDI4R)`, or specifying your local location with
`library(SBDI4R, lib.loc = "C:/pathname/MyLibrary")`.

## Customizing SBDI4R

Various aspects of the SBDI4R package can be customized.

### Caching

SBDI4R can cache most results to local files. This means that if the
same code is run multiple times, the second and subsequent iterations
will be faster. This will also reduce load on the web servers. By
default, this caching is session-based, meaning that the local files are
stored in a temporary directory that is automatically deleted when the R
session is ended. This behaviour can be altered so that caching is
permanent, by setting the caching directory to a non-temporary location.
For example, under Windows, use something like:

```{r}
sbdi_config(cache_directory = file.path("c:","mydata","sbdi_cache")) ## Windows
```

or for Linux:

```{r}
sbdi_config(cache_directory = "~/mydata/sbdi_cache") ## Linux
```

Note that this directory must exist (you need to create it yourself).

All results will be stored in that cache directory and will be used from
one session to the next. They won't be re-downloaded from the server
unless the user specifically deletes those files or changes the caching
setting to "refresh".

If you change the cache_directory to a permanent location, you may wish
to add something like this to your .Rprofile file, so that it happens
automatically each time the SBDI4R package is loaded:

```{r}
setHook(packageEvent("SBDI4R", "onLoad"), 
        function(...) sbdi_config(cache_directory=file.path("~","mydata","sbdi_cache")))
```

Caching can also be turned off entirely by:

```{r}
sbdi_config(caching="off")
```

or set to "refresh", meaning that the cached results will re-downloaded
from the SBDI servers and the cache updated. (This will happen for as
long as caching is set to "refresh" — so you may wish to switch back to
normal "on" caching behaviour once you have updated your cache with the
data you are working on).

### User-agent string

Each request to SBDI servers is accompanied by a "user-agent" string
that identifies the software making the request. This is a standard
behaviour used by web browsers as well. The user-agent identifies the
user requests to SBDI, helping SBDI to adapt and enhance the services
that it provides. By default, the SBDI4R user-agent string is set to
"SBDI4R" plus the SBDI4R version number (e.g. "SBDI4R 1.0").

### E-mail address

Each request to SBDI servers is also accompanied by an "e-mail address"
string that identifies the user making the request. This is a standard
behaviour used by web browsers as well. There is no default for this
field. You can provide your e-mail address as a parameter directly to
each call of the function occurrences(), or you can set it once per
session specifying it in the package configuration:

```{r}
sbdi_config(email="your.valid@emailaddress.com")
```

*NO* other personal identification information is sent. You can see all
configuration settings, including the the user-agent string that is
being used, with the command:

```{r}
sbdi_config()
```

### Debugging

If things aren't working as expected, more detail (particularly about
web requests and caching behaviour) can be obtained by setting the
verbose configuration option:

```{r}
sbdi_config(verbose=TRUE)
```

### Setting the download reason

SBDI requires that you provide a reason when downloading occurrence data
(via the SBDI4R `occurrences()` function). You can provide this as a
parameter directly to each call of `occurrences()`, or you can set it
once per session using:

```{r}
sbdi_config(download_reason_id=your_reason_id)
```

(See `sbdi_reasons()` for valid download reasons, e.g.
download_reason_id=10 for "testing", or 7 for "ecological research", 8
for "systematic research/taxonomy", 3 for "education")

### Other options

If you make a request that returns an empty result set (e.g. an
un-matched name), by default you will simply get an empty data structure
returned to you without any special notification. If you would like to
be warned about empty result sets, you can use:

```{r}
sbdi_config(warn_on_empty=TRUE)
```

### See examples on how to use the package in the next [vignette](https://biodiversitydata-se.github.io/SBDI4R/articles/SBDI4R.html)
