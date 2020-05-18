#' Sites by species
#' 
#' A data.frame is returned as grid cells by species with values in each cell being the number of occurrences of each species. No null (all zero) species should be returned. The coordinates returned are the TOP-LEFT corner of the grid cell.
#'
#' @references Associated web services: \url{https://spatial.bioatlas.se/ws}
#' @references \url{http://www.geoapi.org/3.0/javadoc/org/opengis/referencing/doc-files/WKT.html}
#' 
#' @param taxon string: the identifier to get the species data from the nbn biocache. E.g. "genus:Macropus". See \code{sbdi_fields("occurrence_stored")} for valid field names
#' @param wkt string: Bounding area in Well Known Text (WKT) format. E.g. "POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))".
#' @param gridsize numeric: size of output grid cells in decimal degrees. E.g. 0.1 (=~10km)
#' @param SPdata.frame logical: should the output be returned as a SpatialPointsDataFrame of the sp package?
#' @param verbose logical: show additional progress information? [default is set by sbdi_config()]
#' @return A dataframe or a SpatialPointsDataFrame containing the species by sites data. Columns will include longitude, latitude, and each species present. Values for species are record counts (i.e. number of recorded occurrences of that taxon in each grid cell). The \code{guid} attribute of the data frame gives the guids of the species (columns) as a named character vector
#'
#' @examples
#' \dontrun{
#' # Otters in Scotland based on a 0.1 degree grid
#' ss <- sites_by_species(taxon="genus:Leutra",wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))
#' ",gridsize=0.1,verbose=TRUE)
#' head(ss[,1:6])
#' # equivalent direct POST webservice call:
#' # https://spatial.bioatlas.se/ws/sitesbyspecies?speciesq=genus:Viola&qname=data&
#' #  area=POLYGON((-3%2056,-4%2056,-4%2057,-3%2047,-3%2056))&bs=- https://records.bioatlas.se/ws/
#' #  ws/&movingaveragesize=1&gridsize=0.1&sitesbyspecies=1
#'
#' ## get the guid of the first species (which is the third column of the data frame, since the
#' ##  first two columns are longitude and latitude)
#' attr(ss,"guid")[1]
#' }

# TODO need way to better check input species query. If the query is incorrect, the call will fail with message along the lines of: "Error in sites_by_species(taxon = "gen:Eucalyptus", wkt = "POLYGON((144 -43,148 -43,148 -40,144 -40,144 -43))",  :   Error processing your Sites By Species request. Please try again or if problem persists, contact the Administrator."

## fails: ss <- sites_by_species(taxon="rk_genus:Eucalyptus",wkt="POLYGON((144 -43,148 -43,148 -40,144 -40,144 -43))",verbose=TRUE) 
## fails: ss <- sites_by_species(taxon="scientificNameAuthorship:Maiden",wkt="POLYGON((144 -43,148 -43,148 -40,144 -40,144 -43))",verbose=TRUE)
## fails: ss <- sites_by_species(taxon="parentGuid:http://id.biodiversity.org.au/node/apni/6337078",wkt="POLYGON((144 -43,148 -43,148 -40,144 -40,144 -43))",verbose=TRUE)

#' #@export
sites_by_species <- function(taxon, wkt, gridsize = 0.1, SPdata.frame = FALSE, verbose = sbdi_config()$verbose) {
  
  ALA4R::sites_by_species(taxon, wkt, gridsize, SPdata.frame, verbose )
  
}