#' Fetch a taxon profile given a scientific name or LSID (GUID)
#' 
#' @references Associated SBDI web service: \url{https://api.bioatlas.se/#ws80}
#' 
#' @param scientificname string: scientific name of the taxon of interest (species, genus, family etc) 
#' @param guid string: The Life Science Identifier of the taxon of interest
#' @param verbose logical: show additional progress information? [default is set by sbdi_config()]
#' @return A species profile in the form of a named list, each element of which is generally a data frame. An empty list is returned if no match is found for the supplied name or guid
#' @seealso \code{\link{sbdi_config}}
#' @examples
#' \dontrun{
#'  species_info("Diatoma tenuis")
#'  species_info(guid="NHMSYS0000080188")
#'  species_info("Diatoma",verbose=TRUE)
#' }
#' @export species_info
species_info <- function(scientificname, guid, verbose = sbdi_config()$verbose) {
  
  ALA4R::species_info(scientificname, guid, verbose)
  
}