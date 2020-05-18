#' Species lists
#'
#' Note that this refers to pre-generated lists of species stored on the SBDI servers. The similarly-named but different function \code{\link{specieslist}} provides a different function, namely listing the species matching a query or recorded as present in a search area.
#'
#' @references \url{https://lists.bioatlas.se} and the associated web services at \url{https://lists.bioatlas.se/ws}
#' @param druid string: data resource UID of the list (i.e. the list identifier)
#' @param kvp logical: include key-value pairs? Some lists contain information about the species in the form of key-value pairs
#' @param verbose logical: show additional progress information? 
#'
#' @return data.frame
#'
#' @seealso \code{\link{specieslist}} \code{\link{sbdi_lists}}
#'
#' @examples
#' \dontrun{
#'  all_lists <- sbdi_lists()
#'  ## find the "Field Guide apps species profiles" from Museum Victoria
#'  all_lists[grep("Field Guide",all_lists$listName),]
#'  ## download the vertebrates one
#'  l <- sbdi_list(druid="dr1146")
#' }
#'
#' @export
sbdi_list <- function(druid,kvp=TRUE,verbose=sbdi_config()$verbose){
  
  ALA4R::ala_list(druid,kvp,verbose)
}


#' Find SBDI species lists
#'
#' @references \url{https://lists.bioatlas.se} and the associated web services at \url{https://lists.bioatlas.se/ws}
#' @param guid string: (optional) if provided, return only lists in which this GUID appears
#' @param offset integer: the number of lists to skip. This supports paging
#' @param max integer: the maximum number of lists to return. This supports paging
#' @param verbose logical: show additional progress information? 
#'
#' @return data.frame of list name and other details
#'
#' @seealso \code{\link{sbdi_list}}
#'
#' @examples
#' \dontrun{
#'  ## lists that include the giant African snail Achatina fulica
#'  ##  (which is a notifiable pest species in some states)
#'  l <- sbdi_lists(search_guids("Achatina fulica")$guid)
#' }
#'
#' @export
sbdi_lists <- function(guid,offset=0,max=500,verbose=sbdi_config()$verbose) {
  
 ALA4R::ala_lists(guid,offset,max,verbose) 
  
}