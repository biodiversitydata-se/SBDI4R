#' Returns the name of the cache file associated with the given URL. 
#' 
#' Note that this file may not actually exist, this function just provides the mapping from URL to file name
#' 
#' @references \url{https://api.biodiversitydata.se/}
#' @seealso \code{sbdi_config} for cache settings, particularly the cache directory
#'  
#' @param url string: the URL e.g. \url{https://spatial.biodiversitydata.se/}
#' @return string: the file path and name
#' 
#' @examples
#' \dontrun{
#' sbdi_cache_filename("https://spatial.biodiversitydata.se/")
#' }
#' 
#' @export sbdi_cache_filename
sbdi_cache_filename <- function(url){
  assert_that(is.string(url))
  
  ALA4R::ala_cache_filename(url)
} 
