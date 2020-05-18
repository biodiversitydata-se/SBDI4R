#' Returns the name of the cache file associated with the given URL. Note that this file 
#' may not actually exist, this function just provides the mapping from URL to filename
#' 
#' @references \url{https://api.bioatlas.se/}
#' @seealso \code{sbdi_config} for cache settings, particularly the cache directory
#'  
#' @param url string: the URL
#' @return string: the file path and name
#' 
#' @examples
#' sbdi_cache_filename("https://spatial.bioatlas.se/")
#' 
#' @export sbdi_cache_filename
sbdi_cache_filename <- function(url){
  assert_that(is.string(url))
  
  ALA4R::ala_cache_filename(url)
} 
