#' HTTP GET with caching
#'
#' Convenience wrapper for web GET operations. Caching, setting the user-agent string, and basic checking of the result are handled.
#'
#' @param url string: the url of the page to retrieve
#' @param type string: the expected content type. Either "text" (default), "json", or "filename" (this caches the content directly to a file and returns the filename without attempting to read it in)
#' @param caching string: caching behaviour, by default from sbdi_config()$caching
#' @param verbose string: verbose behaviour, by default from sbdi_config()$verbose
#' @param on_redirect function: passed to check_status_code() 
#' @param on_server_error function: passed to check_status_code()
#' @param on_client_error function: passed to check_status_code()
#' @param encoding encoding
#' @return for type=="text" the content is returned as text. For type=="json", the content is parsed using jsonlite::fromJSON. For "filename", the name of the stored file is returned.
#' @details Depending on the value of caching, the page is either retrieved from the cache or from the url, and stored in the cache if appropriate. The user-agent string is set according to sbdi_config()$user_agent. The returned response (if not from cached file) is also passed to check_status_code().
#' @references \url{https://api.biodiversitydata.se/}
cached_get <- function(url,
                       type="text",
                       caching=sbdi_config()$caching,
                       verbose=sbdi_config()$verbose,
                       on_redirect=NULL,
                       on_client_error=NULL,
                       on_server_error=NULL,
                       encoding=sbdi_config()$text_encoding) {
  
  ALA4R:::cached_get(url,
                     type,
                     caching,
                     verbose,
                     on_redirect,
                     on_client_error,
                     on_server_error,
                     encoding)
  
}