# HTTP POST with caching
#
# Convenience wrapper for web POST operations. Caching, setting the user-agent string, and basic checking of the result are handled.
#
# @param url string: the url of the page to retrieve
# @param body string: the body to POST
# @param type string: the expected content type. Either "text" (default), "json", or "filename" (this caches the content directly to a file and returns the filename without attempting to read it in)
# @param caching string: caching behaviour, by default from sbdi_config()$caching
# @param content_type string: set the Content-Type header to a specific value (needed for e.g. search_names), default is unset
# @param ... additional arguments passed to curlPerform
# @return for type=="text" the content is returned as text. For type=="json", the content is parsed using jsonlite::fromJSON. For "filename", the name of the stored file is returned.
# @details Depending on the value of caching, the page is either retrieved from the cache or from the url, and stored in the cache if appropriate. The user-agent string is set according to sbdi_config()$user_agent. The returned response (if not from cached file) is also passed to check_status_code().
# @references \url{https://api.bioatlas.se/}
# @examples
#
# out = cached_post(url="https://species.bioatlas.se/ws/species/lookup/bulk",body=jsonlite::toJSON(list(names=c("Vulpes vulpes","Lutra"))),type="json")
cached_post <- function(url,body,type="text",caching=sbdi_config()$caching,verbose=sbdi_config()$verbose,content_type,encoding=sbdi_config()$text_encoding,...) {
  ALA4R:::cached_post(url,body,type,caching,verbose,content_type,encoding,...)
}
  