# Check HTTP status code
# 
# Generic check function that checks HTTP status codes coming back from nbn requests.
# 
# @param x string: a status code, or an object of class "response" (from e.g. httr's GET)
# @param on_redirect function: optional function to evaluate in the case of a redirect (3xx) code. By default a warning is issued.
# @param on_client_error function: optional function to evaluate in the case of a client error (4xx) code. By default an error is thrown.
# @param on_server_error function: optional function to evaluate in the case of a server error (5xx) code. By default an error is thrown.
# @param extra_info string: additional diagnostic info that will be shown to the user for 4xx or 5xx codes, where x is not a full response object
# @return integer: simplified status code (0=success (2xx codes), 1=warning (3xx codes))
# @references \url{http://www.w3.org/Protocols/HTTP/HTRESP.html}
# @examples
# \dontrun{
# require(httr)
# out <- GET(url="https://bioatlas.se/")
# check_status_code(out) ## pass the whole response object
# check_status_code(out$headers$status) ## or pass the status code explicitly
# }
check_status_code <- function(x,
                              on_redirect=NULL,
                              on_client_error=NULL,
                              on_server_error=NULL,
                              extra_info="") {
  
  ALA4R:::check_status_code(x,on_redirect,
                            on_client_error,
                            on_server_error,
                            extra_info)
  
}
