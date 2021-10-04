#' Check HTTP status code
#' 
#' Generic check function that checks HTTP status codes coming back from nbn requests.
#' 
#' @param x string: a status code, or an object of class "response" (from e.g. httr's GET)
#' @param on_redirect function: optional function to evaluate in the case of a redirect (3xx) code. By default a warning is issued.
#' @param on_client_error function: optional function to evaluate in the case of a client error (4xx) code. By default an error is thrown.
#' @param on_server_error function: optional function to evaluate in the case of a server error (5xx) code. By default an error is thrown.
#' @param extra_info string: additional diagnostic info that will be shown to the user for 4xx or 5xx codes, where x is not a full response object
#' @return integer: simplified status code (0=success (2xx codes), 1=warning (3xx codes))
#' @references \url{http://www.w3.org/Protocols/HTTP/HTRESP.html}
#' @examples
#' \dontrun{
#' require(httr)
#' out <- GET(url="https://biodiversitydata.se/")
#' check_status_code(out) ## pass the whole response object
#' check_status_code(out$headers$status) ## or pass the status code explicitly
#' }
check_status_code <- function(x,
                              on_redirect=NULL,
                              on_client_error=NULL,
                              on_server_error=NULL,
                              extra_info="") {
  
  # ALA4R:::check_status_code(x,on_redirect,
  #                           on_client_error,
  #                           on_server_error,
  #                           extra_info)
  
  # function (x, on_redirect = NULL, on_client_error = NULL, on_server_error = NULL, 
  #           extra_info = "") 
  # {
    assert_that(is.string(extra_info))
    was_full_response <- FALSE
    if (inherits(x, "response")) {
      was_full_response <- TRUE
      xstatus <- x$headers$status
      if (is.null(xstatus)) {
        xstatus <- as.character(x$status_code)
        was_full_response <- TRUE
      }
      if (is.null(xstatus)) {
        warning("error in http status checking: skipped. ", 
                getOption("ALA4R_server_config")$notify)
        was_full_response <- FALSE
        xstatus <- "200"
      }
    }
    else {
      if (!see_if(is.string(x))) {
        if (!see_if(is.count(x))) {
          stop("expecting either http response object, or status code as\n string or numeric")
        }
        else {
          x <- as.character(x)
        }
      }
      xstatus <- x
    }
    switch(substr(xstatus, 1, 1), `2` = {
      return(0)
    }, `3` = {
      if (!is.null(on_redirect)) {
        assert_that(is.function(on_redirect))
        return(on_redirect(xstatus))
      } else {
        warning("HTTP status code ", xstatus, " received.\nThis may be OK: if there are problems,\n  please notify the package maintainers.")
        return(1)
      }
    }, `4` = {
      if (!is.null(on_client_error)) {
        assert_that(is.function(on_client_error))
        return(on_client_error(xstatus))
      } else {
        diag_msg <- paste0("  Either there was an error with your\n request or in the ", 
                           getOption("ALA4R_server_config")$brand, 
                           " package, or the servers are down. ", 
                           getOption("ALA4R_server_config")$notify)
        if (was_full_response) {
          x <- jsonlite::fromJSON(content(x, type = "text", 
                                          encoding = "UTF-8"))
          if (!is.null(x$message)) {
            diag_msg <- paste(diag_msg, "\nThe error message was:", 
                              x$message, sep = " ")
          }
        } else {
          if (nchar(extra_info) > 0) {
            diag_msg <- paste(diag_msg, "\n  Some additional\n diagnostic information that might\n help:", 
                              extra_info, sep = " ")
          }
        }
        stop("HTTP status code ", xstatus, " received.\n", diag_msg)
      }
    }, `5` = {
      if (!is.null(on_server_error)) {
        assert_that(is.function(on_server_error))
        return(on_server_error(xstatus))
      } else {
        diag_msg <- paste0("  Either there was an error with the request, or the servers may be down\n (try again later). ", 
                           getOption("ALA4R_server_config")$notify)
        if (was_full_response) {
          x <- jsonlite::fromJSON(content(x, type = "text", 
                                          encoding = "UTF-8"))
          if (!is.null(x$message)) {
            diag_msg <- paste(diag_msg, "\nThe error message was:", 
                              x$message, sep = " ")
          }
        } else {
          if (nchar(extra_info) > 0) {
            diag_msg <- paste(diag_msg, "\n  Some additional\n diagnostic information that might\n help:", 
                              extra_info, sep = " ")
          }
        }
        stop("HTTP status code ", xstatus, " received.\n", 
             diag_msg)
      }
    })
    warning("Unexpected HTTP status code ", x, " received.\n  ", 
            getOption("ALA4R_server_config")$notify)
}
