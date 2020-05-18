#' Get or set configuration options that control sbdi4R behaviour
#'
#' @references \url{https://api.bioatlas.se/}
#' @references \url{https://spatial.bioatlas.se/layers/} this will eventually move to the api link
#'
#' Invoking \code{sbdi_config()} with no arguments returns a list with the current values of the options.
#'
#' \code{sbdi_reasons()} returns a data frame with information describing the valid options for \code{download_reason_id}
#'
#' @param \dots Options can be defined using name=value. Valid options are:
#' \itemize{
#'   \item reset: \code{sbdi_config("reset")} will reset the options to their default values
#'   \item caching string: caching can be
#'     "on" (results will be cached, and any cached results will be re-used),
#'     "refresh" (cached results will be refreshed and the new results stored in the cache), or
#'     "off" (no caching, default).
#'   \item cache_directory string: the directory to use for the cache.
#'     By default this is a temporary directory, which means that results will only be cached
#'     within an R session and cleared automatically when the user exits R. The user may wish to set this to a non-temporary directory for
#'     caching across sessions. The directory must exist on the file system.
#'   \item verbose logical: should sbdi4R give verbose output to assist debugging?  (default=FALSE)
#'   \item warn_on_empty logical: should a warning be issued if a request returns an empty result set? (default=FALSE)
#'   \item user_agent string: the user-agent string used with all web requests to the sbdi servers.
#'     Default = "sbdi4R" with version number
#'   \item text_encoding string: text encoding assumed when reading cached files from local disk (default="UTF-8")
#'   \item download_reason_id numeric or string: the "download reason" required by some sbdi services, either as a numeric ID (currently 0--11)
#'   or a string (see \code{sbdi_reasons()} for a list of valid ID codes and names). By default this is NA. Some sbdi services require a valid
#'   download_reason_id code, either specified here or directly to the associated R function.
#' }
#'
#' @return For sbdi_config(), a list of all options. When sbdi_config(...) is
#' called with arguments, nothing is returned but the configuration is set.
#'
#' @examples
#' \dontrun{
#'  sbdi_config()
#'  sbdi_config(caching="off")
#'  sbdi_config()
#'  sbdi_config(download_reason_id=0,verbose=TRUE)
#'  sbdi_config("reset")
#' }
#' @export sbdi_config
sbdi_config <- function(...) {
  # ALA4R::ala_config(...)
  ## get or set options that control ALA4R behaviour
  ## options are stored as a global option with the name defined in ala_option_name
  ala_option_name <- "ALA4R_config"
  user_options <- list(...) ## the options passed by the user
  
  ## default user-agent string
  version_string <- "version unknown"
  suppressWarnings(try(version_string<-utils::packageDescription('SBDI4R')[["Version"]],silent=TRUE)) ## get the SBDI4R version, if we can
  user_agent_string <- paste0("SBDI4R ",version_string)
  
  ## set default options
  default_options <- list(
    caching="on",
    cache_directory=tempdir(),
    user_agent=user_agent_string,
    download_reason_id=NA,
    verbose=FALSE,
    warn_on_empty=FALSE,
    text_encoding="UTF-8",
    email=NULL
  )
  
  ## define allowed options, for those that have restricted values
  allowed_options <- list(caching=c("on","off","refresh"),download_reason_id=c(0:8,10:12))
  ## ideally, the valid download_reason_id values should be populated dynamically from the ala_reasons() function. However if that is called (from here) before the AL4R_config option has been set, then we get infinite recursion. To be addressed later ...
  
  ## has the user asked to reset options to defaults?
  if (identical(tolower(user_options),"reset")) {
    temp <- list(default_options)
    names(temp) <- ala_option_name
    options(temp)
  } else {
    names(user_options) <- tolower(names(user_options))
    ## has the user specified something we don't recognize?
    known_options <- names(default_options)
    unknown <- setdiff(names(user_options),known_options)
    if (length(unknown)>0) {
      stop("unknown ",getOption("ALA4R_server_config")$brand," options: ", str_c(unknown,collapse=", "))
    }
    
    current_options <- getOption(ala_option_name)
    if (is.null(current_options)) {
      ## ALA4R options have not been set yet, so set them to the defaults
      current_options <- default_options
      ## set the global option
      temp <- list(current_options)
      names(temp) <- ala_option_name
      options(temp)
    }
    
    ## convert reason from char to numeric if needed
    if (!is.null(user_options$download_reason_id)) {
      user_options$download_reason_id <- convert_reason(user_options$download_reason_id)
    }
    
    ## override any defaults with user-specified options
    if (length(user_options)>0) {
      for (i in 1:length(user_options)) {
        this_option_name <- names(user_options)[i]
        if (! is.null(allowed_options[[this_option_name]])) {
          ## there are restrictions on the allowed values for this option
          ## could use match.arg here but the output is a bit obscure
          if (! (user_options[[i]] %in% allowed_options[[this_option_name]])) {
            stop("value \"",user_options[[i]],"\" is not a valid choice for ",this_option_name," (should be one of ", str_c(allowed_options[[this_option_name]],collapse=", "),")")
          }
        }
        ## any other specific checks ...
        if (identical(this_option_name,"cache_directory")) {
          if (!see_if(is.notempty.string(user_options[[i]]))) {
            stop("cache_directory should be a string")
          }
          ## strip trailing file separator, if there is one
          user_options[[i]] <- sub("[/\\]+$","",user_options[[i]])
          if (! (file.exists(user_options[[i]]) && file.info(user_options[[i]])$isdir)) {
            ## cache directory does not exist. We could create it, but this is probably better left to the user to manage
            stop("cache directory ",user_options[[i]]," does not exist");
          }
        }
        if (identical(this_option_name,"user_agent")) {
          if (!see_if(is.string(user_options[[i]]))) {
            stop("user_agent should be a string")
          }
        }
        if (identical(this_option_name,"verbose")) {
          if (!see_if(is.flag(user_options[[i]]))) {
            stop("verbose should be TRUE or FALSE")
          }
        }
        if (identical(this_option_name,"warn_on_empty")) {
          if (!see_if(is.flag(user_options[[i]]))) {
            stop("warn_on_empty should be TRUE or FALSE")
          }
        }
        if (identical(this_option_name,"email")) {
          if (!see_if(is.string(user_options[[i]]))) {
            stop("email should be a string")
          }
        }
        
        current_options[this_option_name] <- user_options[[i]]
      }
      ## set the global option
      temp <- list(current_options)
      names(temp) <- ala_option_name
      options(temp)
    } else {
      ## no user options were provided, so user is asking for current options to be returned
      current_options
    }
  }
}

#' @rdname sbdi_config
#' @export ala_config
ala_config <- function(...) {
  sbdi_config( ... )
}


#' @rdname sbdi_config
#' @export
sbdi_reasons <- function() {
  ALA4R::ala_reasons()
}

## internal function, used to define the ALA4R sourceTypeId parameter value, passed by occurrences download and possibly other functions
sbdi_sourcetypeid <- function() {
  ALA4R:::ala_sourcetypeid()
}


convert_reason <- function(reason) {
  ## unexported function to convert string reason to numeric id
 ALA4R:::convert_reason(reason)
}
