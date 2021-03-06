% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbdi_config.R
\name{sbdi_config}
\alias{sbdi_config}
\alias{ala_config}
\alias{sbdi_reasons}
\title{Get or set configuration options that control SBDI4R behaviour}
\usage{
sbdi_config(...)

ala_config(...)

sbdi_reasons()
}
\arguments{
\item{\dots}{Options can be defined using name=value. Valid options are:
\itemize{
  \item reset: \code{sbdi_config("reset")} will reset the options to their default values
  \item caching string: caching can be
    "on" (results will be cached, and any cached results will be re-used),
    "refresh" (cached results will be refreshed and the new results stored in the cache), or
    "off" (no caching, default).
  \item cache_directory string: the directory to use for the cache.
    By default this is a temporary directory, which means that results will only be cached
    within an R session and cleared automatically when the user exits R. The user may wish to set this to a non-temporary directory for
    caching across sessions. The directory must exist on the file system.
  \item verbose logical: should sbdi4R give verbose output to assist debugging?  (default=FALSE)
  \item warn_on_empty logical: should a warning be issued if a request returns an empty result set? (default=FALSE)
  \item user_agent string: the user-agent string used with all web requests to the sbdi servers.
    Default = "sbdi4R" with version number
  \item text_encoding string: text encoding assumed when reading cached files from local disk (default="UTF-8")
  \item download_reason_id numeric or string: the "download reason" required by some sbdi services, either as a numeric ID (currently 0--11)
  or a string (see \code{sbdi_reasons()} for a list of valid ID codes and names). By default this is NA. Some sbdi services require a valid
  download_reason_id code, either specified here or directly to the associated R function.
}}
}
\value{
For sbdi_config(), a list of all options. When sbdi_config(...) is
called with arguments, nothing is returned but the configuration is set.
}
\description{
Get or set configuration options that control SBDI4R behaviour

List valid download reasons
}
\examples{
\dontrun{
 sbdi_config()
 sbdi_config(caching="off")
 sbdi_config()
 sbdi_config(download_reason_id=0,verbose=TRUE)
 sbdi_config("reset")
}
}
\references{
\url{https://api.biodiversitydata.se/}

\url{https://spatial.biodiversitydata.se/layers/} this will eventually move to the api link

Invoking \code{sbdi_config()} with no arguments returns a list with the current values of the options.

\code{sbdi_reasons()} returns a data frame with information describing the valid options for \code{download_reason_id}
}
