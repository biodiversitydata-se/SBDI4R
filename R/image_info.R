#' Fetch information about an image, given its image ID
#'
#' Note that there is currently no web service that provides image information, 
#' and so we are scraping results from pages of the form https://images.bioatlas.se/image/details?imageId=id. 
#' This web scraping may be fragile, and will be replaced by a web-service-based function when one becomes available.
#' 
#' @param id character: IDs of images (e.g. as returned by \code{\link{occurrences}} 
#'  in the imageUrl column). Each ID will be of a format something like "84654e14-dc35-4486-9e7c-40eb2f8d3faa"
#' @param verbose logical: show additional progress information? [default is set by sbdi_config()]
#' @return A data.frame with one row per \code{id}, and at least the columns 
#' imageIdentifier and imageURL
#' @seealso \code{\link{sbdi_config}} \code{\link{occurrences}}
#' @examples
#' \dontrun{
#' ## Using IDs returned from occurrences() function
#' 
#' image(id = c("f31e5f0e-f964-4bc3-b8f3-78f2ad520563",
#'   "39091d30-3016-417f-8860-b7c2a8d37942"))
#' }
#' @export image_info
image_info <- function(id, verbose = sbdi_config()$verbose) {
  
  if (as.character(match.call()[[1]]) == "image_info") {
    warning("image_info() has been renmaed to images() and now can also 
download images. Please use images() instead of image_info()", 
            call. = FALSE)
  }
  
  # suppressWarnings(ALA4R::image_info(id, verbose))
  suppressWarnings(ALA4R::images(id, verbose))
  
}

#' @export
#' @rdname images
images <- image_info

