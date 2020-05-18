#' Retrieve image information using image ids, with an option to download
#' 
#' @references \itemize{
#' \item Associated SBDI web service for images \url{https://images.bioatlas.se/ws}
#' }
#' 
#' @param id character: IDs of images to be downloaded as single string or
#' vector of strings
#' @param download logical: if TRUE download all images and add location to
#' dataframe
#' @param download_path string: (optional) filepath to download images to.
#' If not given and download param is TRUE, will create an images
#' folder
#' @param verbose logical: show additional progress information?
#' [default is set by sbdi_config()]
#' @return Data frame of image results
#' 
#' @examples 
#' \dontrun{
#' ## Retrieve infomation about an image and download
#' images(id="da5fe120-e213-4cd6-9c5f-62346ed2e466", download=TRUE)
#' }
#' @export images

images <- function(id, download=FALSE, download_path,
                   verbose=sbdi_config()$verbose) {
  
  ALA4R::images(id, download, download_path, verbose)
  
}


download_images <- function(data, media_dir, verbose=sbdi_config()$verbose, 
                            sounds = FALSE) {
  
  assert_that(!missing(media_dir))
  
  if(!file.exists(media_dir)) {
    message(sprintf('Media directory does not exist, creating directory %s',
                    media_dir))
    dir.create(media_dir)
  }
  
  ALA4R:::download_images(data, media_dir, verbose, sounds)
  
  
}
