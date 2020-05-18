#' Search for environmental and contextual data layers
#' 
#' @references Associated SBDI web services: \url{https://api.bioatlas.se/#ws11} \url{https://api.bioatlas.se/#ws12} \url{https://api.bioatlas.se/#ws13}
#' @references Descriptions of the spatial layers: \url{https://spatial.bioatlas.se/ws/layers})
#'
#' @param query text string: optional search term against layer metadata. Only layers that include this term in their metadata will be returned.
#' @param type string: either "all" (all possible layers; default), "grids" (gridded environmental layers), or "shapes" (contextual shapefile layers)
#' @param output_format string: controls the print method for the returned object. Either "complete" (the complete data structure is displayed), or "simple" (a simplified version is displayed). Note that the complete data structure exists in both cases: this option only controls what is displayed when the object is printed to the console. The default output format is "simple"
#' @return A data frame of results. The contents (column names) of the data frame will vary depending on the details of the search and the results
#' 
#' @examples
#' \dontrun{
#' search_layers(type="all")
#' search_layers(type="grids",query="park")
#' search_layers(type="shapes",query="park",output_format="simple")
#' }
#' @export
search_layers <- function(query, type = "all", output_format = "simple") {

ALA4R::search_layers(query, type, output_format)

}