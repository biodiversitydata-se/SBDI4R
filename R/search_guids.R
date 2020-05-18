#' Lookup of species identifiers
#' 
#' Provides names, taxonomic classification, and other information for a list of GUIDs.
#'
#' @references The associated SBDI web service: \url{https://api.bioatlas.se/#ws110}
#' 
#' @param guids string: a single GUID or vector of GUIDs
#' @param occurrence_count logical: if TRUE then also return the number of occurrences of each matched GUID.
#' Note that this requires one extra web call for each GUID, and so may be slow.
#' @param output_format string: controls the print method for the returned object. Either "complete" (the complete data structure is displayed), or "simple" (a simplified version is displayed). Note that the complete data structure exists in both cases: this option only controls what is displayed when the object is printed to the console. The default output format is "simple"
#' @return A data frame, which should include one entry (i.e. one data.frame row or one list element) per input GUID. The columns in the data.frame output may vary depending on the results returned by the SBDI server, but should include searchTerm, name, rank, and guid.
#' @examples
#' \dontrun{
#' search_guids(c("5219963",
#' "5219961","this_is_not_a_valid_guid"))
#' }
#' @export
search_guids <- function(guids = c(), occurrence_count = FALSE, output_format = "simple") {
  
  ALA4R::search_guids(guids, occurrence_count, output_format)
  
}