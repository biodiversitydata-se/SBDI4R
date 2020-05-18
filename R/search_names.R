#' Lookup of taxonomic names
#' 
#' Provides GUID, taxonomic classification, and other information for a list of names. 
#' Case-insensitive but otherwise exact matches are used.
#'
#' @references The associated SBDI web service: \url{https://api.bioatlas.se/#ws87}
#' 
#' @param taxa string: a single name or vector of names
#' @param vernacular logical: if TRUE, match on common names as well as scientific names, otherwise match only on scientific names
#' @param guids_only logical: if TRUE, a named list of GUIDs will be returned. Otherwise, a data frame with more comprehensive information for each name will be returned.
#' @param occurrence_count logical: if TRUE (and if \code{guids_only} is FALSE) then also return the number of occurrences of each matched name.
#' Note that this requires one extra web call for each name, and so may be slow. Only applicable if \code{guids_only} is FALSE.
#' @param output_format string: controls the print method for the returned object (only has an effect when \code{guids_only} is FALSE). Either "complete" (the complete data structure is displayed), or "simple" (a simplified version is displayed). Note that the complete data structure exists in both cases: this option only controls what is displayed when the object is printed to the console. The default output format is "simple"
#' @return A data frame of results, or named list of GUIDs if \code{guids_only} is TRUE. The results should include one entry (i.e. one data.frame row or one list element) per input name. The columns in the data.frame output may vary depending on the results returned by the SBDI server, but should include searchTerm, name, rank, and guid.
#' 
#' @examples
#' \dontrun{
#' search_names(c("Diatoma tenuis","Diatoma hyemale var. quadratum",
#' "Macropus"))
#' 
#' search_names(c("Diatoma tenuis","Diatoma hyemale var. quadratum","Macropus"),
#' guids_only=TRUE)
#' 
#' search_names("Diatoma tenuis",vernacular=FALSE)
#' 
#' search_names("Diatoma tenuis",vernacular=TRUE)
#' 
#' ## occurrence counts for matched names
#' search_names(c("Diatoma tenuis","Diatoma hyemale var. quadratum",
#' "Macropus","Thisisnot aname"),occurrence_count=TRUE)
#'
#' ## no occurrence counts because guids_only is TRUE
#' search_names(c("Diatoma tenuis","Diatoma hyemale var. quadratum",
#' "Macropus","Thisisnot aname"),occurrence_count=TRUE,guids_only=TRUE)
#' }
#' @export search_names
search_names <- function(taxa = c(), vernacular = FALSE, guids_only = FALSE,
                         occurrence_count = FALSE, output_format = "simple") {
  
  ALA4R::search_names(taxa, vernacular, guids_only,
                      occurrence_count, output_format)
  
}