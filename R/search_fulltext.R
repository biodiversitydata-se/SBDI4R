#' Full text search
#' 
#' Performs a search across all objects, and selects the closest matches. Generally, the user will provide the search term via the \code{query} parameter, with optional filtering via \code{fq}.
#'
#' @references Associated SBDI web service: \url{https://species.bioatlas.se/ws/}
#'  
#' @param query string: the search term
#' @param fq string: (optional) character string or vector of strings, specifying filters to be applied to the 
#' original query. These are of the form "INDEXEDFIELD:VALUE" e.g. "rk_kingdom:Fungi". 
#' See \code{sbdi_fields("general")} for all the fields that are queryable. 
#' NOTE that fq matches are case-sensitive, but sometimes the entries in the fields are 
#' not consistent in terms of case (e.g. kingdom names "Fungi" and "Plantae" but "ANIMALIA"). 
#' fq matches are ANDed by default (e.g. c("field1:abc","field2:def") will match records that have 
#' field1 value "abc" and field2 value "def"). To obtain OR behaviour, use the form c("field1:abc OR field2:def")
#' @param output_format string: controls the print method for the "data" component of the returned object. Either "complete" (the complete data structure is displayed), or "simple" (a simplified version is displayed). Note that the complete data structure exists in both cases: this option only controls what is displayed when the object is printed to the console. The default output format is "simple"
#' @param start numeric: (optional) (positive integer) start offset for the results
#' @param page_size numeric: (optional) (positive integer) maximum number of records to return. Defaults to the server-side value - currently 10
#' @param sort_by string: (optional) field to sort on. See \code{sbdi_fields("general")} for valid field names
#' @param sort_dir string: (optional) sort direction, either "asc" or "desc"
#' @seealso \code{\link{sbdi_fields}}
#' @return a named list with the components "meta" (search metadata), "facets" (search facet results), and "data" (the search results, in the form of a data frame). The contents (column names) of the data frame will vary depending on the details of the search and the results, but should contain at least the columns \code{guid}, \code{name}, \code{commonName}, \code{rank}, \code{author}, and \code{occurrenceCount}.
#' 
#' @examples
#' \dontrun{
#'  # find information SBDI holds on red kangaroo
#'  search_fulltext("Rabbit")
#'  search_fulltext("Macropus rufus")
#'  search_fulltext("NHMSYS0000080188")
#'
#'  # find genus names like "Oenanthe"
#'  search_fulltext("oenanthe",sort_by="rk_kingdom",fq="rank:genus")
#' }
#' @export search_fulltext
search_fulltext <- function(query, fq="idxtype:TAXON", output_format = "simple", 
                            start, page_size, sort_by, sort_dir) {
  
 ALA4R::search_fulltext(query, fq, output_format, 
                        start, page_size, sort_by, sort_dir)
  
}
