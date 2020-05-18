#' Get list of taxa and their occurrence counts
#'
#' Retrieve a list of taxa matching a search query, within a spatial search area, or both.
#'
#' 
#' @references Associated SBDI web service: \url{https://api.bioatlas.se/#ws78}
#' @references \url{http://www.geoapi.org/3.0/javadoc/org/opengis/referencing/doc-files/WKT.html}
#' 
#' @param taxon string: (optional) query of the form field:value (e.g. "genus:Macropus") or a free text search (e.g. "macropodidae").
#' For reliable results it is recommended to use a specific field where possible (see \code{sbdi_fields("occurrence_indexed")}
#' for valid fields). It is also good practice to quote the taxon name if it contains multiple words, for example
#' \code{taxon="taxon_name:\"Vulpes vulpes\""} (noting, however, that multi-word names are unlikely in the context of a specieslist
#' search, where one would typically be searching for all species within, say, a genus or family)
#' @param wkt string: WKT (well-known text) defining a polygon within which to limit taxon search, e.g. "POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))"
#' @param fq string: character string or vector of strings, specifying filters to be applied to the original query. These are of the form "INDEXEDFIELD:VALUE" e.g. "kingdom:Fungi". See \code{sbdi_fields("occurrence_indexed",as_is=TRUE)} for all the fields that are queryable. NOTE that fq matches are case-sensitive, but sometimes the entries in the fields are not consistent in terms of case (e.g. kingdom names "Fungi" and "Plantae" but "ANIMALIA"). fq matches are ANDed by default (e.g. c("field1:abc","field2:def") will match records that have field1 value "abc" and field2 value "def"). To obtain OR behaviour, use the form c("field1:abc OR field2:def")
#' @return data frame of results, where each row is a taxon, its classification information, and its occurrence count
#' @seealso \code{\link{sbdi_fields}} for occurrence fields that are queryable via the \code{fq} parameter
#' @examples
#' \dontrun{
#' x <- specieslist(taxon="genus:Leuctra",wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))")
#' 
#' x <- specieslist(wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))",fq="rank:species")
#'
#' x <- specieslist(wkt="POLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))",fq="genus:Macropus")
#'
#' x <- specieslist(wkt="POLYGONPOLYGON((-3 56,-4 56,-4 57,-3 57,-3 56))",
#' fq="kingdom:Plantae")
#' ## NOTE that this response might include records with empty or NA kingdom, phylum, or
#' ##  class values, as per the note above.
#' }
#' @export specieslist
specieslist <- function(taxon, wkt, fq) {
  
  ALA4R::specieslist(taxon, wkt, fq)
  
}