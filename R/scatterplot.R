#' Scatterplot of environmental parameters associated with species occurrence data. 
#' Interactive plotly plot of environmental variables and species occurence data
#' 
#' WARNING: This function is under development and its functionality is limited
#' Makes a search using \code{\link{occurrences}}and plots environmental variables 
#' associated with species presence against each other. Alternatively, makes a 
#' search using \code{\link{occurrences}}and plots species occurence data e.g 
#' frequencies against environmental variables. Note that there is a limit of 
#' 500000 records per request when using \code{method="indexed"}. Use the 
#' \code{method="offline"} for larger requests. For small requests, 
#' \code{method="indexed"} likely to be faster. 
#' 
#' @references \itemize{
#' \item Associated Bioatlas web service for record counts: \url{https://api.bioatlas.se/#ws3} # TODO UPDATE
#' \item Associated Bioatlas web service for occurence downloads: \url{https://api.bioatlas.se/#ws4} # TODO UPDATE
#' \item Field definitions: \url{https://docs.google.com/spreadsheet/ccc?key=0AjNtzhUIIHeNdHhtcFVSM09qZ3c3N3ItUnBBc09TbHc}
#' \item WKT reference: \url{http://www.geoapi.org/3.0/javadoc/org/opengis/referencing/doc-files/WKT.html}
#' }
#' 
#' @param taxon string: (optional) query of the form field:value (e.g. "genus:Macropus") 
#' or a free text search (e.g. "macropodidae"). Note that a free-text search is 
#' equivalent to specifying the "text" field (i.e. \code{taxon="Alaba"} is 
#' equivalent to \code{taxon="text:Alaba"}. 
#' The text field is populated with the taxon name along with a handful of other 
#' commonly-used fields, and so just specifying your target taxon (e.g. taxon="Alaba vibex") 
#' will probably work. However, for reliable results it is recommended to use a 
#' specific field where possible (see \code{nbn_fields("occurrence_indexed")} for 
#' valid fields). It is also good practice to quote the taxon name if it contains 
#' multiple words, for example \code{taxon="taxon_name:\"Alaba vibex\""}
#' @param \dots : other options passed to occurences()
#' @return Data frame of occurrence results, with one row per occurrence record. 
#' The columns of the dataframe will depend on the requested fields. The data 
#' frame is plotted with ggplot and output stored as pdf, an interactive ggplot
#' using plotly is displayed on-screen.
#' @seealso \code{\link{sbdi_reasons}} for download reasons; \code{\link{sbdi_config}}
#' @importFrom graphics layout
#' @importFrom plotly plot_ly add_markers
#' @importFrom magrittr %>% 
#' @examples
#' \dontrun{ 
#' scatterplot(taxon="Ectocarpus siliculosus", download_reason_id=10)
#' }
#' @export scatterplot
## TODO: more extensive testing. In particular making the fields = "all" call work. Need to work on design and layout of plotly plot, current design/layout is default.
## TODO FUTURE: adding diversity metrics functionality e.g. calculating species richness and plot against environmental variables.

scatterplot <- function(taxon, ...) {
  
  df <- occurrences(taxon, ...)

  #Plotly - First test with lat vs. long
  plot_ly(df$data, x = ~longitude) %>%
    add_markers(y = ~latitude, name = "lat") %>%
    add_markers(y = ~longitude, name = "long", visible = FALSE) %>%
    graphics::layout(
      title = paste(taxon, "SBDI4R Scatterplot", sep=" "),
      xaxis = list(domain = c(0.1, 1)),
      yaxis = list(title = "lat"),
      showlegend = FALSE,
      updatemenus = list(
        list(
          y = 0.7,
          buttons = list(
            list(method = "update",
                 args = list(list(visible = list(TRUE, FALSE)),
                             list(yaxis = list(title = "lat"))
                             ),
                 label = "lat"),
            
            list(method = "update",
                 args = list(list(visible =  list(FALSE, TRUE)),
                             list(yaxis = list(title = "long"))
                             ),
                 label = "long")
            )
        )
      )
    )
}
