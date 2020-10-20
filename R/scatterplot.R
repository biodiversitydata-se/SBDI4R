#' Scatterplot of environmental parameters associated with species occurrence data. 
#' Interactive plotly plot of environmental variables and species occurrence data
#' 
#' WARNING: This function is under development and its functionality is limited
#' Makes a search using \code{\link{occurrences}}and plots environmental variables 
#' associated with species presence against each other. Alternatively, makes a 
#' search using \code{\link{occurrences}}and plots species occurrence data e.g 
#' frequencies against environmental variables. Note that there is a limit of 
#' 500000 records per request when using \code{method="indexed"}. Use the 
#' \code{method="offline"} for larger requests. For small requests, 
#' \code{method="indexed"} likely to be faster. 
#' 
#' @references \itemize{
#' \item Associated Bioatlas web service for record counts: \url{https://api.bioatlas.se/#ws3} # TODO UPDATE
#' \item Associated Bioatlas web service for occurrence downloads: \url{https://api.bioatlas.se/#ws4} # TODO UPDATE
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
#' @param \dots : other options passed to occurrences()
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
#' @export
## TODO: more extensive testing. In particular making the fields = "all" call work. Need to work on design and layout of plotly plot, current design/layout is default.
## TODO: miss variable containing date to plot on X-axis, current fields only contain year + day?
## TODO FUTURE:changing axis layout to logarithmic scale would be useful 
## TODO FUTURE: adding diversity metrics functionality e.g. calculating species richness and plot against environmental variables.

scatterplot <- function(taxon, ...) {
  
  df <- occurrences(taxon, ...)

  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = TRUE,
    showgrid = TRUE
  )
  
  plot_ly(df$data, x = ~year, y = ~individualCount, type = "scatter", mode = "markers", visible = T) %>%
    plotly::layout(
      title = paste("SBDI4R - Scatterplot",taxon, sep=" = "),
      xaxis = ax,
      yaxis = ax,
      updatemenus = list(
        ## Y-AXIS
        list(
          y = 0.7,
          buttons = list(
            list(method = "restyle",
                 args = list("y", list(df$data$individualCount)),  
                 label = "Y = INDIVIDUAL_COUNT"),
            list(method = "restyle",
                 args = list("y", list(df$data$locality)),  
                 label = "Y = LOCALITY"),
            list(method = "restyle",
                 args = list("y", list(df$data$minimumDepthInMeters)),  
                 label = "Y = MINIMUM_DEPTH (m)"),
            list(method = "restyle",
                 args = list("y", list(df$data$maximumDepthInMeters)),  
                 label = "Y = MAXIMUM_DEPTH (m)"),
            list(method = "restyle",
                 args = list("y", list(df$data$latitude)),  
                 label = "Y = LATITUDE"),
            list(method = "restyle",
                 args = list("y", list(df$data$longitude)),  
                 label = "Y = LONGITUDE"))),
        ## X-AXIS
        list(
          x = 0.7,
          buttons = list(
            list(method = "restyle",
                 args = list("x", list(df$data$year)),  
                 label = "X = YEAR"),
            list(method = "restyle",
                 args = list("x", list(df$data$locality)),  
                 label = "X = LOCALITY"),
            list(method = "restyle",
                 args = list("x", list(df$data$minimumDepthInMeters)),  
                 label = "X = MINIMUM_DEPTH (m)"),
            list(method = "restyle",
                 args = list("x", list(df$data$maximumDepthInMeters)),  
                 label = "X = MAXIMUM_DEPTH (m)"),
            list(method = "restyle",
                 args = list("x", list(df$data$longitude)),  
                 label = "X = LONGITUDE"),
            list(method = "restyle",
                 args = list("x", list(df$data$latitude)),  
                 label = "X = LATITUDE")))
      ))
}
