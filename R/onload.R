.onLoad <- function(libname,pkgname) {
  
  if (pkgname == "SBDI4R") {
  
    ## get default server config from ALA4R package  
    temp <- getOption("ALA4R_server_config")
      
    ## We need to overwrite the server configuration found in ALA4R package with SBDI info and urls
    sbdi_config() 
    version_string <- "version unknown"
    suppressWarnings(try(version_string<-utils::packageDescription('SBDI4R')[["Version"]],silent=TRUE)) ## get the SBDI4R version, if we can
    user_agent_string <- paste0("SBDI4R ",version_string)
    sbdi_config(user_agent=user_agent_string) 
    ## Both APIs are the same (SBDI has recently based their API on ALA)  
    ## Therefore this package is simply a wrapper around ALA4R functions and updates in ALA4R can be
    ## incorporated in SBDI4R
  
    temp$brand <- "SBDI4R" ## the package name that is shown to users in messages and warnings
    temp$support_email <- "support@ala.org.au" ### IS THIS CORRECT?
    temp$max_occurrence_records = 500000
    temp$server_max_url_length = 8150 ## bytes, for Apache with default LimitRequestLine value of 8190, allowing 40 bytes wiggle room. Users will be warned of possible problems when URL exceeds this length
    temp$notify <- "If this problem persists please notify the SBDI4R maintainers by lodging an issue at SBDI4R github repo or emailing 'email@domain.com'" ## the string that will be displayed to users to notify the package maintainers
    temp$reasons_function = "sbdi_reasons" ## the ala_reasons or equivalent function name
    temp$fields_function = "sbdi_fields" ## the sbdi_fields or equivalent function name
    temp$occurrences_function = "occurrences" ## the occurrences or equivalent function name
    temp$config_function = "sbdi_config" ## the ala_config or equivalent function name
    temp$biocache_version = "2.2.3" 
    temp$base_url_spatial = "https://spatial.bioatlas.se/ws/" ## the base url for spatial web services
    temp$base_url_bie = "https://species.bioatlas.se/ws/" ## the base url for BIE web services
    temp$base_url_biocache = "https://records.bioatlas.se/ws/" ## Services for mapping occurrence data, and species breakdowns for geographic areas.
    temp$base_url_biocache_download = "https://records.bioatlas.se/ws/biocache-download/"
    temp$base_url_alaspatial = "https://spatial.bioatlas.se/ws/" ## the base url for older ALA spatial services
    temp$base_url_images = "https://images.bioatlas.se/" ## the base url for the images database. Set to NULL or empty string if not available
    temp$base_url_logger = "https://logger.bioatlas.se/service/logger/" ## the base url for usage logging webservices
    temp$base_url_lists = "https://lists.bioatlas.se/ws/" ## base url for services for creating & editing lists of taxa
    temp$base_url_collections = "https://collections.bioatlas.se/ws/" ## ADDED BY SBDI base url for listing dataresources and Institutions    
    
    ## override any other settings here
    options(ALA4R_server_config = temp)
    

    ##### OTHER WAY
    # if (!"ALA4R_server_config" %in% names(options())) {
    # message("\nNo existing ALA4R server config, using Swedish data sources...\n")
    # options(ALA4R_server_config = server_config)
    # } else {
    #   message("Overwriting existing ALA server config with new...")
    #   options(ALA4R_server_config = server_config)
    # }
  
  }
}
