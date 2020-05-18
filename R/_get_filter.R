



#A messy function to get a vector of filter queries to use in functions asking for fq=. 
#TODO integrate with specieslists, search_layers, and create search_resources
get_filter <- function(type = NULL) {
  
  continue <- TRUE
  res <- c()
  
  while(continue){
    
    if(is.null(type)){
      type <- readline('What type do you want to filter by? Choose: "resource", "specieslist" or "layer". ')
    }
    
    type <- tolower(type)
    
    if(type=="resource"){
      
      # this_url <- SBDI4R:::build_url_from_parts(getOption("ALA4R_server_config")$base_url_collections, 
      #                                                               c("dataResource"))
      
      this_url <- SBDI4R:::build_url_from_parts(getOption("ALA4R_server_config")$base_url_collections, 
                                               c("dataResource"))
      
      
      dataResources<-jsonlite::fromJSON(httr::content(httr::GET(this_url), "text"))
      dataResources <- dataResources[,c(3,1,2)]
      
      print("By what recource do you want to filter? Press enter to get list or answer by uid.")
      r <- readline()
      
      if(r==""){
        print(dataResources)
        print("By what recource do you want to filter? Answer by uid:")
        r <- readline()
      }
      
      if(r %in% dataResources$uid){
        ruid <- r
        uri <- dataResources$uri[dataResources$uid==ruid]
        
        objectResource <- jsonlite::fromJSON(httr::content(httr::GET(uri), "text"))
        objectResource <- objectResource$linkedRecordConsumers[,c(3,1,2)]
        print(objectResource)
        print("By what object do you want to filter? Answer by uid:")
        r <- readline()
        
        if(r %in% objectResource$uid){
          res <- c(res, paste0("institution_uid",":",r))
          
          if(continue("Filter added. Do you want to continue? Type y for yes. ")){
            type <- NULL
          }else{
            continue <- FALSE
          }
          
        }else{
          if(continue("No object with uid found. Do you want to start over? Type y for yes. ")){
            type <- NULL
          }else{
            continue <- FALSE
          }
        }
        
      }else{
        if(continue("No object with uid found. Do you want to start over? Type y for yes. ")){
          type <- NULL
        }else{
          continue <- FALSE
        }
      }
    }else{
      if(continue("Type not correct specified. Do you want to start over? Type y for yes. ")){
        type <- NULL
      }else{
        continue <- FALSE
      }
    }
  }
  
  return(res)
  
}


