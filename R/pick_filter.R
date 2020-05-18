#' List and pick existing filters
#' 
#' Lists and lets the user interactivelly pick a filter from lists.  
#' 
#' @param type string: (optional) type of filter to create query string for. Options are \code{c("resource", "specieslist" , "layer")}. Species lists are waiting development from the API side.
#' @return a string ready to be places in the argument \sQuote{fq} in function \code{\link{occurrences}}. 
#' @seealso \code{\link{occurrences}} for download reasons; \code{\link{sbdi_config}}
#' @examples
#' \dontrun{
#' fq_str <- pick_filter("resource") #"data_resource_uid:dr5"
#' x <- occurrences(taxon="genus:Accipiter", 
#'                   download_reason_id=10,
#'                   fq=fq_str)
#' }
#' @importFrom httr GET content
#' @importFrom jsonlite fromJSON
#' @importFrom utils capture.output
#' @export pick_filter
pick_filter<-function(type = NULL){
  continue <- TRUE
  res <- c()
  
  while(continue){
    
    if(is.null(type)){
      type <- readline("What type do you want to filter by? Choose: 'resource', 'specieslist' or 'layer'. ")
    }
    
    type <- tolower(type)
    
    if(type == "resource"){
      res <- c(res, inst_questionarie())
      continue <- FALSE
    }else if(type == "specieslist"){
      message("not yet implemented")
      continue <- FALSE
    } else if(type == "layer"){  
      res <- layer_questionarie()
      continue <- FALSE
    }else{
      if(continue("'Type' should be one of the following: 'resource', 'specieslist' or 'layer'.\nDo you want to start over? Type 'y' for yes. ")){
        type <- NULL
      }else{
        continue <- FALSE
      }
    }
  }
  
  return(res)
}
 
#' auxiliary function for institutions
inst_questionarie<-function(){
  continue <- TRUE
  res <- c()
  this_url <- SBDI4R:::build_url_from_parts(getOption("ALA4R_server_config")$base_url_collections, 
                                            c("institution.json"))
  
  institutions <- fromJSON(content(GET(this_url), "text"))
  institutions <- institutions[,c("uid","name", "uri")]
  while(continue){
    message("\nWhich institution do you want to get data from? Type the corresponding uid: \n")
    message(paste0(capture.output(institutions[order(institutions[,"uid"]), c("uid","name")]), 
                   collapse = "\n"))
    r <- readline()
   
    if(r %in% institutions$uid){
      ruid <- r
      uri <- institutions$uri[institutions$uid==ruid]
      
      institutionAll <- jsonlite::fromJSON(httr::content(httr::GET(uri), "text"))
      collections <- institutionAll$collections
      dataResource <- institutionAll$linkedRecordProviders
      message("\nBy which collection or data resourcedo you want to filter? Type the corresponding uid or write 'all': \n")
      message("\nCollections: \n")
      message(paste0(capture.output(collections[order(collections$uid),c("uid","name")]), collapse = "\n"))
      message("\nData resources: \n")
      message(paste0(capture.output(dataResource[order(dataResource$uid),c("uid","name")]), collapse = "\n"))
      r <- readline()
      if(r =="all"){
        res <- c(res, paste0("institution_uid",":",ruid))
        continue <- FALSE
      }else{
        if(r %in% c(dataResource$uid, collections$uid)){
          field <- switch( substring(r, 0, 2),
                           "co" = "collection_uid", 
                           "dr" = "data_resource_uid",
                           "dp" = "data_provider_uid")
          res <- c(res, paste0(field,":",r))
          
          if(continue("Filter added. Do you want to continue? Type 'y' for yes. ")){
            type <- NULL
          }else{
            continue <- FALSE
          }
          
        }else{
          if(continue("No resource with that uid was found. Do you want to start over? Type 'y' for yes. ")){
            type <- NULL
          }else{
            continue <- FALSE
          }
        }
      }
      
    }else{
      if(continue("No institution with that uid was found. Do you want to start over? Type 'y' for yes. ")){
        type <- NULL
      }else{
        continue <- FALSE
      }
    }
  }
  
  return(res)
  
}

#' auxiliary function for layers
layer_questionarie<-function(){
  continue <- TRUE
  res <- c()
  layers<-sbdi_fields("layers")
  layers_cl<-layers[grep(pattern = "cl", x = layers$id), 
                    c("id","name", "description", "type",  "shortName")]
  
  while(continue){
    layersDisp <- data.frame("name"=layers_cl$description[order(layers_cl$id)],
                             row.names=sort(layers_cl$id))
    # layersDisp <- data.frame(row.names=layers[,"id"], "name source"=paste(layers$displayname, "/", layers$source))
    ## TODO if the message if longer that 4800 something "bytes" it will be cut in the console. 
    ## Make it so that the tables gets cut and message divided in many messages
    message("\nWhich layer do you want use as filter? Type the corresponding uid: \n")
    message(paste0(capture.output(layersDisp), collapse = "\n"))
    r <- readline()
    
    if(r %in% layers_cl$id){
      lid <- r
      layer_url<-SBDI4R:::build_url_from_parts(getOption("ALA4R_server_config")$base_url_spatial, 
                                                 c("objects", lid))
      objectsLy<-jsonlite::fromJSON(httr::content(httr::GET(layer_url), as = "text", encoding = "UTF-8")) #
      objectsDisp <- data.frame("name"=objectsLy$id[order(objectsLy$pid)],
                               row.names=sort(objectsLy$pid))
      message("\nBy which object in this layer do you want to filter? Type the corresponding 'name': \n")
      message(paste0(capture.output(objectsDisp), collapse = "\n"))
      r <- readline()
      if(r %in% objectsLy$id){
        if(is.na(suppressWarnings(as.numeric(r)))) r <- paste0("%22",r,"%22")
          
        res <- c(res, paste0(lid,":",r))
        
        if(continue("Filter added. Do you want to continue? Type 'y' for yes. ")){
          type <- NULL
        }else{
          continue <- FALSE
        }
        
      }else{
        if(continue("No object with that uid was found. Do you want to start over? Type 'y' for yes. ")){
          type <- NULL
        }else{
          continue <- FALSE
        }
      }
    
    }else{
      if(continue("No layer with that uid was found. Do you want to start over? Type 'y' for yes. ")){
        type <- NULL
      }else{
        continue <- FALSE
      }
    }
  }
  
  return(res)
  
}

 
#' auxiliary function
#' @param msg the message to pass to the function  
  continue <- function(msg){
    print(msg)
    r <- readline()
    r <- tolower(r)
    return(any(c("yes", "y") %in% r))
  }
  