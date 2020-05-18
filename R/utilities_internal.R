## some utility functions used internally within the SBDI4R library: not exported

##----------------------------------------------------------------------------------------------
empty <- function(x) {  ALA4R:::empty(x) }

is.notempty.string <- function(x) {
  
   ALA4R:::is.notempty.string(x)
  
}
##----------------------------------------------------------------------------------------------

## internal function for converting chr data types to numeric or logical
convert_dt <- function(x,test_numeric=TRUE) {
  
  ALA4R:::convert_dt(x,test_numeric)
  
}
##----------------------------------------------------------------------------------------------

clean_string <- function(x) {
  
  ALA4R:::clean_string(x)
 
}

##----------------------------------------------------------------------------------------------

##convert to camel case ... modified from help forum example
## not exported for users: internal SBDI4R use only

tocamel <- function(x, delim = "[^[:alnum:]]", upper = FALSE, sep = "") {
  
  ALA4R:::tocamel(x, delim, upper, sep)
}

##----------------------------------------------------------------------------------------------

## define column names that we will remove from the results because we don't think they will be useful in the SBDI4R context

unwanted_columns <- function(type) {
  ALA4R:::unwanted_columns(type)
}

##----------------------------------------------------------------------------------------------

rename_variables <- function(varnames,type,verbose=sbdi_config()$verbose) {
  
 ALA4R:::rename_variables(varnames,type,verbose) 
}

##--------------------------------------------------------------------------------------------

## construct url path, taking care to remove multiple forward slashes, leading slash
clean_path <- function(...,sep="/") {
  
 ALA4R:::clean_path(...,sep) 
}

##---------------------------------------------------------------------------------------------
## convenience function for building urls
## pass path in one of several ways
##  as single string: build_url_from_parts(base_url,"path/to/thing")
##  as a character vector or list: build_url_from_parts(base_url,c("path","to","thing"))
##  or a combination
build_url_from_parts <- function(base_url,path=NULL,query=list()) {
  
  ALA4R:::build_url_from_parts(base_url,path,query)
}


## wrapper around read.csv but suppressing "incomplete final line" warning
read_csv_quietly <- function(...) {
  ALA4R:::read_csv_quietly(...) 
  
}

replace_nonbreaking_spaces <- function(s) {
  
  ALA4R:::replace_nonbreaking_spaces(s)
}
