#data preprocessing 
#' @param data_id (integer) refers to the id of single cell data
#' @returns 
#' the express matrix of input data, the rows correspond to genes, and coloumns correspond to cells
#' whether the data has been normed by the original author
#' the annotation provided by the original author
get_sc_data_new <- function(data_id){
  
  #human
  
  if( data_id == 1){
    SCEset <- readRDS("Single Cell Data/darmanis.rds")
    m <- assays(SCEset)[[2]]
    rownames(m) <- rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    m <- m[,ann == "neurons"]
    cell_type <- "Human neurons"
  }
  
  if( data_id == 2){
    SCEset <- readRDS("Single Cell Data/yan.rds")
    m <- assays(SCEset)[[2]]
    rownames(m) <- rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    cell_type <- "Human Embryo Devel"
  }
  
   if( data_id == 3){
    SCEset <- readRDS("Single Cell Data/camp1.rds")
    m <- assays(SCEset)[[1]]
    rownames(m) <- rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    cell_type <- "Human Liver"
  }
  
   if(data_id == 4){
    SCEset <- readRDS("Single Cell Data/xin.rds")
    m <- assays(SCEset)[[2]]
    rownames(m) <- rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    cell_type <- "Human Pancreas"
  }

      if(data_id == 5){
    SCEset <- readRDS("Single Cell Data/pollen.rds")
    m <- assays(SCEset)[[2]]
    rownames(m) <- rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    cell_type <- "Human Tissues"
  }
  
  result <- list(m, ann, cell_type)
  return(result)
}

