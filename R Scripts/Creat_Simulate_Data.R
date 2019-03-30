#Simulate from a Multivariate Normal Distribution 
source("LAK.R")
if (!require(MASS)) {
  install.packages("MASS")
  library(MASS)
}

#Create simulate data
#' @param data_id (integer) refers to the id of simulate data
#' @return the simulate data
set_sc_data <- function(data_id){
  if(data_id == 1){
  #data1 k = 3
  Sigma1 <- diag(c(rep(4,100),rep(2,300),rep(3,200))) 
  data <- mvrnorm(n = 15000,c(rep(1,100),rep(2,300),rep(3,200)),Sigma1)
  rownames(data) <- c(1:nrow(data));colnames(data) <- c(1:ncol(data))
  }
  if(data_id == 2){
  #data2 k = 4
  Sigma2 <- diag(c(rep(9,150),rep(9,200),rep(4,100),rep(9,200))) 
  data <- mvrnorm(n = 15000,c(rep(0,150),rep(1,200),rep(2,100),rep(3,200)),Sigma2)
  rownames(data) <- c(1:nrow(data));colnames(data) <- c(1:ncol(data))
  }
  if(data_id == 3){
  #data3 k = 5
  Sigma3 <- diag(c(rep(2,200),rep(1,100),rep(1,100),rep(1,150),rep(2,100))) 
  data <- mvrnorm(n = 15000,c(rep(0,200),rep(1,100),rep(2,100),rep(3,150),rep(4,100)),Sigma3)
  rownames(data) <- c(1:nrow(data));colnames(data) <- c(1:ncol(data))
  }
  if(data_id == 4){
  #data4 k = 6
  Sigma4 <- diag(c(rep(1,100),rep(1,150),rep(1,100),rep(1,60),rep(1,170),rep(1,100))) 
  data <- mvrnorm(n = 15000,c(rep(0,100),rep(1,150),rep(2.5,100),rep(4,60),rep(5,170),rep(6,100)),Sigma4)
  rownames(data) <- c(1:nrow(data));colnames(data) <- c(1:ncol(data))
  }
  if(data_id == 5){
  #data5 k = 7
  Sigma5 <- diag(c(rep(9,100),rep(9,120),rep(6,150),rep(9,70),rep(6,100),rep(9,130),rep(6,200))) 
  data <- mvrnorm(n = 15000,c(rep(0,100),rep(2,120),rep(4,150),rep(6,70),rep(8,100),rep(10,130),rep(12,200)),Sigma5)
  rownames(data) <- c(1:nrow(data));colnames(data) <- c(1:ncol(data))
  }
  return(data)
}
#data1 <- set_sc_data(1)
#data2 <- set_sc_data(2)
#data3 <- set_sc_data(3)
#data4 <- set_sc_data(4)
#data5 <- set_sc_data(5)

#sim_data <- list(data1,data2,data3,data4,data5)
load("RData/Simulate_data.RData")

data1 <- sim_data[[1]]
data2 <- sim_data[[2]]
data3 <- sim_data[[3]]
data4 <- sim_data[[4]]
data5 <- sim_data[[5]]

#data preprocessing 
#' @param data_id (integer) refers to the id of simulate data
#' @returns 
#' the express matrix of input data
#' whether the data has been normed
#' the annotation provided by the original author
get_sc_data <- function(data_id){
  if(data_id == 1){
  #data1 k = 3
  exprs_matrix <- data1
  label <- c(rep(1,100),rep(2,300),rep(3,200))
  normed <- T
  }
  if(data_id == 2){
  #data2 k = 4
  exprs_matrix <- data2
  label <- c(rep(1,150),rep(2,200),rep(3,100),rep(4,200))
  normed <- T
  }
  if(data_id == 3){
  #data3 k = 5
  exprs_matrix <- data3
  label <- c(rep(1,200),rep(2,100),rep(3,100),rep(4,150),rep(5,100))
  normed <- T
  }
  if(data_id == 4){
  #data4 k = 6
  exprs_matrix <- data4
  label <- c(rep(1,100),rep(2,150),rep(3,100),rep(4,60),rep(5,170),rep(6,100))
  normed <- T
  }
  if(data_id == 5){
  #data5 k = 6
  exprs_matrix <- data5
  label <- c(rep(1,100),rep(2,120),rep(3,150),rep(4,70),rep(5,100),rep(6,130),rep(7,200))
  normed <- T
  }
  result <- list(exprs_matrix,normed,label)
  return(result)
}




