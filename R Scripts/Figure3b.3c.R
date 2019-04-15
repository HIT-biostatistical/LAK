setwd("~/LAK/")
source("R Scripts/Creat_Simulate_Data.R")
if (!require(Linnorm)) { # for normalization
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Linnorm", version = "3.8")
  library(Linnorm)
}

if (!require(cluster)) {
  install.packages("cluster")
  library(cluster)
}

if (!require(sparcl)) {
  install.packages("sparcl")
  library(sparcl)
}

if (!require(Rtsne)) {
  install.packages("Rtsne")
  library(Rtsne)
}
if (!require(clues)) {
  install.packages("clues")
  library(clues)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

#k1 <- k_detection(sim_data[[1]])
#k2 <- k_detection(sim_data[[2]])
#k3 <- k_detection(sim_data[[3]])
#k4 <- k_detection(sim_data[[4]])
#k5 <- k_detection(sim_data[[5]])
cluster_numbers <- c(3,4,5,6,7)  #cluster number by LAK
 
#cluster analysis by LAK method, the results were saved in Simulate_ARI_value.RData and Simulate_ARI_value.RData
data1<-sim_data[[1]]
data1_LAK <- LAK(data1,3,normed = T)
data1_ARI_value <- adjustedRand(data1_LAK[[1]]$Cs,get_sc_data(1)[[3]])[[2]]
data1_ARI_value

data5<-sim_data[[5]]
data5_LAK <- LAK(data5,7,normed = T)
data5_ARI_value <- adjustedRand(data5_LAK[[1]]$Cs,get_sc_data(5)[[3]])[[2]]
data5_ARI_value
 
#Simulate_ARI_value <- list(data1_ARI_value,data2_ARI_value,data3_ARI_value,
#                      data4_ARI_value,data5_ARI_value)
#Simulate_matrix_LAK <- list(data1_LAK,data2_LAK,data3_LAK,
#                       data4_LAK,data5_LAK)


load("RData/Simulate_ARI_value.RData")
load("RData/Simulate_matrix_LAK.RData")

data1_LAK <- Simulate_matrix_LAK[[1]]
#data2_LAK <- Simulate_matrix_LAK[[2]]
#data3_LAK <- Simulate_matrix_LAK[[3]]
#data4_LAK <- Simulate_matrix_LAK[[4]]
data5_LAK <- Simulate_matrix_LAK[[5]]
Simulate_ARI_value

# the two-dimensional matrix of ploting clustering result visualization
#' @params 
#' exprs_matrix refer to the express matrix of single cell data
#' LAK_ann refer to whether the data has been normed
#' normed refer to the annotation provided by LAK
#' @return ggplot_matrix
plot_matrix <- function(exprs_matrix,matrix_LAK,normed){
  if(normed == F){  
    pre_filterd_matrix <- exprs_matrix[rowSums(exprs_matrix > 0.5) > 2,]
    normed_matrix <- Linnorm(pre_filterd_matrix)
  }
  if(normed == T){
    pre_filterd_matrix <- exprs_matrix[rowSums(exprs_matrix > 0.5) > 2,]
    normed_matrix <- pre_filterd_matrix
  }
  matrix_with_LAK_genes_normed <- normed_matrix[matrix_LAK[[1]]$ws > 1e-5,]
  after_reduce_dim <- Rtsne(t(matrix_with_LAK_genes_normed),perplexity = 15)
  LAK_ann <- matrix_LAK[[1]]$Cs
  LAK_ann <- as.factor(LAK_ann)
  ggplot_matrix <- data.frame(after_reduce_dim$Y,cluster = LAK_ann)
  return(ggplot_matrix)  
}

Fig.S4a_plot_matrix <- plot_matrix(data1,data1_LAK,normed = T)
#Fig.S4b_plot_matrix <- plot_matrix(data2,data2_LAK,normed = T)
#Fig.S4c_plot_matrix <- plot_matrix(data3,data3_LAK,normed = T)
#Fig.S4d_plot_matrix <- plot_matrix(data4,data4_LAK,normed = T)
Fig.S4e_plot_matrix <- plot_matrix(data5,data5_LAK,normed = T)

load("RData/Figure3b 3c.RData")

cluster_result_visualize <-function(plot_matrix){ 
  ggplot(plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
    geom_point() +
    labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) }

#The visualization of the clustering results
Fig.S4a <- cluster_result_visualize(Fig.S5a_plot_matrix);Fig.S4a
#Fig.S4b <- cluster_result_visualize(Fig.S5b_plot_matrix);Fig.S4b
#Fig.S4c <- cluster_result_visualize(Fig.S5c_plot_matrix);Fig.S4c
#Fig.S4d <- cluster_result_visualize(Fig.Sd_plot_matrix);Fig.S4d
Fig.S4e <- cluster_result_visualize(Fig.S5e_plot_matrix);Fig.S4e

Fig.S4a
Fig.S4e
