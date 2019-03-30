setwd("~/LAK/")
source("LAK.R")
library(Rtsne)
library(Linnorm)
library(sparcl)
library(clues)
library(cluster)
library(ggplot2)
library(factoextra)

Biase <- get_sc_data(4)[[1]]
Biase_LAK <- LAK(Biase,3)
Biase_LAK_ann <- Biase_LAK[[1]]$Cs

#five_matrix_LAK <- list(Patel_LAK,Yan_LAK,Goolam_LAK,Biase_LAK,Treutlein_LAK)
load("RData/Five_matrix_LAK.RData")

#Patel_LAK <- five_matrix_LAK[[1]]
#Yan_LAK <- five_matrix_LAK[[2]]
#Goolam_LAK <- five_matrix_LAK[[3]]
Biase_LAK <- five_matrix_LAK[[4]]
#Treutlein_LAK <- five_matrix_LAK[[5]]

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

#Fig.S2a_plot_matrix <- plot_matrix(Patel,Patel_LAK,normed = T)
#Fig.S2b_plot_matrix <- plot_matrix(Yan,Yan_LAK,normed = F)
#Fig.S2c_plot_matrix <- plot_matrix(Goolam,Goolam_LAK,normed = F)
Fig.S2d_plot_matrix <- plot_matrix(Biase,Biase_LAK,normed = F)
#Fig.S2e_plot_matrix <- plot_matrix(Treutlein,Treutlein_LAK,normed = F)

load("RData/Figure2b.RData")

cluster_result_visualize <-function(plot_matrix){ 
  ggplot(plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) }

#The visualization of the clustering results
#Fig.S2a <- cluster_result_visualize(Fig.S2a_plot_matrix);Fig.S2a
#Fig.S2b <- cluster_result_visualize(Fig.S2b_plot_matrix);Fig.S2b
#Fig.S2c <- cluster_result_visualize(Fig.S2c_plot_matrix);Fig.S2c
Fig.S2d <- cluster_result_visualize(Fig.S2d_plot_matrix);Fig.S2d
#Fig.S2e <- cluster_result_visualize(Fig.S2e_plot_matrix);Fig.S2e

Fig.S2d
