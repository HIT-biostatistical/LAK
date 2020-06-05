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

Camp1 <- get_sc_data(6)[[1]]
Camp1_LAK <- LAK(Camp1,6)
Camp1_LAK_ann <- Camp1_LAK[[1]]$Cs

Li <- get_sc_data(7)[[1]]
Li_LAK <- LAK(Li,7)
Li_LAK_ann <- Li_LAK[[1]]$Cs


#five_matrix_LAK <- list(Patel_LAK,Yan_LAK,Goolam_LAK,Biase_LAK,Treutlein_LAK)
load("RData/Five_matrix_LAK.RData")
load("RData/Camp1_Li_matrix_lak.RData")

#Patel_LAK <- five_matrix_LAK[[1]]
#Yan_LAK <- five_matrix_LAK[[2]]
#Goolam_LAK <- five_matrix_LAK[[3]]
Biase_LAK <- five_matrix_LAK[[4]]
#Treutlein_LAK <- five_matrix_LAK[[5]]
#Camp1_LAK <- camp1_li_matrix_LAK[[1]]
#Li_LAK <- camp1_li_matrix_LAK[[2]]


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
  after_reduce_dim <- Rtsne(t(matrix_with_LAK_genes_normed),perplexity = 15,check_duplicates = F)
  LAK_ann <- matrix_LAK[[1]]$Cs
  LAK_ann <- as.factor(LAK_ann)
  ggplot_matrix <- data.frame(after_reduce_dim$Y,cluster = LAK_ann)
  return(ggplot_matrix)  
}

#Fig.S1a_plot_matrix <- plot_matrix(Patel,Patel_LAK,normed = T)
#Fig.S1b_plot_matrix <- plot_matrix(Yan,Yan_LAK,normed = F)
#Fig.S1c_plot_matrix <- plot_matrix(Goolam,Goolam_LAK,normed = F)
Fig.S1d_plot_matrix <- plot_matrix(Biase,Biase_LAK,normed = F)
#Fig.S1e_plot_matrix <- plot_matrix(Treutlein,Treutlein_LAK,normed = F)
#Fig.S1f_plot_matrix <- plot_matrix(Camp1,Camp1_LAK,normed = F)
#Fig.S1g_plot_matrix <- plot_matrix(Li,Li_LAK,normed = F)

load("RData/Figure2b.RData")

cluster_result_visualize <-function(plot_matrix){ 
  ggplot(plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) }

#The visualization of the clustering results
#Fig.S1a <- cluster_result_visualize(Fig.S1a_plot_matrix);Fig.S1a
#Fig.S1b <- cluster_result_visualize(Fig.S1b_plot_matrix);Fig.S1b
#Fig.S1c <- cluster_result_visualize(Fig.S1c_plot_matrix);Fig.S1c
Fig.S1d <- cluster_result_visualize(Fig.S1d_plot_matrix);Fig.S1d
#Fig.S1e <- cluster_result_visualize(Fig.S1e_plot_matrix);Fig.S1e
#Fig.S1f <- cluster_result_visualize(Fig.S1f_plot_matrix);Fig.S1f
#Fig.S1g <- cluster_result_visualize(Fig.S1g_plot_matrix);Fig.S1g

Fig.2c <- Fig.S1d;Fig.2c
