setwd("~/LAK/")
source("LAK.R")
library(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#library(limma)
library(Linnorm)
library(Rtsne)
library(clues)
library(pheatmap)

load("RData/Five_matrix_LAK.RData")

cluster_numbers<-c(5,7,6,3,6)#The number of clustering we calculate

##Patel:get differentially expressed genes by t-test and W-test
#id = 1
#Patel_k <- cluster_numbers[id]
#Patel <- get_sc_data(id)[[1]] #get preprocessed data
#Patel_normed <- Patel
#Patel_LAK_ann <- five_matrix_LAK[[1]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Patel;k <- Patel_k;normed_matrix <- Patel_normed;LAK_ann <- Patel_LAK_ann #unified variable names
#T_Patel_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#W_Patel_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,F,numsDiffGene=10) #get differentially expressed gene names of ordering by W-test
#rm(id);rm(exprs_matrix);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Yan:get differentially expressed genes by t-test and W-test
#id = 2
#Yan_k <- cluster_numbers[id]
#Yan <- get_sc_data(id)[[1]] #get preprocessed data
#Yan_normed <- Linnorm(Yan) 
#Yan_LAK_ann <- five_matrix_LAK[[2]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Yan;k <- Yan_k;normed_matrix <- Yan_normed;LAK_ann <- Yan_LAK_ann ##unified variable names
#T_Yan_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#W_Yan_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,F,numsDiffGene=10) #get differentially expressed gene names of ordering by W-test
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Goolam:get differentially expressed genes by t-test and W-test
#id = 3
#Goolam_k <- cluster_numbers[id]
#Goolam <- get_sc_data(id)[[1]] #get preprocessed data
#Goolam_normed <- Linnorm(Goolam) 
#Goolam_LAK_ann <- five_matrix_LAK[[3]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Goolam;k <- Goolam_k;normed_matrix <- Goolam_normed;LAK_ann <- Goolam_LAK_ann #unified variable names
#T_Goolam_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#W_Goolam_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,F,numsDiffGene=10) #get differentially expressed gene names of ordering by W-test
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

#Biase:get differentially expressed genes by t-test and W-test
id = 4
Biase_k <- cluster_numbers[id]
Biase <- get_sc_data(id)[[1]] #get preprocessed data
Biase_normed <- Linnorm(Biase) 
Biase_LAK_ann <- five_matrix_LAK[[4]][[1]]$Cs #get cluster annotion of our method
exprs_matrix <- Biase;k <- Biase_k;normed_matrix <- Biase_normed;LAK_ann <- Biase_LAK_ann #unified variable names
T_Biase_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
W_Biase_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,F,numsDiffGene=10) #get differentially expressed gene names of ordering by W-test
rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Treutlein:get differentially expressed genes by t-test and W-test
#id = 5
#Treutlein_k <- cluster_numbers[id]
#Treutlein <- get_sc_data(id)[[1]] #get preprocessed data
#Treutlein_normed <- Linnorm(Treutlein) 
#Treutlein_LAK_ann <- five_matrix_LAK[[5]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Treutlein;k <- Treutlein_k;normed_matrix <- Treutlein_normed;LAK_ann <- Treutlein_LAK_ann #unified variable names
#T_Treutlein_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#W_Treutlein_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,F,numsDiffGene=10) #get differentially expressed gene names of ordering by W-test
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Integrate all data difference expression genes in order to save results
#T_five_diff_gene_names <- list(T_Patel_diff_gene_names,T_Yan_diff_gene_names,T_Goolam_diff_gene_names,
#                             T_Biase_diff_gene_names,T_Treutlein_diff_gene_names) 
#W_five_diff_gene_names <- list(W_Patel_gene_names,W_Yan_diff_gene_names,W_Goolam_diff_gene_names,
#                             W_Biase_diff_gene_names,W_Treutlein_diff_gene_names)


load("RData/T_five_diff_gene_names.RData")
load("RData/W_five_diff_gene_names.RData")

#look up the number of each type
#table(Patel_LAK_ann)
#table(Yan_LAK_ann)
#table(Goolam_LAK_ann)
table(Biase_LAK_ann)
#table(Treutlein_LAK_ann)

##Patel:plot heatmap according to differentially expressed genes by t-test and W-test
#T_Patel_diff_gene_names <- T_five_diff_gene_names[[1]] #get difference expression gene from the results that we save
#W_Patel_diff_gene_names <- W_five_diff_gene_names[[1]] #get difference expression gene from the results that we save
#Fig.S3a <- plot_heatmap(T_Patel_diff_gene_names,1,Patel_k,Patel_normed,Patel_LAK_ann,numsDiffGene=10);Fig.S3a
#Fig.S4a <- plot_heatmap(W_Patel_diff_gene_names,1,Patel_k,Patel_normed,Patel_LAK_ann,numsDiffGene=10);Fig.S4a

##Yan:plot heatmap according to differentially expressed genes by t-test and W-test
#T_Yan_diff_gene_names <- T_five_diff_gene_names[[2]] #get difference expression gene from the results that we save
#W_Yan_diff_gene_names <- W_five_diff_gene_names[[2]] #get difference expression gene from the results that we save
#Fig.S3b <- plot_heatmap(T_Yan_diff_gene_names,2,Yan_k,Yan_normed,Yan_LAK_ann,numsDiffGene=10);Fig.S3b
#Fig.S4b <- plot_heatmap(W_Yan_diff_gene_names,2,Yan_k,Yan_normed,Yan_LAK_ann,numsDiffGene=10);Fig.S4b

##Goolam:plot heatmap according to differentially expressed genes by t-test and W-test
#T_Goolam_diff_gene_names <- T_five_diff_gene_names[[3]] #get difference expression gene from the results that we save
#W_Goolam_diff_gene_names <- W_five_diff_gene_names[[3]] #get difference expression gene from the results that we save
#Fig.S3c <- plot_heatmap(T_Goolam_diff_gene_names,3,Goolam_k,Goolam_normed,Goolam_LAK_ann,numsDiffGene=10);Fig.S3c
#Fig.S4c <- plot_heatmap(W_Goolam_diff_gene_names,3,Goolam_k,Goolam_normed,Goolam_LAK_ann,numsDiffGene=10);Fig.S4c

#Biase:plot heatmap according to differentially expressed genes by t-test and W-test
T_Biase_diff_gene_names <- T_five_diff_gene_names[[4]] #get difference expression gene from the results that we save
W_Biase_diff_gene_names <- W_five_diff_gene_names[[4]] #get difference expression gene from the results that we save
Fig.S3d <- plot_heatmap(T_Biase_diff_gene_names,4,Biase_k,Biase_normed,Biase_LAK_ann,numsDiffGene=10);Fig.S3d
Fig.S4d <- plot_heatmap(W_Biase_diff_gene_names,4,Biase_k,Biase_normed,Biase_LAK_ann,numsDiffGene=10);Fig.S4d

##Treutlein:plot heatmap according to differentially expressed genes by t-test and W-test
#T_Treutlein_diff_gene_names <- T_five_diff_gene_names[[5]] #get difference expression gene from the results that we save
#W_Treutlein_diff_gene_names <- W_five_diff_gene_names[[5]] #get difference expression gene from the results that we save
#Fig.S3e <- plot_heatmap(T_Treutlein_diff_gene_names,5,Treutlein_k,Treutlein_normed,Treutlein_LAK_ann,numsDiffGene=10);Fig.S3e
#Fig.S4e <- plot_heatmap(W_Treutlein_diff_gene_names,5,Treutlein_k,Treutlein_normed,Treutlein_LAK_ann,numsDiffGene=10);Fig.S4e

