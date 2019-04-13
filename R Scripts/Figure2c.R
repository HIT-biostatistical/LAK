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

##Patel:get differentially expressed genes by t-test
#id = 1
#Patel_k <- cluster_numbers[id]
#Patel <- get_sc_data(id)[[1]] #get preprocessed data
#Patel_normed <- Patel
#Patel_LAK_ann <- five_matrix_LAK[[1]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Patel;k <- Patel_k;normed_matrix <- Patel_normed;LAK_ann <- Patel_LAK_ann #unified variable names
#T_Patel_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#Patel_T_Test <- T_matrix(Patel_normed,Patel_LAK_ann)
#rownames(Patel_T_Test) <- row.names(Patel)
#T_Patel_diffgene_Pvalue <- Patel_T_Test[T_Patel_diff_gene_names,]
#rm(id);rm(exprs_matrix);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Yan:get differentially expressed genes by t-test
#id = 2
#Yan_k <- cluster_numbers[id]
#Yan <- get_sc_data(id)[[1]] #get preprocessed data
#Yan_normed <- Linnorm(Yan) 
#Yan_LAK_ann <- five_matrix_LAK[[2]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Yan;k <- Yan_k;normed_matrix <- Yan_normed;LAK_ann <- Yan_LAK_ann ##unified variable names
#T_Yan_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#Yan_T_Test <- T_matrix(Yan_normed,Yan_LAK_ann)
#rownames(Yan_T_Test) <- row.names(Yan)
#T_Yan_diffgene_Pvalue <- Yan_T_Test[T_Yan_diff_gene_names,]
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Goolam:get differentially expressed genes by t-test
#id = 3
#Goolam_k <- cluster_numbers[id]
#Goolam <- get_sc_data(id)[[1]] #get preprocessed data
#Goolam_normed <- Linnorm(Goolam) 
#Goolam_LAK_ann <- five_matrix_LAK[[3]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Goolam;k <- Goolam_k;normed_matrix <- Goolam_normed;LAK_ann <- Goolam_LAK_ann #unified variable names
#T_Goolam_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#Goolam_T_Test <- T_matrix(Goolam_normed,Goolam_LAK_ann)
#rownames(Goolam_T_Test) <- row.names(Goolam)
#T_Goolam_diffgene_Pvalue <- Goolam_T_Test[T_Goolam_diff_gene_names,]
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

#Biase:get differentially expressed genes by t-test
id = 4
Biase_k <- cluster_numbers[id]
Biase <- get_sc_data(id)[[1]] #get preprocessed data
rownames(Biase)[c(which(rownames(Biase)=="ENSMUSG00000088604"),
                  which(rownames(Biase)=="ENSMUSG00000093678"))] <- c("Gm25820","Gm20608")
Biase_normed <- Linnorm(Biase) 
Biase_LAK_ann <- five_matrix_LAK[[4]][[1]]$Cs #get cluster annotion of our method
exprs_matrix <- Biase;k <- Biase_k;normed_matrix <- Biase_normed;LAK_ann <- Biase_LAK_ann #unified variable names
T_Biase_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
Biase_T_Test <- T_matrix(Biase_normed,Biase_LAK_ann)
rownames(Biase_T_Test) <- row.names(Biase)
T_Biase_diffgene_Pvalue <- Biase_T_Test[T_Biase_diff_gene_names,]
rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Treutlein:get differentially expressed genes by t-test
#id = 5
#Treutlein_k <- cluster_numbers[id]
#Treutlein <- get_sc_data(id)[[1]] #get preprocessed data
#Treutlein_normed <- Linnorm(Treutlein) 
#Treutlein_LAK_ann <- five_matrix_LAK[[5]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Treutlein;k <- Treutlein_k;normed_matrix <- Treutlein_normed;LAK_ann <- Treutlein_LAK_ann #unified variable names
#T_Treutlein_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#Treutlein_T_Test <- T_matrix(Treutlein_normed,Treutlein_LAK_ann)
#rownames(Treutlein_T_Test) <- row.names(Treutlein)
#T_Treutlein_diffgene_Pvalue <- Treutlein_T_Test[T_Treutlein_diff_gene_names,]
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

load("RData/T_five_diff_gene_names.RData")
load("RData/W_five_diff_gene_names.RData")

#look up the number of each type
#table(Patel_LAK_ann)
#table(Yan_LAK_ann)
#table(Goolam_LAK_ann)
table(Biase_LAK_ann)
#table(Treutlein_LAK_ann)

##Patel:plot heatmap according to differentially expressed genes by t-test
#Fig.S3a <- plot_heatmap(T_Patel_diff_gene_names,1,Patel_k,Patel_normed,Patel_LAK_ann,numsDiffGene=10);Fig.S3a

##Yan:plot heatmap according to differentially expressed genes by t-test
#Fig.S3b <- plot_heatmap(T_Yan_diff_gene_names,2,Yan_k,Yan_normed,Yan_LAK_ann,numsDiffGene=10);Fig.S3b

##Goolam:plot heatmap according to differentially expressed genes by t-test
#Fig.S3c <- plot_heatmap(T_Goolam_diff_gene_names,3,Goolam_k,Goolam_normed,Goolam_LAK_ann,numsDiffGene=10);Fig.S3c

#Biase:plot heatmap according to differentially expressed genes by t-test
Fig.S3d <- plot_heatmap(T_Biase_diff_gene_names,4,Biase_k,Biase_normed,Biase_LAK_ann,numsDiffGene=10);Fig.S3d

##Treutlein:plot heatmap according to differentially expressed genes by t-test
#Fig.S3e <- plot_heatmap(T_Treutlein_diff_gene_names,5,Treutlein_k,Treutlein_normed,Treutlein_LAK_ann,numsDiffGene=10);Fig.S3e

