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
load("RData/Camp1_Li_matrix_LAK.RData")


cluster_numbers<-c(5,7,6,3,6, 6, 7)#The number of clustering we calculate

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

##Camp1:get differentially expressed genes by t-test
#id = 6
#Camp1_k <- cluster_numbers[id]
#Camp1 <- get_sc_data(id)[[1]] #get preprocessed data
#Camp1_normed <- Linnorm(Camp1) 
#Camp1_LAK_ann <- camp1_li_matrix_LAK[[1]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Camp1;k <- Camp1_k;normed_matrix <- Camp1_normed;LAK_ann <- Camp1_LAK_ann ##unified variable names
#T_Camp1_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#Camp1_T_Test <- T_matrix(Camp1_normed,Camp1_LAK_ann)
#rownames(Camp1_T_Test) <- row.names(Camp1)
#T_Camp1_diffgene_Pvalue <- Camp1_T_Test[T_Camp1_diff_gene_names,]
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

##Li:get differentially expressed genes by t-test
#id = 7
#Li_k <- cluster_numbers[id]
#Li <- get_sc_data(id)[[1]] #get preprocessed data
#Li_normed <- Linnorm(Li) 
#Li_LAK_ann <- camp1_li_matrix_LAK[[2]][[1]]$Cs #get cluster annotion of our method
#exprs_matrix <- Li;k <- Li_k;normed_matrix <- Li_normed;LAK_ann <- Li_LAK_ann ##unified variable names
#T_Li_diff_gene_names <- diff_gene_names_reorder(exprs_matrix,k,normed_matrix,LAK_ann,T,numsDiffGene=10) #get differentially expressed gene names of ordering by T-test
#Li_T_Test <- T_matrix(Li_normed,Li_LAK_ann)
#rownames(Li_T_Test) <- row.names(Li)
#T_Li_diffgene_Pvalue <- Li_T_Test[T_Li_diff_gene_names,]
#rm(id);rm(exprs_matrix);rm(k);rm(normed_matrix);rm(LAK_ann) #delete variable names

load("RData/T_five_diff_gene_names.RData")
load("RData/T_Camp1_Li_diff_gene_names.RData")

#look up the number of each type
#table(Patel_LAK_ann)
#table(Yan_LAK_ann)
#table(Goolam_LAK_ann)
table(Biase_LAK_ann)
#table(Treutlein_LAK_ann)
#table(Camp1_LAK_ann)
#table(Li_LAK_ann)

##Patel:plot heatmap according to differentially expressed genes by t-test
#Fig.S2a <- plot_heatmap(T_Patel_diff_gene_names,1,Patel_k,Patel_normed,Patel_LAK_ann,numsDiffGene=10);Fig.S2a

##Yan:plot heatmap according to differentially expressed genes by t-test
#Fig.S2b <- plot_heatmap(T_Yan_diff_gene_names,2,Yan_k,Yan_normed,Yan_LAK_ann,numsDiffGene=10);Fig.S2b

##Goolam:plot heatmap according to differentially expressed genes by t-test
#Fig.S2c <- plot_heatmap(T_Goolam_diff_gene_names,3,Goolam_k,Goolam_normed,Goolam_LAK_ann,numsDiffGene=10);Fig.S2c

#Biase:plot heatmap according to differentially expressed genes by t-test
Fig.S2d <- plot_heatmap(T_Biase_diff_gene_names,4,Biase_k,Biase_normed,Biase_LAK_ann,numsDiffGene=10);Fig.S2d
Fig.2d <- Fig.S2d;Fig.2d
##Treutlein:plot heatmap according to differentially expressed genes by t-test
#Fig.tre <- plot_heatmap(T_Treutlein_diff_gene_names,5,Treutlein_k,Treutlein_normed,Treutlein_LAK_ann,numsDiffGene=10);Fig.tre

##Camp1:plot heatmap according to differentially expressed genes by t-test
#Fig.S2e <- plot_heatmap(T_Camp1_diff_gene_names,6,Camp1_k,Camp1_normed,Camp1_LAK_ann,numsDiffGene=10);Fig.S2e

##Li:plot heatmap according to differentially expressed genes by t-test
#Fig.S2f <- plot_heatmap(T_Li_diff_gene_names,7,Li_k,Li_normed,Li_LAK_ann,numsDiffGene=10);Fig.S2f

