setwd("~/LAK")
source("LAK.R")
source("R Scripts/get_data.R")
library( plyr )
library(SingleCellExperiment)
library(mclust)

# st <-F 
# for (i in 1:5){
#   sc_new <- get_sc_data_new(i)
#   if(st == T){
#     ES_ite <- sc_new[[1]]
#     if(ncol(ES_ite) > 100){
#       ES_ite <- ES_ite[,sample(ncol(ES_ite),100)]
#     }
#     ES_ite <- as.data.frame(ES_ite)
#     ES_ite$gene <- rownames(ES_ite)
#     cat("ES_ite", i,'\n')
#     ES <- join(ES, ES_ite, type = "inner", by = "gene")
#     cat("join", "\n")
#     cell_types <- c(cell_types,rep(i,ncol(ES_ite)-1))
#   }else{
#     
#     ES <- sc_new[[1]]
#     
#     if(ncol(ES) > 100)  
#     {
#       ES <- ES[,sample(ncol(ES),100)]
#     }
#     
#     ES <- as.data.frame(ES)
#     cell_types <- rep(i,ncol(ES))
#     ES$gene <- rownames(ES)
#     cat("ES", i)
#     st <- T
#   }
# }
# 
# 
# rm(ES_ite,sc_new)
# ES <- subset(ES, select= -gene)
#save(ES, file = "RData")
load("RData/human_data-5.RData")

#LAK_ann <- LAK(ES, num_cluster = 5 ,normed = T)
#after_reduce_dim <- Rtsne(t(ES))

load("RData/Figure_S5.RData")

ES_plot_matrix <- data.frame(after_reduce_dim$Y,cluster = as.factor(LAK_ann[[1]]$Cs))
ggplot(ES_plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


ARIS <- c()
for(k in 1:100){
  LAK_ann <- LAK(ES,num_cluster = 5 , normed = T , s_value = 18.2875)
  ARI <-adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)
  cat(ARI,'\n')
  ARIS <- c(ARIS, ARI)
}

load("RData/match__human_5_tissue_ARI_results.RData")



