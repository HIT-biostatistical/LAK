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

LAK_ann <- LAK(ES, num_cluster = 5 ,normed = T)

ARIS <- c()
for(k in 1:100){
  LAK_ann <- LAK(ES,num_cluster = 5 , normed = T , s_value = 18.2875)
  ARI <-adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)
  cat(ARI,'\n')
  ARIS <- c(ARIS, ARI)
}

#Output the LAK_ann and after_reduce_dim result under ARI=0.940, 0.827, 0.706 
#LAK_ann <- LAK(ES, num_cluster = 5 ,normed = T)
#ari <-adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)
#after_reduce_dim <- Rtsne(t(ES))

load("RData/Figure_S4.RData")

#plot visualization of cluster result under ARI = 0.940
S4a_ES_plot_matrix <- data.frame(S4a_after_reduce_dim$Y,cluster = as.factor(S4a_LAK_ann[[1]]$Cs))
Fig.S4a <- ggplot(S4a_ES_plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Fig.S4a

#plot visualization of cluster result under ARI = 0.827
S4b_ES_plot_matrix <- data.frame(S4b_after_reduce_dim$Y,cluster = as.factor(S4b_LAK_ann[[1]]$Cs))
Fig.S4b <- ggplot(S4b_ES_plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Fig.S4b

#plot visualization of cluster result under ARI = 0.706
S4c_ES_plot_matrix <- data.frame(S4c_after_reduce_dim$Y,cluster = as.factor(S4c_LAK_ann[[1]]$Cs))
Fig.S4c <- ggplot(S4c_ES_plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Fig.S4c

