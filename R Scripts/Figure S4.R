setwd("~/LAK")
source("LAK.R")
source("R Scripts/get_data.R")
library( plyr )
library(SingleCellExperiment)
library(mclust)

st <-F
for (i in 1:5){
  sc_new <- get_sc_data_new(i)
  if(st == T){
    ES_ite <- sc_new[[1]]
    if(ncol(ES_ite) > 100){
      ES_ite <- ES_ite[,sample(ncol(ES_ite),100)]
    }
    ES_ite <- as.data.frame(ES_ite)
    ES_ite$gene <- rownames(ES_ite)
    cat("ES_ite", i,'\n')
    ES <- join(ES, ES_ite, type = "inner", by = "gene")
    cat("join", "\n")
    cell_types <- c(cell_types,rep(i,ncol(ES_ite)-1))
  }else{

    ES <- sc_new[[1]]

    if(ncol(ES) > 100)
    {
      ES <- ES[,sample(ncol(ES),100)]
    }

    ES <- as.data.frame(ES)
    cell_types <- rep(i,ncol(ES))
    ES$gene <- rownames(ES)
    cat("ES", i)
    st <- T
  }
}


rm(ES_ite,sc_new)
ES <- subset(ES, select= -gene)
#load("RData/human_data-5-new.RData")

LAK_ann <- LAK(ES, num_cluster = 5 ,normed = T)

ARIS <- adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)

reco <- c()
records <- list()
n <- 1
for(k in 1:100){
  LAK_ann <- LAK(ES,num_cluster = 5 , normed = T , s_value =  19.03125)
  ARI <-adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)
  cat(ARI,'\n')
  ARIS <- c(ARIS, ARI)
  if(!(ARI %in% reco)){
    reco <- c(reco, ARI)
    records[[n]] <- LAK_ann
    n <- n+1
  }
}

table(ARIS)

LAK_ann <- records[[1]]
adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)
tsne_df <- Rtsne(t(ES),perplexity = 15, dims = 2)
LAK_ann_f <- as.factor(LAK_ann[[1]]$Cs)
plot_matrix <- data.frame(tsne_df$Y,cluster = LAK_ann_f)
ggplot(plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

LAK_ann <- records[[2]]
adjustedRandIndex(LAK_ann[[1]]$Cs, cell_types)
tsne_df <- Rtsne(t(ES),perplexity = 15, dims = 2)
LAK_ann_f <- as.factor(LAK_ann[[1]]$Cs)
plot_matrix <- data.frame(tsne_df$Y,cluster = LAK_ann_f)
ggplot(plot_matrix,aes(x = X1,y = X2,shap = cluster,color = cluster)) +
  geom_point() +
  labs(x = "tSNE1",y = "tSNE2") + theme_classic() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 



#save(ARIS,file = "RData/match__human_5_tissue_ARI_results.RData")
load("RData/match__human_5_tissue_ARI_results.RData")
table(ARIS)
