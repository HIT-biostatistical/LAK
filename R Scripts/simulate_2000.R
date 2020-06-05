
library(splatter)
source("LAK.R")
library(SC3)
library(mclust)

dropout <- function(x){
  length(which(x==0))/(ncol(x)*nrow(x))
}


#params.groups1 <- newSplatParams(batchCells = 2000,  nGenes = 20000,dropout.type="none")
# One small group, one big group
#sim1 <- splatSimulateGroups(params.groups1, group.prob = c(0.3,0.3,0.3, 0.1),
#                            verbose = FALSE)
sim1 <- readRDS("simdata1.rds")
simdata1 <- assays(sim1)$counts
assays(sim1)$logcounts <- log2(as.matrix(simdata1) + 1)
colnames(rowData(sim1))[1]<-"feature_symbol"
dropout(simdata1)
simdata1_ann <- colData(sim1)$Group
ari.lak.simdata1 <- c()
for(i in 1:100){
  simdata1_LAK <- LAK(simdata1, 5,s_error = 0.05, s_value = 19.43789 )
  simdata1_LAK_ann <- simdata1_LAK[[1]]$Cs
  simdata1_LAK_ARI <- adjustedRandIndex(as.integer(as.factor(simdata1_ann)),simdata1_LAK_ann)
  ari.lak.simdata1 <- c(ari.lak.simdata1,simdata1_LAK_ARI)
  cat(i)
}
ari.lak.simdata1
ari.sc3.simdata1 <- c()
for(i in 1:100){
  #simdata1_SCE <- readRDS("simdata1.rds")#14878*1886
  simdata1_sc3 <- sc3(sim1,5, gene_filter = T)
  simdata1_sc3_ann <- colData(simdata1_sc3)$sc3_5_clusters
  simdata1_sc3_ARI <- adjustedRandIndex(as.vector(simdata1_sc3_ann),as.integer(as.factor(simdata1_ann)))
  ari.sc3.simdata1 <- c(ari.sc3.simdata1, simdata1_sc3_ARI)
  cat(i)
  
}
ari.sc3.simdata1
ari_data1 <- data.frame(ari.lak.simdata1,ari.sc3.simdata1)
write.table(ari_data1,file="ari_data1.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=F)
#saveRDS(sim1,file="simdata1.rds")


#params.groups2 <- newSplatParams(batchCells = 2000, nGenes = 20000,dropout.type="experiment",
#                                 dropout.mid=0, dropout.shape = -1)
#sim2 <- splatSimulateGroups(params.groups2, group.prob = c(0.3,0.3,0.3, 0.1),
#                            verbose = FALSE)
sim2 <- readRDS("simdata2.rds")
simdata2 <- assays(sim2)$counts
assays(sim2)$logcounts <- log2(as.matrix(simdata2) + 1)
colnames(rowData(sim2))[1]<-"feature_symbol"
simdata2_ann <- colData(sim2)$Group
ari.lak.simdata2 <- c()
for(i in 1:100){
  simdata2_LAK <- LAK(simdata2, 5,s_error = 0.05, s_value = 17.95195)
  simdata2_LAK_ann <- simdata2_LAK[[1]]$Cs
  simdata2_LAK_ARI <- adjustedRandIndex(as.integer(as.factor(simdata2_ann)),as.vector(simdata2_LAK_ann))
  ari.lak.simdata2 <- c(ari.lak.simdata2,simdata2_LAK_ARI)
  cat(i)
}
ari.lak.simdata2
ari.sc3.simdata2 <- c()
for(i in 1:100){
  #simdata2_SCE <- readRDS("simdata2.rds")#14878*1886
  simdata2_sc3 <- sc3(sim2,5, gene_filter = T)
  simdata2_sc3_ann <- colData(simdata2_sc3)$sc3_5_clusters
  simdata2_sc3_ARI <- adjustedRandIndex(as.vector(simdata2_sc3_ann),as.integer(as.factor(simdata2_ann)))
  ari.sc3.simdata2 <- c(ari.sc3.simdata2, simdata2_sc3_ARI)
  cat(i)
  
}
ari.sc3.simdata2
ari_data2 <- data.frame(ari.lak.simdata2,ari.sc3.simdata2)
write.table(ari_data2,file="ari_data2.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=F)
#saveRDS(sim2,file="simdata2.rds")



#params.groups3 <- newSplatParams(batchCells = 2000, nGenes = 20000,dropout.type="experiment",
#                                 dropout.mid=1, dropout.shape = -1)
#sim3 <- splatSimulateGroups(params.groups3, group.prob = c(0.3,0.3,0.3, 0.1),
#                            verbose = FALSE)
sim3 <- readRDS("simdata3.rds")
simdata3 <- assays(sim3)$counts
assays(sim3)$logcounts <- log2(as.matrix(simdata3) + 1)
colnames(rowData(sim3))[1]<-"feature_symbol"
simdata3_ann <- colData(sim3)$Group
ari.lak.simdata3 <- c()
for(i in 1:100){
  simdata3_LAK <- LAK(simdata3, 5, s_error = 0.05, s_value = 13.99092 )
  simdata3_LAK_ann <- simdata3_LAK[[1]]$Cs
  simdata3_LAK_ARI <- adjustedRandIndex(as.integer(as.factor(simdata3_ann)),as.vector(simdata3_LAK_ann))
  ari.lak.simdata3 <- c(ari.lak.simdata3,simdata3_LAK_ARI)
  cat(i)
}
ari.lak.simdata3
ari.sc3.simdata3 <- c()
for(i in 1:100){
  #simdata3_SCE <- readRDS("simdata3.rds")#14878*1886
  simdata3_sc3 <- sc3(sim3,5, gene_filter = T)
  simdata3_sc3_ann <- colData(simdata3_sc3)$sc3_5_clusters
  simdata3_sc3_ARI <- adjustedRandIndex(as.vector(simdata3_sc3_ann),as.integer(as.factor(simdata3_ann)))
  ari.sc3.simdata3 <- c(ari.sc3.simdata3, simdata3_sc3_ARI)
  cat(i)
  
}
ari.sc3.simdata3
ari_data3 <- data.frame(ari.lak.simdata3,ari.sc3.simdata3)
write.table(ari_data3,file="ari_data3.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=F)
#saveRDS(sim3,file="simdata3.rds")



#params.groups4 <- newSplatParams(batchCells = 2000, nGenes = 20000,dropout.type="experiment",
#                                 dropout.mid=2, dropout.shape = -1)
#sim4 <- splatSimulateGroups(params.groups4, group.prob = c(0.3,0.3,0.3, 0.1),
#                            verbose = FALSE)
sim4 <- readRDS("simdata4.rds")
simdata4 <- assays(sim4)$counts
assays(sim4)$logcounts <- log2(as.matrix(simdata4) + 1)
colnames(rowData(sim4))[1]<-"feature_symbol"
simdata4_ann <- colData(sim4)$Group
ari.lak.simdata4 <- c()
for(i in 1:100){
  simdata4_LAK <- LAK(simdata4, 5,s_error = 0.05, s_value = 11.29992)
  simdata4_LAK_ann <- simdata4_LAK[[1]]$Cs
  simdata4_LAK_ARI <- adjustedRandIndex(as.integer(as.factor(simdata4_ann)),as.vector(simdata4_LAK_ann))
  ari.lak.simdata4 <- c(ari.lak.simdata4,simdata4_LAK_ARI)
  cat(i)
}
ari.lak.simdata4
ari.sc3.simdata4 <- c()
for(i in 1:100){
  #simdata2_SCE <- readRDS("simdata2.rds")#14878*1886
  simdata4_sc3 <- sc3(sim4,5, gene_filter = T)
  simdata4_sc3_ann <- colData(simdata4_sc3)$sc3_5_clusters
  simdata4_sc3_ARI <- adjustedRandIndex(as.vector(simdata4_sc3_ann),as.integer(as.factor(simdata4_ann)))
  ari.sc3.simdata4 <- c(ari.sc3.simdata4, simdata4_sc3_ARI)
  cat(i)
  
}
ari.sc3.simdata4
ari_data4 <- data.frame(ari.lak.simdata4,ari.sc3.simdata4)
write.table(ari_data4,file="ari_data4.csv",quote=FALSE,sep=",",col.names=TRUE,row.names=F)
#saveRDS(sim4,file="simdata4.rds")





dropout(simdata1)
dropout(simdata2)
dropout(simdata3)
dropout(simdata4)

