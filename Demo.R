library(mclust)

setwd("~/LAK")
source("LAK.R")

Biase <-  readRDS("Single Cell Data/biase.rds")
m <- assays(Biase)[[1]][,-(50:56)]
LAK_ann <- LAK(m,3)

Biase_ann <- colData(Biase)$cell_type1[-(50:56)]
Biase_ann_numeric <- c()
id <- names(table(Biase_ann))
for (i  in 1:length(Biase_ann)){
  for (j in 1:length(id)){
    if(Biase_ann[i] == id[j]){
      Biase_ann_numeric <- c(Biase_ann_numeric,j)
      break
    }
  }
}
adjustedRandIndex(LAK_ann[[1]]$Cs, Biase_ann_numeric)
