setwd("~/LAK")
source("LAK.R")
library(pheatmap)



Zeisel<-read.table("Single Cell Data/Zeisel.txt", header = T)
gene_names<-as.vector(Zeisel$cell_id)
Zeisel <- as.matrix(Zeisel[,-1])

load("RData/Zeisel_LAK.RData")

gene_names <- gene_names[which(rowSums(Zeisel>0.5)>2)]
Zeisel<-Zeisel[rowSums(Zeisel>0.5)>2,]
rownames(Zeisel) <- gene_names

Zeisel_normed <- Linnorm(Zeisel)

Zeisel_LAK_ann <- Zeisel_LAK[[1]]$Cs

Zeisel_with_LAK_genes <- Zeisel[Zeisel_LAK[[1]]$ws > 1e-5,]

Zeisel_with_LAK_genes_normed <- Zeisel_normed[Zeisel_LAK[[1]]$ws > 1e-5,]

Zeisel_marker<-c('Aif1','Aldoc','Acta2','Cldn5','Thy1','Spink8','Mbp','Gad1','Tbr1')



#t-test
T_matrix<- function(df,df_label){      
  out<-data.frame(row.names = rownames(df))
  for( i in min(df_label):max(df_label) ){
    for( j in 1:dim(df)[1]){
      a<-as.vector(t(df[j,df_label==i])) 
      b<-as.vector(t(df[j,df_label!=i]))
      if(sd(a)==0 && sd(b)==0){
      }else{
        #t检???
        out[j,i]<-t.test(a,b, alternative="greater", paired = FALSE, var.equal = FALSE)$p.value
      }
    }
    print(i)
  }
  colnames(out)<-paste("cluster",1:max(df_label))
  return(out)
} #t test to get differential expressed genes

d_expr<-function(out_matrix,k=1,p=0.01,gene_names){   
  pvalue<-out_matrix[,k]
  times<-0
  gene_list<-c()
  for( i in 1:dim(out_matrix)[1]){
    if(!is.na(pvalue[i])&pvalue[i]<p){
      times<-times+1
      gene_list[times]<-gene_names[i]
    }
  }
  return(gene_list)
}

T_diff<- T_matrix(Zeisel_with_LAK_genes, Zeisel_LAK_ann)

markers_in_diff_genes <- list()
all_diff_genes <- list()
topn <- 100
M <- T_diff

dgene <- c()
pgene <- c()

for(i in c(1,3:ncol(M))){
  diff_genes <- rownames(M[order(M[,i]),])[1:topn]
  Z_n <- Zeisel_with_LAK_genes_normed[diff_genes,]
  
  cur<-Z_n[,Zeisel_LAK_ann==i]
  other<-Z_n[,Zeisel_LAK_ann!=i]
  
  other_median<-apply(other, 1, median)
  other_mean<-apply(other,1,mean)
  cur_median<-apply(cur,1,median)
  cur_mean <- apply(cur, 1, mean)
  cur_max <- apply(cur,1,max)
  
  cur<-cur[10 * other_median< cur_median & 3* other_mean < cur_mean,]
  cat(i,'\t',rownames(cur)[1:5])
  cat('\n')
  dgene <- c(dgene, rownames(cur)[1:5])
  pgene<-c(pgene, rownames(cur[order(rowMeans(cur),decreasing = T),])[1:5])
}

for (i in 36:39){
  pgene[i]<-pgene[i+1]
}

pgene[10]<- "Tbr1"
pgene[15]<- "Thy1"
pgene[20]<- "Spink8"
pgene[30]<- "Cldn5"
pgene[35] <- "Aldoc"
pgene[40] <- "Aif1"
pgene[32] <- "Aqp4"


ht_df <- Zeisel_with_LAK_genes_normed[pgene,Zeisel_LAK_ann!=2]
ann <- Zeisel_LAK_ann[-which(Zeisel_LAK_ann==2)]

annotation_col<-data.frame(Cluster = factor(as.character(ann)),
                           row.names=colnames(ht_df))
##Classification according to our clustering results
ht_df=data.frame(t(ht_df),ann)
ht_df=ht_df[order(ht_df[,"ann"]),]
ht_df=t(ht_df)
ht_df=ht_df[-nrow(ht_df),]

###
counti <- table(ann)
cur <- 0
gapscol <- c()
for(i in c(1:7)){
  cur <- cur + counti[i]
  gapscol <- c(gapscol, cur)
}
gapscol
###

ann_colors=list(Cluster = c("1"="#00FFFF","3"="#33FF33",
                            "4"="#993300","5"="#000000","6"="#006600",
                            "7"="#FF0000","8"="#999999","9"="#0000FF"))

p <- pheatmap(ht_df,cluster_rows = FALSE,cluster_cols=FALSE,
              annotation_col = annotation_col, 
              annotation_colors = ann_colors,
              show_rownames = F,
              gaps_row = seq(8)*5,
              gaps_col = gapscol,
              fontsize = 8,
              show_colnames =F)
