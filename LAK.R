if (!require(SingleCellExperiment)) { 
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
}

if (!require(Linnorm)) { # for normalization
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Linnorm", version = "3.8")
  library(Linnorm)
}

if (!require(cluster)) {
  install.packages("cluster")
  library(cluster)
}

if (!require(sparcl)) {
  install.packages("sparcl")
  library(sparcl)
}

if (!require(Rtsne)) {
  install.packages("Rtsne")
  library(Rtsne)
}
if (!require(clues)) {
  install.packages("clues")
  library(clues)
}

if (!require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require(pheatmap)) {
  install.packages("pheatmap")
  library(pheatmap)
}

if (!require(MASS)) {
  install.packages("MASS")
  library(MASS)
}

#data preprocessing 
#' @param data_id (integer) refers to the id of single cell data
#' @returns 
#' the express matrix of input data, the rows correspond to genes, and coloumns correspond to cells
#' whether the data has been normed by the original author
#' the annotation provided by the original author
get_sc_data <- function(data_id){
  
  if( data_id == 1){
    SCEset <- readRDS("Single Cell Data/patel.rds")
    m <- data.frame(assays(SCEset)[[1]])
    rownames(m)=rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    normed <- T
  }
  
  if(data_id==2){
    SCEset <- readRDS("Single Cell Data/yan.rds")
    m <- assays(SCEset)[[1]]
    rownames(m)=rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    normed <- F
  }
  
  if(data_id==3){
    SCEset <- readRDS("Single Cell Data/goolam.rds")
    m <- data.frame(assays(SCEset)[[1]])
    rownames(m)=rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    normed <- F
  }
  
  if(data_id==4){
    SCEset <- readRDS("Single Cell Data/biase.rds")
    m <- data.frame(assays(SCEset)[[1]][,-(50:56)])
    rownames(m)=rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1[-(50:56)]
    normed <- F
  }
  
  if(data_id==5){
    SCEset <- readRDS("Single Cell Data/treutlein.rds")
    m <- assays(SCEset)[[1]]
    rownames(m)=rowData(SCEset)$feature_symbol
    ann <- colData(SCEset)$cell_type1
    normed <- F
  }
  
  result <- list(m, normed, ann)
  return(result)
}



get_k_in_clusGap <- function(x, method="firstSEmax", SE.factor = 1)
{
  method <- match.arg(method, choices = eval(formals(maxSE)$method))
  stopifnot((K <- nrow(T <- x$Tab)) >= 1, SE.factor >= 0)
  cat("Clustering Gap statistic [\"clusGap\"] from call:\n", deparse(x$call),
      sprintf("\nB=%d simulated reference sets, k = 1..%d; spaceH0=\"%s\"\n",
              x$B, K, x$spaceH0), sep="")
  nc <- maxSE(f = T[,"gap"], SE.f = T[,"SE.sim"],
              method=method, SE.factor=SE.factor)
  cat(sprintf(" --> Number of clusters (method '%s'%s): %d\n",
              method, if(grepl("SE", method))
                sprintf(", SE.factor=%g",SE.factor) else "", nc))
  invisible(x)
  return(nc)
}



#LAK takes gap statistics as the method to determine the cluster number k
k_detection <- function(exprs_matrix, method = "kmeans", max_K = 20 , B = 100){
  if(method == "kmeans"){
    gap_clust <- clusGap(t(exprs_matrix), kmeans, K.max = max_K, B = B, verbose = interactive())
  }else if(method == "hclust"){
    
    hclusCut <- function(x, k, d.meth = "euclidean", ...)
      list(cluster = cutree(hclust(dist(x, method=d.meth), ...), k=k))
    
    gap_clust <- clusGap(t(exprs_matrix), hclusCut, K.max = max_K, B = B)
  }else
    stop("clustering method should be 'kmeans' or 'hclust'")
  k <- get_k_in_clusGap(gap_clust)
  return(k)
}


parameter_selection <- function(exprs_matrix, num_cluster, cell_times = 2, s_error = 0.1,
                                thre = 0.1, s_start = 5, s_end = 105){
  
  num_cells <- ncol(exprs_matrix)
  iter <- 0
  
  if(num_cells > 50 & num_cells < 10000){
    
    s_mid <- (s_start+s_end)/2
    Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
    num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
    c_t <- num_nz_gene/num_cells
    cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
    
    
    found <- F
    while (s_start <= s_end) {
      
      iter <- iter + 1
      
      if(abs(c_t-cell_times)< thre)
      {
        found <- T
        break()
      }
      
      else if (c_t < cell_times)
      {
        s_start <- s_mid + s_error
        s_mid <- (s_start+s_end)/2
        Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
        num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
        c_t <- num_nz_gene/num_cells
        cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
      }
      
      else
      {
        s_end <- s_mid - s_error
        s_mid <- (s_start+s_end)/2
        Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
        num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
        c_t <- num_nz_gene/num_cells
        cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
      }
    }
    
    if(found == F){
      stop("check parameter 's_start' and 's_end' or use smaller 's_error'")
    }else{
      
      cat("binary search done\n")
      cat("s = ", s_mid, "\n")
      return(Sparse_label)
    }
    
  }else if(num_cells < 50){
    
    s_mid <- (s_start+s_end)/2
    Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
    num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
    cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
    
    
    found <- F
    while (s_start <= s_end) {
      
      iter <- iter + 1
      
      if(abs(num_nz_gene-100)< 5)
      {
        found <- T
        break()
      }
      
      else if (num_nz_gene < 100)
      {
        s_start <- s_mid + s_error
        s_mid <- (s_start+s_end)/2
        Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
        num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
        c_t <- num_nz_gene/num_cells
        cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
      }
      
      else
      {
        s_end <- s_mid - s_error
        s_mid <- (s_start+s_end)/2
        Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
        num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
        c_t <- num_nz_gene/num_cells
        cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
      }
    }
    
    if(found == F){
      stop("check parameter 's_start' and 's_end' or use smaller 's_error'")
    }else{
      cat("binary search done\n")
      cat("s = ", s_mid, "\n")
      return(Sparse_label)
    }
    
  }else{
    
    s_mid <- (s_start+s_end)/2
    Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
    num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
    cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
    
    
    found <- F
    while (s_start <= s_end) {
      
      iter <- iter + 1
      
      if(abs(num_nz_gene-10000)< 500)
      {
        found <- T
        break()
      }
      
      else if (num_nz_gene < 10000)
      {
        s_start <- s_mid + s_error
        s_mid <- (s_start+s_end)/2
        Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
        num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
        c_t <- num_nz_gene/num_cells
        cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
      }
      
      else
      {
        s_end <- s_mid - s_error
        s_mid <- (s_start+s_end)/2
        Sparse_label <- KMeansSparseCluster(t(exprs_matrix), num_cluster, wbounds = s_mid)
        num_nz_gene <- length(which(Sparse_label[[1]]$ws>1e-10))
        c_t <- num_nz_gene/num_cells
        cat("number of non zero genes = ", num_nz_gene, "  s = ", s_mid, "\n")
      }
    }
    
    if(found == F){
      stop("check parameter 's_start' and 's_end' or use smaller 's_error'")
    }else{
      cat("binary search done\n")
      cat("s = ", s_mid, "\n")
      return(Sparse_label)
    }
  }
}


#data preprocessing 
#' @param exprs_matrix (integer) refers to the id of single cell data
#' @param num_cluster (integer) -1 means unknown cluster number and procees Gapclus
#' @param cell_times (integer) control how many genes involved in clustering
#' @param thre (numeric) threshold of stopping the binary search
#' @param s_start (numeric) the start of searching range
#' @param s_end (numeric) the start of searching range
#' @param normed (logical) indicate whether the data has been normlized, if not, do Linnorm
#' @param s_value (numeric) -1 means selecting the SparseKmeans parameter s automatically 
#' @param B (integer) iteration times in Gapclus, can result huge time consumption if very large
#' @returns LAK_result
#' the express matrix of input data, the rows correspond to genes, and coloumns correspond to cells
#' whether the data has been normed by the original author
#' the annotation provided by the original author
LAK <- function(exprs_matrix, num_cluster = -1, cell_times = 2,  s_error = 0.1,
                thre = 0.1,  s_start = 5, s_end = 100, normed = F, s_value = -1 , 
                B = 10)
{
   
  if(normed == F){
    
    cat("normlization.....")
    pre_filterd_matrix<-exprs_matrix[rowSums(exprs_matrix>0.5)>2,]
    normed_matrix <- Linnorm(pre_filterd_matrix)
    cat("done" , "\n")
    
  }else normed_matrix <- exprs_matrix
  
  if(num_cluster == -1){
    cat("detection of cluster number k.....\n")
    num_cluster <- k_detection(normed_matrix, B = B)
    cat("done", '\t', "k = ", num_cluster, '\n')
  }
  
  
  if(s_value == -1){ 
    cat("parameter selection....",'\n')
    LAK_result <- parameter_selection(exprs_matrix = normed_matrix, 
                                         num_cluster = num_cluster, 
                                         cell_times,  s_error,
                                         thre, s_start, s_end)
    cat("parameter selection done" , "\n")
  }else{
    LAK_result <- KMeansSparseCluster(t(normed_matrix), num_cluster, wbounds = s_value)
  }
  
  return(LAK_result)
}


#T-test by each row of the data matrix
#' @params
#' exprs_matrix refers to normalized data matrix
#' LAK_ann refers to the annotion of our cluster by our cluster method
#' @return P-values for each row
T_matrix <- function(exprs_matrix,LAK_ann){
  out <- data.frame(row.names = rownames(df))
  for( i in min(LAK_ann):max(LAK_ann) ){ #label is our cluster label by LAK
    for( j in 1:dim(exprs_matrix)[1]){
      a <- as.vector(t(exprs_matrix[j,LAK_ann == i])) #extract column of label i in the j row 
      b <- as.vector(t(exprs_matrix[j,LAK_ann != i])) #extract column of label not i in the j row
      if(sd(a) == 0 && sd(b) == 0){
      }else{
        out[j,i] <- t.test(a,b, alternative="greater", paired = FALSE, var.equal = FALSE)$p.value #W-test for each row
      }
    }
    print(i)
  }
  colnames(out) <- paste("cluster",1:max(LAK_ann))
  return(out)
}
#W-test by each row of the data matrix
#' @params
#' exprs_matrix refers to normalized data matrix
#' LAK_ann refers to the annotion of our cluster by our cluster method
#' @return P-values for each row
W_matrix <- function(exprs_matrix,df_ann){ 
  out <- data.frame(row.names = rownames(df)) 
  for( i in min(LAK_ann):max(LAK_ann) ){  #label is our cluster label by LAK
    for( j in 1:dim(exprs_matrix)[1]){
      a <- as.vector(t(exprs_matrix[j,LAK_ann == i])) #extract column of label i in the j row 
      b <- as.vector(t(exprs_matrix[j,LAK_ann != i])) #extract column of label not i in the j row
      if(sd(a) == 0 && sd(b) == 0){
      }else{
        out[j,i] <- wilcox.test(a,b, alternative="greater")$p.value #W-test for each row
      }
    }
  }
  colnames(out) <- paste("cluster",1:max(LAK_ann))
  return(out)
} 
#select the genetic names under P<0.001
#' @params
#' out_matrix refers to the p-values by T-test or W-test
#' k=1 refers to extract the first column of out_matrix
#' p=0.001 refers to p-value equal 0.001(generally as a threshold for hypothesis testing)
#' gene_names refers to all gene names of single cell data
#' @return a list of selected gene names
d_expr <- function(out_matrix,k = 1,p = 0.001,gene_names){
  pvalue <- out_matrix[,k]
  times <- 0
  gene_list <- c()
  for( i in 1:dim(out_matrix)[1]){
    if(!is.na(pvalue[i]) & pvalue[i] < p){
      times <- times+1
      gene_list[times] <- gene_names[i]
    }
  }
  return(gene_list)
}
#select the genetic names
#' @params
#' diff_pvalue_matrix refers to the p-values by T-test or W-test
#' topn=numsDiffGene(inter) in order to selecting for numsDiffGene(inter) different expression genes each type
#' p=1e-3 refers to p-value equal 1e-3(generally as a threshold for hypothesis testing)
#' test refers to T-test or W-test
#' @return differential expression gene names
get_diff_exprs_gene_names <- function(diff_pvalue_matrix,topn = numsDiffGene,p_value = 1e-3,test){
  pgene <- c()
  if(test == T){    #control the data selecting conditions by T-test
  max_med <- c(-0.4,0.5,0.5,0.5,0.5)
  max_mea <- c(-0.4,0.5,0.5,0.5,0.5)
  min_med <- c(3,1.5,1,3,0.5)
  }
  if(test == F){    #control the data selecting conditions by W-test
  max_med <- c(-0.4,0.5,0.5,0.5,0.4)
  max_mea <- c(-0.4,0.5,0.5,0.5,0.5)
  min_med <- c(3,1.5,1,3,0.5)
  }
  max_med <- max_med[id]
  max_mea <- max_mea[id]
  min_med <- min_med[id]
  for (i in 1:k){
    df <- normed_matrix
    igene <- d_expr(diff_pvalue_matrix,i,p_value,gene_names = rownames(exprs_matrix)) #select the genetic names under P<0.01
    df <- df[igene,]
    cur <- df[,LAK_ann == i]  #extrect column of label i
    other <- df[,LAK_ann != i] #other column
    other_median <- apply(other, 1, median) 
    other_mean <- apply(other,1,mean)
    cur_median <- apply(cur,1,median)
    #selecting the differentially expressed genes satisfying these three conditions
    cur <- cur[other_median < max_med & other_mean < max_mea & cur_median > min_med,] 
    #integrate the selected differentially expressed genes
    pgene <- c(pgene, rownames(cur[order(rowSums(cur),decreasing = T),])[1:topn])
  }
  return(pgene)
}
#ordering the differential expression gene
#' @params
#' gene refers to the differentially expressed gene names of selection
#' clusterOrder refers to specified sequence
#' numsDiffGene in order to selecting for the quantity of different expression genes each type
#' @return the differential expression gene names of ordering
get_reOrderGene <- function(gene,clusterOrder,numsDiffGene){
  reOrderGene <- c()
  for (i in clusterOrder){
    reOrderGene <- c(reOrderGene,gene[((i-1)*numsDiffGene+1):(i*numsDiffGene)])
  }
  return(reOrderGene)
}
#selecting the differential expression gene and ordering the gene names
#' @params
#' exprs_matrix refers to single cell data matrix after preprocessing
#' k refers to the cluster number by our method
#' normed_matrix refers to normalized data matrix
#' LAK_ann refers to the annotion of our cluster by our cluster method
#' test refers to T-test or W-test
#' numsDiffGene=10 in order to selecting for 10 different expression genes each type
#' @return the differential expression gene names of ordering
diff_gene_names_reorder <- function(exprs_matrix,k,normed_matrix,LAK_ann,test,numsDiffGene=10){
  
  test <- test
  if(test == T){
  diff_gene_pvalue <- T_matrix(normed_matrix, LAK_ann) #T-test
  }
   if(test == F){
  diff_gene_pvalue <- W_matrix(normed_matrix, LAK_ann) #W-test
  }
  #selecting the differential expression gene names
  get_pgene <- get_diff_exprs_gene_names(diff_gene_pvalue,topn = numsDiffGene,p_value =1e-2,test)
  #specify the cluster order
  patelClusterOrder <- c(1,2,3,4,5)
  yanClusterOrder <- c(1,2,3,4,5,6,7)
  goolamClusterOrder <- c(1,2,3,4,5,6)
  biaseClusterOrder <- c(1,2,3)
  TreutleinClusterOrder <- c(1,2,3,4,5,6)
  ClusterOrder <- list(patelClusterOrder,yanClusterOrder,goolamClusterOrder,
                     biaseClusterOrder,TreutleinClusterOrder) #integrate the cluster orders
  reOrderGene <- get_reOrderGene(get_pgene,ClusterOrder[[id]],numsDiffGene=numsDiffGene)
  return(reOrderGene)
}
#ploting heatmap according to the differentially expressed genes
#' @params
#' diff_expr_genenames refers to the differential expression gene names
#' k refers to the cluster number by our method
#' normed_matrix refers to normalized data matrix
#' LAK_ann refers to the annotion of our cluster by our cluster method
#' numsDiffGene=10 refers to every 10 genes a gap by ploting heatmap 
#' @return the result of pheatmap
plot_heatmap <- function(diff_expr_genenames,id,k,normed_matrix,LAK_ann,numsDiffGene){
  k <- k
  normed_matrix <- normed_matrix
  LAK_ann <- LAK_ann
  diff_gene_normed_matrix <- normed_matrix[diff_expr_genenames,] #extracting the data of all the differentially expressed genes
  annotation_col <- data.frame(Cluster = factor(as.character(LAK_ann)),row.names=colnames(diff_gene_normed_matrix)) #annotation col in ploting heatmap
  rownames(annotation_col) <- colnames(diff_gene_normed_matrix)
  #the following four steps are for classification according to our clustering results
  diff_gene_normed_matrix1 <- data.frame(t(diff_gene_normed_matrix),LAK_ann) 
  diff_gene_normed_matrix2 <- diff_gene_normed_matrix1[order(diff_gene_normed_matrix1[,"LAK_ann"]),] #ordering cells according to LAK annotation
  colorder_diff_gene_normed_matrix <- t(diff_gene_normed_matrix2)[-nrow(t(diff_gene_normed_matrix2)),]  #Delete the row of the label
  ann_colors <- list(Cluster = c("1" = "#00FFFF","2" = "#FF6699","3" = "#33FF33",
                              "4" = "#993300","5" = "#000000","6" = "#006600",
                              "7" = "#FF0000","8" = "#999999","9" = "#0000FF",
                              "10" = "#006699","11" = "#FF00FF","12" = "#FFFF00")[1:k])
  #table(LAK_ann) #it can view all kinds of sample size in order to constructing the gaps_col
  patel_gaps_col <- c(120,199,287,365)
  yan_gaps_col <- c(24,28,44,58,60,85)
  goolam_gaps_col <- c(46,62,74,90,100)
  biase_gaps_col <- c(9,29)
  treutlein_gaps_col <- c(41,54,58,62,70)
  gapscol <- list(patel_gaps_col,yan_gaps_col,goolam_gaps_col,
                biase_gaps_col,treutlein_gaps_col) #integrate them in irder to extracting it easily at function pheatmap
  heatmap <- pheatmap(colorder_diff_gene_normed_matrix,cluster_rows = FALSE,cluster_cols = FALSE,annotation_col = annotation_col, 
          annotation_colors = ann_colors,
          gaps_row = seq(1,k-1)*numsDiffGene,
          gaps_col = gapscol[[id]],
          show_colnames = F)
  return(heatmap)
}


