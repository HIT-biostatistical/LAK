library(ggplot2)

#set path
setwd("~/LAK")
load("RData/Figure4.RData")
cbbPalette <- c("#33FF33","#FF0000","#0000FF","#000000","#999999","#FFFF00",
                "#00FFFF", "#FF00FF", "#993300", "#006699", "#FF9966", "#006600",
                "#AAAAAA")
num_cells <- c(430, 90, 124, 49, 80)
num_genes <- c(5948,20214, 41480,25737, 23271)
authors  <- c("Patel", "Yan", "Goolam", "Biase", "Treutlein")


ARI <- c()
cell_multi <- c()
nz_gene_proportion <- c()
author <- c()
s <- seq(1.2,30,0.2)

for( i in 1:5){
  ARI <- c(ARI, as.vector(aris[[i]]))
  cell_multi <- c(cell_multi, as.vector(non_zero_weights[[i]])/num_cells[i])
  nz_gene_proportion <- c(nz_gene_proportion, as.vector(non_zero_weights[[i]])/num_genes[i])
  
  author <- c(author , rep(authors[i],145))
}

figure3 <- data.frame(ARI, m= cell_multi, nz_gene_proportion, author, s)

a <- ggplot(figure3, aes(m,ARI)) +
  geom_point(aes(color = author))+
  scale_colour_manual(values=cbbPalette)+ theme_bw()+
  #scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+
  theme(title=element_text(size=15),axis.title=element_text(size=15),
        legend.text = element_text(size=15),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
a

b <- ggplot(figure3, aes(s, m)) + 
  geom_point(aes(color = author))+
  scale_colour_manual(values=cbbPalette)+ theme_bw()+
  #scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+
  theme(title=element_text(size=15),axis.title=element_text(size=15),
        legend.text = element_text(size=15),
        panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
b
