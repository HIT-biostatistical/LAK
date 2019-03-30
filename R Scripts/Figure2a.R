library(ggplot2)
library(mclust)

setwd("~/LAK/")
source("LAK.R")
load("RData/Figure2a.RData")

LAK_goolam <- ARIS_result$Goolam
LAK_biase <- ARIS_result$Biase
LAK_t <- ARIS_result$Treutlein
LAK_y <- ARIS_result$Yan

datagroup<-factor(c(rep("goolam",401),rep("biase",401),
                    rep("treulein",401),rep("yan",401)))

hari<-c(LAK_goolam,ari.linnorm.goolam,ari.pcaReduce.goolam,ari.sc3.goolam,ari.sincera.goolam,
        LAK_biase,ari.linnorm.biase,ari.pcaReduce.biase,ari.sc3.biase,ari.sincera.biase,
        LAK_t,ari.linnorm.treu,ari.pcaReduce.treu,ari.sc3.treu,ari.sincera.treu,
        LAK_y,ari.linnorm.yan,ari.pcaReduce.yan,ari.sc3.yan,ari.sincera.yan) 
methodgroup<-factor(rep(c(rep("LAK",100),rep("Tsne+hclust",100),
                          rep("pcaReduce",100),rep("SC3",100),"Sincera"),4))    
hARI<-data.frame(hari,methodgroup,datagroup)
hARI1 <- within(hARI,{
  methodgroup <- factor(methodgroup,levels = c("LAK", "Tsne+hclust","pcaReduce","SC3","Sincera"))
})                

ggplot(data=hARI1,aes(x=methodgroup,y=hari))+
  geom_bar(fun.y=median,stat="summary",aes(fill=methodgroup))+
  geom_point(position="jitter",size=1)+
  labs(title="comparison",x="Method",y="ARI")+
  guides(fill=guide_legend(title="Method"))+
  theme(plot.title = element_text(size=20,hjust=0.5),
        axis.title.x=element_text(size=14,hjust=0.5),
        axis.title.y=element_text(size=14,hjust=0.5),
        # axis.text.x=element_text(size=6,angle=30,vjust=0.5),
        axis.text.x=element_blank(),
        panel.background = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"))+
  facet_grid(.~datagroup)







