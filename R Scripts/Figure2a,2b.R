library(ggplot2)
library(mclust)
library("cowplot")
library(splatter)
setwd("C:/Users/12441/Desktop/LAK CODE/")
source("LAK.R")
load("Figure2a.RData")



#save(ARIS_result,ari.linnorm.biase,ari.pcaReduce.biase,ari.sc3.biase,ari.sincera.biase,
#     ari.linnorm.treu,ari.pcaReduce.treu,ari.sc3.treu,ari.sincera.treu,
#     ari.linnorm.yan,ari.pcaReduce.yan,ari.sc3.yan,ari.sincera.yan,
#     ari.linnorm.goolam,ari.pcaReduce.goolam,ari.sc3.goolam,ari.sincera.goolam,
#     ari.camp1.tshc,ari.pcaReduce.camp1,ari.sc3.camp1,ari.sincera.camp1,
#     ari.tshc.li,ari.pcaReduce.li,ari.sc3.li,ari.sincera.li,
#     ari.tshc.patel,ari.pcaReduce.patel,ari.sc3.patel,ari.sincera.patel,
#     ari.tshc.zeisel,ari.pcaReduce.zeisel,ari.sc3.zeisel,ari.sincera.zeisel,
#     file="Figure2a.RData")
LAK_patel <- ARIS_result$Patel
LAK_zeisel <- ARIS_result$Zeisel
LAK_goolam <- ARIS_result$Goolam
LAK_biase <- ARIS_result$Biase
LAK_t <- ARIS_result$Treutlein
LAK_y <- ARIS_result$Yan
LAK_camp1 <- ARIS_result$Camp1
LAK_li <- ARIS_result$Li


datagroup<-factor(c(rep("Biase",401), rep("Treulein",401),rep("Yan",401)))
datagroup1<-factor(c(rep("Goolam",401),rep("Camp1",401),rep("Li",401)))

hari<-c(LAK_biase,ari.linnorm.biase,ari.pcaReduce.biase,ari.sc3.biase,ari.sincera.biase,
        LAK_t,ari.linnorm.treu,ari.pcaReduce.treu,ari.sc3.treu,ari.sincera.treu,
        LAK_y,ari.linnorm.yan,ari.pcaReduce.yan,ari.sc3.yan,ari.sincera.yan) 
hari1<-c(LAK_goolam,ari.linnorm.goolam,ari.pcaReduce.goolam,ari.sc3.goolam,ari.sincera.goolam,
         LAK_camp1,ari.camp1.tshc,ari.pcaReduce.camp1,ari.sc3.camp1,ari.sincera.camp1,
        LAK_li,ari.tshc.li,ari.pcaReduce.li,ari.sc3.li,ari.sincera.li) 
methodgroup<-factor(rep(c(rep("LAK",100),rep("Tsne+hclust",100),
                          rep("pcaReduce",100),rep("SC3",100),"Sincera"),3))    
hARI<-data.frame(hari,methodgroup,datagroup)
hARI_1<-data.frame(hari1,methodgroup,datagroup1)
hARI1 <- within(hARI,{
  methodgroup <- factor(methodgroup,levels = c("LAK", "Tsne+hclust","pcaReduce","SC3","Sincera"))
})
hARI1_1 <- within(hARI_1,{
  methodgroup <- factor(methodgroup,levels = c("LAK", "Tsne+hclust","pcaReduce","SC3","Sincera"))
})

fig2a_1 <- ggplot(data=hARI1,aes(x=methodgroup,y=hari))+
  geom_bar(fun.y=median,stat="summary",aes(fill=methodgroup))+
  geom_point(position="jitter",size=1)+
  labs(title="Comparison",x="Method",y="ARI")+
  guides(fill=guide_legend(title="Method"))+
  theme(plot.title = element_text(size=20,hjust=0.5),
        axis.title.x=element_text(size=14,hjust=0.5),
        axis.title.y=element_text(size=14,hjust=0.5),
        # axis.text.x=element_text(size=6,angle=30,vjust=0.5),
        axis.text.x=element_blank(),
        panel.background = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"))+
  facet_grid(.~datagroup)
fig2a_2 <- ggplot(data=hARI1_1,aes(x=methodgroup,y=hari1))+
  geom_bar(fun.y=median,stat="summary",aes(fill=methodgroup))+
  geom_point(position="jitter",size=1)+
  labs(x="Method",y="ARI")+
  guides(fill=guide_legend(title="Method"))+
  theme(plot.title = element_text(size=20,hjust=0.5),
        axis.title.x=element_text(size=14,hjust=0.5),
        axis.title.y=element_text(size=14,hjust=0.5),
        # axis.text.x=element_text(size=6,angle=30,vjust=0.5),
        axis.text.x=element_blank(),
        panel.background = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"))+
  facet_grid(.~datagroup1)

fig2b<-plot_grid(fig2a_1,fig2a_2,nrow=2)
fig2b

dropout <- function(x){
  length(which(x==0))/(ncol(x)*nrow(x))
}
patel_dropout <- dropout(Patel)
yan_dropout <- dropout(yan)
goolam_dropout <- dropout(Goolam)
biase_dropout <- dropout(Biase)
treutlein_dropout <- dropout(Treutlein)
camp1_dropout <- dropout(camp1)
li_dropout <- dropout(li)
zeisel_dropout <- dropout(zeisel)
dropout_rate<- data.frame(dropout_rate=c(patel_dropout,yan_dropout,goolam_dropout,biase_dropout,
                                        treutlein_dropout,camp1_dropout,li_dropout,zeisel_dropout),
                         row.names=c("patel","yan","goolam","biase","treutlein","camp1","li","zeisel"))


stable_biase <- c(max(table(LAK_biase))/100,max(table(ari.linnorm.biase))/100,
                max(table(ari.pcaReduce.biase))/100,max(table(ari.sc3.biase)/100))
stable_treu <- c(max(table(LAK_t))/100,max(table(ari.linnorm.treu))/100,
                max(table(ari.pcaReduce.treu))/100,max(table(ari.sc3.treu)/100))
stable_yan <- c(max(table(LAK_y))/100,max(table(ari.linnorm.yan))/100,
                max(table(ari.pcaReduce.yan))/100,max(table(ari.sc3.yan)/100))
stable_camp1 <- c(max(table(LAK_camp1))/100,max(table(ari.tshc.camp1))/100,
                max(table(ari.pcaReduce.camp1))/100,max(table(ari.sc3.camp1)/100))
stable_goolam <- c(max(table(LAK_goolam))/100,max(table(ari.linnorm.goolam))/100,
                max(table(ari.pcaReduce.goolam))/100,max(table(ari.sc3.goolam)/100))
stable_li <- c(max(table(LAK_li))/100,max(table(ari.tshc.li))/100,
                  max(table(ari.pcaReduce.li))/100,max(table(ari.sc3.li)/100))
stable_patel <- c(max(table(LAK_patel))/100,max(table(ari.tshc.patel))/100,
               max(table(ari.pcaReduce.patel))/100,max(table(ari.sc3.patel)/100))
stable_zeisel <- c(max(table(LAK_zeisel))/100,max(table(ari.tshc.zeisel))/100,
                  max(table(ari.pcaReduce.zeisel))/100,max(table(ari.sc3.zeisel)/100))

stable_data <- rbind(biase=stable_biase,treu=stable_treu,yan=stable_yan,
                          camp1=stable_camp1,goolam=stable_goolam,li=stable_li,patel = stable_patel,
                     zeisel=stable_zeisel)
colnames(stable_data) <- c("LAK","Tsne_hclust","pcaReduce","SC3")
#save(stable_data,dropout_rate,
#     file = "Figure2b.RData")
load("Figure2b.RData")
