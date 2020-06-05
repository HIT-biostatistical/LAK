
library(splatter)
source("LAK.R")
library(SC3)
library(mclust)

dropout <- function(x){
  length(which(x==0))/(ncol(x)*nrow(x))
}
sim1 <- readRDS("simdata1_1000.rds")
simdata1 <- assays(sim1)$counts
dropout(simdata1)


sim2 <- readRDS("simdata2_1000.rds")
simdata2 <- assays(sim2)$counts
dropout(simdata2)

sim3 <- readRDS("simdata3_1000.rds")
simdata3 <- assays(sim3)$counts
dropout(simdata3)

sim4 <- readRDS("simdata4_1000.rds")
simdata4 <- assays(sim4)$counts
dropout(simdata4)




setwd("~")
ari.simdata.1000 <- read.csv("ari_simdata_1000.csv", header = T)
ari.simdata.2000 <- read.csv("ari_simdata_2000.csv", header = T)
stable_data1_1000 <- c(max(table(ari.simdata.1000$ari.lak.simdata1))/100, 
                      max(table(ari.simdata.1000$ari.sc3.simdata1))/100)
stable_data2_1000 <- c(max(table(ari.simdata.1000$ari.lak.simdata2))/100, 
                       max(table(ari.simdata.1000$ari.sc3.simdata2))/100)
stable_data3_1000 <- c(max(table(ari.simdata.1000$ari.lak.simdata3))/100, 
                       max(table(ari.simdata.1000$ari.sc3.simdata3))/100)
stable_data4_1000 <- c(max(table(ari.simdata.1000$ari.lak.simdata4))/100, 
                       max(table(ari.simdata.1000$ari.sc3.simdata4))/100)

stable_data_1000 <- rbind(stable_data1_1000, stable_data2_1000, stable_data3_1000, stable_data4_1000)

stable_data1_2000 <- c(max(table(ari.simdata.2000$ari.lak.simdata1))/100, 
                       max(table(ari.simdata.2000$ari.sc3.simdata1))/100)
stable_data2_2000 <- c(max(table(ari.simdata.2000$ari.lak.simdata2))/100, 
                       max(table(ari.simdata.2000$ari.sc3.simdata2))/100)
stable_data3_2000 <- c(max(table(ari.simdata.2000$ari.lak.simdata3))/100, 
                       max(table(ari.simdata.2000$ari.sc3.simdata3))/100)
stable_data4_2000 <- c(max(table(ari.simdata.2000$ari.lak.simdata4))/100, 
                       max(table(ari.simdata.2000$ari.sc3.simdata4))/100)

stable_data_2000 <- rbind(stable_data1_2000, stable_data2_2000, stable_data3_2000, stable_data4_2000)

drop_sim_1000 <- c(0.5916628, 0.7312638, 0.8066918, 0.8776254)
drop_sim_2000 <- c(0.5932867, 0.7324426,0.8073914, 0.8780138)


plot_sim1000 <- data.frame(Data = rep(c("Simulated data 1", "Simulated data 2", "Simulated data 3",
                                        "Simulated data 4"), 2),
                           Methods = c(rep("LAK", 4), rep("SC3", 4)), Zero_ratios = drop_sim_1000,
                           Stable_ratios = c(stable_data_1000[,1], stable_data_1000[,2]))
plot_sim2000 <- data.frame(Data = rep(c("Simulated data 1", "Simulated data 2", "Simulated data 3",
                                        "Simulated data 4"), 2),
                           Methods = c(rep("LAK", 4), rep("SC3", 4)), Zero_ratios = drop_sim_2000,
                           Stable_ratios = c(stable_data_2000[,1], stable_data_2000[,2]))
p1 <- ggplot(plot_sim1000, aes(Data, Methods))+
  geom_point(aes(size = Zero_ratios,color = Stable_ratios))+# 
  theme_classic()+
  scale_color_gradient(low="green",high = "red")+
  labs(size="Zero ratios of datasets",color="Cluster stability", x="Simulated data of 1000 cells",y="Methods",title="" 
  )  +#+  theme(axis.text.x = element_text(angle = 45, hjust = 1.1, vjust = 1.1))
  guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2)
  )
p2 <- ggplot(plot_sim2000, aes(Data, Methods))+
  geom_point(aes(color = Stable_ratios,size = Zero_ratios))+# 
  theme_classic()+
  scale_color_gradient(low="green",high = "red")+
  labs(size="Zero ratios of datasets",color="Cluster stability",x="Simulated data of 2000 cells",y="Methods",title="" 
  )+
  guides(
    colour = guide_colourbar(order = 1),
    size = guide_legend(order = 2)
  )
library("cowplot")
plot_grid(p1,p2, ncol=2,labels = c("b","c"))



