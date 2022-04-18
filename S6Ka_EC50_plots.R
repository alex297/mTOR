rm(list=ls())

# this code plots S6Ka and EC50 based on the output files from
# S6Ka_EC50_calculate.R
# and calculate_FC.R

s6k <- read.csv("S6Ka_max_EC50_all.csv", header = T,stringsAsFactors = F)
#s6k <- read.csv("S6Ka_max_EC50_all_AD.csv", header = T,stringsAsFactors = F)
s6k$ligand <- factor(s6k$ligand,levels = c("ins","pdgf","nrg"))
cells <- unique(s6k$cell)
N <- length(cells)
s6k$cell <- factor(s6k$cell,levels =s6k$cell[1:N])
cel <- c("EC art","EC cap","EC ven","SMC vas","SMC art","peric st","peric ECM", "astro hpc","astro ctx","microglia","OD","OPC","neuron")

library(ggplot2)
# S6Ka plot; Fig.3B
levels(s6k$cell) <- cel
ggplot(s6k,aes(x=cell, y = S6K_max,fill = ligand)) + geom_bar(stat = "identity",position = "stack",alpha = 0.6) +
  labs(x="",y= "S6Ka max, a.u.", size = (12))+
  theme_bw() +  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("green","royalblue","grey30"),labels = c("insulin", "pdgf","nrg"))+
  #theme_classic()+
  theme(axis.text.x = element_text(vjust=0.4),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),) +
  theme(legend.title = element_blank())+ labs(x = "") +
  theme(axis.line = element_line(colour = "grey30"),panel.border = element_blank()) #+
theme(legend.position = "none")
#jpeg("Fig3B.jpeg", width = 4, height =2.8, units = 'in', res = 600)
#dev.off()   

# EC50 plot; Fig.3C
ggplot(s6k) +  geom_line(aes(x = cell,y = T05,color = ligand,group=ligand),size=0.3) + 
  geom_point(aes(x = cell,y = T05,color = ligand,group=ligand),size=2)+
  labs(x="",y= "EC50, nM", size = (12))+ scale_x_discrete(labels=cel)+
  theme_bw() +  theme(axis.text.x=element_text(angle=90)) +#ylim(0.01, 10)+
  scale_color_manual(values=c("green","royalblue","grey30"),labels = c("insulin", "pdgf","nrg"))+
  #theme_classic()+
  theme(axis.text.x = element_text(vjust=0.4),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),) +
  theme(legend.title = element_blank())+ labs(x = "") +
  theme(axis.line = element_line(colour = "grey30"),panel.border = element_blank()) +
  scale_y_continuous(limits = c(0.01, 10),trans='log10')  #+theme(legend.position = "none")
#jpeg("Fig3C.jpeg", width = 4, height =2.8, units = 'in', res = 600)
#dev.off()   

# For Figs.3E,F (sensitivity to changes in PTEN, TSC) fold changes (FC) are calculated in calculate_FC.R
FC <- read.csv("fold_changes_0422.csv", header = T,stringsAsFactors = F)
cells <- FC$cell
cel <- c("EC art","EC cap","EC ven","SMC vas","SMC art","peric st","peric ecm", "astro hpc","astro ctx","microglia","OD","OPC","neuron")
FC$cell <- factor(FC$cell)
FC$cell <- ordered(FC$cell, levels = c("EC_artery","EC_capillary","EC_venous","SMC_vasc","SMC_artery","pericytes_st","pericytes_ECM",
                                       "astro_Hpc","astro_ctx","microglia","oligodegro","OPC","neuron"))
FC_ins <- FC[,c(1,5)] # for PTEN
FC_ins <- FC[,c(1,8)] # for TSC
colnames(FC_ins)[2] <- "FC"
FC_ins$ligand <- "insulin"
FC_pdgf <- FC[,c(1,6)] 
FC_pdgf <- FC[,c(1,9)] 
colnames(FC_pdgf)[2] <- "FC"
FC_pdgf$ligand <- "pdgf"
FC_nrg <- FC[,c(1,7)] 
FC_nrg <- FC[,c(1,10)] 
colnames(FC_nrg)[2] <- "FC"
FC_nrg$ligand <- "nrg"
FC_ <- rbind(FC_ins,FC_pdgf,FC_nrg)
library("ggplot2")
levels(FC_$cell) <- cel
FC_$ligand <- factor(FC_$ligand)
FC_$ligand <- ordered(FC_$ligand,levels = c("insulin","pdgf","nrg"))
ggplot(FC_,aes(x=cell, y = FC,fill = ligand)) + geom_bar(color="grey50",stat = "identity",position = "stack",alpha = 0.3) +
  labs(x="",y= "sensitivity to changes in PTEN", size = (12))+
  theme_bw() +  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("green","royalblue","grey30"),labels = c("insulin", "pdgf","nrg"))+
  #theme_classic()+
  theme(axis.text.x = element_text(vjust=0.4),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),) +
  theme(legend.title = element_blank())+ labs(x = "") + scale_y_continuous(breaks=seq(0,8,2)) +
  theme(axis.line = element_line(colour = "grey30"),panel.border = element_blank()) 

#jpeg("Fig3F.jpeg", width = 4, height =2.9, units = 'in', res = 600)
#jpeg("Fig3E.jpeg", width = 4, height =2.9, units = 'in', res = 600)
#dev.off()   

