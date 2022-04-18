rm(list=ls())
#BiocManager::install("")
# this code plots S6Ka and 1/EC50 based on the output files from
# S6Ka_EC50_calculate.R

library(ggplot2)
# Fig.4
s6k_CT <- read.csv("S6Ka_max_EC50_all.csv", header = T,stringsAsFactors = F)
s6k_CT$group <- "CT"
s6k_AD <- read.csv("S6Ka_max_EC50_all_AD.csv", header = T,stringsAsFactors = F)
s6k_AD$group <- "AD"
s6k <- rbind(s6k_CT,s6k_AD)
rownames(s6k) <- NULL
s6k$group <- factor(s6k$group, levels= c("CT","AD"))
s6k$ligand <- factor(s6k$ligand,levels = c("ins","pdgf","nrg"))
s6k$cell <- factor(s6k$cell,levels =s6k$cell[1:13])
cel <- c("EC art","EC cap","EC ven","SMC vas","SMC art","peric st","peric ECM", "astro hpc","astro ctx","microglia","OD","OPC","neuron")

# S6Ka, Fig.4A
lab1 <- "S6Ka max, a.u."
ggplot(s6k,aes(x = cell, y = S6K_max)) + 
  geom_line(aes(linetype = group, color = ligand,group=interaction(group, ligand))) +
  geom_point(aes(color = ligand,group=ligand),size=1)+
  theme_bw() +  theme(axis.text.x=element_text(angle=90)) +
  scale_colour_manual(values=c("green","royalblue1","grey30"),labels = c("insulin","pdgf","nrg"))+
  theme(legend.title = element_blank())+ labs(x = "")+ labs(y = lab1)+ theme(axis.text.x = element_text(vjust=0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  scale_x_discrete(labels=cel)
#jpeg("Fig4A.jpeg", width = 4, height =2.8, units = 'in', res = 600)
#dev.off()

# potency, Fig.4B
lab1 <- "1/EC50, 1/nM"
ggplot(s6k,aes(x = cell, y = EC50_reciprocal)) + 
  geom_line(aes(linetype = group, color = ligand,group=interaction(group, ligand))) +
  geom_point(aes(color = ligand,group=ligand),size=1)+
  theme_bw() +  theme(axis.text.x=element_text(angle=90)) +
  scale_colour_manual(values=c("green","royalblue1","grey30"),labels = c("insulin","pdgf","nrg"))+
  theme(legend.title = element_blank())+ labs(x = "")+ labs(y = lab1)+ theme(axis.text.x = element_text(vjust=0.4)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey"))+
  scale_x_discrete(labels=cel) #+
  geom_errorbar(aes(ymin=recipr_from, ymax=recipr_to,linetype = group, color = ligand,group=interaction(group, ligand)), width=.3)+ # ,position=position_dodge(.9)
#jpeg("Fig4B.jpeg", width = 4, height =2.8, units = 'in', res = 600)
#dev.off()



