rm(list=ls())

s6k_CT <- read.csv("S6Ka_max_EC50_all.csv", header = T,stringsAsFactors = F)
s6k_AD <- read.csv("S6Ka_max_EC50_all_AD.csv", header = T,stringsAsFactors = F)
S6K_AD_to_CT <- 100*(s6k_AD$S6K_max-s6k_CT$S6K_max)/(s6k_CT$S6K_max)
pot_AD_to_CT <- 100*(s6k_AD$EC50_reciprocal-s6k_CT$EC50_reciprocal)/(s6k_CT$EC50_reciprocal)
s6k_CT_ <- s6k_CT
s6k_CT <- s6k_CT[,c(1,2,4)]
s6k_CT$group <- "S6Ka"
colnames(s6k_CT)[2] <- "value"
s6k_CT_ <- s6k_CT_[,c(1,5,4)]
s6k_CT_$group <- "Potency"
colnames(s6k_CT_)[2] <- "value"
s6k_CT_$value <- s6k_CT_$value/130
s6k_CT__ <- s6k_CT_
s6k_CT_$cell <- paste0(s6k_CT_$cell,"_1")
s6k_CT__$value <- NA
s6k_CT__$group <- "empty"
s6k_CT__$cell <- paste0(s6k_CT__$cell,"_2")
s6k_ct <- rbind(s6k_CT,s6k_CT_,s6k_CT__)
s6k_ct$cell <- factor(s6k_ct$cell)
s6k_ct$cell <- ordered(s6k_ct$cell, levels = c("EC_artery","EC_artery_1","EC_artery_2","EC_capillary","EC_capillary_1","EC_capillary_2","EC_venous","EC_venous_1","EC_venous_2",
                                               "SMC_vasc","SMC_vasc_1","SMC_vasc_2","SMC_artery","SMC_artery_1","SMC_artery_2","pericytes_st","pericytes_st_1","pericytes_st_2",
                                               "pericytes_ECM","pericytes_ECM_1","pericytes_ECM_2","astro_Hpc","astro_Hpc_1","astro_Hpc_2","astro_ctx","astro_ctx_1","astro_ctx_2",
                                               "microglia","microglia_1","microglia_2","oligodegro","oligodegro_1","oligodegro_2","OPC","OPC_1","OPC_2","neuron","neuron_1","neuron_2"))

cel <- c("","EC art","","","EC cap","","","EC ven","","","SMC vas","","","SMC art","","","peric st","",
         "","peric ECM","", "","astro hpc","","","astro ctx","","","microglia","","","OD","","","OPC","","","neuron","")
library("ggplot2")
s6k_ct$ligand <- factor(s6k_ct$ligand)
s6k_ct$ligand <- ordered(s6k_ct$ligand,levels=c("ins","pdgf","nrg"))
s6k_ct$ID <- "bars"
s6k_line <- s6k_ct
s6k_line$ID <- "line_S6K_ratio"
rownames(s6k_line) <- NULL
s6k_line <- s6k_line[40:78,]
s6k_line$value <- S6K_AD_to_CT
s6k_line1 <- s6k_ct
s6k_line1$ID <- "line_Pot_ratio"
rownames(s6k_line1) <- NULL
s6k_line1 <- s6k_line1[40:78,]
s6k_line1$value <- pot_AD_to_CT
s6k_ct_ <- rbind(s6k_ct,s6k_line,s6k_line1)
coeff1 <- 700
coeff2 <- 50
lab2 <- expression(paste("S6Ka"["AD"]," / S6Ka"["CT"]))
lab3 <- expression(paste("Potency"["AD"]," / Potency"["CT"]))
rownames(s6k_ct_) <- NULL
s6k_ct_$value[127] <- 151 # change to 302 on final Fig, after axis break in Inskape
s6k_ct_$value[166] <- 125.5 # change to 251 on final Fig, after axis break in Inskape

# to determine scale for lines check values in s6k_line, s6k_line1
ggplot(s6k_ct_,aes(x=cell, y = value,fill = ligand)) + geom_bar(data=s6k_ct_[s6k_ct_$ID=="bars",],aes(alpha = as.factor(group)),stat = "identity",position = "stack") +
  labs(x="",y= "", size = (12))+guides(alpha = guide_legend(override.aes = list(size = 2)),fill = guide_legend(override.aes = list(size = 2))) +
  theme_bw() +  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=c("green","royalblue","grey30"),labels = c("insulin", "pdgf","nrg"))+
  scale_alpha_manual(values = c(0.5,1,0),labels = c("Potency/130","S6Ka", "")) + # guide='none'
   theme(axis.text.x = element_text(vjust=0.4),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),) +
  theme(legend.title = element_blank())+ labs(x = "") +
  theme(axis.line = element_line(colour = "grey30"),panel.border = element_blank()) +
  scale_x_discrete(labels=cel) +
  geom_line(data=s6k_ct_[s6k_ct_$ID=="line_S6K_ratio",],aes(x = as.numeric(cell), y = (value+coeff2)/coeff1, color = ligand),size=0.5) + 
  geom_point(data=s6k_ct_[s6k_ct_$ID=="line_S6K_ratio",],aes(x = as.numeric(cell), y = (value+coeff2)/coeff1, color = ligand),size=1.5) +
  scale_colour_manual(values=c("green","royalblue","grey30"),labels = c( "","",lab2)) +
  geom_line(data=s6k_ct_[s6k_ct_$ID=="line_Pot_ratio",],linetype = "dashed",aes(x = as.numeric(cell), y = (value+coeff2)/coeff1, color = ligand),size=0.5) + 
  geom_point(data=s6k_ct_[s6k_ct_$ID=="line_Pot_ratio",],shape=17,aes(x = as.numeric(cell), y = (value+coeff2)/coeff1, color = ligand),size=1.5) +
  scale_y_continuous(sec.axis = sec_axis(~.*coeff1-coeff2,name="change in AD, % (S6Ka, Potency)"),limits = c(0, 0.35))+labs(y= "S6Ka; Potency", size = (12))+
  scale_colour_manual(values=c("green","royalblue","grey30"),labels = c( "","AD change in S6Ka","AD change in Potency")) +
  geom_hline(yintercept=0.072, color ="red",size=0.5)
#svg("mTOR_summary.svg", width = 6, height =4)
#dev.off()   


