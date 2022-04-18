rm(list=ls())

# plots use S6Ka values for 1000 cells with random abundances
# S6Ka was calculated using random_cells.R, which output files
# mod_abundances_random_1221.csv; S6Ka_refmod_random_1221.csv

  s6k <- read.csv("S6Ka_refmod_random_1221.csv", header = T, stringsAsFactors = FALSE)
  ini <- read.csv("mod_abundances_random_1221.csv", header = T, stringsAsFactors = FALSE)
  IR <- as.numeric(unname(ini[1,-1]))
  IRS <- as.numeric(unname(ini[3,-1]))
  PTEN <- as.numeric(unname(ini[2,-1]))
  PKB <- as.numeric(unname(ini[4,-1]))
  TSC <- as.numeric(unname(ini[5,-1]))
  S6K <- as.numeric(unname(ini[8,-1]))
  PDGFR <- as.numeric(unname(ini[9,-1]))
  NRGR <- as.numeric(unname(ini[10,-1]))
  s6k_ins <- as.numeric(unname(s6k[1,-1]))
  s6k_pdgf <- as.numeric(unname(s6k[2,-1]))
  s6k_nrg <- as.numeric(unname(s6k[3,-1]))
  
brain_ini <- read.csv("snRNAseq_grouped_upquant_0422.csv", header = T, stringsAsFactors = FALSE)
ini_ <- brain_ini[c(8:13,4,2),-11] # initial conditions for "IR","IRS","PTEN","AKT","TSC","S6K","PDGFR", NRGR
rownames(ini_) <- NULL
cells <- colnames(ini_)[-1]
cells
cells <- c("EC","EC","EC","peric/SMC","peric/SMC","peric/SMC","peric/SMC","astrocytes","astrocytes","OD/OPC","OD/OPC","neuron")
s6k_ <- read.csv("S6Ka_max_EC50_all.csv", header = T,stringsAsFactors = F)
s6k_ <- s6k_[-c(10,23,36),]
IR_ <- as.numeric(unname(ini_[1,-1]))
IRS_ <- as.numeric(unname(ini_[2,-1]))
PTEN_ <- as.numeric(unname(ini_[3,-1]))
PKB_ <- as.numeric(unname(ini_[4,-1]))
TSC_ <- as.numeric(unname(ini_[5,-1]))
S6K_ <- as.numeric(unname(ini_[6,-1]))
PDGFR_ <- as.numeric(unname(ini_[7,-1]))
NRGR_ <- as.numeric(unname(ini_[8,-1]))
S6Ka_max_ins <- s6k_[1:12,2]
S6Ka_max_pdgf <- s6k_[13:24,2]
S6Ka_max_nrg <- s6k_[25:36,2]

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(length(unique(cells)))

library(ggplot2)
library(plyr)

# Fig. 2A
# random cells only; coloured by ligand

mult <- PKB*S6K/TSC/PTEN
R1 <- round(cor(log10(mult),log10(s6k_ins)),2) #0.87
R2 <- round(cor(log10(mult),log10(s6k_pdgf)),2) #0.94
R3 <- round(cor(log10(mult),log10(s6k_nrg*2)),2) #0.94
lab1 <- expression(paste("log"[10]," PKB*S6K/TSC/PTEN"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(rep(mult,3)),s6ka=log10(c(s6k_ins,s6k_pdgf,s6k_nrg*2)),
                  group=c(rep("insulin",1000),rep("PDGF",1000),rep("NRG",1000)), shap=c(rep(19,1000),rep(19,1000),rep(1,1000)),ID="random")
dat_all <- rbind(dat)
dat_all$group <- factor(dat_all$group,levels = c("insulin","PDGF","NRG"))
cols <- c("green","royalblue","grey30")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.title = element_blank(),legend.position = "top",legend.justification='left') +
  annotate(x=0.4, y=-0.06,label="R=0.87", geom="text", size=4,color="green",fontface="bold") +
  annotate(x=0.8, y=0.2,label="R=0.94", geom="text", size=4,color="royalblue",fontface="bold") +
  annotate(x=1.6, y=0.4,label="R=0.94", geom="text", size=4,color="grey30",fontface="bold") +
  geom_point(data=dat_all[dat_all$ID=="brain",],size=2) + 
  scale_color_manual(values = cols) 

#jpeg("Fig2A.jpeg", width = 3.5, height =4, units = 'in', res = 600)
#dev.off()


# Fig.2B
# expression of mTOR components abundances_0422_control.csv was downloaded from https://twc-stanford.shinyapps.io/human_bbb/
expr <- read.csv("abundances_0422_control.csv", header = T, stringsAsFactors = FALSE)
expr <- expr[c(13,22,27),]
N <- (dim(expr)[2]-1)/2
expr_max <- apply(expr[,2:(N+1)], 1, function(c) max(c))
expr[1,-1] <- expr[1,-1]/expr_max[1]
expr[2,-1] <- expr[2,-1]/expr_max[2]
expr[3,-1] <- expr[3,-1]/expr_max[3]
x = 1:N # x axis width

par(mar = c(7.1, 4.1, 4.1, 2.1))    
plot(x,expr[1,2:(N+1)], pch=16, lwd=2,cex=0.7,xaxt='n',bty="n", ylab="expression", xlab='',type = "b", lty = 1,ylim = c(0,1.3)) # , xlim=c(0, 16),ylim = c(0,2)
min <- as.numeric(expr[1,2:(N+1)]+expr[1,(N+2):(2*N+1)])
max <- as.numeric(expr[1,2:(N+1)]-expr[1,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3)
lines(x,expr[2,2:(N+1)],lwd=2,col="red",pch=16, cex=0.7, type = "b")
min <- as.numeric(expr[2,2:(N+1)]+expr[2,(N+2):(2*N+1)])
max <- as.numeric(expr[2,2:(N+1)]-expr[2,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="red")
lines(x,expr[3,2:(N+1)],lwd=2,col="blue",pch=16, cex=0.7, type = "b")
min <- as.numeric(expr[3,2:(N+1)]+expr[3,(N+2):(2*N+1)])
max <- as.numeric(expr[3,2:(N+1)]-expr[3,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="blue")
colo <- c("black","red","blue")
symb <- rep(16,dim(expr)[1])
legend("top",ncol=3, bty="n",legend=expr$gene, col=colo, pch=symb, cex=1)
axis(1, at = 1:N, labels = colnames(expr)[2:(N+1)], las=3,cex.axis=1,tck=-0.01)
#jpeg("Fig2B.jpeg", width = 5, height =4, units = 'in', res = 600)
#dev.off()



# insulin stimulation; random cells and brain cells
# Fig. 2C

#mult <- IR
#mult <- IRS
#mult <- PTEN
#mult <- PKB
#mult <- TSC
#mult <- S6K
mult <- PKB*S6K/TSC/PTEN
#mult <- PKB*S6K/PTEN
#mult <- PKB/PTEN
#mult <- IR*PKB*S6K/PTEN
#mult <- IRS*PKB*S6K/PTEN
#mult <- IR*PKB*S6K/TSC/PTEN
#mult <- IRS*PKB*S6K/TSC/PTEN
#mult <- IR*IRS*PKB*S6K/TSC/PTEN
#mult <- IR*PKB*S6K/TSC
#mult <- PKB*S6K/TSC
#mult_ <- IR_
#mult_ <- IRS_
#mult_ <- PTEN_
#mult_ <- PKB_
#mult_ <- TSC_
#mult_ <- S6K_
mult_ <- PKB_*S6K_/TSC_/PTEN_
#mult_ <- PKB_*S6K_/PTEN_
#mult_ <- PKB_/PTEN_
#mult_ <- IR_*PKB_*S6K_/PTEN_
#mult_ <- IRS_*PKB_*S6K_/PTEN_
#mult_ <- IR_*PKB_*S6K_/TSC_
#mult_ <- PKB_*S6K_/TSC_
#mult_ <- IR_*PKB_*S6K_/TSC_/PTEN_
#mult_ <- IRS_*PKB_*S6K_/TSC_/PTEN_
#mult_ <- IR_*IRS_*PKB_*S6K_/TSC_/PTEN_

R1 <- round(cor(log10(mult),log10(s6k_ins)),2) 
lab1 <- expression(paste("log"[10]," PKB*S6K/TSC/PTEN"))
lab2 <- expression(paste("log"[10]," S6Ka max"))

dat <- data.frame(x=log10(rep(mult,1)),s6ka=log10(c(s6k_ins)),group=c(rep("random",1000)), shap=c(rep(19,1000)),ID="random")
dat1 <- data.frame(x=log10(rep(mult_,1)),s6ka=log10(c(S6Ka_max_ins)),group=cells, shap=c(rep(19,length(cells))),ID="brain")
R1_ <- round(cor(dat1$x,dat1$s6ka),2) #
dat2 <- dat1
dat2$group <- "peric/SMC"
dat2$ID <- "brain1"
dat_all <- rbind(dat1,dat,dat2)
cols <- c(color_list,"grey")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = group, colour=ID, fill = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",],colour="grey70")+ #scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],aes(colour=ID),colour="grey70",method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ #xlim()+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=11),legend.text=element_text(size=11),
        axis.text.y = element_text(size=11),legend.title = element_blank()) +
  annotate(x=0., y=-0.0,label="R=0.87", geom="text", size=4,fontface="bold",color="grey70") +
  annotate(x=0.5, y=-0.7,label="R=0.92", geom="text", size=4,color="magenta",fontface="bold") +
  geom_smooth(data=dat_all[dat_all$ID=="brain1",],colour="magenta",method = "lm",alpha=0.2,show.legend = FALSE)+
  geom_point(data=dat_all[dat_all$ID=="brain",],colour="grey30",size=3) + labs(title="insulin stimulus") +
  scale_shape_manual(values = c(rep(22,length(cells)),21)) + scale_fill_manual(values = cols) 

#jpeg("Fig2C.jpeg", width = 4, height =4, units = 'in', res = 600)
#dev.off()


# NRG
# Fig. 2D

#mult <- NRGR
#mult <- PTEN
#mult <- PKB*S6K/TSC/PTEN
#mult <- NRGR*PKB*S6K/PTEN
#mult <- NRGR*S6K/TSC/PTEN
#mult <- NRGR*S6K/PTEN
#mult <- NRGR/PTEN
#mult <- NRGR*S6K
#mult <- NRGR*PKB*S6K/TSC
mult <- NRGR*PKB*S6K/TSC/PTEN
#mult_ <- NRGR_
#mult_ <- PTEN_
#mult_ <- PKB_*S6K_/TSC_/PTEN_
#mult_ <- NRGR_*PKB_*S6K_/PTEN_
#mult_ <- NRGR_*S6K_/TSC_/PTEN_
#mult_ <- NRGR_*S6K_/PTEN_
#mult_ <- NRGR_/PTEN_
#mult_ <- NRGR_*S6K_
mult_ <- NRGR_*PKB_*S6K_/TSC_/PTEN_
#mult_ <- NRGR_*PKB_*S6K_/TSC_

R1 <- round(cor(log10(mult),log10(s6k_nrg)),2)
lab1 <- expression(paste("log"[10]," NRGR*PKB*S6K/TSC/PTEN"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(rep(mult,1)),s6ka=log10(c(s6k_nrg)),group=c(rep("random",1000)), shap=c(rep(19,1000)),ID="random")
dat1 <- data.frame(x=log10(rep(mult_,1)),s6ka=log10(c(S6Ka_max_nrg)),group=cells, shap=c(rep(19,length(cells))),ID="brain")
R1_ <- round(cor(dat1$x,dat1$s6ka),2) #
dat2 <- dat1
dat2$group <- "peric/SMC"
dat2$ID <- "brain1"
dat_all <- rbind(dat1,dat,dat2)
cols <- c(color_list,"grey")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = group, colour=ID, fill = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",],colour="grey70")+ #scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],aes(colour=ID),colour="grey70",method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ #xlim()+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=11),legend.text=element_text(size=11),
        axis.text.y = element_text(size=11),legend.title = element_blank()) +
  annotate(x=1, y=-0.25,label="R=0.94", geom="text", size=4,fontface="bold",color="grey70") +
  annotate(x=0.5, y=-0.7,label="R=0.96", geom="text", size=4,color="magenta",fontface="bold") +
  geom_smooth(data=dat_all[dat_all$ID=="brain1",],colour="magenta",method = "lm",alpha=0.2,show.legend = FALSE)+
  geom_point(data=dat_all[dat_all$ID=="brain",],colour="grey30",size=3) + labs(title="NRG stimulus") +
  scale_fill_manual(values = cols) + scale_shape_manual(values = c(rep(22,length(cells)),21))

#jpeg("Fig3D.jpeg", width = 4, height =4, units = 'in', res = 600)
#dev.off()


# plot for IRS1 for insulin only (Fig. S4)
mult <- IRS
mult_ <- IRS_
R1 <- round(cor(log10(mult),log10(s6k_ins)),2) #0.12
lab1 <- expression(paste("log"[10]," IRS1"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(rep(mult,1)),s6ka=log10(c(s6k_ins)),group=c(rep("random",1000)), shap=c(rep(19,1000)),ID="random")
dat_all <- dat
cols <- "grey"
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = group, fill = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",],colour="grey70")+ #scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],colour="grey70",method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ #xlim()+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=11),legend.text=element_text(size=11),
        axis.text.y = element_text(size=11),legend.title = element_blank()) +
  annotate(x=0.4, y=0.07,label="R=0.12", geom="text", size=4,fontface="bold",color="grey70") +
  geom_point(data=dat_all[dat_all$ID=="brain",],colour="grey30",size=3) + labs(title="insulin stimulus") +
  scale_fill_manual(values = cols) + scale_shape_manual(values = c(rep(22,16),21))


# PDGF
# for chosen cells with non-zero PDGFR (PDGFR >= 1% of max) 
cells_chosen <- cells
cells <- cells_chosen[c(4:7,11,12)]
cells[5] <- "OPC"
PTEN_chosen <- PTEN_
PTEN_ <- PTEN_chosen[c(4:7,11,12)]
PKB_chosen <- PKB_
PKB_ <- PKB_chosen[c(4:7,11,12)]
TSC_chosen <- TSC_
TSC_ <- TSC_chosen[c(4:7,11,12)]
S6K_chosen <- S6K_
S6K_ <- S6K_chosen[c(4:7,11,12)]
PDGFR_chosen <- PDGFR_
PDGFR_ <- PDGFR_chosen[c(4:7,11,12)]
S6Ka_max_pdgf_chosen <- S6Ka_max_pdgf
S6Ka_max_pdgf <- S6Ka_max_pdgf_chosen[c(4:7,11,12)]
color_list_chosen <- color_list
color_list <- color_list_chosen[c(3:5)]

#mult <- PDGFR
#mult <- PTEN
#mult <- PKB
#mult <- TSC
#mult <- S6K
#mult <- PKB*S6K/TSC/PTEN
#mult <- PKB*S6K/TSC
#mult <- PKB*S6K
#mult <- S6K/TSC
#mult <- PKB/TSC
#mult <- PDGFR*PKB*S6K/TSC
#mult <- PDGFR*S6K/TSC/PTEN
#mult <- PDGFR*S6K/TSC
#mult <- PDGFR*S6K/PTEN
mult <- PDGFR*PKB*S6K/TSC/PTEN
#mult_ <- PDGFR_
#mult_ <- PTEN_
#mult_ <- PKB_
#mult_ <- TSC_
#mult_ <- PKB_*S6K_/TSC_/PTEN_
#mult_ <- PKB_*S6K_/TSC_
#mult_ <- PKB_*S6K_
#mult_ <- S6K_/TSC_
#mult_ <- PKB_/TSC_
#mult_ <- S6K_
#mult_ <- PDGFR_*PKB_*S6K_/TSC_
#mult_ <- PDGFR_*S6K_/TSC_/PTEN_
#mult_ <- PDGFR_*S6K_/TSC_
#mult_ <- PDGFR_*S6K_/PTEN_
mult_ <- PDGFR_*PKB_*S6K_/TSC_/PTEN_

R1 <- round(cor(log10(mult),log10(s6k_pdgf)),2)    
lab1 <- expression(paste("log"[10]," PDGFR*PKB*S6K/TSC/PTEN"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(rep(mult,1)),s6ka=log10(c(s6k_pdgf)),group=c(rep("random",1000)), shap=c(rep(19,1000)),ID="random")
dat1 <- data.frame(x=log10(rep(mult_,1)),s6ka=log10(c(S6Ka_max_pdgf)),group=cells, shap=c(rep(19,length(cells))),ID="brain")
R1_ <- round(cor(dat1$x,dat1$s6ka),2) #0.78/0.79/0.91
dat2 <- dat1
dat2$group <- "peric/SMC"
dat2$ID <- "brain1"
dat_all <- rbind(dat1,dat,dat2)
cols <- c(color_list,"grey")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = group,colour=ID, fill = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",],colour="grey70")+ #scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],aes(colour=ID),colour="grey70",method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ #xlim()+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=11),legend.text=element_text(size=11),
        axis.text.y = element_text(size=11),legend.title = element_blank()) +
  annotate(x=0.4, y=-0.04,label="R=0.84", geom="text", size=4,fontface="bold",color="grey70") +
  annotate(x=0.6, y=-0.6,label="R=0.85", geom="text", size=4,color="magenta",fontface="bold") +
  geom_smooth(data=dat_all[dat_all$ID=="brain1",],colour="magenta",method = "lm",alpha=0.3,show.legend = FALSE)+
  geom_point(data=dat_all[dat_all$ID=="brain",],colour="grey30",size=3) + labs(title="PDGF stimulus") +
  scale_fill_manual(values = cols) + scale_shape_manual(values = c(rep(22,length(cells)),21))
#jpeg("FigS4.jpeg", width = 4, height =4, units = 'in', res = 600)
#dev.off()


# For Fig. S4
# plot for receptors for all stimulus
mult1 <- IR
mult2 <- PDGFR
mult3 <- NRGR
R1 <- round(cor(log10(mult1),log10(s6k_ins)),2) #0.08
R2 <- round(cor(log10(mult2),log10(s6k_pdgf)),2) #-0.03
R3 <- round(cor(log10(mult3),log10(s6k_nrg*2)),2) #0.2
lab1 <- expression(paste("log"[10]," R (IR | PDGFR | NRGR)"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult1,mult2,mult3)),s6ka=log10(c(s6k_ins,s6k_pdgf,s6k_nrg*2)),
                  group=c(rep("insulin",1000),rep("PDGF",1000),rep("NRG",1000)), shap=c(rep(19,1000),rep(16,1000),rep(1,1000)),ID="random")
dat1 <- data.frame(x=log10(rep(mult_,3)),s6ka=log10(c(S6Ka_max_ins,S6Ka_max_pdgf,S6Ka_max_nrg*2)),
                   group=c(rep("insulin_brain",16),rep("PDGF_brain",16),rep("NRG_brain",16)), shap=c(rep(19,16),rep(19,16),rep(19,16)),ID="brain")
dat_all <- rbind(dat,dat1)
dat_all$group <- factor(dat_all$group,levels = c("insulin","PDGF","NRG","insulin_brain","PDGF_brain","NRG_brain"))
cols <- c("black","blue","brown","green","yellow","cyan")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.title = element_blank(),legend.position = "top",legend.justification='left') +
  annotate(x=0.4, y=0.1,label="R=0.08", geom="text", size=4,fontface="bold") +
  annotate(x=0.7, y=-0.3,label="R=-0.03", geom="text", size=4,color="blue",fontface="bold") +
  annotate(x=0.94, y=-0.7,label="R=0.2", geom="text", size=4,color="brown",fontface="bold") +
  scale_color_manual(values = cols) 

# Fig.S4 common variables in 3 pathway (insulin, PDGF, NRG)
#mult <- PTEN
#mult <- PKB
mult <- TSC
#mult <- S6K
#mult_ <- PTEN_
#mult_ <- PKB_
mult_ <- TSC_
#mult_ <- S6K_
R1 <- round(cor(log10(mult),log10(s6k_ins)),2)
R2 <- round(cor(log10(mult),log10(s6k_pdgf)),2)
R3 <- round(cor(log10(mult),log10(s6k_nrg*2)),2)
lab1 <- expression(paste("log"[10]," TSC"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(rep(mult,3)),s6ka=log10(c(s6k_ins,s6k_pdgf,s6k_nrg*2)),
                  group=c(rep("insulin",1000),rep("PDGF",1000),rep("NRG",1000)), shap=c(rep(19,1000),rep(16,1000),rep(1,1000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels = c("insulin","PDGF","NRG"))
cols <- c("black","blue","brown")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.title = element_blank(),legend.position = "top",legend.justification='left') +
  annotate(x=0.3, y=0.1,label="R=-0.33", geom="text", size=4,fontface="bold") +
  annotate(x=0.6, y=-0.1,label="R=-0.43", geom="text", size=4,color="blue",fontface="bold") +
  annotate(x=0.8, y=-0.4,label="R=-0.47", geom="text", size=4,color="brown",fontface="bold") +
  scale_color_manual(values = cols) 

# plot for combination of several functions of abundances, Fig. S4
# for insulin 
mult1 <- IR*IRS/PTEN/TSC*PKB*S6K
mult2 <- IR*IRS/PTEN/TSC*S6K
mult3 <- IRS/PTEN/TSC*PKB*S6K
mult4 <- IR/PTEN/TSC*PKB*S6K
R1 <- round(cor(log10(mult1),log10(s6k_ins)),2) 
R2 <- round(cor(log10(mult2),log10(s6k_ins)),2) 
R3 <- round(cor(log10(mult3),log10(s6k_ins)),2)
R4 <- round(cor(log10(mult4),log10(s6k_ins)),2)
lab1 <- expression(paste("log"[10]," function of abundances"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult1,mult2,mult3,mult4)),s6ka=log10(rep(s6k_ins,4)),
                  group=c(rep("IR*IRS/PTEN/TSC*PKB*S6K",1000),rep("IR*IRS/PTEN/TSC*S6K",1000),rep("IRS/PTEN/TSC*PKB*S6K",1000),rep("IR/PTEN/TSC*PKB*S6K",1000)), shap=c(rep(16,4000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels = c("IR*IRS/PTEN/TSC*PKB*S6K","IR*IRS/PTEN/TSC*S6K","IRS/PTEN/TSC*PKB*S6K","IR/PTEN/TSC*PKB*S6K"))
cols <- c("black","grey50","blue","red")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.position = c(1.0,0.2),legend.justification='right',legend.background=element_blank()) +
  annotate(x=0.6, y=-0.04,label="R=0.8", geom="text", size=4,fontface="bold") + labs(color = "                  insulin stimulus")+
  annotate(x=2.4, y=-0.5,label="R=0.71", geom="text", size=4,fontface="bold",color="grey50") + labs(fill = "Dose (mg)")+
  annotate(x=2.4, y=-0.3,label="R=0.83", geom="text", size=4,fontface="bold",color="blue") + labs(fill = "Dose (mg)")+
  annotate(x=1.7, y=0.0,label="R=0.83", geom="text", size=4,fontface="bold",color="red") + 
  scale_color_manual(values = cols) 

mult1 <- PKB*S6K/TSC
mult2 <- PKB*S6K/PTEN
mult3 <- PKB*S6K/PTEN/TSC
R1 <- round(cor(log10(mult1),log10(s6k_ins)),2) 
R2 <- round(cor(log10(mult2),log10(s6k_ins)),2) 
R3 <- round(cor(log10(mult3),log10(s6k_ins)),2)
lab1 <- expression(paste("log"[10]," function of abundances"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult1,mult2,mult3)),s6ka=log10(rep(s6k_ins,3)),
                  group=c(rep("PKB*S6K/TSC",1000),rep("PKB*S6K/PTEN",1000),rep("PKB*S6K/PTEN/TSC",1000)), shap=c(rep(16,3000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels = c("PKB*S6K/TSC","PKB*S6K/PTEN","PKB*S6K/PTEN/TSC"))
cols <- c("black","blue","red")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.position = c(1.0,0.2),legend.justification='right',legend.background=element_blank()) +
  annotate(x=0.5, y=-0.,label="R=0.85", geom="text", size=4,fontface="bold") + labs(color = "        insulin stimulus")+
  annotate(x=1.7, y=0.0,label="R=0.81", geom="text", size=4,fontface="bold",color="blue") + 
  annotate(x=2., y=-0.3,label="R=0.87", geom="text", size=4,fontface="bold",color="red") + 
  scale_color_manual(values = cols) 

mult1 <- IRS*S6K/PTEN/TSC
mult2 <- S6K/PTEN/TSC
mult3 <- S6K/PTEN
mult4 <- S6K/TSC
R1 <- round(cor(log10(mult1),log10(s6k_ins)),2)
R2 <- round(cor(log10(mult2),log10(s6k_ins)),2)
R3 <- round(cor(log10(mult3),log10(s6k_ins)),2)
R4 <- round(cor(log10(mult4),log10(s6k_ins)),2)
lab1 <- expression(paste("log"[10]," function of abundances"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult1,mult2,mult3,mult4)),s6ka=log10(rep(s6k_ins,4)),
                  group=c(rep("IRS*S6K/PTEN/TSC",1000),rep("S6K/PTEN/TSC",1000),rep("S6K/PTEN",1000),rep("S6K/TSC",1000)), shap=c(rep(16,1000),rep(16,1000),rep(16,1000),rep(16,1000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels = c("IRS*S6K/PTEN/TSC","S6K/PTEN/TSC","S6K/PTEN","S6K/TSC"))
cols <- c("black","blue","grey50","red")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+#geom_point() +  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.position = c(1.02,0.2),legend.justification='right',legend.background=element_blank()) +
  annotate(x=2, y=-0.5,label="R=0.74", geom="text", size=4,fontface="bold") + labs(color = "insulin stimulus")+
  annotate(x=2, y=-0.01,label="R=0.8", geom="text", size=4,fontface="bold",color="blue") + labs(fill = "Dose (mg)")+
  annotate(x=1.5, y=-0.14,label="R=0.74", geom="text", size=4,fontface="bold",color="grey50") + 
  annotate(x=1., y=-0.01,label="R=0.8", geom="text", size=4,fontface="bold",color="red") + 
  scale_color_manual(values = cols) 

# for PDGF Fig. S4
mult1 <- PDGFR*PKB*S6K/PTEN/TSC
mult2 <- PKB*S6K/PTEN/TSC
mult3 <- PDGFR*S6K/PTEN/TSC
mult4 <- PDGFR*PKB*S6K/PTEN
R1 <- round(cor(log10(mult1),log10(s6k_pdgf)),2) 
R2 <- round(cor(log10(mult2),log10(s6k_pdgf)),2) 
R3 <- round(cor(log10(mult3),log10(s6k_pdgf)),2) 
R4 <- round(cor(log10(mult4),log10(s6k_pdgf)),2) 
lab1 <- expression(paste("log"[10]," multiplication of abundances"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult3,mult4,mult1,mult2)),s6ka=log10(rep(s6k_pdgf,4)),
                  group=c(rep("PDGFR*S6K/PTEN/TSC",1000),rep("PDGFR*PKB*S6K/PTEN",1000),rep("PDGFR*PKB*S6K/PTEN/TSC",1000),rep("PKB*S6K/PTEN/TSC",1000)), shap=c(rep(16,1000),rep(16,1000),rep(16,1000),rep(16,1000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels = c("PDGFR*PKB*S6K/PTEN/TSC","PKB*S6K/PTEN/TSC","PDGFR*S6K/PTEN/TSC","PDGFR*PKB*S6K/PTEN"))
cols <- c("blue","red","black","grey50")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.position = c(1.02,0.26),legend.justification='right',legend.background=element_blank()) +
  annotate(x=2.2, y=-0.3,label="R=0.71", geom="text", size=4,fontface="bold") + labs(color = "             PDGF stimulus")+
  annotate(x=2.3, y=-0.01,label="R=0.84", geom="text", size=4,fontface="bold",color="blue") + labs(fill = "Dose (mg)")+
  annotate(x=2.4, y=-0.15,label="R=0.7", geom="text", size=4,fontface="bold",color="grey50") + 
  annotate(x=1.7, y=0.03,label="R=0.94", geom="text", size=4,fontface="bold",color="red") + 
  geom_point(data=dat_all[dat_all$ID=="brain",],size=2) + 
  scale_color_manual(values = cols) 

mult1 <- PDGFR*PKB*S6K/TSC
mult2 <- PKB/PTEN/TSC
mult3 <- PKB*S6K/PTEN
mult4 <- PKB*S6K/TSC
R1 <- round(cor(log10(mult1),log10(s6k_pdgf)),2) 
R2 <- round(cor(log10(mult2),log10(s6k_pdgf)),2) 
R3 <- round(cor(log10(mult3),log10(s6k_pdgf)),2)
R4 <- round(cor(log10(mult4),log10(s6k_pdgf)),2)
lab1 <- expression(paste("log"[10]," multiplication of abundances"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult1,mult2,mult3,mult4)),s6ka=log10(rep(s6k_pdgf,4)),
                  group=c(rep("PDGFR*PKB*S6K/TSC",1000),rep("PKB/PTEN/TSC",1000),rep("PKB*S6K/PTEN",1000),rep("PKB*S6K/TSC",1000)), shap=c(rep(16,1000),rep(16,1000),rep(16,1000),rep(16,1000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels = c("PDGFR*PKB*S6K/TSC","PKB/PTEN/TSC","PKB*S6K/PTEN","PKB*S6K/TSC"))
cols <- c("black","grey50","blue","red")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.position = c(1.02,0.26),legend.justification='right',legend.background=element_blank()) +
  annotate(x=2.3, y=-0.01,label="R=0.81", geom="text", size=4,fontface="bold") + labs(color = "             PDGF stimulus")+
  annotate(x=2.2, y=-0.3,label="R=0.66", geom="text", size=4,fontface="bold",color="grey50") +
  annotate(x=2.4, y=-0.15,label="R=0.82", geom="text", size=4,fontface="bold",color="blue") + 
  annotate(x=1.4, y=0.09,label="R=0.92", geom="text", size=4,fontface="bold",color="red") + 
  geom_point(data=dat_all[dat_all$ID=="brain",],size=2) + 
  scale_color_manual(values = cols) 

# for NRG, Fig. S4
mult1 <- NRGR*PKB*S6K/PTEN/TSC
mult2 <- PKB*S6K/PTEN/TSC
mult3 <- PKB*S6K/PTEN
mult4 <- PKB*S6K/TSC
R1 <- round(cor(log10(mult1),log10(s6k_nrg)),2) 
R2 <- round(cor(log10(mult2),log10(s6k_nrg)),2) 
R3 <- round(cor(log10(mult3),log10(s6k_nrg)),2)
R4 <- round(cor(log10(mult4),log10(s6k_nrg)),2) 
lab1 <- expression(paste("log"[10]," function of abundances"))
lab2 <- expression(paste("log"[10]," S6Ka max"))
dat <- data.frame(x=log10(c(mult1,mult2,mult3,mult4)),s6ka=log10(rep(s6k_nrg,4)),
                  group=c(rep("PKB*S6K/PTEN",1000),rep("PKB*S6K/TSC",1000),rep("NRGR*PKB*S6K/PTEN/TSC",1000),rep("PKB*S6K/PTEN/TSC",1000)), shap=c(rep(16,1000),rep(16,1000),rep(16,1000),rep(16,1000)),ID="random")
dat_all <- dat
dat_all$group <- factor(dat_all$group,levels =  c("PKB*S6K/PTEN","PKB*S6K/TSC","NRGR*PKB*S6K/PTEN/TSC","PKB*S6K/PTEN/TSC"))
cols <- c("black","grey50","blue","red")
names(cols) <- levels(dat_all$group)
ggplot(dat_all,aes(x = x,y=s6ka,shape = shap, color = group)) + 
  geom_point(data=dat_all[dat_all$ID=="random",])+ scale_shape_identity()+  
  geom_smooth(data=dat_all[dat_all$ID=="random",],method = "lm",show.legend = FALSE)+
  theme_classic() + labs(x=lab1, y = lab2)+ theme(
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=11),legend.text=element_text(size=11),
    axis.text.y = element_text(size=11),legend.position = c(1.02,0.2),legend.justification='right',legend.background=element_blank()) +
  annotate(x=2.3, y=-0.01,label="R=0.8", geom="text", size=4,fontface="bold") + labs(color = "                  NRG stimulus")+
  annotate(x=2.1, y=-0.4,label="R=0.91", geom="text", size=4,fontface="bold",color="grey50") +
  annotate(x=1.5, y=-0,label="R=0.94", geom="text", size=4,fontface="bold",color="blue") + 
  annotate(x=1.3, y=0.2,label="R=0.94", geom="text", size=4,fontface="bold",color="red") + 
  scale_color_manual(values = cols) 

#jpeg("FigS4.jpeg", width = 4, height =4, units = 'in', res = 600)
#dev.off()

