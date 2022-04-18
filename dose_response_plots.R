rm(list=ls())

# this code plots S6Ka dependence on ligand concentrations, based on output files from ode runs

#These files should be loaded from subdirectory /outputs while running this code, for example by setwd() to this subdirectory

# The output files were generated using
# mTOR_output_ligand_range_ins.R
# mTOR_output_ligand_range_PDGF.R
# mTOR_output_ligand_range_NRG.R

insul <- 10
ins <- 10
while (ins > 0.01) {
  ins <- ins/1.04
  insul <- c(insul,ins)
}

# plot dose dependence for PDGF, Fig. 3A

ins_ref <- read.csv("pdgf_ref_EC_artery.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_EC_artery.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_EC_artery.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(PDGF=insul,EC_artery=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_EC_capillary.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_EC_capillary.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_EC_capillary.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,EC_capillary=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_EC_venous.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_EC_venous.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_EC_venous.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,EC_venous=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_SMC_vasc.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_SMC_vasc.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_SMC_vasc.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,SMC_vasc=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_SMC_artery.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_SMC_artery.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_SMC_artery.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,SMC_artery=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_pericytes_st.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_pericytes_st.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_pericytes_st.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,pericytes_st=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_pericytes_ECM.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_pericytes_ECM.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_pericytes_ECM.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,pericytes_ECM=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_astro_Hpc.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_astro_Hpc.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_astro_Hpc.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,astro_Hpc=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_astro_ctx.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_astro_ctx.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_astro_ctx.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,astro_ctx=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_microglia.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_microglia.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_microglia.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,microglia=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_oligodegro.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_oligodegro.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_oligodegro.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,oligodegro=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_OPC.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_OPC.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_OPC.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,OPC=as.numeric(t(ins_ref[21,-1])))
ins_ref <- read.csv("pdgf_ref_neuron.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("nrg_ref_neuron.csv", header = T, stringsAsFactors = FALSE)
#ins_ref <- read.csv("ins_ref_neuron.csv", header = T, stringsAsFactors = FALSE)
S6K_pdgf <- data.frame(S6K_pdgf,neuron=as.numeric(t(ins_ref[21,-1])))
L <- dim(S6K_pdgf)[1]
background <- S6K_pdgf[L,-1]
for (i in 1:(L-1)) {
  background <- rbind(background,S6K_pdgf[L,-1])
}
S6K_pdgf_ <- S6K_pdgf[,-1]-background
S6K_pdgf[,-1] <- S6K_pdgf_

#par(mfrow=c(1,3)) 
plot(insul,S6K_pdgf$EC_artery, bty="n",pch=16,lwd = 2,cex.lab=1.4,cex.axis=1.3,main = "",lty="solid",log = 'x',col="darkgreen",ylim=c(0, 0.15),xaxt="n",ylab='S6Ka, a.u.', xlab='PDGF, nM',type = "l") 
#plot(insul,S6K_pdgf$EC_artery, bty="n",pch=16,lwd = 2,cex.lab=1.4,cex.axis=1.3,main = "",lty="solid",log = 'x',col="darkgreen",ylim=c(0, 0.06),xaxt="n",ylab='S6Ka, a.u.', xlab='NRG, nM',type = "l") 
#plot(insul,S6K_pdgf$EC_artery, bty="n",pch=16,lwd = 2,cex.lab=1.4,cex.axis=1.3,main = "",lty="solid",log = 'x',col="darkgreen",ylim=c(0, 0.15),xaxt="n",ylab='S6Ka, a.u.', xlab='insulin, nM',type = "l") 
lines(insul,S6K_pdgf$EC_capillary, pch=16,lwd = 2,cex=0.1,lty="solid",col="green",xaxt="n",type = "l")
lines(insul,S6K_pdgf$EC_venous, pch=16,lwd = 2,cex=0.1,lty="solid",col="green3",xaxt="n",type = "l")
lines(insul,S6K_pdgf$SMC_vasc, pch=16,lwd = 2,cex=0.1,lty="solid",col="deepskyblue3",xaxt="n",type = "l")
lines(insul,S6K_pdgf$SMC_artery, pch=16,lwd = 2,cex=0.1,lty="solid",col="deepskyblue4",xaxt="n",type = "l")
lines(insul,S6K_pdgf$pericytes_st, pch=16,lwd = 2,cex=0.1,lty="solid",col="royalblue1",xaxt="n",type = "l")
lines(insul,S6K_pdgf$pericytes_ECM, pch=16,lwd = 2,cex=0.1,lty="solid",col="royalblue2",xaxt="n",type = "l")
lines(insul,S6K_pdgf$astro_Hpc, pch=16,lwd = 2,cex=0.1,lty="solid",col="orange",xaxt="n",type = "l")
lines(insul,S6K_pdgf$astro_ctx, pch=16,lwd = 2,cex=0.1,lty="solid",col="orange3",xaxt="n",type = "l")
lines(insul,S6K_pdgf$microglia, pch=16,lwd = 2,cex=0.1,lty="solid",col="mediumorchid1",xaxt="n",type = "l")
lines(insul,S6K_pdgf$OPC, pch=16,lwd = 2,cex=0.1,lty="solid",col="cyan",xaxt="n",type = "l")
lines(insul,S6K_pdgf$oligodegro, pch=16,lwd = 2,cex=0.1,lty="solid",col="grey",xaxt="n",type = "l")
lines(insul,S6K_pdgf$neuron, pch=16,lwd = 2,cex=0.1,lty="solid",col="brown3",xaxt="n",type = "l")
axis(1, at = c(0.0001,0.001,0.01,0.1,1,10),cex.axis=1.4, las=1)

cel <- c("EC art","EC cap","EC ven","SMC vas","SMC art","peric st","peric ecm", "astro hpc","astro ctx","microglia","OD","OPC","neuron")
legend("top",cex = 1,bty="n",ncol=7, legend=cel, col=c("darkgreen","green","green3","deepskyblue3","deepskyblue4","royalblue1","royalblue2","orange","orange3","mediumorchid1","cyan","grey","brown3"),lty=rep("solid",length(cel)),lwd = 2)

#svg("Fig3A.svg", width = 12, height =4)
#jpeg("Fig3A.jpeg", width = 5, height =5.5, units = 'in', res = 600)
#dev.off()


