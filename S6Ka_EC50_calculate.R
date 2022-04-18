rm(list=ls())

# this code calculates S6Ka and EC50 based on the output files from
# ode runs of mTOR simulations for wide range of ligand concentrations 
# these output files are stored in subdirectory /outputs and have names like
# ins_ref_neuron.csv, ins_AD_ref_neuron.csv, et.cet. 

#These files should be loaded while running this code, for example by setwd() to subdirectory with these files as shown below

# The output files were generated using
# mTOR_output_ligand_range_ins.R
# mTOR_output_ligand_range_PDGF.R
# mTOR_output_ligand_range_NRG.R
# and for AD:
# mTOR_output_ligand_range_ins_AD.R
# mTOR_output_ligand_range_PDGF_AD.R
# mTOR_output_ligand_range_NRG_AD.R

ligand <- 10   
ins <- 10
while (ins > 0.01) {
  ins <- ins/1.04
  ligand <- c(ligand,ins)
}
brain_ini <- read.csv("snRNAseq_grouped_upquant_0422.csv", header = T, stringsAsFactors = FALSE)
#brain_ini <- read.csv("snRNAseq_grouped_upquant_AD_0422.csv", header = T, stringsAsFactors = FALSE)
cells <- colnames(brain_ini)[-1]
S6K_ins <- data.frame(cell=cells,S6K_max=0,T05=0,ligand="ins")
S6K_pdgf <- data.frame(cell=cells,S6K_max=0,T05=0,ligand="pdgf")
S6K_nrg <- data.frame(cell=cells,S6K_max=0,T05=0,ligand="nrg")
N <- length(cells)

# to loaded the files required below use, for windows platform:
path <- getwd()
setwd(paste0(path,"/outputs"))

for(j in 1:N) {
  cell <- cells[j]
  name <- paste0("ins_ref_",cell,".csv")
  #name <- paste0("ins_AD_ref_",cell,".csv")
  res <- read.csv(name, header = T, stringsAsFactors = FALSE)
  s6k <- as.numeric(unname(res[21,-1]))
  S6K_ins$S6K_max[j] <- max(s6k)
  s6k <- s6k-min(s6k)
  if (max(s6k)>0) {
    s6k <- s6k/max(s6k)
    for (i in 1:length(s6k)) {
      if (s6k[i]>0.5) {
        i_max <- i
      }
    }
  } else {
    s6k <- 0
    i_max <- 1
  }
  S6K_ins$T05[j] <- ligand[i_max]
  name <- paste0("pdgf_ref_",cell,".csv")
  #name <- paste0("pdgf_AD_ref_",cell,".csv")
  res <- read.csv(name, header = T, stringsAsFactors = FALSE)
  s6k <- as.numeric(unname(res[21,-1]))
  S6K_pdgf$S6K_max[j] <- max(s6k)
  s6k <- s6k-min(s6k)
  if (max(s6k)>0) {
    s6k <- s6k/max(s6k)
    for (i in 1:length(s6k)) {
      if (s6k[i]>0.5) {
        i_max <- i
      }
    }
  } else {
    s6k <- 0
    i_max <- 1
  }
  S6K_pdgf$T05[j] <- ligand[i_max]
  name <- paste0("nrg_ref_",cell,".csv")
  #name <- paste0("nrg_AD_ref_",cell,".csv")
  res <- read.csv(name, header = T, stringsAsFactors = FALSE)
  s6k <- as.numeric(unname(res[21,-1]))
  S6K_nrg$S6K_max[j] <- max(s6k)
  s6k <- s6k-min(s6k)
  if (max(s6k)>0) {
    s6k <- s6k/max(s6k)
    for (i in 1:length(s6k)) {
      if (s6k[i]>0.5) {
        i_max <- i
      }
    }
  } else {
    s6k <- 0
    i_max <- 1
  }
  S6K_nrg$T05[j] <- ligand[i_max]
}

S6K_nrg$S6K_max <- S6K_nrg$S6K_max*2
S6K_nrg$T05 <- S6K_nrg$T05/5
s6k <-rbind(S6K_ins,S6K_pdgf,S6K_nrg)
s6k$cell <- factor(s6k$cell, levels = s6k$cell[1:N])
s6k <- data.frame(s6k,EC50_reciprocal=1/s6k$T05)

# change working directory back before writing:
setwd(path)

#write.csv(as.matrix(s6k), file = "S6Ka_max_EC50_all.csv", row.names=F)  # with scaling of NRG
#write.csv(as.matrix(s6k), file = "S6Ka_max_EC50_all_AD.csv", row.names=F)  
