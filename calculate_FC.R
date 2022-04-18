# calculate sensitivity for PTEN/TSC for maximum of S6Ka

# this code plots used output files from ode runs, which look like:
# ins_PTEN05_neuron.csv; ins_PTEN2_neuron.csv; ins_TSC05_neuron.csv; ins_TSC2_neuron.csv, et.cet

#These files should be loaded from subdirectory /outputs while running this code, for example by setwd() to this subdirectory

# The output files were generated using
# mTOR_output_ligand_range_ins.R
# mTOR_output_ligand_range_PDGF.R
# mTOR_output_ligand_range_NRG.R

ligand <- 10   
ins <- 10
while (ins > 0.01) {
  ins <- ins/1.04
  ligand <- c(ligand,ins)
}
brain_ini <- read.csv("snRNAseq_grouped_upquant_0422.csv", header = T, stringsAsFactors = FALSE)
cells <- colnames(brain_ini)[-1]
N <- length(cells)

ins_max <- NULL
nrg_max <- NULL
pdgf_max <- NULL
fc_PTEN_ins <- NULL
fc_PTEN_nrg <- NULL
fc_PTEN_pdgf <- NULL
fc_TSC_ins <- NULL
fc_TSC_nrg <- NULL
fc_TSC_pdgf <- NULL
for(j in 1:N) {
  cell <- cells[j]
  name <- paste0("ins_PTEN05_",cell,".csv")
  res_PTEN05 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  name <- paste0("ins_PTEN2_",cell,".csv")
  res_PTEN2 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  s6k <- as.numeric(unname(res_PTEN05[21,-1]))
  i <- which(s6k==max(s6k))
  ins_max[j] <- ligand[i]
  fc_PTEN_ins[j] <- (res_PTEN05[21,i+1]-res_PTEN05[21,179])/(res_PTEN2[21,i+1]-res_PTEN2[21,179])
  name <- paste0("ins_TSC05_",cell,".csv")
  res_TSC05 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  name <- paste0("ins_TSC2_",cell,".csv")
  res_TSC2 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  fc_TSC_ins[j] <- (res_TSC05[21,i+1]-res_TSC05[21,179])/(res_TSC2[21,i+1]-res_TSC2[21,179])
  name <- paste0("nrg_PTEN05_",cell,".csv")
  res_PTEN05 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  name <- paste0("nrg_PTEN2_",cell,".csv")
  res_PTEN2 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  s6k <- as.numeric(unname(res_PTEN05[21,-1]))
  i <- which(s6k==max(s6k))
  nrg_max[j] <- ligand[i]
  fc_PTEN_nrg[j] <- (res_PTEN05[21,i+1]-res_PTEN05[21,179])/(res_PTEN2[21,i+1]-res_PTEN2[21,179])
  name <- paste0("nrg_TSC05_",cell,".csv")
  res_TSC05 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  name <- paste0("nrg_TSC2_",cell,".csv")
  res_TSC2 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  fc_TSC_nrg[j] <- (res_TSC05[21,i+1]-res_TSC05[21,179])/(res_TSC2[21,i+1]-res_TSC2[21,179])
  name <- paste0("pdgf_PTEN05_",cell,".csv")
  res_PTEN05 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  name <- paste0("pdgf_PTEN2_",cell,".csv")
  res_PTEN2 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  s6k <- as.numeric(unname(res_PTEN05[21,-1]))
  i <- which(s6k==max(s6k))
  pdgf_max[j] <- ligand[i]
  fc_PTEN_pdgf[j] <- (res_PTEN05[21,i+1]-res_PTEN05[21,179])/(res_PTEN2[21,i+1]-res_PTEN2[21,179])
  name <- paste0("pdgf_TSC05_",cell,".csv")
  res_TSC05 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  name <- paste0("pdgf_TSC2_",cell,".csv")
  res_TSC2 <- read.csv(name, header = T, stringsAsFactors = FALSE)
  fc_TSC_pdgf[j] <- (res_TSC05[21,i+1]-res_TSC05[21,179])/(res_TSC2[21,i+1]-res_TSC2[21,179])
}
FC <- data.frame(cell=cells,ins_max=ins_max,pdgf_max=pdgf_max,nrg_max=nrg_max,FC_PTEN_ins=fc_PTEN_ins,FC_PTEN_pdgf=fc_PTEN_pdgf,
                 FC_PTEN_nrg=fc_PTEN_nrg,FC_TSC_ins=fc_TSC_ins,FC_TSC_pdgf=fc_TSC_pdgf,FC_TSC_nrg=fc_TSC_nrg)

#write.csv(as.matrix(FC), file = "fold_changes_0422.csv", row.names=F)
