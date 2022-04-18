rm(list=ls())
#BiocManager::install("")

expr <- read.csv("abundances_0422_control.csv", header = T, stringsAsFactors = FALSE)
expr1 <- read.csv("abundances_0422_AD.csv", header = T, stringsAsFactors = FALSE)

N <- (dim(expr)[2]-1)/2
x = 1:N # x axis width
expr_norm_all <- NULL

group <- c("EGF","BTC","NRG1","NRG2","NRG3","NRG4","HBEGF","TGFA")
gene <- "NRG"
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2

par(mar = c(7.1, 4.1, 4.1, 2.1))    # Fig. S3
plot(x,expr_[1,2:(N+1)], pch=16, cex=0.7,xaxt='n',bty="n", ylab=gene, xlab='',type = "b", lty = 1,ylim = c(0,8)) # NRG
min <- as.numeric(expr_[1,2:(N+1)]+expr_[1,(N+2):(2*N+1)])
max <- as.numeric(expr_[1,2:(N+1)]-expr_[1,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3)
lines(x,expr_[2,2:(N+1)],col="red",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[2,2:(N+1)]+expr_[2,(N+2):(2*N+1)])
max <- as.numeric(expr_[2,2:(N+1)]-expr_[2,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="red")
lines(x,expr_[3,2:(N+1)],col="blue",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[3,2:(N+1)]+expr_[3,(N+2):(2*N+1)])
max <- as.numeric(expr_[3,2:(N+1)]-expr_[3,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="blue")
lines(x,expr_[4,2:(N+1)],col="green",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[4,2:(N+1)]+expr_[4,(N+2):(2*N+1)])
max <- as.numeric(expr_[4,2:(N+1)]-expr_[4,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="green")
lines(x,expr_[5,2:(N+1)],col="purple",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[5,2:(N+1)]+expr_[5,(N+2):(2*N+1)])
max <- as.numeric(expr_[5,2:(N+1)]-expr_[5,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="purple")
lines(x,expr_[6,2:(N+1)],col="magenta",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[6,2:(N+1)]+expr_[6,(N+2):(2*N+1)])
max <- as.numeric(expr_[6,2:(N+1)]-expr_[6,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="magenta")
lines(x,expr_[7,2:(N+1)],col="grey",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[7,2:(N+1)]+expr_[7,(N+2):(2*N+1)])
max <- as.numeric(expr_[7,2:(N+1)]-expr_[7,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="grey")
lines(x,expr_[8,2:(N+1)],col="orange",pch=16, cex=0.7, ylab="",type = "b")
min <- as.numeric(expr_[8,2:(N+1)]+expr_[8,(N+2):(2*N+1)])
max <- as.numeric(expr_[8,2:(N+1)]-expr_[8,(N+2):(2*N+1)])
arrows(x, min, x, max, length=0.05, angle=90, code=3,col="orange")
#colo <- c("black")
#colo <- c("black","red")
#colo <- c("black","red","blue")
#colo <- c("black","red","blue","green")
#colo <- c("black","red","blue","green","purple","magenta")
colo <- c("black","red","blue","green","purple","magenta","grey","orange")
#colo <- c("blue","green","purple","magenta")
symb <- rep(16,dim(expr_)[1])
legend("topleft",ncol=2, bty="n",legend=expr_$gene,col=colo, pch=symb, cex=1)
axis(1, at = 1:N, labels = colnames(expr)[2:(N+1)], las=3,cex.axis=1,tck=-0.01)

expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("CAV1","EGFR","ERBB","ERBB1","ERBB2","ERBB3","ERBB4","FSHR")
gene <- c("ERBB")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2

expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

par(mar = c(7.1, 4.1, 4.1, 2.1))    # plor normalized Fig.S6
plot(x,expr_norm, pch=16, cex=0.7,xaxt='n',bty="n", ylab=gene, xlab='',type = "b", lty = 1,ylim = c(0,5))
min <- as.numeric(expr_norm-expr_err)
max <- as.numeric(expr_norm+expr_err)
arrows(x, min, x, max, length=0.05, angle=90, code=3)
axis(1, at = 1:N, labels = colnames(expr)[2:(N+1)], las=3,cex.axis=1,tck=-0.01)

group <- c("PDGFA","PDGFB","PDGFC","PDGFD")
gene <- c("PDGF")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("PDGFRA","PDGFRB")
gene <- c("PDGFR")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

cc <- data.frame(gene="f",t(rep(0,26)))
colnames(cc) <- colnames(expr_norm_all)
expr_norm_all <- rbind(expr_norm_all,cc)
cc <- data.frame(gene="ff",t(rep(0,26)))
colnames(cc) <- colnames(expr_norm_all)
expr_norm_all <- rbind(expr_norm_all,cc)

group <- c("IGF1","IGF2")
gene <- c("IGF")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("INSR","IGF1R","IGF2R")
gene <- c("IR")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("IRS1","IRS2")
gene <- c("IRS")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- "PTEN"
gene <- c("PTEN")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("AKT1","AKT2","AKT3")
gene <- c("AKT")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("TSC1","TSC2")
gene <- c("TSC")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)

group <- c("RPS6KB1","RPS6KB2")
gene <- c("S6K")
expr_ <- expr[expr$gene %in% group,]
expr_ <- expr_[order(match(expr_$gene,group)),]
expr_sum <- colSums(expr_[,2:(N+1)])
expr_err <- apply(expr_[,c((N+2):(2*N+1))], 2, function(c) max(c))
expr1_ <- expr1[expr1$gene %in% group,]
expr1_ <- expr1_[order(match(expr1_$gene,group)),]
expr1_sum <- colSums(expr1_[,2:(N+1)])
expr_sum_tot <- c(expr_sum,expr1_sum)
norm_coef <- quantile(expr_sum_tot[expr_sum_tot>0])[4]
expr_norm <- expr_sum/norm_coef
expr_err <- expr_err/norm_coef
expr_norm[unname(which(expr_norm==0))] <- min(expr_norm[expr_norm>0])/2
expr_norm_all_0 <- data.frame(gene=gene,t(expr_norm),t(expr_err))
expr_norm_all <- rbind(expr_norm_all, expr_norm_all_0)
#write.csv(as.matrix(expr_norm_all), file = "snRNAseq_grouped_upquant_0422.csv", row.names=F)

# for AD expression exchange lines 5,6
expr1 <- read.csv("abundances_0422_control.csv", header = T, stringsAsFactors = FALSE)
expr <- read.csv("abundances_0422_AD.csv", header = T, stringsAsFactors = FALSE)
# and repeat all above steps
#write.csv(as.matrix(expr_norm_all), file = "snRNAseq_grouped_upquant_0422_AD.csv", row.names=F)

