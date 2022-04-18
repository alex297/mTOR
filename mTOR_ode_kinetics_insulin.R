rm(list=ls())
#BiocManager::install("deSolve")
library(deSolve)
#ode
model_0 <- function (t, x, params) {
  ## extract variables
  IRm <- x[1]
  IRins <- x[2]
  IRmYP <- x[3]
  IRiYP <- x[4]
  IRi <- x[5]
  IRS1 <- x[6]
  IRS1YP <- x[7]
  IRS1YPS307 <- x[8]
  IRS1S307 <- x[9]
  X <- x[10]
  Xp <- x[11]
  PKB <- x[12]
  PKBT308 <- x[13]
  PKBS473 <- x[14]
  PKBT308S473 <- x[15]
  mTORC1 <- x[16]
  mTORC1a <- x[17]
  mTORC2 <- x[18]
  mTORC2a <- x[19]
  S6K <- x[20]
  S6KT389 <- x[21]
  S6 <- x[22]
  S6S235 <- x[23]
  r <- x[24]
  c1 <- x[25]
  c2 <- x[26]
  L <- x[27]
  c2e <- x[28] # endcyted c2
  ## extract parameters
  insulin <- params["insulin"]
  k3_0 <- params["k3_0"]
  k1a <- params["k1a"]
  k1c <- params["k1c"]
  k1d <- params["k1d"]
  k1f <- params["k1f"]
  k1g <- params["k1g"]
  k1r <- params["k1r"]
  k2_0 <- params["k2_0"]
  k2a <- params["k2a"]
  k2b <- params["k2b"]
  k2c <- params["k2c"]
  k2d <- params["k2d"]
  k2f <- params["k2f"]
  k2g <- params["k2g"]
  k3a <- params["k3a"]
  k3b <- params["k3b"]
  k4a <- params["k4a"]
  k4b <- params["k4b"]
  k4c <- params["k4c"]
  k4f <- params["k4f"]
  k4h <- params["k4h"]
  k5b <- params["k5b"]
  k5c <- params["k5c"]
  k5d <- params["k5d"]
  k9b1 <- params["k9b1"]
  k9b2 <- params["k9b2"]
  k9f1 <- params["k9f1"]
  k9f2 <- params["k9f2"]
  km9 <- params["km9"]
  TSCtot <- params["TSCtot"]
  k5a1 <- params["k5a1"]
  kiTSC <- params["kiTSC"]
  k5 <- params["k5"]
  kmx <- params["kmx"]
  kx <- params["kx"]
  KDL <- params["KDL"]
  kf <- params["kf"]
  ke <- params["ke"]
  kLmax <- params["kLmax"]
  KML <- params["KML"]
  kPI3K <- params["kPI3K"]
  aPI3K <- params["aPI3K"]
  k4d <- params["k4d"]
  k5e <- params["k5e"]
  PTEN <- params["PTEN"]
  kPTEN <- params["kPTEN"]
  k4c_0 <- params["k4c_0"]
  k4_0 <- params["k4_0"]
  R0 <- params["R0"]
  kr <- KDL*kf
  PI3K <- 0.5*(1+kPI3K+2*aPI3K*R0*c2-sqrt((1+kPI3K+2*aPI3K*R0*c2)^2-8*aPI3K*R0*c2))
  
  ## model equations
  dIRmdt <- k1r*IRi-k1a*IRm*insulin+k1g*IRmYP
  dIRinsdt <- k1a*IRm*insulin-k1c*IRins
  dIRmYPdt <- k1c*IRins-k1d*IRmYP-k1g*IRmYP
  dIRiYPdt <- k1d*IRmYP-k1f*IRiYP*Xp
  dIRidt <- k1f*IRiYP*Xp-k1r*IRi
  dIRS1dt <- k2b*IRS1YP+k2g*IRS1S307-k2_0*IRS1-k2a*IRS1*IRiYP
  dIRS1YPdt <- k2a*IRS1*IRiYP+k2d*IRS1YPS307-k2b*IRS1YP-k2c*IRS1YP*mTORC1a
  dIRS1YPS307dt <- k2c*IRS1YP*mTORC1a-k2d*IRS1YPS307-k2f*IRS1YPS307
  dIRS1S307dt <- k2_0*IRS1+k2f*IRS1YPS307-k2g*IRS1S307
  dXdt <- k3b*Xp-k3a*X*(k3_0+IRS1YPS307)
  dXpdt <- k3a*X*(k3_0+IRS1YPS307)-k3b*Xp
  dPKBdt <- k4b*PKBT308+k4h*PKBS473-k4a*PKB*IRS1YP/(PTEN+kPTEN)-k4d*PKB*PI3K/(PTEN+kPTEN)-k4_0*PKB
  dPKBT308dt <- k4a*PKB*IRS1YP/(PTEN+kPTEN)+k4d*PKB*PI3K/(PTEN+kPTEN)-k4b*PKBT308-k4c*PKBT308*mTORC2a+k4_0*PKB-k4c_0*PKBT308
  dPKBS473dt <- k4f*PKBT308S473-k4h*PKBS473-k4a*PKBS473*IRS1YP/(PTEN+kPTEN)-k4d*PKBS473*PI3K/(PTEN+kPTEN)-k4_0*PKBS473
  dPKBT308S473dt <- k4c*PKBT308*mTORC2a+k4a*PKBS473*IRS1YP/(PTEN+kPTEN)+k4d*PKBS473*PI3K/(PTEN+kPTEN)-k4f*PKBT308S473+k4c_0*PKBT308+k4_0*PKBS473
  dmTORC1dt <- k5b*mTORC1a-k5*mTORC1*(PKBT308S473+k5a1*PKBT308)/(TSCtot+kiTSC)
  #dmTORC1dt <- k5b*mTORC1a-0*k5*mTORC1*(PKBT308S473+k5a1*PKBT308)/(TSCtot+kiTSC)
  dmTORC1adt <- k5*mTORC1*(PKBT308S473+k5a1*PKBT308)/(TSCtot+kiTSC)-k5b*mTORC1a
  #dmTORC1adt <- 0*k5*mTORC1*(PKBT308S473+k5a1*PKBT308)/(TSCtot+kiTSC)-k5b*mTORC1a
  dmTORC2dt <- k5d*mTORC2a-k5c*mTORC2*IRiYP/(PTEN+kPTEN)-k5e*mTORC2*PI3K/(PTEN+kPTEN)
  dmTORC2adt <- k5c*mTORC2*IRiYP/(PTEN+kPTEN)+k5e*mTORC2*PI3K/(PTEN+kPTEN)-k5d*mTORC2a
  dS6Kdt <- k9b1*S6KT389-k9f1*S6K*mTORC1a/(km9+mTORC1a)
  dS6KT389dt <- k9f1*S6K*mTORC1a/(km9+mTORC1a)-k9b1*S6KT389
  dS6dt <- k9b2*S6S235-k9f2*S6*S6KT389
  dS6S235dt <- k9f2*S6*S6KT389-k9b2*S6S235
  drdt <- kr*c1-kf*L*r+kmx*c2
  dc1dt <- kf*L*r-kr*c1+kmx*c2-2*kx*R0*c1*c1
  dc2dt <- kx*R0*c1*c1-(kmx+ke)*c2
  dLdt <- R0*(kr*c1-kf*L*r+kmx*c2)-kLmax*L/(1+L/KML)
  dc2edt <- ke*c2
  ## combine results into a single vector
  dxdt <- c(dIRmdt,dIRinsdt,dIRmYPdt,dIRiYPdt,dIRidt,dIRS1dt,dIRS1YPdt,dIRS1YPS307dt,dIRS1S307dt,dXdt,dXpdt,
            dPKBdt,dPKBT308dt,dPKBS473dt,dPKBT308S473dt,dmTORC1dt,dmTORC1adt,dmTORC2dt,dmTORC2adt,dS6Kdt,dS6KT389dt,dS6dt,dS6S235dt,
            drdt,dc1dt,dc2dt,dLdt,dc2edt)
  ## return result as a list
  list(dxdt)
}

# for normal conditions
TSC <- 1
# for TSC-ko (on Fig. S1) use TSC <- 0.001
# for rapamicin treatment (on Fig. S1) use commented eqs in the model with 0 rates of mTORC1 activation

# first simulate background activation in absence of insulin
parms <- c(insulin=0,k1a=0.633,k1c=0.877,k1d=31,k1f=184000,k1g=1940,k1r=0.547,k2_0=0.0423,k2a=323,k2b=3420,
           k2c=576000,k2d=281,k2f=2.91,k2g=0.267,k3_0=0.01,k3a=0.0001,k3b=0.0988,k4a=579000,k4b=34.8,k4c=446,k4f=30,k4h=0.536,k5b=24.8,
           k5c=2,k5d=1.06,k9b1=0.0444,k9b2=31,k9f1=0.9,k9f2=333,km9=58.7,TSCtot=TSC,k5a1=0.03,kiTSC=0.001,k5=40,
           kmx=0.07,kx=30,KDL=1.5,kf=1,ke=0.2,kLmax=0.011,KML=1.66,kPI3K=0.3,aPI3K=8000,k4d=30,k5e=0.2,PTEN=1,kPTEN=1e-5,
           R0=0.01,k4_0=4,k4c_0=0.3) 

times <- seq(from=0,to=600,by=1/6)
xstart <- c(IRm=1,IRins=0,IRmYP=0,IRiYP=0,IRi=0,IRS1=1,IRS1YP=0,IRS1YPS307=0,IRS1S307=0,X=1,Xp=0,PKB=1,PKBT308=0,PKBS473=0,PKBT308S473=0,
            TORC1=1,TORC1a=0,TORC2=1,TORC2a=0,S6K=1,S6KT389=0,S6=1,S6S235=0,r=1,c1=0,c2=0,L=0,c2e=0)

library(tidyverse)
ode(
  func=model_0,
  y=xstart,
  times=times,
  parms=parms
) %>%
  as.data.frame() -> out

xstart <- as.numeric(as.vector(out[3601,-1]))

# simulate addition of 10 nM of insulin (for Fig. S2)
parms <- c(insulin=10,k1a=0.633,k1c=0.877,k1d=31,k1f=184000,k1g=1940,k1r=0.547,k2_0=0.0423,k2a=323,k2b=3420,
           k2c=576000,k2d=281,k2f=2.91,k2g=0.267,k3_0=0.01,k3a=0.0001,k3b=0.0988,k4a=579000,k4b=34.8,k4c=446,k4f=30,k4h=0.536,k5b=24.8,
           k5c=2,k5d=1.06,k9b1=0.0444,k9b2=31,k9f1=0.9,k9f2=333,km9=58.7,TSCtot=TSC,k5a1=0.03,kiTSC=0.001,k5=40,
           kmx=0,kx=30,KDL=1.5,kf=1,ke=0.2,kLmax=0.011,KML=1.66,kPI3K=0.3,aPI3K=8000,k4d=100,k5e=0.3,PTEN=1,kPTEN=1e-5,
           R0=0.01,k4_0=4,k4c_0=0.3)

times <- seq(from=0,to=120,by=1/6)
ode(
  func=model_0,
  y=xstart,
  times=times,
  parms=parms    # method = "ode45"
  ) %>%
  as.data.frame() -> out

# Fig. S1 - normal conditions (insulin = 1 nM)
plot(out[,1],(out[,15]+out[,16]), pch=16,col="blue",cex=0.1,lty="solid", ylim=c(0, 1),main = "insulin = 1 nM",xlim=c(0,30),ylab='PKB_S473, a.u.', xlab='time, min',type = "b") #xlim=c(0, 16),
#lines(out[,1],(out[,15]+out[,16]), pch=16,cex=0.1,lty="solid", col="green",type = "b") # add simulation of TSC-ko
#lines(out[,1],(out[,15]+out[,16]), pch=16,cex=0.1,lty="solid", col="red",type = "b") # add simulation of rapamicin treatment
#legend("top",ncol=2, legend=c("+rap","-rap"),col=c("red","blue"), pch=c(16,16), cex=1)

plot(out[,1],(out[,14]+out[,16]), pch=16,col="blue",cex=0.1,lty="solid", ylim=c(0, 1),xlim=c(0,30),main = "insulin = 1 nM",ylab='PKB_T308, a.u.', xlab='time, min',type = "b") #xlim=c(0, 16),
#lines(out[,1],(out[,14]+out[,16]), pch=16,cex=0.1,lty="solid", col="green",type = "b") # add simulation of TSC-ko
#lines(out[,1],(out[,14]+out[,16]), pch=16,cex=0.1,lty="solid", col="red",type = "b") # add simulation of rapamicin treatment
#legend("top",ncol=2, legend=c("+rap","-rap"),col=c("red","blue"), pch=c(16,16), cex=1)

plot(out[,1],out[,18], pch=16,cex=0.1,lty="solid",col="blue", ylab='mTORC1a, a.u.',main = "insulin = 1 nM",ylim=c(0, 1.1),xlim=c(0,30), xlab='time, min',type = "b") #xlim=c(0, 16),
#lines(out[,1],out[,18], pch=16,cex=0.1,lty="solid", type = "b",col="green") # add simulation of TSC-ko
#lines(out[,1],out[,18], pch=16,cex=0.1,lty="solid", type = "b",col="red") # add simulation of rapamicin treatment
#legend(14,0.9,ncol=1, legend=c("TSC=1","TSC=0.001,-rep","TSC=1(0.001),+rap"),col=c("blue","green","red"), pch=c(16,16,16), cex=1)

plot(out[,1],out[,20], pch=16,cex=0.1,lty="solid",col="blue", ylab='mTORC2a, a.u.',main = "insulin = 1 nM", xlab='time, min',xlim=c(0,30),ylim=c(0,0.01),type = "b") 
#lines(out[,1],out[,20], pch=16,cex=0.1,lty="solid", col="green",type = "b") # add simulation of TSC-ko
#lines(out[,1],out[,20], pch=16,cex=0.1,lty="solid", col="red",type = "b") # add simulation of rapamicin treatment

# Fig. S2
x_PKB_T308 <- c(0,10,30)
y_PKB_T308 <- c(1,2.03,2.77)
ky_PKB_T308 <- 4.36
ky_PKB_S473 <- 1.6
plot(out[,1],(out[,14]+out[,16])*ky_PKB_T308, pch=16,col="blue",cex=0.1,lty="solid",ylim=c(0, 4),xlim=c(0,30), ylab='PKB_T308, a.u.', xlab='time, min',type = "b") #xlim=c(0, 16),
points(x_PKB_T308,y_PKB_T308,pch=16, cex=0.5,type = "b",col="red")

x_PKB_S473 <- c(0,0.65,0.8,1.73,2.85,5,7,10,15,20,30)
y_PKB_S473 <- c(0.11,0.22,0.68,1.12,1.25,1.37,1,1.35,1.27,1.29,1.28)
plot(out[,1],(out[,15]+out[,16])*ky_PKB_S473, pch=16,col="blue",cex=0.1,lty="solid", ylim=c(0, 1.5),xlim=c(0,30),ylab='PKB_S473, a.u.', xlab='time, min',type = "b") #xlim=c(0, 16),
points(x_PKB_S473,y_PKB_S473,pch=16, cex=0.5,type = "b",col="red")

x_IR <- c(0,0.65,0.81,1.14,2.1,3,5,7,10,15,20,30)
y_IR <- c(0.32,2.24,1.92,1.8,1.48,1.28,1.42,1.2,1.32,1.12,1.14,1)
ky_IR <- 525
plot(out[,1],(out[,4]+out[,5])*ky_IR, pch=16,col="blue",cex=0.1,lty="solid", ylab='IRp, a.u.', xlab='time, min',type = "b", xlim=c(0, 30)) #xlim=c(0, 16),
points(x_IR,y_IR,pch=16, cex=0.5,type = "b",col="red")

x_IRS <- c(0,0.7,1,1.2,2.1,3,5,7,10,15,20,30)
y_IRS <- c(0.93,3.95,3.71,3.27,2.99,2.59,2.43,2.23,2.43,2.2,1.84,2.34)
ky_IRS1_YP <- 58
ky_IRS1_S307 <- 6.7 
plot(out[,1],(out[,8]+out[,9])*ky_IRS1_YP, pch=16,col="blue",ylim=c(0, 6),cex=0.1,lty="solid", ylab='IRS_YP, a.u.', xlab='time, min',type = "b", xlim=c(0, 30)) #xlim=c(0, 16),
points(x_IRS,y_IRS,pch=16, cex=0.5,type = "b",col="red")

x_IRS_S307 <- c(0,0.65,0.8,0.96,2,3,5,7,10,15,20,30)
y_IRS_S307 <- c(1.79,2.37,2.59,3.17,4.11,2.77,3.39,3.75,3.57,2.81,3.79,2.99)
plot(out[,1],(out[,9]+out[,10])*ky_IRS1_S307, pch=16,col="blue",cex=0.1,lty="solid", ylim=c(0, 6),xlim=c(0, 30), ylab='IRS_S307, a.u.', xlab='time, min',type = "b") #xlim=c(0, 16),
points(x_IRS_S307,y_IRS_S307,pch=16, cex=0.5,type = "b",col="red")
