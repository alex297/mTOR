rm(list=ls())
#BiocManager::install("deSolve")
library(deSolve)
#?ode
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

TSC <- 1
parms <- c(insulin=0,k1a=0.633,k1c=0.877,k1d=31,k1f=184000,k1g=1940,k1r=0.547,k2_0=0.0423,k2a=323,k2b=3420,
           k2c=576000,k2d=281,k2f=2.91,k2g=0.267,k3_0=0.01,k3a=0.0001,k3b=0.0988,k4a=579000,k4b=34.8,k4c=446,k4f=30,k4h=0.536,k5b=24.8,
           k5c=2,k5d=1.06,k9b1=0.0444,k9b2=31,k9f1=0.9,k9f2=333,km9=58.7,TSCtot=TSC,k5a1=0.03,kiTSC=0.001,k5=40,
           kmx=0.07,kx=30,KDL=9,kf=1,ke=0.06,kLmax=0.011,KML=1.66,kPI3K=100,aPI3K=8000,k4d=30,k5e=0.2,PTEN=1,kPTEN=1e-5,
           R0=0.01,k4_0=4,k4c_0=0.3)

L0 <- 10 # choose 10; 1; 0.1 nM NRG for Fig. S2C

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
xstart[27] <- L0
times <- seq(from=0,to=120,by=1/6)
ode(
  func=model_0,
  y=xstart,
  times=times,
  parms=parms    # method = "ode45"
  ) %>%
  as.data.frame() -> out

ky_PKB_S473 <- 3
x_AKT <- c(0,60,115,300,600,1800)/60 
y_AKT_01 <- c(0,3.3,5.9,7.3,4.8,3.5)/100*ky_PKB_S473 # NRG=0.1 nM
y_AKT_1 <- c(0,11.7,18.3,15.9,11.9,8.6)/100*ky_PKB_S473 # NRG=1 nM
y_AKT_10 <- c(0,18.1,39.9,33.5,29.1,24.4)/100*ky_PKB_S473 # NRG=10 nM

# Fig. S2 C
plot(out[,1],out[,15]+out[,16]-out[1,15]-out[1,16], lwd = 2,lty="solid",col="grey70",ylab='PKB-Ser473, a.u.', xlab='time, min',type = "l", ylim=c(0, 1.5),xlim=c(0, 30)) # for L0=10
lines(out[,1],out[,15]+out[,16]-out[1,15]-out[1,16], lwd = 2,lty="solid",col="grey50",type = "l") # for L0=1
lines(out[,1],out[,15]+out[,16]-out[1,15]-out[1,16], lwd = 2,lty="solid",col="grey0",type = "l") # for L0=0.1
points(x_AKT,y_AKT_10,pch=16, type = "p", lwd = 2,col="grey70")
points(x_AKT,y_AKT_1,pch=16, type = "p", lwd = 2,col="grey50")
points(x_AKT,y_AKT_01,pch=16, type = "p", lwd = 2,col="grey0")
legend("top",ncol=3, legend=c("   0.1 nM","1 nM","10 nM"), bty = "n",
       col=c("grey0","grey50","grey70"), pch=c(16,16,16), cex=1)

