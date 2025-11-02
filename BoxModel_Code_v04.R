#####################################################################################################
# Box model of Otago coast
# Author: Sebastiaan van de Velde - University of Otago
# contact: sebastiaan.vandevelde@otago.ac.nz
#####################################################################################################

#=====================================================================================================
# Information/versions
#=====================================================================================================

# Model units:
# micromol
# cm
# year

# Model consists of two boxes representing the Otago coastal waters and seafloor. The atmosphere is upper boundary,
# and lateral water exchange on the non landward sides
# transport is via air-sea gas transfer, particle settling, and 
# a fixed fraction of mineralisation is allocated to anaerobic processes (= pyrite precipitation)

# Model state variables are: dissolved inorganic carbon (DIC)
#                            total alkalinity (TA)
#                            calcium (Ca)
#                            calcium carbonate (CaCO3)

# V01: 14/10/2024 SVDV


#=====================================================================================================
# load packages
#=====================================================================================================

library(pracma)
library(ReacTran)
library(marelac)
library(AquaEnv)
library(seacarb)

#source("../../Models/pHModule/SolveSaphePackage.R")

load("Sens1DLookUpTable_processed.Rdata")

#=====================================================================================================
# Aux. functions
#=====================================================================================================

#----------------------------------------------------------------------------------------------------
# Numerical limit for reactions
#----------------------------------------------------------------------------------------------------

FSAT <- function(C,K,n) (C/K)^n/((C/K)^n + 1)

#-----------------------------------------------------------------------------
# Function: initialize.params 
#-----------------------------------------------------------------------------

initialize.params <- function(PL=NULL){

  if (is.null(PL)) PL <- list()
  PL$N.var <- 5
  
  if (is.null(PL$C.lim)) PL$C.lim <- 1E-3 # [?mol cm-3] numerical limit for reactions
  PL$CFF <- 1e-3*1e4/365.25 # umol cm-2 yr-1 to mmol m-2 d-1
  
  # Geometry parameters
  
  if (is.null(PL$h.coast))    PL$h.coast <- 30*1e2      # [cm] Average water height coastal zone
  if (is.null(PL$A.coast))    PL$A.coast <- 33*125*1e10 # [cm2] 30km to shelf break, 125 km from Balclutha to Shag point
  PL$V.coast <- PL$A.coast*PL$h.coast                   # [cm3] 
  
  if (is.null(PL$h.seafloor)) PL$h.seafloor <- 5. # [cm] height of 'active' zone in seafloor   
  PL$V.seafloor <- PL$A.coast*PL$h.seafloor       # [cm3] 
  
  # Transport parameters
  
  if (is.null(PL$Q.clutha))   PL$Q.clutha   <- 600*1e6*3600*24*365.25 # [cm-3 yr-1] average discharge Clutha (https://web.archive.org/web/20050416222951/http://www.mfe.govt.nz/publications/ser/ser1997/html/chapter7.6.html)
  if (is.null(PL$Q.taieri))   PL$Q.taieri   <- 600*1e6*3600*24*365.25 # [cm-3 yr-1] average discharge Clutha (https://web.archive.org/web/20050416222951/http://www.mfe.govt.nz/publications/ser/ser1997/html/chapter7.6.html)
  if (is.null(PL$Q.offshore)) PL$Q.offshore <- PL$V.coast/(60/365.25) # [cm-3 yr-1] rate of water exchange between coastal and offshore waters based on residence time of 30-90 days (Liu et al., GRL, 2019)
  if (is.null(PL$v.settling)) PL$v.settling <- 1e2*365.25             # [cm yr-1] 1 m d-1 in coastal waters
  if (is.null(PL$v.burial))   PL$v.burial   <- 0.1                    # [cm yr-1] need to re-estimate
  
  if (is.null(PL$k.transfer))  PL$k.transfer <- 5.*1e2*365.25 # [cm yr-1] gas transfer velocity - interested over long-term, so assume rapid equilibration (Raymond and Cole, 2001, https://doi.org/10.2307/1352954)
  
  # Environmental parameters
  
  if (is.null(PL$S))   PL$S   <- 34.5 # [-] Salinity
  if (is.null(PL$TC))  PL$TC  <- 12. # [degC] average temperature coast
  if (is.null(PL$P))   PL$P   <- 1.  # [bar] pressure
  
  if (is.null(PL$sed_dens)) PL$sed_dens <- 2.6 # [g cm-3]
  if (is.null(PL$por))      PL$por      <- 0.4 # [-] porosity
  #PL$tort <- 1 - 2*log(PL$por)
  
  if (is.null(PL$constant.pCO2)) PL$constant.pCO2 <- T   # [logical] run with fixed atmospheric CO2?
  if (is.null(PL$pCO2atm))       PL$pCO2atm       <- 285 # [ppm] atmospheric CO2 concentration
  
  if (is.null(PL$TA.ocean))    PL$TA.ocean    <- 2.3*sw_dens(S=PL$S,t=PL$TC)*1e-3  # [umol cm-3] TA concentration ocean water - 2.3 umol kg-1 (Munida)
  if (is.null(PL$Ca.ocean))    PL$Ca.ocean    <- 10.0                              # [umol cm-3] Ca concentration ocean water
  if (is.null(PL$CaCO3.ocean)) PL$CaCO3.ocean <- 0.0                               # [umol cm-3] CaCO3 concentration ocean water

  if (is.null(PL$constant.clutha)) PL$constant.clutha <- TRUE # [logical] run with fixed riverine input?
  if (is.null(PL$TA.clutha))       PL$TA.clutha       <- 0.60 # [umol cm-3] TA concentration clutha (Kim et al., 1996, Marine and Freshwater Research)
  if (is.null(PL$Ca.clutha))       PL$Ca.clutha       <- 0.25 # [umol cm-3] Ca concentration clutha
  if (is.null(PL$CaCO3.clutha))    PL$CaCO3.clutha    <- 0.0  # [umol cm-3] CaCO3 concentration clutha

  if (is.null(PL$constant.taieri)) PL$constant.taieri <- TRUE # [logical] run with fixed riverine input?
  if (is.null(PL$TA.taieri))       PL$TA.taieri       <- 0.50 # [umol cm-3] TA concentration taieri (Kim et al., 2001, Marine and Freshwater Research)
  if (is.null(PL$Ca.taieri))       PL$Ca.taieri       <- 0.14 # [umol cm-3] Ca concentration taieri
  if (is.null(PL$CaCO3.taieri))    PL$CaCO3.taieri    <- 0.0  # [umol cm-3] CaCO3 concentration taieri
  
  if (is.null(PL$constant.calc)) PL$constant.calc  <- T                                  # [logical] run with fixed calcification?
  if (is.null(PL$calcification)) PL$calcification  <- 1e4*1e6/100*1e6/(23e3*2e3*4.5*1e6) # [umol cm-3 yr-1] calcification rate
  # take calcification rate from otago harbour and express per volume
  # 10,000 tonnes per year = 1e4*1e6 g / 100 g/mol 
  # Volume Otago harbour = 23km*2km*4.5m // 23e3*2e3*4.5 m3 
  if (is.null(PL$calrel)) PL$calrel <- 1 # integer [1/2] - 1, linear relation between calcification and omega (Andersson et al., Geology)
                                         #                 2, exponential relation between calcification and omega (Andersson et al., Geology)
  if (is.null(PL$Cal_ref)) PL$Cal_ref <- 1. # reference calcification
    
  if (is.null(PL$fixed.diss)) PL$fixed.diss  <- T        # [logical] run with fixed seafloor dissolution?
  if (is.null(PL$Cmin.sed))   PL$Cmin.sed    <- 81.17762 # [umol cm-3 yr-1] organi carbon mineralization rate
  if (is.null(PL$PRT.sed))    PL$PRT.sed     <- 5./24.   # [d] porewater residence time 
  if (is.null(PL$SeafloorDissolution)) PL$SeafloorDissolution <- 0.0 # fixed seafloor dissolution rate
  
  PL$Ksp.Cal  <- aquaenv(S=PL$S,t=PL$TC,Pa=PL$P)$Ksp_calcite   # [M^2] Calcite solubility product
  PL$Ksp.Ara  <- aquaenv(S=PL$S,t=PL$TC,Pa=PL$P)$Ksp_aragonite # [M^2] Aragonite solubility product
  
  return(PL)
}

#-----------------------------------------------------------------------------
# Function: set.inicons 
#-----------------------------------------------------------------------------

set.inicons <- function(PL){
  
  DIC.0 <- carb(flag=24,var1=285,var2=PL$TA.ocean/sw_dens(S=PL$S,t=PL$TC),S=PL$S,T=PL$TC,warn="n")$DIC*sw_dens(S=PL$S,t=PL$TC)
  
  inicons <- c(DIC.0,PL$TA.ocean,PL$Ca.ocean,PL$CaCO3.ocean,0.0)
  
  return(inicons)
}

#-----------------------------------------------------------------------------
# Function: transient interpolation function
# function that returns the dissolution rate in the seabed, based on
# mineralization rate, residence time, and saturation state
#-----------------------------------------------------------------------------

transient.int.function <- function(t,trend){
  
  value <- rep(NA,length(t))
  
  for (i in 1:length(t)){
    #print(t[i]) 
    NN <- which(abs(t[i]-trend$t)== min(abs(t[i]-trend$t)))
    #print(NN)
    if (length(NN)>1) {
      value[i] <- mean(trend$c[NN])
    } else {
      #if ((trend$t[NN]==t[i]) | (NN == 1) | (NN == nrow(trend))) { 
      if ((trend$t[NN]==t[i])) { 
        value[i] <- trend$c[NN]
      } else{
        if (trend$t[NN]>t[i]){
          value[i] <- ((trend$t[NN] - t[i])*trend$c[NN-1] + (t[i] - trend$t[NN-1])*trend$c[NN])/(trend$t[NN]-trend$t[NN-1])
        } 
        if (trend$t[NN]<t[i]) {
          value[i] <- ((trend$t[NN+1] - t[i])*trend$c[NN] + (t[i] - trend$t[NN])*trend$c[NN+1])/(trend$t[NN+1]-trend$t[NN])
        }
      }
    }
    
  }
  
  return(value)
}

#-----------------------------------------------------------------------------
# Function: Seafloor.Dissolution
# function that returns the dissolution rate in the seabed, based on
# PIC (in wt%) and saturation state
# returns CaCO3 dissolution in umol cm-3 yr-1 
#-----------------------------------------------------------------------------

Seafloor.Dissolution <- function(PIC,Ohm){
  
  Temp.diss.0  <- 0
  Temp.diss.01 <- 74.39	 - 53.8*Ohm  + 9.73*Ohm*Ohm
  Temp.diss.02 <- 107.95 - 74.33*Ohm + 12.77*Ohm*Ohm
  Temp.diss.03 <- 121.9  - 79.8*Ohm  + 12.93*Ohm*Ohm
  Temp.diss.04 <- 127.83 - 79.67*Ohm + 12.02*Ohm*Ohm
  Temp.diss.05 <- 131.22 - 79.07*Ohm + 11.38*Ohm*Ohm
  Temp.diss.06 <- 133.49 - 78.59*Ohm + 10.98*Ohm*Ohm
  Temp.diss.07 <- 135.17 - 78.28*Ohm + 10.74*Ohm*Ohm
  Temp.diss.08 <- 136.48 - 78.11*Ohm + 10.6*Ohm*Ohm
  Temp.diss.09 <- 137.54 - 78.02*Ohm + 10.52*Ohm*Ohm
  Temp.diss.1  <- 138.43 - 77.99*Ohm + 10.47*Ohm*Ohm
  Temp.diss.5  <- 146.35 - 78.47*Ohm + 10.47*Ohm*Ohm
  Temp.diss.10 <- 147.77 - 78.59*Ohm + 10.49*Ohm*Ohm
  Temp.diss.20 <- 148.64 - 78.66*Ohm + 10.5*Ohm*Ohm

  
  CaCO3.diss <- approx(x=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,5.,10.,20.),
                       y=c(Temp.diss.0,Temp.diss.01,Temp.diss.02,Temp.diss.03, 
                           Temp.diss.04,Temp.diss.05,Temp.diss.06,Temp.diss.07,
                           Temp.diss.08,Temp.diss.09,Temp.diss.1,Temp.diss.5, 
                           Temp.diss.10,Temp.diss.20),
                       xout=PIC,rule=2)$y
  
  if (Ohm >= 2.55) CaCO3.diss <- 0.0
    
  return(CaCO3.diss)
}

#-----------------------------------------------------------------------------
# Function: check.mass.balance 
#-----------------------------------------------------------------------------

check.mass.balance <- function(out,PL){
   
  DIC.deficit   <- (out$f_DIC_clutha   + out$f_DIC_taieri   - out$f_DIC_lateral   - out$f_DIC_exchange_out   + out$f_DIC_exchange_in   + out$f_SWI_DIC + out$f_airsea_DIC)/PL$V.coast   -    out$WaterColumnCalc
  Ca.deficit    <- (out$f_Ca_clutha    + out$f_Ca_taieri    - out$f_Ca_lateral    - out$f_Ca_exchange_out    + out$f_Ca_exchange_in    + out$f_SWI_Ca)/PL$V.coast    -    out$WaterColumnCalc
  TA.deficit    <- (out$f_TA_clutha    + out$f_TA_taieri    - out$f_TA_lateral    - out$f_TA_exchange_out    + out$f_TA_exchange_in    + out$f_SWI_TA)/PL$V.coast    - 2.*out$WaterColumnCalc
  CaCO3.deficit <- (out$f_CaCO3_clutha + out$f_CaCO3_taieri - out$f_CaCO3_lateral - out$f_CaCO3_exchange_out + out$f_CaCO3_exchange_in - out$f_SWI_CaCO3)/PL$V.coast +    out$WaterColumnCalc
  
  CaCO3.seafloor.deficit <- (out$f_SWI_CaCO3 - out$f_sed_CaCO3)/PL$V.seafloor - out$SeafloorDiss
  
  return(list(DIC.deficit=DIC.deficit,Ca.deficit=Ca.deficit,TA.deficit=TA.deficit,CaCO3.deficit=CaCO3.deficit,CaCO3.seafloor.deficit=CaCO3.seafloor.deficit))
}

#-----------------------------------------------------------------------------
# Post-processing functions: plot.transient, ExtractForPlot 
#-----------------------------------------------------------------------------

plot.transient <- function(out,PL,data=NULL){
  
  par(mfrow=c(1,1))
  plot(x=time_seq,y=out[,"2"]/sw_dens(S=PL$S,t=PL$TC)*1e6,
       xlab="time (year)", ylab=expression("[A"['T']*"] ("*mu*"mol kg"^'-1'*")"),
       type="l",ylim=range(c(out[,"2"]/sw_dens(S=PL$S,t=PL$TC)*1e6,as.numeric(data$AT)),na.rm=T )
       )
  if(!is.null(data)){
     date.numeric <- as.numeric(unlist(strsplit(as.character(Station.DataMatch$S9$time),split="-"))[seq(1,length.out=length(unlist(strsplit(as.character(Station.DataMatch$S9$time),split="-")))/3,by=3)]) +
                    as.numeric(unlist(strsplit(as.character(Station.DataMatch$S9$time),split="-"))[seq(2,length.out=length(unlist(strsplit(as.character(Station.DataMatch$S9$time),split="-")))/3,by=3)])/12 + 
                    as.numeric(unlist(strsplit(as.character(Station.DataMatch$S9$time),split="-"))[seq(3,length.out=length(unlist(strsplit(as.character(Station.DataMatch$S9$time),split="-")))/3,by=3)])/30
    points(x=date.numeric,y=as.numeric(data$AT),pch=1,col="red")
    
  } 
  legend("bottomright",bty='n',legend=c("model","data"),lty=c(1,NA),pch=c(NA,1),col=c("black","red"))
  
  #plot(x=time_seq,y=out[,"SeafloorDiss"]*PL$V.seafloor,
  #     xlab="time (year)", ylab=expression("R_CaCO"['3']*"] ("*mu*"mol cm"^'-3'*")"),
  #     type="l",ylim=range(c(out[,"SeafloorDiss"]*PL$V.seafloor,out[,"WaterColumnCalc"]*PL$V.coast)))
  #lines(x=time_seq,y=out[,"WaterColumnCalc"]*PL$V.coast,lty=2,col="red")
  #legend("bottomright",bty='n',legend=c("Dissolution","Production"),lty=c(1,2),col=c("black","red"))
  
  
  #plot(x=time_seq,y=out[,"5"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por),
  #     xlab="time (year)", ylab=expression("[CaCO"['3']*" (wt%)"),type="l")
  
  #plot.new()
  
}

ExtractForPlot <- function(xvar,yvar,zvar,input){
  
  x <- sort(unique(input[,xvar]))
  y <- sort(unique(input[,yvar]))
  
  extract <- matrix(nrow=length(x),ncol=length(y))
  
  for (i in 1:length(x)){
    for (j in 1:length(y)){
      
      sel          <- (input[,xvar] == x[i]) &  (input[,yvar] == y[j])
      extract[i,j] <- input[sel,zvar]
      
    }
  }
  
  return(list(x=x,y=y,z=extract))
  
}

#====================================================================================================
# Model function
#====================================================================================================

OtagoCoast.box <- function(t, state, PL){
  with(as.list(c(PL)),{
    
    DIC.coast      <- state[1]
    TA.coast       <- state[2]
    Ca.coast       <- state[3]
    CaCO3.coast    <- state[4]
    CaCO3.seafloor <- state[5]
    
    #print(DIC.coast )
    #print(TA.coast )
    #print(Ca.coast  )
    #print(CaCO3.coast  )
    #print(CaCO3.seafloor  )
    
    # Boundary conditions
    if (constant.pCO2) { 
      pCO2.t <- pCO2atm
    } else { 
      pCO2.t <- approx(x=pCO2atm$t,y=pCO2atm$c,xout=t,rule=2)$y#transient.int.function(t,trend=pCO2atm)
    }
    #print(t)
    #print(pCO2.t)
    
    if (constant.clutha) { 
      Q.c     <- Q.clutha
      TA.c    <- TA.clutha
      Ca.c    <- Ca.clutha
      CaCO3.c <- CaCO3.clutha
    } else { 
      Q.c     <- approx(x=Q$t,y=Q$Clutha,xout=t,rule=2)$y*1e6*3600*24*365.25#transient.int.function(t,trend=Q.clutha)
      TA.c    <- TA.clutha#transient.int.function(t,trend=TA.clutha)
      Ca.c    <- Ca.clutha#transient.int.function(t,trend=Ca.clutha)
      CaCO3.c <- CaCO3.clutha#transient.int.function(t,trend=CaCO3.clutha)
    }
    if (constant.taieri) { 
      Q.t     <- Q.taieri
      TA.t    <- TA.taieri
      Ca.t    <- Ca.taieri
      CaCO3.t <- CaCO3.taieri
    } else { 
      Q.t     <- approx(x=Q$t,y=Q$Taieri+Q$Waipori,xout=t,rule=2)$y*1e6*3600*24*365.25#transient.int.function(t,trend=Q.taieri)
      TA.t    <- TA.taieri#transient.int.function(t,trend=TA.taieri)
      Ca.t    <- Ca.taieri#transient.int.function(t,trend=Ca.taieri)
      CaCO3.t <- CaCO3.taieri#transient.int.function(t,trend=CaCO3.taieri)
    }
    if (constant.clutha & constant.taieri){
      Q.o <- Q.offshore
    } else {
      Q.o <- Q.offshore - Q.clutha - Q.taieri
    }
    
    # Derived variables
    DIC.c     <- carb(flag=24,var1=pCO2.t,var2=TA.c/sw_dens(S=0,t=TC),S=0,T=TC,warn="n")$DIC*sw_dens(S=0,t=TC)
    DIC.t     <- carb(flag=24,var1=pCO2.t,var2=TA.t/sw_dens(S=0,t=TC),S=0,T=TC,warn="n")$DIC*sw_dens(S=0,t=TC)
    #DIC.ocean  <- carb(flag=24,var1=pCO2.t-20.,var2=TA.ocean/sw_dens(S=35.,t=TC),S=35.,T=TC,warn="n")$DIC*sw_dens(S=35.,t=TC)
    DIC.ocean  <- carb(flag=24,var1=pCO2.t,var2=TA.ocean/sw_dens(S=35.,t=TC),S=35.,T=TC,warn="n")$DIC*sw_dens(S=35.,t=TC)
    
    CO2_coast  <- carb(flag=15,var1=TA.coast/sw_dens(S=S,t=TC),var2=DIC.coast/sw_dens(S=S,t=TC),S=S,T=TC,warn="n")$CO2
    #Carbspec   <- carb(flag=24,var1=(pCO2.t-40.),var2=TA.coast/sw_dens(S=S,t=TC),S=S,T=TC,warn="n")
    Carbspec   <- carb(flag=24,var1=(pCO2.t),var2=TA.coast/sw_dens(S=S,t=TC),S=S,T=TC,warn="n")
    CO2_eq     <- Carbspec$CO2
    Ohm.Ca     <- Carbspec$OmegaCalcite
    Ohm.Ar     <- Carbspec$OmegaAragonite
    pCO2.coast <- Carbspec$pCO2
    
    ## Reactions
    # Water column precipitation
    if (constant.calc) { 
      WaterColumnCalc <- (calcification)*V.coast*(Ohm.Ca/Ksp.Cal>1)*FSAT(C=TA.coast,K=1.e-3,n=5)
    } else {
      WaterColumnCalc <- ((21.3*Ohm.Ca + 12)/Cal_ref*calcification * (calrel==1) + 
                         ((228*(1-exp(-Ohm.Ca/0.69)) - 128)/Cal_ref*calcification) * (calrel==2) )*V.coast*(Ohm.Ca/Ksp.Cal>1)
    }
    # Seafloor dissolution (returns CaCO3 dissolution in umol cm-3 yr-1)
    PIC.seafloor <- CaCO3.seafloor*1e-6*molweight(c("Ca"))[[1]]/sed_dens*1e2/(1.-por) # convert CaCO3 from umol cm-3 to wt% C
    if (fixed.diss) { 
      SeafloorDiss <- SeafloorDissolution
    } else {
      #SeafloorDiss <- Seafloor.Dissolution(PIC=PIC.seafloor,Ohm=Ohm.Ar,LookUpTable=SeafloorLookUpTable)
      SeafloorDiss <- Seafloor.Dissolution(PIC=PIC.seafloor,Ohm=Ohm.Ar)
    }
    
    ## Transport
    # Riverine delivery [umol yr-1]
    f_DIC_clutha   <- Q.c * DIC.c
    f_TA_clutha    <- Q.c * TA.c
    f_Ca_clutha    <- Q.c * Ca.c
    f_CaCO3_clutha <- Q.c * CaCO3.c
    
    f_DIC_taieri   <- Q.t * DIC.t
    f_TA_taieri    <- Q.t * TA.t
    f_Ca_taieri    <- Q.t * Ca.t
    f_CaCO3_taieri <- Q.t * CaCO3.t
    
    # lateral transport [umol yr-1]
    f_DIC_lateral   <- (Q.c + Q.t) * DIC.coast
    f_TA_lateral    <- (Q.c + Q.t) * TA.coast
    f_Ca_lateral    <- (Q.c + Q.t) * Ca.coast
    f_CaCO3_lateral <- (Q.c + Q.t) * CaCO3.coast

    # shelf-slope exchange [umol yr-1]
    f_DIC_exchange_out   <- Q.o * DIC.coast
    f_TA_exchange_out    <- Q.o * TA.coast
    f_Ca_exchange_out    <- Q.o * Ca.coast
    f_CaCO3_exchange_out <- Q.o * CaCO3.coast
    
    f_DIC_exchange_in   <- Q.o * DIC.ocean
    f_TA_exchange_in    <- Q.o * TA.ocean
    f_Ca_exchange_in    <- Q.o * Ca.ocean
    f_CaCO3_exchange_in <- Q.o * CaCO3.ocean
    
    # Seafloor-water exchange [umol yr-1]
    f_SWI_DIC   <-    SeafloorDiss * V.seafloor
    f_SWI_TA    <- 2.*SeafloorDiss * V.seafloor
    f_SWI_Ca    <-    SeafloorDiss * V.seafloor
    f_SWI_CaCO3 <- v.settling * CaCO3.coast * A.coast

    # Sediment burial [umol yr-1]
    f_sed_CaCO3 <- v.burial * CaCO3.seafloor * A.coast
    
    # Air-sea gas exchange [umol yr-1]
    f_airsea_DIC <- A.coast*k.transfer*(CO2_eq - CO2_coast)  
    
    ## Assemble differential equations
    ddt.DIC.coast    <- (f_DIC_clutha   + f_DIC_taieri   - f_DIC_lateral   + f_DIC_exchange_in   - f_DIC_exchange_out   -   WaterColumnCalc + f_SWI_DIC + f_airsea_DIC)/V.coast
    ddt.TA.coast     <- (f_TA_clutha    + f_TA_taieri    - f_TA_lateral    + f_TA_exchange_in    - f_TA_exchange_out    - 2*WaterColumnCalc + f_SWI_TA)/V.coast
    ddt.Ca.coast     <- (f_Ca_clutha    + f_Ca_taieri    - f_Ca_lateral    + f_Ca_exchange_in    - f_Ca_exchange_out    -   WaterColumnCalc + f_SWI_Ca)/V.coast
    ddt.CaCO3.coast  <- (f_CaCO3_clutha + f_CaCO3_taieri - f_CaCO3_lateral + f_CaCO3_exchange_in - f_CaCO3_exchange_out +   WaterColumnCalc - f_SWI_CaCO3)/V.coast
    
    ddt.CaCO3.seafloor  <- (f_SWI_CaCO3 - f_sed_CaCO3 - SeafloorDiss*V.seafloor)/V.seafloor
    
    return(list(
      # Total rates of change
      c(ddt.DIC.coast, ddt.TA.coast, ddt.Ca.coast, ddt.CaCO3.coast, ddt.CaCO3.seafloor), 
      pH = Carbspec$pH, Ohm.Ca = Ohm.Ca, Ohm.Ar = Ohm.Ar, pCO2.coast = pCO2.coast,
      
      # Reactions
      "SeafloorDiss" = SeafloorDiss, "WaterColumnCalc" = WaterColumnCalc/V.coast, 

      # Fluxes
      f_DIC_clutha  = f_DIC_clutha,  f_TA_clutha  = f_TA_clutha,  f_Ca_clutha  = f_Ca_clutha,  f_CaCO3_clutha  = f_CaCO3_clutha,
      f_DIC_taieri  = f_DIC_taieri,  f_TA_taieri  = f_TA_taieri,  f_Ca_taieri  = f_Ca_taieri,  f_CaCO3_taieri  = f_CaCO3_taieri,
      f_DIC_lateral = f_DIC_lateral, f_TA_lateral = f_TA_lateral, f_Ca_lateral = f_Ca_lateral, f_CaCO3_lateral = f_CaCO3_lateral,
      
      f_DIC_exchange_out = f_DIC_exchange_out, f_TA_exchange_out = f_TA_exchange_out, f_Ca_exchange_out = f_Ca_exchange_out, f_CaCO3_exchange_out = f_CaCO3_exchange_out,
      f_DIC_exchange_in  = f_DIC_exchange_in,  f_TA_exchange_in  = f_TA_exchange_in,  f_Ca_exchange_in  = f_Ca_exchange_in,  f_CaCO3_exchange_in  = f_CaCO3_exchange_in,
      
      f_SWI_DIC   = f_SWI_DIC, f_SWI_TA = f_SWI_TA, f_SWI_Ca = f_SWI_Ca, f_SWI_CaCO3 = f_SWI_CaCO3,
      f_sed_CaCO3 = f_sed_CaCO3, f_airsea_DIC = f_airsea_DIC  
    )
    )
    
  })}

