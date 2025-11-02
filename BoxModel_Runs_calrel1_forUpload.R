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
# source code
#=====================================================================================================

source("BoxModel_Code_v04.R")

#=====================================================================================================
# load data
#=====================================================================================================

load("../Data/Munida2403/MunidaData.Rdata")

#====================================================================================================
# Model validation and uncertainties
#====================================================================================================
#----------------------------------------------------------------------------------------------------
# Get historical pCO2 data
#----------------------------------------------------------------------------------------------------

# pCO2 is extracted from Baring Head https://data.niwa.co.nz/collections/climate
# column names of the dataframe 'pCO2.trend' are "t" and"c"
# t should be in numerical year (e.g. 2012.1 or 2024.5)
# C is in ppm/uatm

#----------------------------------------------------------------------------------------------------
# Get historical River discharge data
#----------------------------------------------------------------------------------------------------

# river discharge is downloaded from https://envdata.orc.govt.nz/AQWebPortal and is in m3/s

#----------------------------------------------------------------------------------------------------
# Data/estimates available
#----------------------------------------------------------------------------------------------------
# Coastal water residence time: 30-90 days (Liu et al., GRL, 2019)
# Coastal water residence time: could be as low as 14 days (https://doi.org/10.1016/j.csr.2024.105248)
# Sediment inorganic carbon content: 14-40 wt% Carbonate (Bostock et al., 2018) - Carter et al., 1986
# Coastal water TA concentration: 2.17 - 2.30 umol kg-1 (Munida)
# Calcification rate: 100-1000 g CaCO3 m-2 yr-1 
#----------------------------------------------------------------------------------------------------

# set initial conditions
PL <- list()
PL <- initialize.params(PL=PL)

set.seed(1)
iterations <- 10000

Validation.output <- matrix(data=rep(NA,22*iterations),
                            ncol=22,nrow=iterations,
                            dimnames=list(1:iterations,
                                          c("Q.exchange","WaterCalc","v.settling","v.burial",
                                            "TA.clutha","TA.taieri","CaCO3.clutha","CaCO3.taieri",
                                            "Cmin.sed","PRT.sed",
                                            "SeafloorDiss","DIC","TA","Ca","CaCO3.water","CaCO3.seafloor","pH","pCO2","Ohm.Ca",
                                            "Ohm.Ara","f_SWI_TA","f_SWI_DIC"))
)


# shelf exchange [cm-3 yr-1] Liu et al., 2019, Johnsen et al., 2014
ss          <- rnorm(n=iterations,mean=30,sd=10.)
ss[ss<(10)] <- NA
ss[ss>(90)] <- NA
Validation.output[,"Q.exchange"] <- PL$V.coast/(ss/365.25)
# calcification rate [umol cm-3 yr-1] Carter, 1986
ss                   <- rnorm(n=iterations,mean=1150,sd=500.)
ss[ss<0]             <- NA
ss[ss>(1150+3*500.)] <- NA
Validation.output[,"WaterCalc"] <- ss/100*1e6*(PL$A.coast*1e-4)/PL$V.coast
# settling velocity [cm yr-1] De Kruyf/McCandliss
ss         <- rlnorm(n=iterations,meanlog=-1.,sdlog=log(2.5))*1e-2 #cm s-1
ss[ss<0.]  <- NA
ss[ss>10.] <- NA
Validation.output[,"v.settling"] <- ss*3600*24*365.25              
# burial velocity [cm yr-1] sand accumulation from Carter (1.39 Mton yr-1), 2.6 g cm-3 sed density, surface area of study area
ss                    <- rnorm(n=iterations,mean=0.01,sd=0.005)
ss[ss<0]              <- NA
ss[ss>(0.01+3*0.005)] <- NA
Validation.output[,"v.burial"] <- ss  
# TA clutha [umol cm-3] (Kim et al., 1996, Marine and Freshwater Research)
ss                 <- rnorm(n=iterations,mean=0.6,sd=0.3)
ss[ss<0]           <- NA
ss[ss>(0.6+3*0.3)] <- NA
Validation.output[,"TA.clutha"] <- ss # [umol cm-3] TA concentration clutha (Kim et al., 1996, Marine and Freshwater Research)
# TA Taieri [umol cm-3] (Kim et al., 2001, Marine and Freshwater Research)
ss                  <- rnorm(n=iterations,mean=0.5,sd=0.25)
ss[ss<0]            <- NA
ss[ss>(0.5+3*0.25)] <- NA
Validation.output[,"TA.taieri"] <- ss # [umol cm-3] TA concentration taieri (Kim et al., 2001, Marine and Freshwater Research)

# CaCO3 clutha [umol cm-3] 3.14 Mt y‐1 (Carter)
# 0.26 wt% [0.24-0.28] (Muller et al.)
ss                   <- rnorm(n=iterations,mean=0.26,sd=0.02)
ss[ss<(0.26-3*0.02)] <- NA
ss[ss>(0.26+3*0.02)] <- NA
Validation.output[,"CaCO3.clutha"] <- ss/100*3.14*1e6*1e6/molweight(c("CaCO3"))/(mean(Q$Clutha,na.rm=T)*1e6*3600*24*365.25)*1e6 # [umol cm-3] CaCO3 concentration clutha (derived from Carter)
# CaCO3 taieri [umol cm-3] 0.6 Mt y‐1 (Carter)
Validation.output[,"CaCO3.taieri"] <- ss/100*0.6*1e6*1e6/molweight(c("CaCO3"))/(mean(Q$Taieri+Q$Waipori,na.rm=T)*1e6*3600*24*365.25)*1e6 # [umol cm-3] CaCO3 concentration clutha (derived from Carter)
# Cmin [umol cm-3 yr-1] (Burdige et al., 2007, Chem Rev) Cmin for 0-200 m
ss                <- rnorm(n=iterations,mean=30,sd=20)
ss[ss<0]          <- NA
ss[ss>(30+3*20)]  <- NA
Validation.output[,"Cmin.sed"] <- ss*1e3*1e-4/10*365.25
# PRT [d] 
# Lunstrum & Berelson (2022): 80 to 1130 L m-2 d-1
# Reimers et al (2004): 300 to 1,300 L m-2 d-1
# Precht and Huettel (2004; Santos et al. (2012); McGinnis et al. (2014): ~100 L m-2 d-1
# McLachlan (1989): 8,000 L m-2 d-1 (high-energy environment)
ss                <- rnorm(n=iterations,mean=150,sd=75)
ss[ss<0]          <- NA
ss[ss>(150+3*75)] <- NA
Validation.output[,"PRT.sed"] <- (ss/24)/24.


Validation.output <- Validation.output[complete.cases(Validation.output[ , 1:10]),]

par(mfrow=c(2,5))
for (i in 1:10) hist(Validation.output[,i],main=colnames(Validation.output)[i],xlab=colnames(Validation.output)[i])


for (i in 1:nrow(Validation.output)){

    PL <- list()
    # initialize default parameter list
    PL   <- initialize.params(PL=PL)
    # change parameters of choice
    PL$pCO2atm       <- approx(x=pCO2.trend$t,y=pCO2.trend$c,xout=1998.)$y
    PL$Q.clutha      <- mean(Q$Clutha,na.rm=T)*1e6*3600*24*365.25#approx(x=Q$t,y=Q$Clutha,xout=1998.,rule=2)$y*1e6*3600*24*365.25
    PL$Q.taieri      <- mean(Q$Taieri+Q$Waipori,na.rm=T)*1e6*3600*24*365.25#approx(x=Q$t,y=Q$Taieri+Q$Waipori,xout=1998,rule=2)$y*1e6*3600*24*365.25
    
    PL$CaCO3.clutha  <- Validation.output[i,"CaCO3.clutha"]
    PL$CaCO3.taieri  <- Validation.output[i,"CaCO3.taieri"]
    PL$Q.offshore    <- Validation.output[i,"Q.exchange"] - (PL$Q.clutha + PL$Q.taieri)
    PL$calcification <- Validation.output[i,"WaterCalc"]
    PL$v.settling    <- Validation.output[i,"v.settling"]
    PL$v.burial      <- Validation.output[i,"v.burial"]
    PL$TA.clutha     <- Validation.output[i,"TA.clutha"]
    PL$TA.taieri     <- Validation.output[i,"TA.taieri"]
    PL$Cmin.sed      <- Validation.output[i,"Cmin.sed"]
    PL$PRT.sed       <- Validation.output[i,"PRT.sed"]
    
    PL$fixed.diss <- FALSE
    
    # update parameter list
    PL   <- initialize.params(PL=PL)
    
    # set initial conditions
    yini <- set.inicons(PL=PL)
    
    # solve model
    out <- steady(y = yini, func = OtagoCoast.box, parms = PL, method="stode",positive=T)
    head(out)
    
    Validation.output[i,"SeafloorDiss"]   <- out$SeafloorDiss
    Validation.output[i,"DIC"]            <- out$y[1]
    Validation.output[i,"TA"]             <- out$y[2]
    Validation.output[i,"Ca"]             <- out$y[3]
    Validation.output[i,"CaCO3.water"]    <- out$y[4]
    Validation.output[i,"CaCO3.seafloor"] <- out$y[5]
    Validation.output[i,"pH"]             <- out$pH
    Validation.output[i,"pCO2"]           <- out$pCO2
    Validation.output[i,"Ohm.Ca"]         <- out$Ohm.Ca[[1]]
    Validation.output[i,"Ohm.Ara"]        <- out$Ohm.Ar[[1]]
    
    Validation.output[i,"f_SWI_TA"]       <- out$f_SWI_TA[[1]]
    Validation.output[i,"f_SWI_DIC"]      <- out$f_SWI_DIC[[1]]
    
    print(paste("PIC = ",out$y[5]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por), "wt%"))
    print(paste("[TA]=",out$y[2]))
    print(paste("pH=",out$pH))
    
}

save(Validation.output,file="MonteCarloValidation_250506.Rdata")


selection <- ((Validation.output[,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por)) <40.) & 
             ((Validation.output[,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por)) >5.) & 
             ((Validation.output[,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6) < 2275.) & 
             ((Validation.output[,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6) > 2232.) &
             (Validation.output[,"pCO2"] > 275.) & (Validation.output[,"pCO2"] < 375.)
  

x11(width=12,height=6)
par(mfrow=c(2,5))

hist(PL$V.coast/(Validation.output[selection,"Q.exchange"]+PL$Q.clutha+PL$Q.taieri)*365.25,xlab="Residence time (days)",main="")

hist(Validation.output[selection,"WaterCalc"]*100/1e6/(PL$A.coast*1e-4)*PL$V.coast,xlab="Calcification (g m-2 yr-1)",main="")
abline(v=100.,col="black",lty=2)
abline(v=1000.,col="black",lty=2)
hist(Validation.output[selection,"TA.clutha"]*1000,xlab="TA Clutha (umol L-1)",main="")
abline(v=600,col="black",lty=2)
hist(Validation.output[selection,"TA.taieri"]*1000,xlab="TA Taieri (umol L-1)",main="")
abline(v=500,col="black",lty=2)

hist((Validation.output[selection,"PRT.sed"])*24,xlab="Porewater residence time (hours)",main="")
hist((Validation.output[selection,"Cmin.sed"]),xlab="Sediment mineralization rate (umol cm-3 yr-1)",main="")


hist(Validation.output[selection,"f_SWI_TA"]/(PL$A.coast/2)*1e-3*1e4/365.25,xlab="F_TA (mmol m-2 d-1)",main="")
hist(Validation.output[selection,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por),xlab="CaCO3 (wt%)",main="")
abline(v=14.,col="black",lty=2)
hist(Validation.output[selection,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6,xlab="[AT] (umol kg-1)",main="")
abline(v=2262.3,col="black",lty=2)
#2232, 2275, 2262.3
hist(Validation.output[selection,"Ohm.Ara"],xlab="Aragonite saturation state",main="")


savePlot(filename="Spinup_250506.png",type="png")

#====================================================================================================
# Transient pCO2 Munida specific
#====================================================================================================
#----------------------------------------------------------------------------------------------------
# load parameter sets 
#----------------------------------------------------------------------------------------------------

PL <- list()
PL <- initialize.params(PL=PL)

load("MonteCarloValidation_250506.Rdata")

selection <- ((Validation.output[,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por)) <40.) & 
  ((Validation.output[,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por)) >5.) & 
  ((Validation.output[,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6) < 2275.) & 
  ((Validation.output[,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6) > 2232.) &
  (Validation.output[,"pCO2"] > 275.) & (Validation.output[,"pCO2"] < 350.)


filtered.parametersets <- Validation.output[selection,]

#----------------------------------------------------------------------------------------------------
# Get historical pCO2 data New Zealand
#----------------------------------------------------------------------------------------------------

# pCO2 is extracted from Baring Head https://data.niwa.co.nz/collections/climate
# column names of the dataframe 'pCO2.trend' are "t" and"c"
# t should be in numerical year (e.g. 2012.1 or 2024.5)
# C is in ppm/uatm

#----------------------------------------------------------------------------------------------------
# Get historical River discharge data
#----------------------------------------------------------------------------------------------------

# river discharge is downloaded from https://envdata.orc.govt.nz/AQWebPortal and is in m3/s

#----------------------------------------------------------------------------------------------------
# Run transient simulations
#----------------------------------------------------------------------------------------------------

TransientRun <- list()

# initialize default parameter list

PL <- list()
PL <- initialize.params(PL=PL)

# change parameters of choice

for (i in 1:nrow(filtered.parametersets)){
  
  PL$pCO2atm         <- approx(x=pCO2.trend$t,y=pCO2.trend$c,xout=1998.)$y
  PL$constant.pCO2   <- TRUE
  PL$Q.clutha        <- mean(Q$Clutha)*1e6*3600*24*365.25
  PL$Q.taieri        <- mean(Q$Taieri+Q$Waipori,na.rm=T)*1e6*3600*24*365.25
  PL$constant.clutha <- TRUE
  PL$constant.taieri <- TRUE
  
  PL$fixed.diss      <- FALSE
  
  PL$CaCO3.clutha  <- filtered.parametersets[i,"CaCO3.clutha"]
  PL$CaCO3.taieri  <- filtered.parametersets[i,"CaCO3.taieri"]
  PL$Q.offshore    <- filtered.parametersets[i,"Q.exchange"] - (PL$Q.clutha + PL$Q.taieri)
  PL$calcification <- filtered.parametersets[i,"WaterCalc"]
  PL$v.settling    <- filtered.parametersets[i,"v.settling"]
  PL$v.burial      <- filtered.parametersets[i,"v.burial"]
  PL$TA.clutha     <- filtered.parametersets[i,"TA.clutha"]
  PL$TA.taieri     <- filtered.parametersets[i,"TA.taieri"]
  PL$Cmin.sed      <- filtered.parametersets[i,"Cmin.sed"]
  PL$PRT.sed       <- filtered.parametersets[i,"PRT.sed"]
  
  
  # update parameter list
  
  PL   <- initialize.params(PL=PL)
  
  # set initial conditions
  
  yini <- set.inicons(PL=PL)
  
  # solve model
  out <- steady(y = yini, func = OtagoCoast.box, parms = PL, method="stode",positive=T)
  head(out)
  print(check.mass.balance(out,PL))
  
  # Turn on transient pCO2
  
  PL$constant.pCO2   <- FALSE
  PL$pCO2atm         <- pCO2.trend # 0 = januari 1972
  PL$constant.clutha <- TRUE
  PL$constant.taieri <- TRUE
  PL$Q               <- Q
  PL$Q.offshore      <- filtered.parametersets[i,"Q.exchange"] - (PL$Q.clutha + PL$Q.taieri)
  PL$constant.calc   <- FALSE
  PL$calrel          <- 1
  PL$Cal_ref         <- (21.3*out$Ohm.Ca[[1]] + 12)
  
  # update parameter list
  
  PL   <- initialize.params(PL=PL)
  
  # set initial conditions
  
  yini <- out$y
  
  # Time sequence
  
  t_start  <- 1998   # year
  t_stop   <- 2030   # year
  n_output <- (t_stop-t_start)*12 + 1 # number of output times: per month
  
  time_seq <- seq(from = t_start, to = t_stop, length = n_output) # years
  
  # solve model
  transient.out <- ode(y = yini, times = time_seq, func = OtagoCoast.box, parms = PL, method="vode",maxsteps=20000)
  head(transient.out)
  tail(transient.out)
  
  # plot output
  
  plot.transient(out=transient.out,PL,data=Station.DataMatch$S9)
  
  # save output
  
  TransientRun[[i]] <- list()
  TransientRun[[i]]$output <- transient.out
  TransientRun[[i]]$PL     <- PL
  
}

save(TransientRun,file="OtagoCoastBox_TransientRun_250506.Rdata")

#====================================================================================================
# Transient pCO2 Munida specific - no sediment response
#====================================================================================================
#----------------------------------------------------------------------------------------------------
# load parameter sets 
#----------------------------------------------------------------------------------------------------

# set initial conditions
PL <- list()
PL <- initialize.params(PL=PL)

load("MonteCarloValidation_250506.Rdata")

selection <- ((Validation.output[,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por)) <40.) & 
  ((Validation.output[,"CaCO3.seafloor"]*1e-6*molweight(c("CaCO3"))[[1]]/PL$sed_dens*1e2/(1.-PL$por)) >5.) & 
  ((Validation.output[,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6) < 2275.) & 
  ((Validation.output[,"TA"]/sw_dens(S=PL$S,t=PL$TC)*1e6) > 2232.) &
  (Validation.output[,"pCO2"] > 275.) & (Validation.output[,"pCO2"] < 375.)


filtered.parametersets <- Validation.output[selection,]

#----------------------------------------------------------------------------------------------------
# Get historical pCO2 data New Zealand
#----------------------------------------------------------------------------------------------------

# pCO2 is extracted from Baring Head https://data.niwa.co.nz/collections/climate
# column names of the dataframe 'pCO2.trend' are "t" and"c"
# t should be in numerical year (e.g. 2012.1 or 2024.5)
# C is in ppm/uatm

#----------------------------------------------------------------------------------------------------
# Get historical River discharge data
#----------------------------------------------------------------------------------------------------

# river discharge is downloaded from https://envdata.orc.govt.nz/AQWebPortal and is in m3/s

#----------------------------------------------------------------------------------------------------
# Run transient simulations
#----------------------------------------------------------------------------------------------------

TransientRun <- list()

# initialize default parameter list

PL <- list()
PL <- initialize.params(PL=PL)

# change parameters of choice

for (i in 1:nrow(filtered.parametersets)){
  
  PL$pCO2atm         <- approx(x=pCO2.trend$t,y=pCO2.trend$c,xout=1998.)$y
  PL$constant.pCO2   <- TRUE
  PL$Q.clutha        <- mean(Q$Clutha)*1e6*3600*24*365.25
  PL$Q.taieri        <- mean(Q$Taieri+Q$Waipori,na.rm=T)*1e6*3600*24*365.25
  PL$constant.clutha <- TRUE
  PL$constant.taieri <- TRUE
  
  PL$fixed.diss      <- TRUE
  PL$SeafloorDissolution <- filtered.parametersets[i,"SeafloorDiss"]
  
  PL$CaCO3.clutha  <- filtered.parametersets[i,"CaCO3.clutha"]
  PL$CaCO3.taieri  <- filtered.parametersets[i,"CaCO3.taieri"]
  PL$Q.offshore    <- filtered.parametersets[i,"Q.exchange"] - (PL$Q.clutha + PL$Q.taieri)
  PL$calcification <- filtered.parametersets[i,"WaterCalc"]
  PL$v.settling    <- filtered.parametersets[i,"v.settling"]
  PL$v.burial      <- filtered.parametersets[i,"v.burial"]
  PL$TA.clutha     <- filtered.parametersets[i,"TA.clutha"]
  PL$TA.taieri     <- filtered.parametersets[i,"TA.taieri"]
  PL$Cmin.sed      <- filtered.parametersets[i,"Cmin.sed"]
  PL$PRT.sed       <- filtered.parametersets[i,"PRT.sed"]
  
  
  # update parameter list
  
  PL   <- initialize.params(PL=PL)
  
  # set initial conditions
  
  yini <- set.inicons(PL=PL)
  
  # solve model
  out <- steady(y = yini, func = OtagoCoast.box, parms = PL, method="stode",positive=T)
  head(out)
  print(check.mass.balance(out,PL))
  
  # Turn on transient pCO2
  
  PL$constant.pCO2   <- FALSE
  PL$pCO2atm         <- pCO2.trend # 0 = januari 1972
  PL$constant.clutha <- TRUE
  PL$constant.taieri <- TRUE
  PL$Q               <- Q
  PL$Q.offshore      <- filtered.parametersets[i,"Q.exchange"] - (PL$Q.clutha + PL$Q.taieri)
  PL$constant.calc   <- FALSE
  PL$calrel          <- 1
  PL$Cal_ref         <- (21.3*out$Ohm.Ca[[1]] + 12)
  
  # update parameter list
  
  PL   <- initialize.params(PL=PL)
  
  # set initial conditions
  
  yini <- out$y
  
  # Time sequence
  
  t_start  <- 1998   # year
  t_stop   <- 2030   # year
  n_output <- (t_stop-t_start)*12 + 1 # number of output times: per month
  
  time_seq <- seq(from = t_start, to = t_stop, length = n_output) # years
  
  # solve model
  transient.out <- ode(y = yini, times = time_seq, func = OtagoCoast.box, parms = PL, method="lsoda")
  head(transient.out)
  tail(transient.out)
  
  # plot output
  
  plot.transient(out=transient.out,PL,data=Station.DataMatch$S9)
  
  # save output
  
  TransientRun[[i]] <- list()
  TransientRun[[i]]$output <- transient.out
  TransientRun[[i]]$PL     <- PL
  
}

save(TransientRun,file="OtagoCoastBox_TransientRun_NoDynamicSeafloor_250506.Rdata")
