# Analyses for Figure 3
# January 2019
# Isla Myers-Smith

packrat::init()

# Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(viridis)

# Libraries for Bayesian phenology models
library(rjags)
library(R2jags)

# Customised ggplot2 theme function ----
theme_QHI <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = -45, hjust = -0.05),
          axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size=20, vjust=1, hjust=0),
          axis.ticks.length=unit(.25, "cm"),
          legend.text = element_text(size=16, face="italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))
}

# Phenology analyses ----

# Phenology
qiphen <- read.csv("data/qiki_phen_with_before_2017.csv", stringsAsFactors = F)
# Change column name Spp to SPP to match with script (column name was changes in 2 Nov. 2017 data update)
colnames(qiphen)[which(colnames(qiphen)=="Spp")] <- "SPP"

# P2 Model
# Keep only observations for which our variable of interest ("P2") was observed
qi <- qiphen[!is.na(qiphen$P2),]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an event was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qi[,c("P2_before","P2")])

# Make sure that P2 is always a larger number than P2_before, or else JAGS will be unhappy (NA's not accepted)
which(qi$P2_before>qi$P2)
qi[is.na(qi$P2_before),]

# Can set P2_before values that don't meet this criteria to 10 days before the "P2" date
# qi$P2_before[qi$P2_before > qi$P2 | is.na(qi$P2_before)] <- qi$P2[qi$P2_before > qi$P2 | is.na(qi$P2_before)] - 10

# Subset to SALARC
qi <- subset(qi, SPP == "SALARC")

# PREPARE DATA FOR JAGS IMPORT

qi <- qi[order(qi$Year, qi$Plot.ID),]
qi$PlotNum <- as.numeric(factor(qi$Plot.ID, levels=unique(qi$Plot.ID)))

jags.dat <- list(
  lim = cbind(qi$P2_before,qi$P2), # ANNE REMOVED -1 (subtracted from before date)
  y1s = ifelse(qi$P2_before < qi$P2,1,0), #column of 1's indicates iterval censoring
  year = qi$Year - min(qi$Year) + 1, #year must start with 1
  nyear = length(unique(qi$Year)),
  plot = qi$PlotNum,
  n = nrow(qi),
  nplot = length(unique(qi$PlotNum)),
  xhat = c(sort(unique(qi$Year - min(qi$Year) + 1))),
  nxhat = length(unique(qi$Year - min(qi$Year) + 1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      # PRIORS
      a ~ dunif(0, 365)
      b ~ dnorm(0, 0.001)
      for (i in 1:nplot){
      aplot[i] ~ dnorm(0, tau.plot)
      }
      sigma.plot ~ dunif(0, 100)
      tau.plot <- 1/(sigma.plot*sigma.plot)
      sigma.year ~ dunif(0, 100)
      tau.year <- 1/(sigma.year*sigma.year)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
      #LIKELIHOOD
      for (j in 1:nyear){
      a.year[j]~dnorm(muyear[j], tau.year)
      muyear[j] <- a+b*j
      }
      for (i in 1:n){
      y1s[i]~dinterval(t[i],lim[i,])
      t[i]~dnorm(pred[i],tau)
      pred[i] <- aplot[plot[i]]+a.year[year[i]]
      }
      # DERIVED QUANTITIES
      for (j in 1:nxhat){
      phat[j] <- a + b*xhat[j]
      }
      }
      ","models/int_cens_P2.jags")

# INITIAL VALUES ("t" is mandatory)
inits <- function() list(aplot = rnorm(jags.dat$nplot,0,2), sigma.plot = runif(1,0,1), sigma = runif(1,0,1), t = as.vector(apply(jags.dat$lim, 1, mean)))

# PARAMETERS TO MONITOR
params <- c("a", "b", "phat","aplot", "sigma", "sigma.plot", "sigma.year") ## ADD PARAMETERS TO MONITOR

modoutP2 <- jags(jags.dat, inits, params, model.file="models/int_cens_P2.jags", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, DIC=FALSE, working.directory=NULL, progress.bar = "text")

#plot(modoutP2) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP2$BUGSoutput$summary[,c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')])
hist(coefs$Rhat) # all should be below 1.1
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH RESULTS
coefs.predP2 <- coefs[coefs$Type=="phat",]

# Add species and year information
coefs.predP2$SPP <- qi$SPP[1]
coefs.predP2$YearC <- c(1:17)
coefs.predP2$Year <- coefs.predP2$YearC + min(qi$Year) - 1

# Safe coefs for reference later
coefsP2 <- coefs

# Safe subset of phenology data for later use
qiP2 <- qi

# # 5% and 95% quantiles
# P2quants <- rbind(round(quantile(modoutP2$BUGSoutput$sims.list$a, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$b, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$sigma.year, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$sigma.plot, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$sigma, probs = c(0.05, 0.95)), 2))
# rownames(P2quants) <- c("a", "b", "sigma.year", "sigma.plot", "sigma")
# P2quants

# P3 Model
# Keep only observations for which our variable of interest ("P3") was observed
qi <- qiphen[!is.na(qiphen$P3),]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an even was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qi[,c("P3_before","P3")])

# Make sure that P3 is always a larger number than P3_before, or else JAGS will be unhappy (NA's not accepted)
which(qi$P3_before>qi$P3)
qi[is.na(qi$P3_before),]

# PREPARE DATA FOR JAGS IMPORT
qi$PlotNum <- as.numeric(as.factor(qi$Plot.ID))
qi$SpeciesNum <- as.numeric(as.factor(qi$SPP))
qi <- qi[order(qi$SPP, qi$Year, qi$PlotNum),]

jags.dat <- list(
  lim = cbind(qi$P3_before,qi$P3),
  y1s = ifelse(qi$P3_before<qi$P3,1,0), #column of 1's indicates iterval censoring
  year = qi$Year - min(qi$Year) + 1, #year must start with 1 (this becomes necessary in the more complex models)
  nyear = length(unique(qi$Year)),
  plot = qi$PlotNum,
  species = qi$SpeciesNum,
  n = nrow(qi),
  nspp = length(unique(qi$SpeciesNum)),
  nplot = length(unique(qi$PlotNum)),
  xhat = matrix(c(sort(unique(qi$Year - min(qi$Year) + 1))), nrow = 3,ncol = length(unique(qi$Year - min(qi$Year) + 1)),byrow = T),
  nxhat = length(unique(qi$Year - min(qi$Year) + 1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      # PRIORS
      for (i in 1:nspp){
      aspp[i] ~ dunif(0, 365)
      bspp[i] ~ dnorm(0, 0.001)
      }
      for (i in 1:nplot){
      aplot[i] ~ dnorm(0, tau.plot)
      }
      sigma.plot ~ dunif(0, 100)
      tau.plot <- 1/(sigma.plot*sigma.plot)
      sigma.year ~ dunif(0, 100)
      tau.year <- 1/(sigma.year*sigma.year)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
      #LIKELIHOOD
      for (j in 1:nspp){
      for (k in 1:nyear){
      a.sppyear[j,k]~dnorm(musppyear[j,k], tau.year)
      musppyear[j,k] <- aspp[j]+bspp[j]*k
      }
      }
      for (i in 1:n){
      y1s[i]~dinterval(t[i],lim[i,])
      t[i]~dnorm(pred[i],tau)
      pred[i] <- aplot[plot[i]]+a.sppyear[species[i],year[i]]
      }
      # DERIVED QUANTITIES
      for (i in 1:nspp){
      for (j in 1:nxhat){
      phat[i,j] <- aspp[i] + bspp[i]*xhat[i,j]
      }
      }
      diff1to2 <- aspp[1] - aspp[2]
      diff1to3 <- aspp[1] - aspp[3]
      diff2to3 <- aspp[2] - aspp[3]
      }
      ","models/int_cens_P3.jags")

# INITIAL VALUES ("t" is mandatory)
inits <- function() list(aplot = rnorm(jags.dat$nplot,0,2), sigma.plot = runif(1,0,1), sigma = runif(1,0,1), t = as.vector(apply(jags.dat$lim, 1, mean)))

# PARAMETERS TO MONITOR
params <- c("aspp", "bspp", "phat","aplot", "sigma", "sigma.plot", "sigma.year", "diff1to2", "diff1to3", "diff2to3") ## ADD PARAMETERS TO MONITOR

modoutP3 <- jags(jags.dat,inits, params, model.file="models/int_cens_P3.jags", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, DIC=FALSE, working.directory=NULL, progress.bar = "text")

#plot(modoutP3) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP3$BUGSoutput$summary[,c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')])
hist(coefs$Rhat)
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH RESULTS
coefs.predP3 <- coefs[coefs$Type=="phat",] 

# Add species and year information
coefs.predP3$SpeciesNum <- rep(c(1:3), times = length(unique(jags.dat$year)))
coefs.predP3$SPP <- qi$SPP[match(coefs.predP3$SpeciesNum, qi$SpeciesNum)]
coefs.predP3$YearC <- rep(c(1:17), each = jags.dat$nspp)
coefs.predP3$Year <- coefs.predP3$YearC + min(qi$Year) - 1
# remove prediction for ERIVAG 2001 as no data is available for this year
coefs.predP3 <- coefs.predP3[!(coefs.predP3$Year == 2001 & coefs.predP3$SPP == "ERIVAG"),]

# Save coefs and subset of phenology data for later use
qiP3 <- qi
coefsP3 <- coefs

# 5% and 95% quantiles
# P3quants <- rbind(round(quantile(modoutP3$BUGSoutput$sims.list$aspp[,1], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$bspp[,1], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$aspp[,2], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$bspp[,2], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$aspp[,3], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$bspp[,3], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$sigma.year, probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$sigma.plot, probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$sigma, probs = c(0.05, 0.95)), 2))
# rownames(P3quants) <- c("a1", "b1", "a2", "b2", "a3", "b3", "sigma.year", "sigma.plot", "sigma")
# P3quants

# P5 Model
# Keep only observations for which our variable of interest ("P5") was observed
qi <- qiphen[!is.na(qiphen$P5),]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an even was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qi[,c("P5_before","P5")])

# Make sure that P5 is always a larger number than P5_before, or else JAGS will be unhappy (NA's not accepted)
which(qi$P5_before>qi$P5)
qi[is.na(qi$P5_before),]

# Subset to SALARC
qi <- subset(qi, SPP == "SALARC")

# Thereis only one data point for 2001, remove that year from the dataset
qi <- qi[qi$Year > 2001,]

# PREPARE DATA FOR JAGS IMPORT
#qi$Plot.Num <- qi$Plot.ID
#qi$Plot.ID <- recode(qi$Plot.ID, F1 = "1", F2 = "2", F3 = "3", F4 = "4", F5 = "5", F6 = "6", F7 = "7", F8 = "8", F9 = "9", F10 = "10", M1 = "11", M2 = "12", M3 = "13", M4 = "14", M5 = "15", M6 = "16", M7 = "17", M8 = "18", M9 = "19", M10 = "20") #order doesn't matter here
#qi$PlotNum <- as.numeric(as.character(qi$Plot.ID))
qi$PlotNum <- as.numeric(as.factor(qi$Plot.ID))
qi <- qi[order(qi$Year, qi$PlotNum),]

jags.dat <- list(
  lim = cbind(qi$P5_before,qi$P5), #changed P5_after to P5 as this is actually the end of the window
  y1s = ifelse(qi$P5_before<qi$P5,1,0), #column of 1's indicates iterval censoring
  year = qi$Year - min(qi$Year) + 1, #year must start with 1 (this becomes necessary in the more complex models)
  nyear = length(unique(qi$Year)),
  plot = qi$PlotNum,
  n = nrow(qi),
  nplot = length(unique(qi$PlotNum)),
  xhat = c(sort(unique(qi$Year - min(qi$Year) + 1))),
  nxhat = length(unique(qi$Year - min(qi$Year) + 1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      # PRIORS
      a ~ dunif(0, 365)
      b ~ dnorm(0, 0.001)
      for (i in 1:nplot){
      aplot[i] ~ dnorm(0, tau.plot)
      }
      sigma.plot ~ dunif(0, 100)
      tau.plot <- 1/(sigma.plot*sigma.plot)
      sigma.year ~ dunif(0, 100)
      tau.year <- 1/(sigma.year*sigma.year)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
      #LIKELIHOOD
      for (j in 1:nyear){
      a.year[j]~dnorm(muyear[j], tau.year)
      muyear[j] <- a+b*j
      }
      for (i in 1:n){
      y1s[i]~dinterval(t[i],lim[i,])
      t[i]~dnorm(pred[i],tau)
      pred[i] <- aplot[plot[i]]+a.year[year[i]]
      }
      # DERIVED QUANTITIES
      for (j in 1:nxhat){
      phat[j] <- a + b*xhat[j]
      }
      }
      ","models/int_cens_P5.jags")

# INITIAL VALUES
inits <- function() list(aplot = rnorm(jags.dat$nplot,0,2), sigma.plot = runif(1,0,1), sigma = runif(1,0,1), t = as.vector(apply(jags.dat$lim, 1, mean)))

# PARAMETERS TO MONITOR
params <- c("a", "b", "phat","aplot", "sigma", "sigma.plot", "sigma.year") ## ADD PARAMETERS TO MONITOR

modoutP5 <- jags(jags.dat, inits, params, model.file="models/int_cens_P5.jags", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, DIC=FALSE, working.directory=NULL, progress.bar = "text")

#plot(modoutP5) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP5$BUGSoutput$summary[,c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')])
hist(coefs$Rhat)
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH RESULTS
# Graph the trends over time for each species. Include both the raw data and the modeled predictions.
coefs.predP5 <- coefs[coefs$Type=="phat",] #change this name depending on what you called it

# Add species and year information
coefs.predP5$SPP <- qi$SPP[1]
coefs.predP5$YearC <- c(1:16) #change this if using all years!
coefs.predP5$Year <- coefs.predP5$YearC + min(qi$Year) - 1

#Save subset of data and coefs for later use
qiP5 <- qi
coefsP5 <- coefs

# 5% and 95% quantiles
# P5quants <- rbind(round(quantile(modoutP5$BUGSoutput$sims.list$a, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$b, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$sigma.year, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$sigma.plot, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$sigma, probs = c(0.05, 0.95)), 2))
# rownames(P5quants) <- c("a", "b", "sigma.year", "sigma.plot", "sigma")
# P5quants

# Growing Season Length

# DIFFERENCE between P2 and P5 for Salix arctica only
# length of time between green-up and senescence (Salix only)
# look at P6 instead? (date of last leaf turning yellow instead of first)
# NB: There is VIRTUALLY NO data for senescense in 2001 (1 data point) so the model below only runs from 2002 onwards
# Could also model this DIRECTLY - i.e. P5-P2 per plot per year - tried this but it doesn't work with censoring
# Need to add in right- and left- censoring

# Keep only observations for which our variables of interest ("P2" and "P5") were observed
qiphenP2 <- qiphen[!is.na(qiphen$P2) & qiphen$SPP=="SALARC",]
qiphenP5 <- qiphen[!is.na(qiphen$P5) & qiphen$SPP=="SALARC",]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an even was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qiphenP2[,c("P2_before","P2")])

# Make sure that P2_after is always a larger number than P2_before, or else JAGS will be unhappy (NA's not accepted)
which(qiphenP2$P2_before>qiphenP2$P2)
which(qiphenP5$P5_before>qiphenP5$P5)
qiphenP2[is.na(qiphenP2$P2_before) & !is.na(qiphenP2$P2),]
qiphenP5[is.na(qiphenP5$P5_before) & !is.na(qiphenP5$P5),]

# As there is virtually no data for senescense in 2001, remove the year from the datasets
qiphenP2gs <- qiphenP2 #[qiphenP2$Year > 2001,] # no longer removing 2001 for P2 because we want it to match with other P2 analysis!
qiphenP5gs <- qiphenP5[qiphenP5$Year > 2001,]

# PREPARE DATA FOR JAGS IMPORT

#qiphenP2gs$Plot.Num <- qiphenP2gs$Plot.ID
#qiphenP2gs$Plot.ID <- recode(qiphenP2gs$Plot.ID, F1 = "1", F2 = "2", F3 = "3", F4 = "4", F5 = "5", F6 = "6", F7 = "7", F8 = "8", F9 = "9", F10 = "10", M1 = "11", M2 = "12", M3 = "13", M4 = "14", M5 = "15", M6 = "16", M7 = "17", M8 = "18", M9 = "19", M10 = "20")
#qiphenP2gs$PlotNum <- as.numeric(as.character(qiphenP2gs$Plot.ID))
qiphenP2gs <- qiphenP2gs[order(qiphenP2gs$Year, qiphenP2gs$Plot.ID),]
qiphenP2gs$PlotNum <- as.numeric(factor(qiphenP2gs$Plot.ID, levels = unique(qiphenP2gs$Plot.ID)))

#qiphenP5gs$Plot.Num <- qiphenP5gs$Plot.ID
#qiphenP5gs$Plot.ID <- recode(qiphenP5gs$Plot.ID, F1 = "1", F2 = "2", F3 = "3", F4 = "4", F5 = "5", F6 = "6", F7 = "7", F8 = "8", F9 = "9", F10 = "10", M1 = "11", M2 = "12", M3 = "13", M4 = "14", M5 = "15", M6 = "16", M7 = "17", M8 = "18", M9 = "19", M10 = "20")
#qiphenP5gs$PlotNum <- as.numeric(as.character(qiphenP5gs$Plot.ID))
qiphenP5gs <- qiphenP5gs[order(qiphenP5gs$Year, qiphenP5gs$Plot.ID),]
qiphenP5gs$PlotNum <- as.numeric(factor(qiphenP5gs$Plot.ID, levels = unique(qiphenP5gs$Plot.ID)))

jags.dat <- list(
  limP2=cbind(qiphenP2gs$P2_before,qiphenP2gs$P2),
  limP5=cbind(qiphenP5gs$P5_before,qiphenP5gs$P5),
  y1sP2=rep(1,length(qiphenP2gs$P2_before)), #column of 1's indicates iterval censoring
  y1sP5=rep(1,length(qiphenP5gs$P5_before)),
  yearP2=qiphenP2gs$Year-min(qiphenP2gs$Year)+1, #year must start with 1
  yearP5=qiphenP5gs$Year-min(qiphenP5gs$Year)+1,
  nyearP2=length(unique(qiphenP2gs$Year)),
  nyearP5=length(unique(qiphenP5gs$Year)),
  plotP2=qiphenP2gs$PlotNum,
  plotP5=qiphenP5gs$PlotNum,
  nP2=nrow(qiphenP2gs),
  nP5=nrow(qiphenP5gs),
  nplotP2=length(unique(qiphenP2gs$PlotNum)),
  nplotP5=length(unique(qiphenP5gs$PlotNum)),
  xhatP2=c(sort(unique(qiphenP2gs$Year-min(qiphenP2gs$Year)+1))),
  xhatP5=c(sort(unique(qiphenP5gs$Year-min(qiphenP5gs$Year)+1))),
  nxhatP2=length(unique(qiphenP2gs$Year-min(qiphenP2gs$Year)+1)),
  nxhatP5=length(unique(qiphenP5gs$Year-min(qiphenP5gs$Year)+1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      
      #PRIORS P2
      aP2~dunif(0,365)
      bP2~dnorm(0,0.001)
      #Priors for plot random effect
      for (i in 1:nplotP2){
      aplotP2[i]~dnorm(0,tau.plotP2)
      }
      sigma.plotP2~dunif(0,100)
      tau.plotP2 <- 1/(sigma.plotP2*sigma.plotP2)
      sigmaP2~dunif(0,100)
      tauP2 <- 1/(sigmaP2*sigmaP2)
      sigma.yearP2~dunif(0,100)
      tau.yearP2 <- 1/(sigma.yearP2*sigma.yearP2)
      
      #PRIORS P5
      aP5~dunif(0,365)
      bP5~dnorm(0,0.001)
      #Priors for plot random effect
      for (i in 1:nplotP5){
      aplotP5[i]~dnorm(0,tau.plotP5)
      }
      sigma.plotP5~dunif(0,100)
      tau.plotP5 <- 1/(sigma.plotP5*sigma.plotP5)
      sigmaP5~dunif(0,100)
      tauP5 <- 1/(sigmaP5*sigmaP5)
      sigma.yearP5~dunif(0,100)
      tau.yearP5 <- 1/(sigma.yearP5*sigma.yearP5)
      
      #LIKELIHOOD P2
      for (k in 1:nyearP2){
      #this adds in the random variation
      a.yearP2[k]~dnorm(muyearP2[k], tau.yearP2)
      #this is the predicted mean effect per year
      muyearP2[k] <- aP2+bP2*k
      }
      for (i in 1:nP2){
      y1sP2[i]~dinterval(tP2[i],limP2[i,])
      #y1s is a column of 1's (because t[i] always falls between lim[i,1] and lim[i,2])
      tP2[i]~dnorm(predP2[i],tauP2)
      predP2[i] <- aplotP2[plotP2[i]]+a.yearP2[yearP2[i]]
      }
      
      #LIKELIHOOD P5
      for (k in 1:nyearP5){
      #this adds in the random variation
      a.yearP5[k]~dnorm(muyearP5[k], tau.yearP5)
      #this is the predicted mean effect per year
      muyearP5[k] <- aP5+bP5*k
      }
      for (i in 1:nP5){
      y1sP5[i]~dinterval(tP5[i],limP5[i,])
      #y1s is a column of 1's (because t[i] always falls between lim[i,1] and lim[i,2])
      tP5[i]~dnorm(predP5[i],tauP5)
      predP5[i] <- aplotP5[plotP5[i]]+a.yearP5[yearP5[i]]
      }
      # DERIVED QUANTITIES
      for (j in 1:nxhatP2){
      yhatP2[j] <- aP2+bP2*xhatP2[j]
      }
      for (j in 1:nxhatP5){
      yhatP5[j] <- aP5+bP5*xhatP5[j]
      }
      slopeDiff <- bP5-bP2
      for (j in 1:nxhatP5){
      yhatDiff[j] <- yhatP5[j]-yhatP2[j+1] #add 1 to yhatP2 because it starts one year earlier
      }
      }
      ","models/growing_season_pheno.jags")

# INITIAL VALUES

inits <- function() list(aplotP2=rnorm(jags.dat$nplotP2,0,2), tP2=as.vector(apply(jags.dat$limP2,1,mean)), sigma.plotP2=runif(1,0,1), sigma.yearP2=runif(1,0,1), sigmaP2=runif(1,0,1), aplotP5=rnorm(jags.dat$nplotP2,0,2), tP5=as.vector(apply(jags.dat$limP5,1,mean)), sigma.plotP5=runif(1,0,1), sigma.yearP5=runif(1,0,1), sigmaP5=runif(1,0,1))

# PARAMETERS TO MONITOR

params <- c("aplotP2","aplotP5","sigma.plotP2","sigma.plotP5","aP2","aP5","sigmaP2","sigmaP5","bP2","bP5","yhatP2","yhatP5","sigma.yearP2","sigma.yearP5","slopeDiff","yhatDiff")

modoutP2P5 <- jags(jags.dat,inits, params, model.file="models/growing_season_pheno.jags", n.chains=3,n.iter=30000,n.burnin=15000, n.thin=2, DIC=FALSE, working.directory=NULL, progress.bar = "text")

#plot(modoutP2P5) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP2P5$BUGSoutput$summary)
hist(coefs$Rhat)
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH DIFFERENCE BETWEEN P2 and P5
coefs.diff <- coefs[coefs$Type=="yhatDiff",]
coefs.diff$YearC <- 1:16
coefs.diff$Year <- coefs.diff$YearC+min(qiphenP2gs$Year)-1+1 # +1 b/c started at year 2

# Save coefs for later reference
coefsP2P5.all <- coefs
coefsP2P5 <- coefs.diff

# Get probability that the change in growing season length is greater than 0
slopeDiff <- modoutP2P5$BUGSoutput$sims.list$slopeDiff
pnorm(0, mean(slopeDiff), sd(slopeDiff), lower.tail=F)
#probability that a normally distributed random number with this mean and sd will be greater than 0

# # 5% and 95% quantiles
# P2P5quants <- round(quantile(modoutP2P5$BUGSoutput$sims.list$slopeDiff, probs = c(0.05, 0.95)), 2)
# P2P5quants

# Phenology figures ----

phenology_plot <- ggplot() +
  geom_point(data = qiP2, aes(x = Year, y = as.vector(apply(cbind(qiP2$P2_before, qiP2$P2), 1, mean))), 
             colour = "#007052", alpha = 0.5, size = 3) +
  geom_ribbon(data = coefs.predP2, aes(x = Year, ymin = `2.5%`, ymax = `97.5%`), fill = "#007052",
              alpha = 0.5, show.legend=FALSE) +
  geom_line(data = coefs.predP2, aes(x = Year, y = mean), colour = "#007052", size = 1, show.legend=FALSE, linetype = "dashed") +
  geom_point(data = qiphen[qiphen$SPP=="SALARC" & !is.na(qiphen$P5),], 
             aes(x = Year, y = as.vector(apply(cbind(qiphen$P5before[qiphen$SPP=="SALARC" & !is.na(qiphen$P5)],
                                                     qiphen$P5[qiphen$SPP=="SALARC" & !is.na(qiphen$P5)]), 1, mean))), colour = "#e8a600", alpha = 0.5, size = 3) +
  geom_ribbon(data = coefs.predP5, aes(x = Year, ymin = `2.5%`, ymax = `97.5%`), fill = "#e8a600", alpha = 0.5, 
              show.legend = FALSE) +
  geom_line(data = coefs.predP5, aes(x = Year, y = mean), colour = "#e8a600", size = 1, show.legend = FALSE, linetype = "dashed") +
  scale_x_continuous(breaks = c(2001, 2006, 2011, 2017)) +
  scale_y_continuous(breaks = c(140, 180, 220, 260, 300)) +
  coord_cartesian(ylim = c(130, 310)) +
  labs(x = " ", y = "Day of Year\n") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"))

ggsave("plots/Qikiqtarkuk_phenology.png", plot = phenology_plot, width = 7, height = 5, units = "in", dpi = 300)

# Satellite phenology analyses ----

greenup <- load("data/ITEX_greenup.RData")
qikiq_greenup <- greenup.all %>% filter(site == "i36")
kanger_greenup <- greenup.all %>% filter(site == "i20")

# Make data.frames for plotting
greenup.greenup <- qikiq_greenup %>% dplyr::select(site, year, greenup.date.05, greenup.date.50, greenup.date.95)
greenup.max <- qikiq_greenup %>% dplyr::select(site, year, ndvi.date.max)
greenup.scen <- qikiq_greenup %>% dplyr::select(site, year, senescence.date.05, senescence.date.50, senescence.date.95)
greenup.GSL <- qikiq_greenup %>% dplyr::select(site, year, greenup.date.50, senescence.date.50, gs.length.50)
# greenup.modVal <- qikiq_greenup %>% dplyr::select(site, year, modVal) %>% group_by(site, year) %>% mutate(DOY = list(modVal[[1]][,1]), modVal = list(modVal[[1]][,2])) %>% unnest()
# greenup.modVal <-  merge(qikiq_greenup, greenup.modVal)

GSstart <- greenup.greenup %>% summarise(mean(greenup.date.05))
GSend <- greenup.scen %>% summarise(mean(senescence.date.95))

# Growing season length - by site - greenup.date.05 and senescence.date.05
phenex_plot <- ggplot(greenup.GSL) +
  geom_point(aes(x = year, y = greenup.date.50), colour = "#007052", alpha = 0.5, size = 3) +
  geom_point(aes(x = year, y = senescence.date.50), colour = "#e8a600", alpha = 0.5, size = 3) +
  geom_smooth(method=lm, aes(x = year, y = greenup.date.50), colour = "#007052", fill = "#007052", alpha = 0.5, show.legend=F, linetype = "dotted") + 
  geom_smooth(method=lm, aes(x = year, y = senescence.date.50), colour = "#e8a600", fill = "#e8a600", alpha = 0.5, show.legend=F, linetype = "dotted") + 
  ylab("MODISv6 NDVI green up dates\n") +
  xlab("") +
  scale_x_continuous(breaks = c(2001, 2006, 2011, 2017)) +
  scale_y_continuous(breaks = c(140, 180, 220, 260, 300)) +
  coord_cartesian(ylim = c(130, 310)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"))

ggsave("plots/Qikiqtarkuk_phenex.png", plot = phenex_plot, width = 7, height = 5, units = "in", dpi = 300)

# Figure of all satellite data

load("data/MODIS6_ShrubHub_ITEX.RData")
MODIS <- MODIS %>% filter(id == "i36") %>% mutate(ndvi = (NDVI/10000)) %>% mutate(ndvi2 = replace(ndvi, ndvi <  0, 0)) %>% dplyr::select(DayOfYear, year, ndvi, ndvi2)
MODIS2 <- MODIS 

MODIS_zero <- read.csv("data/MODIS_zero_winter.csv")
MODIS_zero <- MODIS_zero %>% dplyr::select(DayOfYear, year, ndvi, ndvi2)
MODIS_zero$year <- as.character(MODIS_zero$year)
MODIS_zero$ndvi <- as.numeric(as.character(MODIS_zero$ndvi))
MODIS_zero$ndvi2 <- as.numeric(as.character(MODIS_zero$ndvi2))

MODIS3 <- rbind(MODIS2, MODIS_zero)
MODIS3 <- MODIS3 %>% na.omit()

phenex_plot2 <- ggplot(MODIS3) + 
  geom_smooth(aes(x=DayOfYear, y=ndvi2, colour = year), alpha = 0.5, size = 0.5, se = F, show.legend=T) +
  geom_rect(mapping=aes(xmin=100, xmax=300, ymin=-0.1, ymax=0.01), color="white", fill="white") +
  geom_point(aes(x=DayOfYear, y=ndvi2, colour = year), alpha = 0.5, size = 2) +
  geom_vline(xintercept = 160, colour = "#007052", size = 1, linetype = "dashed") +
  geom_vline(xintercept = 220, colour = "#e8a600", size = 1, linetype = "dashed") +
  geom_vline(xintercept = 152, colour = "#007052", size = 1, linetype = "dotted") +
  geom_vline(xintercept = 254, colour = "#e8a600", size = 1, linetype = "dotted") +
  ylab("NDVI\n") +
  xlab("Day of Year") +
  scale_color_viridis(discrete=T) +
  #annotate("text", label = ("A. MODISv6 Greenup Dates"), x = 2004, y = 200, size = 5) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.1, 1)) +
  scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350), limits = c(100, 300)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.25, "cm"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.title=element_blank())

ggsave("plots/Qikiqtarkuk_phenex2.png", plot = phenex_plot2, width = 8, height = 5, units = "in", dpi = 300)
