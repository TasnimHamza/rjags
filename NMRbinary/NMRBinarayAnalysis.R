#*********************************************************************************
#             Load the libraries needed             
#*********************************************************************************
library(meta)
library(metafor)
library(netmeta)
library(readxl)

library(devtools)
install_github("esm-ispm-unibe-ch/NMAJags",force = TRUE)
library(NMAJags)
library(R2jags)
#*********************************************************************************
#             Network Meta-Analysis           
#*********************************************************************************


## Data
Depression <- read_excel("Depression.xlsx")
Depression <- as.data.frame(Depression)
str(Depression)
ls(package:NMAJags)

## 1. Create a random variable called year of randomisation
year <- sample(1900:2017,nrow(Depression),replace = T) ## could we use the year_pub but we should centrelize it first
str(Depression)
Depression$year2000 <- NA
Depression$year2000 <- unlist(sapply(unique(Depression$studyID), 
  function(x){Depression$year2000[Depression$studyID==x] <- rep(sample(1970:2017,1,replace = F),sum(Depression$studyID==x))}))
##!!!!!!!
##!!!!!!! what you created above is a variable with a different year per arm. this is impossible, a study should have the same year!!!!
## !!!! create as an exercise a random variable that has the same year for the arms of the same study using the apply functiona

## 2. Read the model from NMRbinary.R file
source('NMRbinary.R')

## 3.Run network meta-regression for response versus year of randomisation

# a. transform the data into a list suitable for JAGS analysis
Depression$year2000 <- Depression$year2000-2000
# Depression$study_year_2001 <- Depression$study_year-2001
# Depression1 <- Depression[complete.cases(Depression),]

NMRdataBinary=make.jagsNMA.data(studyid=studyID,t=drug_name,r=Responders,n=Ntotal,data=Depression,othervar = year2000,type="binary",reference = "Agomelatine")
str(NMRdataBinary)

##!!!!!!!
# because of the thinning you have very few samples left, increase the n.iter and the burn in and use also parallele jags to increase speed
# b.run Jags and create a jags object
NMRinJAGSBin<- jags(data = NMRdataBinary, inits = NULL,
                    parameters.to.save = c("ORref","tau",'LOR'), n.chains = 2, n.iter = 10000,
                    n.burnin = 5000,DIC=F,n.thin=10,
                    model.file = modelNMRBinary)

NMRinJAGSBinP <- jags.parallel(data = NMRdataBinary, inits = NULL,
              parameters.to.save = c("ORref","tau",'beta'), n.chains = 2, n.iter = 100000,
              n.burnin = 10000,DIC=F,n.thin=10,
              model.file = modelNMRBinary)
# These are our results
print(NMRinJAGSBinP)
save(NMRinJAGSBinP,file="NMAinJAGSall.RData",envir = .GlobalEnv)
load("NMAinJAGSall.RData")
############################
# CODA: Convergence analysis
############################

# 1. convergence to stationarity: Trace plot

traceplot(NMRinJAGSBinP,varname="tau" ) 
##!!!!!!! see in the traceplot there are few samples and early bad mixing
##+++++++++++++ Now, after increase n.iter and n.burnin the traceplot looks better 

traceplot(NMRinJAGSBinP,varname="ORref" )
traceplot(NMRinJAGSBinP,varname="LOR" )

  #Traceplot for the ORref and OR dont look so nice, probaly bad mixing.

# 2. convergence to ergodic average:n.eff; effective sample size for estimating the mean
n.eff <- NMRinJAGSBinP$BUGSoutput$summary[,'n.eff']

 ## small n.eff -> large ACF (autocorrelation) -> Bad mixing

# 3. Gelman diag plot and Rhat: If the chains converge then they should be overlapped 
    # after some iteration and so sqrt(Rhat) almost 1. 
Rhat <- NMRinJAGSBinP$BUGSoutput$summary[,'Rhat'] 
sqrt(Rhat) < 1.01 ## So the chains are not converges to the same point, no convegence for the chains

#forestplot
y=as.vector(NMRinJAGSBinP$BUGSoutput$mean$ORref)
ES <- log(y)
seES=as.vector(NMRinJAGSBinP$BUGSoutput$sd$ORref)
m1 <- metagen(ES[-17],seES[-17])
forest(m1)
## i get weired results! ====> 
##!!!!!!! you need to plot only the logOR against the reference, otherwise you try to plot ALL pairwise ORs which is impossible!
###!!!!! also, your seES refer to the logOR, so you cannot mix them with OR - you need to set ES=logOR and seES=se(logOR)

##+++++++++++++++ After modifying the code accordingly, for the drug 17 I get an incredible value SE, so I remove it from the plot.
##+++++++++++++++ Do you have any idea why this could happen? or just by chance!
ES[17]
seES[17]

## I am trying to see how the data look like for drug17 vs the reference drug
Depression$Responders[Depression$drug_name=="Vortioxetine"]  ## n.events for using Vortioxetine drug.
Depression$studyID[Depression$drug_name=="Vortioxetine"]     ## studyID for the drug17.
Depression$Ntotal[Depression$drug_name=="Vortioxetine"]      ## Ntotal for the drug17.
Depression$drug_name[Depression$studyID==516]      
## relative treatment effect for drug17 is computed indirectly, so it was difficult to check why we got this result!

# the reference drug
Depression$Responders[Depression$drug_name=="Agomelatine"]  ## n.events for using Vortioxetine drug.
Depression$studyID[Depression$drug_name=="Agomelatine"]     ## studyID for the reference drug.
Depression$Ntotal[Depression$drug_name=="Agomelatine"]      ## Ntotal for the reference drug.

## I still dont now why I get huge SE!

## estimate of common tau
NMRinJAGSBin$BUGSoutput$mean$tau
m1$tau

##!!!!!!! you can compare the tau from the meta-regression to the tau from the meta-analysis to see the importance of the covariate
### you need to monitor the coefficients, plot them and interpret them 

# plot
traceplot(NMRinJAGSBinP,varname="beta" )

# Interpretation
# beta is the change in the relative treatment effect (logOR) due to varying relative publication year 
     ## (related to the year 2000).



