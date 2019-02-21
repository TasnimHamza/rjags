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

## 2. Read the model from NMRbinary.R file
source('NMRbinary.R')

## 3.Run network meta-regression for response versus year of randomisation

# a. transform the data into a list suitable for JAGS analysis
Depression$year_2005 <- year-2005
# Depression$study_year_2001 <- Depression$study_year-2001
# Depression1 <- Depression[complete.cases(Depression),]

NMRdataBinary=make.jagsNMA.data(studyid=studyID,t=drug_name,r=Responders,n=Ntotal,data=Depression,othervar = year_2005,type="binary",reference = "Agomelatine")
str(NMRdataBinary)

# b.run Jags and create a jags object
NMRinJAGSBin<- jags(data = NMRdataBinary, inits = NULL,
                    parameters.to.save = c("ORref","tau",'OR'), n.chains = 2, n.iter = 10000,
                    n.burnin = 5000,DIC=F,n.thin=10,
                    model.file = modelNMRBinary)

# These are our results
print(NMRinJAGSBin)
save(NMRinJAGSBin,file="NMAinJAGSall.RData",envir = .GlobalEnv)

############################
# CODA: Convergence analysis
############################

# 1. convergence to stationarity: Trace plot
traceplot(NMRinJAGSBin,varname="tau" )
traceplot(NMRinJAGSBin,varname="ORref" )
traceplot(NMRinJAGSBin,varname="OR" )

  #Traceplot for the ORref and OR dont look so nice, probaly bad mixing.

# 2. convergence to ergodic average:n.eff; effective sample size for estimating the mean
NMRinJAGSBin$BUGSoutput$summary[,'n.eff']

 ## small n.eff -> large ACF (autocorrelation) -> Bad mixing

# 3. Gelman diag plot and Rhat: If the chains converge then they should be overlapped 
    # after some iteration and so sqrt(Rhat) almost 1. 
Rhat <- NMRinJAGSBin$BUGSoutput$summary[,'Rhat'] 
sqrt(Rhat) < 1.01 ## So the chains are not converges to the same point, no convegence for the chains

#forestplot
y=as.vector(NMRinJAGSBin$BUGSoutput$mean$OR)
seES=as.vector(NMRinJAGSBin$BUGSoutput$sd$OR)
m1 <- metagen(y,seES)
forest(m1)
## i get weired results!

## estimate of common tau
NMRinJAGSBin$BUGSoutput$mean$tau


