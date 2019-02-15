library(R2jags)

#######################
# Data
#######################
mydata <- list(ns=3,y=cbind(c(2,3,4),c(5,6,7)),v=cbind(c(1.2,1.3,1.4),c(.4,1.5,1.6)),
               n=cbind(c(12,13,14),c(14,15,16)))
mydata$s.pooled <- with(mydata,sqrt(((n[,1]-1)*v[,1]+(n[,2]-1)*v[,2])/(n[,1]+n[,2]-2)))
#mydata$v.md <- with(mydata,(n[,1]+n[,2])*s.pooled/(n[,1]*n[,2]))


#######################
# Model 1: for MD
#######################
PMAcontsMD <- function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    y[i,1] ~ dnorm(mu[i,1],w[i,1]) #likelihood in one arm (why do we assume likelihood for each arm not for each )
    y[i,2] ~ dnorm(mu[i,2],w[i,2]) #likelihood in the other arm

    ## Parameterisation
    mu[i,1] <- u[i]
    mu[i,2] <- u[i] +theta[i]
    theta[i] ~dnorm(MD,prec)
    
    w[i,1]<- 1/v[i,1]
    w[i,2]<- 1/v[i,2]
    # MD <- mu[i,1]-mu[i,2]
    # SMD<- (y[i,1]-y[i,2])/s.pooled[i]

  }
  
  #prior distributions
  for (i in 1:ns) {u[i] ~ dnorm(0,.01)}
  tau ~ dunif(0,1)   #dnorm(0,100)%_%T(0,)                                 
  prec<- 1/pow(tau,2)
  MD ~ dnorm(0,100)

}


#######################
# run the model
#######################

PMAinJAGS.MD<- jags(mydata,inits = NULL,parameters.to.save = c('MD','tau'), n.chains = 2, n.iter = 10000,
                 n.burnin = 1000, DIC=F, model.file = PMAcontsMD)

#results
print(PMAinJAGS.MD)

#check chain mixing
traceplot(PMAinJAGS.MD)






#######################
# Model 2: for SMD
#######################
PMAcontsSMD <- function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    y[i,1] ~ dnorm(mu[i,1],w[i,1]) #likelihood in one arm (why do we assume likelihood for each arm not for each )
    y[i,2] ~ dnorm(mu[i,2],w[i,2]) #likelihood in the other arm
    
    ## Parameterisation
    mu[i,1] <- u[i]
    mu[i,2] <- u[i] +(delta[i]*s.pooled[i])
    
    delta[i]~dnorm(SMD,prec)
    w[i,1]<- 1/v[i,1]
    w[i,2]<- 1/v[i,2]
    # MD <- mu[i,1]-mu[i,2]
    # SMD<- (y[i,1]-y[i,2])/s.pooled[i]
    
  }
  
  #prior distributions
  for (i in 1:ns) {u[i] ~ dnorm(0,.01)}
  tau ~ dunif(0,1)   #dnorm(0,100)%_%T(0,)                                 
  prec<- 1/pow(tau,2)
  SMD ~ dnorm(0,100)
  
  ## Outcome of interest
}


#######################
# run the model
#######################

PMAinJAGS.SMD<- jags(mydata,inits = NULL,parameters.to.save = c('SMD','tau'), n.chains = 2, n.iter = 10000,
                 n.burnin = 1000, DIC=F, model.file = PMAcontsSMD)

#results
print(PMAinJAGS.SMD)

#check chain mixing
traceplot(PMAinJAGS.SMD)


