library(R2jags)
help(jags)


#######################
#this is the data
#######################
mydata = list(ns=3,r=cbind(c(2,3,4),c(4,5,6)),n=cbind(c(12,13,14),c(14,15,16)))

#######################
#then make the model
#######################
PMAbinary=function() {
  
  for(i in 1:ns) { 
    
    #likelihood
    r[i,1] ~ dbin(p[i,1],n[i,1])#likelihood in one arm
    r[i,2] ~ dbin(p[i,2],n[i,2])#likelihood in the other arm
    
    #parametrisation          
    logit(p[i,1])<- u[i]
    logit(p[i,2])<- u[i] + theta[i] # OR should I use logit(p[i,2])<-lpi2 and  lpi2 ~ sum(u[i] + theta[i])
    theta[i] ~ dnorm(mean,prec)
    }
  
    #prior distributions
      for (i in 1:ns) {u[i] ~ dnorm(0,.01)}
      tau ~ dunif(0,1)   #dnorm(0,100)%_%T(0,)                                 
      prec<- 1/pow(tau,2)
      mean ~ dnorm(0,100)
}#end of model

#######################
# initial values
#######################

initialval = NULL
#initialval = list(list(tau=0.2,mean=0.3))

#######################
# run the model
#######################

PMAinJAGS<- jags(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 2, n.iter = 10000, n.burnin = 1000, DIC=F, model.file = PMAbinary)

#results
print(PMAinJAGS)

#check chain mixing
traceplot(PMAinJAGS)

#what else is there
names(PMAinJAGS)
names(PMAinJAGS$BUGSoutput)
PMAinJAGS$BUGSoutput$summary


#///////////////////////////

#######################
# parallelising
#######################

start_time <- Sys.time()
PMAinJAGS<- jags(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 4, n.iter = 100000, n.burnin = 1000, DIC=F, model.file = PMAbinary)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
PMAinJAGS<- jags.parallel(mydata,initialval,parameters.to.save = c("tau","mean"), n.chains = 4, n.iter = 100000, n.burnin = 1000, DIC=F, model.file = PMAbinary)
end_time <- Sys.time()
end_time - start_time




