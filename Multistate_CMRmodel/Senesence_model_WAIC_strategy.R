####### Model flamant Nimble ###
rm(list=ls())

setwd("//hpcm.cluster.calcul.hpc/Shared-2/ROQUES/Flamingo/Senescence_3states/WAIC/")

library(nimble)

load('CH_95FR_winterAD')
load('strategy_95FR_based_winterAD')

CH<-CH_FR




nind<-dim(CH)[1]
K<-dim(CH)[2]
ngroup<-length(unique(strategy))

get.first <- function(x) min(which(x!=0)) 
f <- apply(CH, 1, get.first)

CH[CH==0]<-4 ## recoded non seen

age<-matrix(NA,dim(CH)[1], dim(CH)[2])
## age
for ( i in 1:nind){
  for (t in f[i]:K){
    
    age[i,t]<-t-f[i]+1
    
  }
}


maxage<-max(age, na.rm=T)

ageeff<-age/100

# --------------------------------------------------------------

# Let's define the model. To do so, some notation first:

# OBSERVATIONS

# 0 = Non detected at the colony
# 1 = detected at the colony
# 2 = detected incubating at the colony
# 3 = detected with a chick at the colony

# STATES

# 1 = alive non-breeder
# 2 = alive failed breeder at the egg stage 
# 3 = alive successful breeder with one fledged chick
# 4 = Dead


# p : detection prob.
# -------------------------------------
# pNB : detection prob. of non-breeder
# pB : detection prob. of breeder




# phi: survival prob.
# -------------------------------------
# phiNB : survival prob. of non-breeder
# phiB  : survival prob. of breeder


# b : breeding prob.
# -------------------------------------
# 1-bB : breeding prob. of non-breeder
# bB  : bredding prob. of breeder


# g : successful breeding prob.
# -------------------------------------
# 1-gB: chick hatching prob. of non-breeder
# gB  : chick hatching prob. of breeder


code <- nimbleCode({
  
  
  for (t in 1:K){
    
    # phiB[t]<-mean.phiB
    # phiNB[t]<-mean.phiNB
    
    # gB[t]<-mean.gB
    # gNB[t]<-mean.gNB
    # 
    pB[t]~dunif(0,1)
    pNB[t]~dunif(0,1)
    
    # bB[t]<-mean.bB
    # bNB[t]<-mean.bNB
  }
  
  
  for (i in 1:nind){
    for (t in f[i]:K){
      logit(phi[i,t])<-mu[strategy[i]]+beta[strategy[i]]*ageeff[i,t]+beta2[strategy[i]]*ageeff[i,t]^2+eps[i]
      
      logit(B[i,t])<-muB[strategy[i]]+betaB[strategy[i]]*ageeff[i,t]+beta2B[strategy[i]]*ageeff[i,t]^2+epsB[i]
      
      logit(gB[i,t])<-mugB[strategy[i]]+betagB[strategy[i]]*ageeff[i,t]+beta2gB[strategy[i]]*ageeff[i,t]^2+epsgB[i]
    }
  }
  
  
  # for (i in 1:nind){
  #   for (t in f[i]:f[i]+2){
  #     phi[i,t]<-1
  #     
  #     B[i,t]<-0
  #     
  #     gB[i,t]<-0
  #   }
  # }
  
  
  for (i in 1:nind){
    eps[i]~dnorm(0, sd=sd)
    epsB[i]~dnorm(0, sd=sdB)
    epsgB[i]~dnorm(0, sd=sdgB)
  }
  for (g in 1:ngroup){
    
    mu[g]~dnorm(0, sd=10)
    beta[g]~dnorm(0, sd=10)
    beta2[g]~dnorm(0, sd=10)
    
    muB[g]~dnorm(0, sd=10)
    betaB[g]~dnorm(0, sd=10)
    beta2B[g]~dnorm(0, sd=10)
    
    mugB[g]~dnorm(0, sd=10)
    betagB[g]~dnorm(0, sd=10)
    beta2gB[g]~dnorm(0, sd=10)
  }
  
  sd~dunif(0,10)
  sdB~dunif(0,10)
  sdgB~dunif(0,10)
  
  # mean.phiB~dunif(0,1)
  # mean.phiNB~dunif(0,1)
  
  # mean.bB~dunif(0,1)
  # mean.bNB~dunif(0,1)
  # 
  # mean.gB~dunif(0,1)
  # mean.gNB~dunif(0,1)
  # 
  
  ### transition matrix ##
  
  for (i in 1:nind){
    for (t in f[i]:K){
      
      #step survive
      
      ps1[1,i,t,1]<-phi[i,t]
      ps1[1,i,t,2]<-0
      ps1[1,i,t,3]<-0
      ps1[1,i,t,4]<-1-phi[i,t]
      
      ps1[2,i,t,1]<-0
      ps1[2,i,t,2]<-phi[i,t]
      ps1[2,i,t,3]<-0
      ps1[2,i,t,4]<-1-phi[i,t]
      
      ps1[3,i,t,1]<-0
      ps1[3,i,t,2]<-0
      ps1[3,i,t,3]<-phi[i,t]
      ps1[3,i,t,4]<-1-phi[i,t]
      
      ps1[4,i,t,1]<-0
      ps1[4,i,t,2]<-0
      ps1[4,i,t,3]<-0
      ps1[4,i,t,4]<-1
      
      #step breeding
      
      ps2[1,i,t,1]<-1-B[i,t]
      ps2[1,i,t,2]<-B[i,t]
      ps2[1,i,t,3]<-0
     
      
      ps2[2,i,t,1]<-1-B[i,t]
      ps2[2,i,t,2]<-B[i,t]
      ps2[2,i,t,3]<-0

      
      ps2[3,i,t,1]<-1-B[i,t]
      ps2[3,i,t,2]<-B[i,t]
      ps2[3,i,t,3]<-0

      
      ps2[4,i,t,1]<-0
      ps2[4,i,t,2]<-0
      ps2[4,i,t,3]<-1

      
      # step breeding success
      
      ps3[1,i,t,1]<-1
      ps3[1,i,t,2]<-0
      ps3[1,i,t,3]<-0
      ps3[1,i,t,4]<-0
      
      ps3[2,i,t,1]<-0
      ps3[2,i,t,2]<-1-gB[i,t]
      ps3[2,i,t,3]<-gB[i,t]
      ps3[2,i,t,4]<-0
      
      ps3[3,i,t,1]<-0
      ps3[3,i,t,2]<-0
      ps3[3,i,t,3]<-0
      ps3[3,i,t,4]<-1
     
      
      ##observation
      
      po[1,i,t,1]<-pNB[t]
      po[1,i,t,2]<-0
      po[1,i,t,3]<-0
      po[1,i,t,4]<-1-pNB[t]
      
      po[2,i,t,1]<-0
      po[2,i,t,2]<-pB[t]
      po[2,i,t,3]<-0
      po[2,i,t,4]<-1-pB[t]
      
      po[3,i,t,1]<-0
      po[3,i,t,2]<-0#pB[t]*delta
      po[3,i,t,3]<-pB[t]#*(1-delta)
      po[3,i,t,4]<-1-pB[t]
      
      po[4,i,t,1]<-0
      po[4,i,t,2]<-0
      po[4,i,t,3]<-0
      po[4,i,t,4]<-1
      
    }
  }
  
  # form the matrix product
  #--------------------------
  for (i in 1:nind){
    for (t in f[i]:K) {
      
      ps[1:4,i,t,1:4] <- ps1[1:4,i,t,1:4] %*% ps2[1:4,i,t,1:3] %*% ps3[1:3,i,t,1:4]
    }
  } 
  
  
  for (i in 1:nind)  # for each ind
  {
    
    # First capture occasion
    
    z[i,f[i]] <- Y[i,f[i]]
    #----------------------------------
    
    #subsequent occasions
    for (t in (f[i]+1):K)  # loop over time
      
    {
      
      # draw states at t given states at t-1
      z[i,t] ~ dcat(ps[z[i,t-1],i,t-1,1:4])
      
      # draw observations at t given states at t
      Y[i,t] ~ dcat(po[z[i,t],i,t-1,1:4])
      
    }
    
  }
  
  
  
})


known.state.ms <- function(ms, notseen){ # notseen: label for 
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,]))) 
    state[i,m] <- NA
  }
  return(state) }


ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}


nind <- dim(CH)[1]
K <- dim(CH)[2]
z <- known.state.ms(CH,4)


constants <- list(nind = nind, K = K, f=f, strategy=strategy, ngroup=ngroup, age=age, maxage=maxage, ageeff=ageeff)#delta=0.08

data <- list(Y = CH, z=z)

inits <- list(z = ms.init.z(CH,f),
              pB=runif(K,0,1), pNB=runif(K,0,1),
              mu=rnorm(3,0,10),
              beta=rnorm(3,0,10),
              beta2=rnorm(3,0,10),
              muB=rnorm(3,0,10),
              betaB=rnorm(3,0,10),
              beta2B=rnorm(3,0,10),
              mugB=rnorm(3,0,10),
              betagB=rnorm(3,0,10),
              beta2gB=rnorm(3,0,10),
              eps=rnorm(nind,0,0.1),
              sd=runif(1,0,10), 
              epsB=rnorm(nind,0,0.1),
              sdB=runif(1,0,10), 
              epsgB=rnorm(nind,0,0.1),
              sdgB=runif(1,0,10))

Rmodel <- nimbleModel(code, constants, data, inits, calculate = F)

#Rmodel$calculate()

params <- c('eps','sd','beta','beta2', 'mu',
            'epsB','sdB','betaB','beta2B', 'muB',
            'epsgB','sdgB','betagB','beta2gB', 'mugB',
            'pB', 'pNB', 'z')

conf <- configureMCMC(Rmodel, monitors = params, enableWAIC = T)

#conf$printSamplers()
#conf$printMonitors()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

if(FALSE) {
  compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
  Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
}

ni <- 60000
nt <- 4
nb <- 30000
nc <- 2

set.seed(0)

samples <- runMCMC(Cmcmc, niter = ni, thin = nt, nburnin = nb, nchains = nc, WAIC = TRUE)
save(samples, file="95_breedFR_Senescence_strategyWAIC.Rdata")