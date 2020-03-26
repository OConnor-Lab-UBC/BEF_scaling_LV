library(synchrony)

BEF_simulation <- function(type = "spatial", env_gamma = 2){
  #fixed parameters
  species <- 100
  sigma <- 0.25 #environmental niche breadth
  r_max <- 5 #max growth rate
  alpha <- matrix(runif(n = species*species, min = 0, max = 0.25), species, species) #per capita interspecific competition matrix
  diag(alpha) <- 1 #per capita intraspecific competition
  z <- seq(-0.2, 1.2, length = species) #environmental optima
  ext_thresh <- 0.05 #extinction threshold
  reseed_value <- 0.035
  
  if(type == "spatial"){
    Tmax <- 150 #length of simulation
    patches <- 80
    #environment
    env.df<-data.frame(env = phase.partnered(n = patches, gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1], patch = 1:patches)
    #species responses to environmental variation
    r <- matrix(r_max*exp(-((rep(z, each = patches)-env.df$env)/(2*sigma))^2), patches, species)
  }
  
  
  if(type == "temporal"){
    env_time <- 80
    Tmax <- env_time*2 #length of simulation
    patches <- 1
    #environment
    env.hold<-phase.partnered(n = env_time, gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1]
    env.df <- data.frame(env = c(env.hold,rev(env.hold)), time = 1:Tmax)
  }
  
  if(type == "spatial"){
    Nsave_all <- array(NA,dim=c(patches,species,species))
  }
  if(type == "temporal"){
    Nsave_all <- array(NA,dim=c(Tmax,species,species))
  }
  pb <- txtProgressBar(min = 0, max = species, style = 3)
  for(j in 1:species){
    N<-matrix(0,ncol=species,nrow=patches) #empty matrix of abundances
    seed.sp <- sample(1:species,size = j,replace = F) #choose which species to include
    N[,seed.sp]<-1 #seed initial abundances
    
    #community dynamics
    Nsave<-array(NA,dim=c(patches,species,Tmax))
    Esave<-matrix(NA,nrow = Tmax, ncol = patches)
    
    for(i in 1:Tmax){
      if(type == "temporal"){
        r <- r_max*exp(-((z-env.df$env[i])/(2*sigma))^2)
      } 
      if (type == "spatial") {
        Esave[i,]<-env.df$env
      }
      N[,seed.sp][N[,seed.sp]==0] <- reseed_value
      N <- N*r/(1+N%*%alpha) 
      N[N<ext_thresh] <- 0 #apply extinction threshold
      Nsave[,,i]<-N
    }
    if(type == "spatial"){
      Nsave_all[,,j] <- N
    }
    if(type == "temporal"){
      Esave[,1]<-env.df$env
      Nsave_all[,,j] <- t(Nsave[1,,])
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)
  return(list(Nsave_all,Esave))
}