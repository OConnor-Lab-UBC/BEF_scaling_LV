library(synchrony)

BEF_simulation <-  function(type = "spatial", env_gamma = 2, rep = 1, r_max = 5, sigma = 0.25, alpha_max = 0.25, max_scale = 80, alpha_ii_sd = 0){
  set.seed(rep)
  
  common_env <- phase.partnered(n = max_scale, gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1]
  repeat{
    burn_in_env <- phase.partnered(n = max_scale+1, gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1]
    burn_in_env_stnd <- c(burn_in_env - (burn_in_env[max_scale+1] - common_env[1]))[-(max_scale+1)]

    if(env_gamma == 0) {
      temporal_env <- c(burn_in_env[-c(max_scale+1)],common_env)
    }
    if(max(burn_in_env_stnd) < 1.3 & min(burn_in_env_stnd) > - 0.3) {break}
  }
  temporal_env <- c(burn_in_env_stnd,common_env)
  plot(temporal_env, type = "l")
  abline(v = max_scale, lty = 2)
  
  if(type == "spatial"){
    Tmax <- 150 #length of simulation
    patches <- max_scale
    #environment
    env.df <- data.frame(env = common_env, patch = 1:patches)
    #species responses to environmental variation
  }
  
  if(type == "temporal"){
    #Tmax <- (env_time+15)*2 #length of simulation
    Tmax <- length(temporal_env)
    patches <- 1
    #environment
    #env.df <- data.frame(env = c(rev(common_env),common_env), time = 1:Tmax)
    env.df <- data.frame(env = temporal_env, time = 1:Tmax)
  }
  
  #fixed parameters
  species <- 100
  alpha <- matrix(runif(n = species*species, min = 0, max = alpha_max), species, species) #per capita interspecific competition matrix
  diag(alpha) <- rnorm(n = species, mean = 1, sd = alpha_ii_sd) #per capita intraspecific competition
  z <- seq(-0.2, 1.2, length = species) #environmental optima
  ext_thresh <- 0.05 #extinction threshold
  reseed_value <- 0.035
  
  if(type == "spatial"){
    Nsave_all <- array(NA,dim=c(patches,species,species))
    r <- matrix(r_max*exp(-((rep(z, each = patches)-env.df$env)/(2*sigma))^2), patches, species)
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
