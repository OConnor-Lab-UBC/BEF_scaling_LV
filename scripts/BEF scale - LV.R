library(tidyverse)
library(broom)
library(ggExtra)
library(igraph)
library(vegan)
library(synchrony)
library(viridis)
library(minpack.lm)
library(cowplot)
#library(drc)

meta_dyn_model_loss_gain<-function(disp = 0.01, type = "competitive", spatial_env = TRUE, temporal_env = FALSE, loss=TRUE, global_disp=FALSE, env_gamma = 1){
  patches<-80
  species<-100
  interval<-150
  Tmax<-interval*species
  
  Env_perform<-function(env,z,zmax=z+5,sig_p){
    Tmat<-matrix(rep(env,each=species),species,length(env)) 
    wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
    #wT2<-1-((Tmat-z)/(z-zmax))^2
    #wT[Tmat>=z]<-wT2[Tmat>=z]
    wT[wT<0]<-0
    wT<-wT-1
    return(wT)
  }
  
  #environmental trait
  z<-seq(-0.5,1.5,length=species)#runif(n = species, min = -0.5,max = 1.5)
  sig_p<-0.05 #rate of performance decay as z differs from local environment
  envPeriod<-5000 # temporal period of environmental sin wave
  
  n <- patches # number of points you want on the unit circle
  landscape <- 1+(data.frame(t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n)))))+1)/2*999
  names(landscape)<-c("x","y")
  
  distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
  
  distance_mat<-1*(round(distance_mat1)==round(unique(c(distance_mat1))[order(unique(c(distance_mat1)))[2]]))
  diag(distance_mat)<-0
  connections<-distance_mat
  
  graph<-as.undirected(graph.adjacency(distance_mat))
  graph<-set.vertex.attribute(graph,"x coordinate",value=landscape$x)
  graph<-set.vertex.attribute(graph,"y coordinate",value=landscape$y)
  graph<-set.edge.attribute(graph,"weight",value=distance_mat1[cbind(as.numeric(get.edgelist(graph)[,1]),  as.numeric(get.edgelist(graph)[,2]))])
  
  disp_mat <- connections*0.5
  
  if(spatial_env == TRUE){
    env.df<-data.frame(env = phase.partnered(n = 80,gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1], patch = 1:patches)
    plot(env.df$env, type ='l')
    
    
    #env <- 0.5 * (sin(seq(0, spatial_period*pi, length.out = patches) + (2 * pi * 1) / envPeriod) + 1) # local environmental conditions at time t
    #plot(env, type="l")
    
  } else {
    env.df<-data.frame(env = 0.5, patch = 1:patches)
  }
  A<-t(Env_perform(env.df$env,z,sig_p = sig_p)*2000)
  matplot(A, type = 'l')
  
  r<-0.3 #intrinsic rate of increase
  
  bdiag1<- -.2
  if(type == "neutral"){
    B=matrix(-0.2,nrow=species,ncol=species)
    N<-sapply(X = 1:species,FUN = function(X){
      hold<-rep(0,patches)
      hold[sample(1:patches,size = 5,replace=F)]<-5
      return(hold)})
  }
  
  if(type == "competitive"){
    B<-matrix(runif(species*patches)*-0.05,nrow=species,ncol=species)#matrix(runif(species*patches)*-0.15,nrow=species,ncol=species)
    N<-matrix(1,ncol=species,nrow=patches)
  }
  
  if(type == "priority"){
    B<-matrix(rnorm(species*species,mean = bdiag1+0.01,sd=0.005),nrow=species,ncol=species)
    N<-sapply(X = 1:species,FUN = function(X){
      hold<-rep(0,patches)
      hold[sample(1:patches,size = 5,replace=F)]<-5
      return(hold)})
  }
  diag(B)<-bdiag1
  
  colonizeOrder<-sample(1:species,size = species,replace = F)
  
  if(loss==TRUE){
    N<-matrix(1,ncol=species,nrow=patches)
  } else{
    N<-matrix(0,ncol=species,nrow=patches)
    N[,colonizeOrder[1]]<-1
  }
  
  Nsave<-array(NA,dim=c(patches,species,species))
  Esave<-matrix(NA,nrow = Tmax, ncol = patches)
  
  colV<-seq(from = interval,to=Tmax,by=interval)
  
  pb <- txtProgressBar(min = 0, max = Tmax, style = 3)
  for(i in 1:Tmax){
    if(sum(i==colV)==1){
      if(loss==TRUE){
        N[,colonizeOrder[which(colV==i)+1]]<-0
      } else{
        N[,colonizeOrder[which(colV==i)+1]]<-1
      }
      Nsave[,,which(colV==i)]<-N
    }
    # if(temporal_env==TRUE){
    #   if(spatial_env==TRUE){
    #     env <- 0.5 * (sin(seq(0, 2*pi, length.out = patches) + (2 * pi * i) / envPeriod) + 1) # local environmental conditions at time t
    #   } else{
    #     env <- 0.5 * (sin(rep(0,patches) + (2 * pi * i) / envPeriod) + 1) # local environmental conditions at time t
    #     
    #   }
    #   
    #   A<-t(Env_perform(env,z,sig_p = sig_p)*2000)
    # }
    
    Esave[i,]<-env.df$env
    
    
    
    
    
    Nt<-N*exp(r+N%*%B+A)
    Nt<-Nt+disp*disp_mat%*%Nt-Nt*disp
    
    #Nt<-Nt+rnorm(n = patches*species,mean = 0,sd=Nt*0.01)*(Nt>0)
    Nt[Nt<0.05]<-0
    N<-Nt
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(Nsave,Esave))  
}

meta_dyn_model_sp_pool<-function(disp = 0.01, type = "competitive", spatial_env = TRUE, temporal_env = FALSE, env_gamma = 1){
  patches<-80
  species<-100
  interval<-150
  Tmax<-interval
  
  Env_perform<-function(env,z,zmax=z+5,sig_p){
    Tmat<-matrix(rep(env,each=species),species,length(env)) 
    wT<-exp(-((Tmat-z)/2*rep(sig_p,each=patches))^2)
    #wT2<-1-((Tmat-z)/(z-zmax))^2
    #wT[Tmat>=z]<-wT2[Tmat>=z]
    wT[wT<0]<-0
    wT<-wT-1
    return(wT)
  }
  
  #environmental trait
  z<-seq(-0.2,1.2,length=species)#runif(n = species, min = -0.5,max = 1.5)
  sig_p<-0.05 #rate of performance decay as z differs from local environment
  envPeriod<-5000 # temporal period of environmental sin wave
  
  n <- patches # number of points you want on the unit circle
  landscape <- 1+(data.frame(t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n)))))+1)/2*999
  names(landscape)<-c("x","y")
  
  distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
  
  distance_mat<-1*(round(distance_mat1)==round(unique(c(distance_mat1))[order(unique(c(distance_mat1)))[2]]))
  diag(distance_mat)<-0
  connections<-distance_mat
  
  graph<-as.undirected(graph.adjacency(distance_mat))
  graph<-set.vertex.attribute(graph,"x coordinate",value=landscape$x)
  graph<-set.vertex.attribute(graph,"y coordinate",value=landscape$y)
  graph<-set.edge.attribute(graph,"weight",value=distance_mat1[cbind(as.numeric(get.edgelist(graph)[,1]),  as.numeric(get.edgelist(graph)[,2]))])
  
  disp_mat <- connections*0.5
  
  if(spatial_env == TRUE){
    env.df<-data.frame(env = phase.partnered(n = 80,gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1], patch = 1:patches)
    #plot(env.df$env, type ='l')
    
    
    #env <- 0.5 * (sin(seq(0, spatial_period*pi, length.out = patches) + (2 * pi * 1) / envPeriod) + 1) # local environmental conditions at time t
    #plot(env, type="l")
    
  } else {
    env.df<-data.frame(env = 0.5, patch = 1:patches)
  }
  A<-t(Env_perform(env.df$env,z,sig_p = sig_p)*2000)
  #matplot(A, type = 'l')
  
  r<-0.3 #intrinsic rate of increase
  
  bdiag1<- -.2
  if(type == "neutral"){
    B=matrix(-0.2,nrow=species,ncol=species)
    N<-sapply(X = 1:species,FUN = function(X){
      hold<-rep(0,patches)
      hold[sample(1:patches,size = 5,replace=F)]<-5
      return(hold)})
  }
  
  if(type == "competitive"){
    B<-matrix(runif(species*patches)*-0.05,nrow=species,ncol=species)#matrix(runif(species*patches)*-0.15,nrow=species,ncol=species)
    N<-matrix(1,ncol=species,nrow=patches)
  }
  
  if(type == "priority"){
    B<-matrix(rnorm(species*species,mean = bdiag1+0.01,sd=0.005),nrow=species,ncol=species)
    N<-sapply(X = 1:species,FUN = function(X){
      hold<-rep(0,patches)
      hold[sample(1:patches,size = 5,replace=F)]<-5
      return(hold)})
  }
  diag(B)<-bdiag1
  
  
  Nsave<-array(NA,dim=c(patches,species,species))
  Esave<-matrix(NA,nrow = Tmax, ncol = patches)
  
  pb <- txtProgressBar(min = 0, max = species, style = 3)
  for(j in 1:species){
    N<-matrix(0,ncol=species,nrow=patches)
    N[,sample(1:species,size = j,replace = F)]<-1
    
    for(i in 1:Tmax){
      Esave[i,]<-env.df$env
      
      Nt<-N*exp(r+N%*%B+A)
      Nt<-Nt+disp*disp_mat%*%Nt-Nt*disp
      
      #Nt<-Nt+rnorm(n = patches*species,mean = 0,sd=Nt*0.01)*(Nt>0)
      Nt[Nt<0.05]<-0
      N<-Nt
      
    }
    Nsave[,,j]<-N
    setTxtProgressBar(pb, j)
  }
  close(pb)
  return(list(Nsave,Esave))  
}

all.data.df<-coef.df<-data.frame()

for(gamma in c(-1,0,2)){
  for(r in 1:20){
    print(paste("rep = ",r,"; gamma = ", gamma, sep = ""))
    #hold<-meta_dyn_model_loss_gain(spatial_env = TRUE, env_gamma = gamma,disp = 0)
    hold<-meta_dyn_model_sp_pool(spatial_env = TRUE, env_gamma = gamma, disp = 0.001)
    #hold<-meta_dyn_model_loss_gain(spatial_env = TRUE, env_gamma = gamma, disp = 0.001)
    
    hold.data.run<-data.frame()
    coef.df_run<-data.frame()
    for(k in 1:80){
      hold.data<-data.frame()
      for(j in 1:100){
        if(k == 1){
          div<-renyi(hold[[1]][1:k,,j],scales = 0:1,hill = TRUE)
        } else{
        div<-renyi(colSums(hold[[1]][1:k,,j]),scales = 0:1,hill = TRUE)}
        bmass<-sum(hold[[1]][1:k,,j])
        hold.data<-bind_rows(hold.data,data.frame(SR = div[1], D = div[2],bmass = bmass, scale = k, sp.pool = j))
        hold.data$gamma<-gamma
        hold.data$rep<-r
        }
      BEF.nl<-coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2]
      BEF.pl<-coef(lm(log(bmass)~log(SR), data = hold.data[hold.data$SR>0,]))[2]
      plot(bmass/k ~SR, data = hold.data, pch=19, main = k, xlim=c(0,100), ylim = c(0,12))
      lines(predict(nlsLM(formula = bmass/k ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)),newdata =data.frame(SR =1:max(hold.data$SR))), col = 2)
      abline(v = BEF.nl, col = 3)
      
      
      coef.df_run<-bind_rows(coef.df_run, data.frame(BEF.nl = BEF.nl, BEF.pl = BEF.pl, scale = k, gamma = gamma, rep = r))
      all.data.df<-bind_rows(all.data.df,hold.data)
      hold.data.run<-bind_rows(hold.data.run, hold.data)
    }
    
    coef.df_run<-left_join(coef.df_run,hold.data.run %>% 
                             group_by(sp.pool) %>% 
                             mutate(alpha_SR = min(SR)) %>% 
                             mutate(beta = SR/alpha_SR), by = c("scale", "gamma", "rep"))
    
    coef.df<-bind_rows(coef.df, coef.df_run)
  }
}

save(coef.df, file = "./data/LV_data.RData")

means<-coef.df %>% 
  group_by(gamma, scale) %>% 
  summarise(lower = quantile(BEF.nl, probs = 0.25), upper = quantile(BEF.nl, probs = 0.75), BEF.nl = mean(BEF.nl))

Fig.a<- means%>% 
  ggplot(aes(x=scale,y=BEF.nl, color = gamma, group = gamma, fill = gamma))+
  geom_hline(yintercept = filter(means, scale == 1)$BEF.nl, color = c("dodgerblue", "grey", "red"), lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "noise")+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("spatial scale (# of local patches)")+
  ylab("BEF half saturation richness")

Fig.b<-coef.df %>% 
  group_by(gamma, scale) %>% 
  summarise(lower = quantile(beta, probs = 0.25,na.rm = T), upper = quantile(beta, probs = 0.75,na.rm=T), beta = median(beta,na.rm = T)) %>% 
  ggplot(aes(x=scale,y=beta, color = gamma, group = gamma, fill = gamma))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "noise")+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("spatial scale (# of local patches)")+
  ylab("spatial beta diversity")

plot_grid(Fig.a, Fig.b, labels = c("a", "b"))
ggsave("./figures/Spatial BEF scaling.png", height = 4, width = 9)
