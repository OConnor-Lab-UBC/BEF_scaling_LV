library(synchrony)
library(vegan)
library(viridis)
library(tidyverse)
library(minpack.lm)
library(cowplot)

meta_dyn_model_sp_pool_time<-function(disp = 0.1, spatial_env = TRUE, temporal_env = TRUE, space_gamma = 0, time_gamma = 2){
  patches<-80
  species<-100
  interval<-150
  Tmax<-interval
  
  type= "competitive"
  
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
  
  n <- patches # number of points you want on the unit circle
  landscape <- 1+(data.frame(t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n)))))+1)/2*999
  names(landscape)<-c("x","y")
  
  distance_mat1<-as.matrix(dist(landscape,method = "euclidean",diag = T,upper=T))
  
  distance_mat<-1*(round(distance_mat1)==round(unique(c(distance_mat1))[order(unique(c(distance_mat1)))[2]]))
  diag(distance_mat)<-0
  connections<-distance_mat
  
  disp_mat <- connections*0.5
  
  if(spatial_env == TRUE){
    env.df<-data.frame(env = phase.partnered(n = 80,gamma = space_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1], patch = 1:patches)
  } else {
    env.df<-data.frame(env = 0.5, patch = 1:patches)
  }
  
  if(temporal_env == TRUE){
    temp.env<-phase.partnered(n = Tmax,gamma = time_gamma, mu = 0, sigma = 0.25)$timeseries[,1]
  } else {
    temp.env<-rep(0, Tmax)
  }
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
  
  
  Nsave<-array(NA,dim=c(patches,species,Tmax,species))
  Esave<-matrix(NA,nrow = Tmax, ncol = patches)
  
  pb <- txtProgressBar(min = 0, max = species, style = 3)
  for(j in 1:species){
    N<-matrix(0,ncol=species,nrow=patches)
    seed.sp<-sample(1:species,size = j,replace = F)
    N[,seed.sp]<-1
    
    for(i in 1:Tmax){
      env_t<-env.df$env+temp.env[i]
      Esave[i,]<-env_t
      
      A<-t(Env_perform(env_t,z,sig_p = sig_p)*2000)
      N[,seed.sp]<-N[,seed.sp]+0.049
      
      Nt<-N*exp(r+N%*%B+A)
      Nt<-Nt+disp*disp_mat%*%Nt-Nt*disp
      
      #Nt<-Nt+rnorm(n = patches*species,mean = 0,sd=Nt*0.01)*(Nt>0)
      Nt[Nt<0.05]<-0
      N<-Nt
      Nsave[,,i,j]<-N
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)
  return(list(Nsave,Esave))  
}

coef.df<-data.frame()
for(r in 1:20){
  for(t_g in c(-1,0,2)){
    print(paste("rep = ",r,"; time gamma = ", t_g, sep = ""))
    
    hold<- meta_dyn_model_sp_pool_time(disp = 0.01,spatial_env = T,temporal_env = T,space_gamma = 0,time_gamma = t_g)
    
    hold.data.run<-data.frame()
    coef.df_run<-data.frame()
    for(i in 1:100){ #time window
      hold.data<-data.frame()
      for(j in 1:100){ #species pool
        if(i == 1){
          div<-renyi(t(hold[[1]][40,,50:(i+50-1),j]),scales = 0:1,hill = TRUE)
        } else{
          div<-renyi(colSums(t(hold[[1]][40,,50:(i+50-1),j])),scales = 0:1,hill = TRUE)}
        bmass<-sum(hold[[1]][40,,50:(i+50-1),j])
        hold.data<-bind_rows(hold.data,data.frame(SR = div[1], D = div[2],bmass = bmass, t_scale = i, sp.pool = j))
        hold.data$time_gamma<-t_g
        hold.data$space_gamma<-space_gamma
        hold.data$rep<-r
      }
      
      hold.data.run<-bind_rows(hold.data.run, hold.data)
      
      BEF.nl<-coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2]
      BEF.pl<-coef(lm(log(bmass)~log(SR), data = hold.data[hold.data$SR>0,]))[2]
      plot(bmass/i ~SR, data = hold.data, pch=19, main = i, xlim=c(0,70), ylim = c(0,12))
      lines(predict(nlsLM(formula = bmass/i ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)),newdata =data.frame(SR =1:max(hold.data$SR))), col = 2)
      abline(v = BEF.nl, col = 3)
      
      
      coef.df_run<-bind_rows(coef.df_run, data.frame(BEF.nl = BEF.nl, BEF.pl = BEF.pl, t_scale = i, time_gamma = t_g, space_gamma = space_gamma, rep = r))
    }
    coef.df_run<-left_join(coef.df_run,hold.data.run %>% 
      group_by(sp.pool) %>% 
      mutate(alpha_SR = min(SR)) %>% 
      mutate(beta = SR/alpha_SR), by = c("t_scale", "time_gamma", "space_gamma", "rep"))
    
    
    coef.df<-bind_rows(coef.df, coef.df_run)
  }
}


means<-coef.df %>% 
  group_by(time_gamma, t_scale) %>% 
  summarise(lower = quantile(BEF.nl, probs = 0.25), upper = quantile(BEF.nl, probs = 0.75), BEF.nl = mean(BEF.nl))

Fig.a<- means%>% 
ggplot(aes(x=t_scale,y=BEF.nl, color = time_gamma, group = time_gamma, fill = time_gamma))+
  geom_hline(yintercept = filter(means, t_scale == 1)$BEF.nl, color = c("dodgerblue", "grey", "red"), lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "noise")+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab("BEF half saturation richness")

Fig.b<-coef.df %>% 
  group_by(time_gamma, t_scale) %>% 
  summarise(lower = quantile(beta, probs = 0.25,na.rm = T), upper = quantile(beta, probs = 0.75,na.rm=T), beta = median(beta,na.rm = T)) %>% 
  ggplot(aes(x=t_scale,y=beta, color = time_gamma, group = time_gamma, fill = time_gamma))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "noise")+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab("temporal beta diversity")

plot_grid(Fig.a, Fig.b, labels = c("a", "b"))
ggsave("./figures/Temporal BEF scaling.png", height = 4, width = 9)
