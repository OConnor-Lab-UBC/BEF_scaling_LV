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

meta_dyn_model_sp_pool<-function(disp = 0, type = "competitive", spatial_env = TRUE, temporal_env = FALSE, env_gamma = 1){
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
    env.df<-data.frame(env = phase.partnered(n = patches,gamma = env_gamma, mu = 0.5, sigma = 0.25)$timeseries[,1], patch = 1:patches)
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
    seed.sp <- sample(1:species,size = j,replace = F)
    N[,seed.sp]<-1
    
    for(i in 1:Tmax){
      Esave[i,]<-env.df$env
      
      #N[,seed.sp][N[,seed.sp]==0] <- 0.049
      
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
hold.raw.data <-data.frame()
for(r in 1:20){
  for(gamma in c(0,1,2)){
    print(paste("rep = ",r,"; gamma = ", gamma, sep = ""))
    #hold<-meta_dyn_model_loss_gain(spatial_env = TRUE, env_gamma = gamma,disp = 0)
    hold<-meta_dyn_model_sp_pool(spatial_env = TRUE, env_gamma = gamma, disp = 0)
    #hold<-meta_dyn_model_loss_gain(spatial_env = TRUE, env_gamma = gamma, disp = 0.001)
    
    hold.data.run<-data.frame()
    coef.df_run<-data.frame()
    for(k in c(1:30,32,34,36,38,40,45,50,55,60,65,70,75,80)){
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
      #plot(bmass/k ~SR, data = hold.data, pch=19, main = k, xlim=c(0,100), ylim = c(0,12))
      hold.raw.data <- rbind(hold.raw.data, data.frame(rep = r, scale = k, SR = hold.data$SR, bmass = hold.data$bmass, gamma = gamma))
      #lines(predict(nlsLM(formula = bmass/k ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)),newdata =data.frame(SR =1:max(hold.data$SR))), col = 2)
      #abline(v = BEF.nl, col = 3)
      
      
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

save(coef.df, hold.raw.data, file = "./data/LV_data.RData")


b.df <- data.frame()
for(i in 1:4){
  b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
           data = filter(hold.raw.data, rep == 2, gamma == 0, scale == i), 
           start = c(a = 20, b = 10)))[2], 
           scale = i))
}

A<- hold.raw.data %>% 
  filter(rep == 2, gamma == 0) %>%
  filter(scale %in% c(1,2,3,4)) %>% 
  ggplot(aes(x = SR, y = bmass/scale, fill = factor(scale), group = scale))+
  scale_fill_viridis_d(end = 0.8, option = "B", name = "spatial\nscale")+
  ylab("average biomass")+
  xlab("species richness")+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  geom_vline(data = b.df, aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  theme(legend.justification=c(1,0), legend.position=c(1,0))

B<- ggplot(b.df, aes(x = scale, y = b))+
  scale_fill_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  ylab("BEF half saturation richness")+
  xlab("scale")+
  geom_line()+
  geom_point(pch = 21, size = 4, aes(fill = factor(scale)))
plot_grid(A, B)
ggsave("./figures/spatial_raw_BEF_example.pdf", height = 4.5, width = 9)

hold.raw.data %>% 
  filter(rep == 1, gamma == 0) %>%
  filter(scale %in% c(1)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = factor(scale), group = scale))+
  scale_fill_viridis_d(end = 0.8, option = "B", name = "spatial\nscale")+
  ylab("community biomass")+
  xlab("species richness")+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  geom_vline(data = filter(b.df,scale == 1), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)
ggsave("./figures/spatial_raw_BEF_example_1.pdf", height = 5, width = 6)

b.df <- data.frame()
for(i in c(1,2,3,5,10,20,30,40,60,80)){
  for(g_sel in c(0,1,2)){
  hold.data = filter(hold.raw.data, rep == 1, gamma == g_sel, scale == i)
  b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
                                                data = hold.data, 
                                                start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2], 
                                 scale = i,
                                 gamma = g_sel))
}}

hold.raw.data %>% 
  filter(rep == 2) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = factor(scale), group = scale))+
  scale_fill_viridis_d(end = 0.8, option = "B", name = "spatial\nscale")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  geom_vline(data = b.df, aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  facet_wrap(~gamma)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)
ggsave("./figures/spatial_raw_BEF.png", height = 4*1.5, width = 10*1.5)


means<-coef.df %>% 
  group_by(gamma, scale) %>% 
  summarise(lower = quantile(BEF.nl, probs = 0.25), upper = quantile(BEF.nl, probs = 0.75), BEF.nl = median(BEF.nl))

Fig.a<- means%>% 
  ggplot(aes(x=scale,y=BEF.nl, color = gamma, group = gamma, fill = gamma))+
  geom_hline(yintercept = filter(means, scale == 1)$BEF.nl, color = c("grey","#ED7F64","red"), lty = 2)+
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

plot_grid(Fig.a, Fig.b, labels = c("a)", "b)"))
ggsave("./figures/Spatial BEF scaling.png", height = 4, width = 9)
