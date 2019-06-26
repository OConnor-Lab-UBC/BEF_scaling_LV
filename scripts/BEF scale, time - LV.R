library(synchrony)
library(vegan)
library(viridis)
library(tidyverse)
library(minpack.lm)
library(cowplot)

meta_dyn_model_sp_pool_time<-function(time_gamma = 2, r = 0.3){
  patches<-1
  species<-100
  interval<-1000
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
  sig_p<- 0.1#0.05 #rate of performance decay as z differs from local environment
  
  temp.env<-phase.partnered(n = Tmax,gamma = time_gamma, mu = 0, sigma = 0.25)$timeseries[,1]
  
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
  
  
  Nsave<-array(NA,dim=c(species,Tmax,species))
  Esave<-matrix(NA,nrow = Tmax, ncol = patches)
  
  pb <- txtProgressBar(min = 0, max = species, style = 3)
  for(j in 1:species){
    N<-matrix(0,ncol=species,nrow=patches)
    seed.sp<-sample(1:species,size = j,replace = F)
    N[,seed.sp]<-1
    
    for(i in 1:Tmax){
      env_t<-temp.env[i]+0.5
      Esave[i,]<-env_t
      
      A<-t(Env_perform(env_t,z,sig_p = sig_p)*2000)
      #if(i %in% seq(10, Tmax, by = 5)){
      N[,seed.sp]<-N[,seed.sp]+0.049
      #}
      
      Nt<-N*exp(r+N%*%B+A)
      
      #Nt<-Nt+rnorm(n = patches*species,mean = 0,sd=Nt*0.01)*(Nt>0)
      Nt[Nt<0.05]<-0
      N<-Nt
      Nsave[,i,j]<-N
    }
    setTxtProgressBar(pb, j)
  }
  close(pb)
  return(list(Nsave,Esave))  
}

coef.df<-data.frame()
hold.raw.data <-data.frame()
for(r in 1:20){
  for(t_g in c(0,1,2)){
    print(paste("rep = ",r,"; time gamma = ", t_g, sep = ""))
    
    hold<- meta_dyn_model_sp_pool_time(time_gamma = t_g)
    
    hold.data.run<-data.frame()
    coef.df_run<-data.frame()
    for(i in 1:500){ #time window
      hold.data<-data.frame()
      for(j in 1:100){ #species pool
        if(i == 1){
          div<-renyi(t(hold[[1]][,500:(i+500-1),j]),scales = 0:1,hill = TRUE)
        } else{
          div<-renyi(colSums(t(hold[[1]][,500:(i+500-1),j])),scales = 0:1,hill = TRUE)}
        bmass<-sum(hold[[1]][,50:(i+50-1),j])
        hold.data<-bind_rows(hold.data,data.frame(SR = div[1], D = div[2],bmass = bmass, t_scale = i, sp.pool = j))
        hold.data$time_gamma<-t_g
        hold.data$rep<-r
      }
      
      hold.data.run<-bind_rows(hold.data.run, hold.data)
      
      BEF.nl<-coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2]
      plot(bmass/i ~SR, data = hold.data, pch=19, main = i, xlim=c(0,70), ylim = c(0,12))
      lines(predict(nlsLM(formula = bmass/i ~ (a * SR)/ (SR + b), data = hold.data,start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)),newdata =data.frame(SR =1:max(hold.data$SR))), col = 2)
      abline(v = BEF.nl, col = 3)
      hold.raw.data <- rbind(hold.raw.data, data.frame(rep = r, t_scale = i, SR = hold.data$SR, bmass = hold.data$bmass, gamma = t_g))
      
      coef.df_run<-bind_rows(coef.df_run, data.frame(BEF.nl = BEF.nl, t_scale = i, time_gamma = t_g, rep = r))
    }
    coef.df_run<-left_join(coef.df_run,hold.data.run %>% 
                             group_by(sp.pool) %>% 
                             mutate(alpha_SR = min(SR)) %>% 
                             mutate(beta = SR/alpha_SR), by = c("t_scale", "time_gamma", "rep"))
    
    
    coef.df<-bind_rows(coef.df, coef.df_run)
  }
}

save(coef.df, hold.raw.data, file = "./data/LV_time.RData")

hold.raw.data %>% 
  filter(rep == 6) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass, color = factor(t_scale), group = t_scale))+
  geom_point()+
  scale_color_viridis_d()+
  facet_wrap(~gamma)+
  ylab("average biomass")#+
# geom_smooth(method = "nls", formula = bmass/scale ~ a * SR / (b + SR), start = c(a = 20, b = 10))
ggsave("./figures/temporal_raw_BEF.pdf", height = 4, width = 10)

b.df <- data.frame()
for(i in c(1,2,3,5,10,20,40,100,200,500)){
  for(g_sel in c(0,1,2)){
  hold.data = filter(hold.raw.data, rep == 1, gamma == g_sel, t_scale == i)
  b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
                                                data = hold.data, 
                                                start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2], 
                t_scale = i,
                gamma = g_sel))
}}

hold.raw.data %>% 
  filter(rep == 1) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = factor(t_scale), group = t_scale))+
  scale_fill_viridis_d(end = 0.8, option = "B", name = "spatial\nscale")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  #geom_smooth(method = "nls", 
  #            formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  geom_vline(data = b.df, aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  facet_wrap(~gamma)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)
ggsave("./figures/temporal_raw_BEF.png", height = 4*1.5, width = 10*1.5)


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

plot_grid(Fig.a, Fig.b, labels = c("a)", "b)"))
ggsave("./figures/Temporal BEF scaling.png", height = 4, width = 9)
