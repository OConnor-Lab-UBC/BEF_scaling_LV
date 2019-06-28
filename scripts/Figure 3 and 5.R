library(tidyverse)
library(broom)
library(ggExtra)
library(igraph)
library(vegan)
library(synchrony)
library(viridis)
library(minpack.lm)
library(cowplot)

load("./data/LV_data.RData")

means<-coef.df %>% 
  group_by(gamma, scale) %>% 
  summarise(lower = quantile(BEF.nl, probs = 0.25), upper = quantile(BEF.nl, probs = 0.75), BEF.nl = median(BEF.nl))

Fig5.a<- means%>% 
  ggplot(aes(x=scale,y=BEF.nl, color = gamma, group = gamma, fill = gamma))+
  geom_hline(yintercept = filter(means, scale == 1)$BEF.nl, color = c("grey","#ED7F64","red"), lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "gamma")+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("spatial scale (# of local patches)")+
  ylab(expression(paste("BEF half saturation richness (", italic(b[i]),")"), sep = ""))+
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  coord_cartesian(ylim = c(3,68))

Fig3.a<-coef.df %>% 
  group_by(gamma, scale) %>% 
  summarise(lower = quantile(beta, probs = 0.25,na.rm = T), upper = quantile(beta, probs = 0.75,na.rm=T), beta = median(beta,na.rm = T)) %>% 
  ggplot(aes(x=scale,y=beta, color = gamma, group = gamma, fill = gamma))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "gamma")+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("spatial scale (# of local patches)")+
  ylab(expression(paste("spatial ", beta, " diversity", sep = "")))+
  theme(legend.justification=c(0,1), legend.position=c(0,1))+
  ylim(c(1,3.8))

load("./data/LV_time.RData")

means<-coef.df %>% 
  group_by(time_gamma, t_scale) %>% 
  summarise(lower = quantile(BEF.nl, probs = 0.25), upper = quantile(BEF.nl, probs = 0.75), BEF.nl = median(BEF.nl))

Fig5.b<- means%>% 
  ggplot(aes(x=t_scale,y=BEF.nl, color = time_gamma, group = time_gamma, fill = time_gamma))+
  geom_hline(yintercept = filter(means, t_scale == 1)$BEF.nl, color = c("grey","#ED7F64","red"), lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "noise", guide = FALSE)+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab(expression(paste("BEF half saturation richness (", italic(b[i]),")"), sep = ""))+
  coord_cartesian(ylim = c(3,68))


Fig3.b<-coef.df %>% 
  group_by(time_gamma, t_scale) %>% 
  summarise(lower = quantile(beta, probs = 0.25,na.rm = T), upper = quantile(beta, probs = 0.75,na.rm=T), beta = median(beta,na.rm = T)) %>% 
  ggplot(aes(x=t_scale,y=beta, color = time_gamma, group = time_gamma, fill = time_gamma))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line()+
  scale_color_gradient2(low = "dodgerblue",mid = "grey",high = "red", name = "noise", guide = FALSE)+
  scale_fill_gradient2(low = "dodgerblue",mid = "grey",high = "red", guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab(expression(paste("temporal ", beta, " diversity", sep = "")))+
  ylim(c(1,3.8))

plot_grid(Fig3.a, Fig3.b, labels = c("a)", "b)"))
ggsave("./figures/Beta_diversity.png", height = 4, width = 9) 

plot_grid(Fig5.a, Fig5.b, labels = c("a)", "b)"))
ggsave("./figures/BEF_scaling.png", height = 4, width = 9) 
