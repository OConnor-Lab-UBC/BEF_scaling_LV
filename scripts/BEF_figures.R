library(tidyverse)
library(cowplot)
library(vegan)
library(minpack.lm)

source("./scripts/BEF_scale_functions.R")

colV <- c("grey40", "deeppink1", "red3")

#Figure 2#####
sp_0 <- BEF_simulation(type = "spatial" ,env_gamma = 0)
env_sp0<-(decostand(c(sp_0[[2]][1,],-0.2,1.2),method = "range")[1:80])*100

sp_1 <- BEF_simulation(type = "spatial" ,env_gamma = 1)
env_sp1 <- (decostand(c(sp_1[[2]][1,],-0.2,1.2),method = "range")[1:80])*100

sp_2 <- BEF_simulation(type = "spatial" ,env_gamma = 2)
env_sp2 <- (decostand(c(sp_2[[2]][1,],-0.2,1.2),method = "range")[1:80])*100

tp_0 <- BEF_simulation(type = "temporal" ,env_gamma = 0)
env_tp0 <- (decostand(c(tp_0[[2]][81:160],-0.2,1.2),method = "range")[1:80])*100

tp_1 <- BEF_simulation(type = "temporal" ,env_gamma = 1)
env_tp1 <- (decostand(c(tp_1[[2]][81:160],-0.2,1.2),method = "range")[1:80])*100

tp_2 <- BEF_simulation(type = "temporal" ,env_gamma = 2)
env_tp2 <- (decostand(c(tp_2[[2]][81:160],-0.2,1.2),method = "range")[1:80])*100

Fig2a<- data.frame(N = c(sp_0[[1]][,,100]), species = rep(1:100, each = 80), patch = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = patch, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F, breaks = c(0.05,1,1.4))+
  geom_line(data = data.frame(species = env_sp0, patch = 1:80), aes(size = NULL),color = "red", size = 0.5, alpha = 1)+
  ggtitle("gamma = 0")+
  theme_classic()+
  ylim(0,100)

Fig2b<- data.frame(N = c(sp_1[[1]][,,100]), species = rep(1:100, each = 80), patch = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = patch, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F, breaks = c(0.05,1,1.4))+
  geom_line(data = data.frame(species = env_sp1, patch = 1:80), aes(size = NULL),color = "red", size = 0.5, alpha = 1)+
  ggtitle("gamma = 1")+
  theme_classic()+
  ylim(0,100)

Fig2c<- data.frame(N = c(sp_2[[1]][,,100]), species = rep(1:100, each = 80), patch = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = patch, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F, breaks = c(0.05,1,1.4))+
  geom_line(data = data.frame(species = env_sp2, patch = 1:80), aes(size = NULL),color = "red", size = 0.5, alpha = 1)+
  ggtitle("gamma = 2")+
  theme_classic()+
  ylim(0,100)

Fig2d<- data.frame(N = c(tp_0[[1]][81:160,,100]), species = rep(1:100, each = 80), time = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = time, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F, breaks = c(0.05,1,1.4))+
  geom_line(data = data.frame(species = env_tp0, time = 1:80), aes(size = NULL),color = "red", size = 0.5, alpha = 1)+
  ggtitle("gamma = 0")+
  theme_classic()+
  ylim(0,100)

Fig2e<- data.frame(N = c(tp_1[[1]][81:160,,100]), species = rep(1:100, each = 80), time = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = time, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F, breaks = c(0.05,1,1.4))+
  geom_line(data = data.frame(species = env_tp1, time = 1:80), aes(size = NULL),color = "red", size = 0.5, alpha = 1)+
  ggtitle("gamma = 1")+
  theme_classic()+
  ylim(0,100)

Fig2f<- data.frame(N = c(tp_2[[1]][81:160,,100]), species = rep(1:100, each = 80), time = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = time, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F, breaks = c(0.05,1,1.4))+
  geom_line(data = data.frame(species = env_tp2, time = 1:80), aes(size = NULL),color = "red", size = 0.5, alpha = 1)+
  ggtitle("gamma = 2")+
  theme_classic()+
  ylim(0,100)

plot_grid(Fig2a, Fig2b, Fig2c, Fig2d, Fig2e, Fig2f, nrow = 2, labels = "AUTO")
ggsave("./figures/Fig.2.png", width = 10*1.2, height = 6*1.2)
ggsave("./figures/Fig.2.pdf", width = 10*1.2, height = 6*1.2)

#run simulations####
reps <- 100
hold.raw.data <-data.frame()
coef.df <- data.frame()
for(r in 1:reps){
  for(gamma in c(0,1,2)){
    print(paste("rep = ",r,"; gamma = ", gamma, sep = ""))
    spatial_run <- BEF_simulation(type = "spatial", env_gamma = gamma, rep = r)
    spatial_com <- spatial_run[[1]]
    temporal_run <- BEF_simulation(type = "temporal", env_gamma = gamma, rep = r)
    temporal_com <- temporal_run[[1]][81:160,,]
    
    hold.data.run<-data.frame()
    start_id <- 1
    for(k in c(1:30,32,34,36,38,40,45,50,55,60,65,70,75,80)){
      time_seq <- 80:(80-k+1)
      hold.data<-data.frame()
      for(j in 1:100){
        if(k == 1){
          div_space<-renyi(spatial_com[time_seq,,j],scales = 0:1,hill = TRUE)
          div_time<-renyi(temporal_com[time_seq,,j],scales = 0:1,hill = TRUE)
          bmass_time<-sum(temporal_com[time_seq,,j])
          bmass_space<-sum(spatial_com[time_seq,,j])
        } else{
          div_space<-renyi(colSums(spatial_com[time_seq,,j]),scales = 0:1,hill = TRUE)
          div_time<-renyi(colSums(temporal_com[time_seq,,j]),scales = 0:1,hill = TRUE)
          bmass_space<-mean(rowSums(spatial_com[time_seq,,j]))
          bmass_time<-mean(rowSums(temporal_com[time_seq,,j]))
        }

        hold.data<-bind_rows(hold.data, data.frame(SR = c(div_space[1], div_time[1]), D = c(div_space[2], div_time[2]), bmass = c(bmass_space, bmass_time), scale = k, sp.pool = j, scenario = c("space", "time")))
        hold.data$gamma<-gamma
        hold.data$rep<-r
      }
      coef.hold_space <- coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = filter(hold.data, scenario == "space"),start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))
      coef.hold_time <- coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = filter(hold.data, scenario == "time"), start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))
      coef.hold.df <- bind_rows(coef.hold_space,coef.hold_time)
      coef.hold.df$scenario = c("space", "time")
      coef.hold.df$rep <- r
      coef.hold.df$scale <- k
      coef.hold.df$gamma <- gamma
      coef.df<-bind_rows(coef.df, coef.hold.df)
      
      hold.raw.data <- rbind(hold.raw.data, hold.data)
    }
  }
}
#save(hold.raw.data, coef.df, file ="./data/BEF_scale_simulation_date.RData")
load("./data/BEF_scale_simulation_date.RData")

#Figure 3####
diversity.df <- hold.raw.data %>% 
  filter(sp.pool == 100) %>% 
  group_by(scenario, gamma, scale) %>% 
  summarise(lower = quantile(SR, probs = 0.25), upper = quantile(SR, probs = 0.75), richness = mean(SR))


Fig.3a<- diversity.df %>% 
  filter(scenario == "space") %>% 
  ggplot(aes(x=scale,y=richness, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")))+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("spatial scale (# of local patches)")+
  ylab("species richness")+
  theme(legend.justification=c(1,0), legend.position=c(1,0.001))+
  coord_cartesian(ylim = c(min(diversity.df$lower), max(diversity.df$upper)))

Fig.3b<- diversity.df %>% 
  filter(scenario == "time") %>% 
  ggplot(aes(x=scale,y= richness, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")), guide = F)+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab("species richness")+
  coord_cartesian(ylim = c(min(diversity.df$lower), max(diversity.df$upper)))

biomass.df <- hold.raw.data %>% 
  filter(sp.pool == 100) %>% 
  group_by(scenario, gamma, scale) %>% 
  summarise(lower = quantile(bmass, probs = 0.25), upper = quantile(bmass, probs = 0.75), bmass = mean(bmass))

Fig.3c<- biomass.df %>% 
  filter(scenario == "space") %>% 
  ggplot(aes(x=scale,y=bmass, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")), guide = F)+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("spatial scale (# of local patches)")+
  ylab("average biomass per patch")+
  theme(legend.justification=c(1,0), legend.position=c(1,0.001))+
  coord_cartesian(ylim = c(min(biomass.df$lower), max(biomass.df$upper)))

Fig.3d<- biomass.df %>% 
  filter(scenario == "time") %>% 
  ggplot(aes(x=scale,y= bmass, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")), guide = F)+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab("average biomass per time step")+
  theme(legend.justification=c(1,0), legend.position=c(1,0.001))+
  coord_cartesian(ylim = c(min(biomass.df$lower), max(biomass.df$upper)))

plot_grid(Fig.3a, Fig.3b, Fig.3c, Fig.3d,labels = "AUTO")
ggsave("./figures/Fig.3.png", height = 8, width = 9) 

#Figure 4####
sum.coef <- coef.df %>% 
  group_by(scenario, scale, gamma) %>% 
  summarise(lower = quantile(b, probs = 0.25), upper = quantile(b, probs = 0.75), b = mean(b)) %>% 
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80))

Fig.4a <- hold.raw.data %>%
  filter(gamma == 0, scenario == "space") %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  group_by(scale, sp.pool) %>% 
  summarise(SR = mean(SR), bmass = mean(bmass)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "scale", trans = "log10", breaks = c(1,3,10,30,80))+
  ylab("mean biomass per patch")+
  xlab("species richness")+
  geom_vline(data = filter(sum.coef, gamma == 0, scenario=="space"), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 0")+
  theme_classic()+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))+
  theme(legend.justification=c(1,0), legend.position=c(1,0.001))

Fig.4b <- hold.raw.data %>%
  filter(gamma == 1, scenario == "space") %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  group_by(scale, sp.pool) %>% 
  summarise(SR = mean(SR), bmass = mean(bmass)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("mean biomass per patch")+
  xlab("species richness")+
  geom_vline(data = filter(sum.coef, gamma == 1, scenario=="space"), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 1")+
  theme_classic()+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

Fig.4c <- hold.raw.data %>%
  filter(gamma == 2, scenario == "space") %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  group_by(scale, sp.pool) %>% 
  summarise(SR = mean(SR), bmass = mean(bmass)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("mean biomass per patch")+
  xlab("species richness")+
  geom_vline(data = filter(sum.coef, gamma == 2, scenario=="space"), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 2")+
  theme_classic()+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

Fig.4d <- hold.raw.data %>%
  filter(gamma == 0, scenario == "time") %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  group_by(scale, sp.pool) %>% 
  summarise(SR = mean(SR), bmass = mean(bmass)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("mean biomass per time step")+
  xlab("species richness")+
  geom_vline(data = filter(sum.coef, gamma == 0, scenario=="time"), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  theme_classic()+
  ggtitle("gamma = 0")+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

Fig.4e <- hold.raw.data %>%
  filter(gamma == 1, scenario == "time") %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  group_by(scale, sp.pool) %>% 
  summarise(SR = mean(SR), bmass = mean(bmass)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("mean biomass per time step")+
  xlab("species richness")+
  geom_vline(data = filter(sum.coef, gamma == 0, scenario=="time"), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  theme_classic()+
  ggtitle("gamma = 1")+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

Fig.4f <- hold.raw.data %>%
  filter(gamma == 2, scenario == "time") %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  group_by(scale, sp.pool) %>% 
  summarise(SR = mean(SR), bmass = mean(bmass)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("mean biomass per time step")+
  xlab("species richness")+
  geom_vline(data = filter(sum.coef, gamma == 2, scenario=="time"), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  theme_classic()+
  ggtitle("gamma = 2")+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

plot_grid(Fig.4a, Fig.4b, Fig.4c, Fig.4d, Fig.4e, Fig.4f, labels = "AUTO",nrow = 2)
ggsave("./figures/Fig.4.png", height = 10*0.7, width = 16*0.7)

#Figure 5####
sum.coef <- coef.df %>% 
  group_by(scenario, scale, gamma) %>% 
  summarise(lower = quantile(b, probs = 0.25), upper = quantile(b, probs = 0.75), b = mean(b)) 

Fig.5a <- sum.coef %>%
  filter(scenario == "space") %>% 
  ggplot(aes(x = scale,y = b, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  geom_hline(yintercept = filter(sum.coef, scale == 1, scenario == "space")$b, color = colV, lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")))+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("spatial scale (# of patches)")+
  ylab(expression(paste("BEF half saturation richness (", italic(b[i]),")"), sep = ""))+
  coord_cartesian(ylim = c(min(sum.coef$lower), max(sum.coef$upper)))+
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.background = element_rect(fill = NA, color = NA))

Fig.5b <- sum.coef %>%
  filter(scenario == "time") %>% 
  ggplot(aes(x = scale,y = b, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  geom_hline(yintercept = filter(sum.coef, scale == 1, scenario == "time")$b, color = colV, lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = "gamma", guide = F)+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("temporal scale (# of time steps)")+
  ylab(expression(paste("BEF half saturation richness (", italic(b[i]),")"), sep = ""))+
  coord_cartesian(ylim = c(min(sum.coef$lower), max(sum.coef$upper)))

plot_grid(Fig.5a, Fig.5b, labels = "AUTO")
ggsave("./figures/Fig.5.png", height = 4, width = 9) 
