library(tidyverse)
reps <- 20
hold.raw.data <-data.frame()
coef.df <- data.frame()
sensitivity_trials <- c("r_max" , "sigma")#c("r_max", "alpha_max", "sigma", "max_scale")
for(r in 1:reps){
  for(gamma in c(0,1,2)){
    for(sens_type in sensitivity_trials){
      r_max <- 5
      alpha_max <- 0.25
      sigma <- 0.25
      max_scale <- 80
      for(sens_2 in c("low", "high")){
        if(sens_type == "r_max") {
          if(sens_2 == "low"){
            r_max <- 3
          } else {
            r_max <- 10
          }
        } 
        if(sens_type == "alpha_max") {
          if(sens_2 == "low"){
            alpha_max <- 0.1
          } else {
            alpha_max <- 0.5
          }
        } 
        if(sens_type == "sigma") {
          if(sens_2 == "low"){
            sigma <- 0.2
          } else {
            sigma <- 0.3
          }
        } 
        if(sens_type == "max_scale") {
          if(sens_2 == "low"){
            max_scale <- 40
          } else {
            max_scale <- 160
          }
        } 
        print(paste("rep = ",r,"; gamma = ", gamma, ", ", sens_type, ", ", sens_2, sep = ""))
        spatial_run <- BEF_simulation(type = "spatial", env_gamma = gamma, rep = r, 
                                      r_max = r_max, 
                                      alpha_max = alpha_max, 
                                      sigma = sigma, 
                                      max_scale = max_scale)
        spatial_com <- spatial_run[[1]]
        temporal_run <- BEF_simulation(type = "temporal", env_gamma = gamma, rep = r,
                                       r_max = r_max, 
                                       alpha_max = alpha_max, 
                                       sigma = sigma, 
                                       max_scale = max_scale)
        temporal_com <- temporal_run[[1]][(max_scale+1):(max_scale*2),,]
        
        hold.data.run<-data.frame()
        start_id <- 1
        for(k in c(1:30,32,34,36,38,seq(40,max_scale, by = 5))){
          time_seq <- max_scale:(max_scale-k+1)
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
            hold.data$sens_type <- sens_type
            hold.data$sens_2 <- sens_2
          }
          coef.hold_space <- coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = filter(hold.data, scenario == "space"),start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))
          coef.hold_time <- coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), data = filter(hold.data, scenario == "time"), start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))
          coef.hold.df <- bind_rows(coef.hold_space,coef.hold_time)
          coef.hold.df$scenario = c("space", "time")
          coef.hold.df$rep <- r
          coef.hold.df$scale <- k
          coef.hold.df$gamma <- gamma
          coef.hold.df$sens_type <- sens_type
          coef.hold.df$sens_2 <- sens_2
          coef.df <-bind_rows(coef.df, coef.hold.df)
          
          hold.raw.data <- rbind(hold.raw.data, hold.data)
        }
      }
    }
  }
}


save(coef.df, file = "./data/sensitivity.RData")

sum.coef <- coef.df %>% 
  group_by(sens_type, sens_2, scenario, scale, gamma) %>% 
  summarise(lower = quantile(b, probs = 0.25), upper = quantile(b, probs = 0.75), b = mean(b))

colV <- c("grey40", "deeppink1", "red3")

sum.coef %>%
  filter(sens_2 == "low") %>% 
  ggplot(aes(x = scale,y = b, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  facet_grid(sens_type ~ scenario, scales = "free")+
  #geom_hline(yintercept = filter(sum.coef, scale == 1, scenario == "space")$b, color = colV, lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")))+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("spatial scale (# of patches)")+
  ylab(expression(paste("BEF half saturation richness (", italic(b[i]),")"), sep = ""))+
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.background = element_rect(fill = NA, color = NA))
ggsave("./figures/Fig_S1.pdf", height = 8, width = 8)


sum.coef %>%
  filter(sens_2 == "high") %>% 
  ggplot(aes(x = scale,y = b, color = factor(gamma), group = gamma, fill = factor(gamma)))+
  facet_grid(sens_type ~ scenario, scales = "free")+
  #geom_hline(yintercept = filter(sum.coef, scale == 1, scenario == "space")$b, color = colV, lty = 2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, col = NA)+
  geom_line(size = 1)+
  scale_color_manual(values = colV, name = expression(paste("env. ", gamma, sep = "")))+
  scale_fill_manual(values = colV, guide = F)+
  theme_classic()+
  xlab("spatial scale (# of patches)")+
  ylab(expression(paste("BEF half saturation richness (", italic(b[i]),")"), sep = ""))+
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.background = element_rect(fill = NA, color = NA))
ggsave("./figures/Fig_S2.pdf", height = 8, width = 8)

#intraspecific variation in alpha_ii sensitivity analysis####
reps <- 20
hold.raw.data <-data.frame()
coef.df <- data.frame()
for(r in 1:reps){
  for(gamma in c(0,1,2)){
    print(paste("rep = ",r,"; gamma = ", gamma, sep = ""))
    spatial_run <- BEF_simulation(type = "spatial", env_gamma = gamma, rep = r, alpha_ii_sd = 0.25)
    spatial_com <- spatial_run[[1]]
    temporal_run <- BEF_simulation(type = "temporal", env_gamma = gamma, rep = r, alpha_ii_sd = 0.25)
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

sum.coef <- coef.df %>% 
  group_by(scenario, scale, gamma) %>% 
  summarise(lower = quantile(b, probs = 0.25), upper = quantile(b, probs = 0.75), b = mean(b)) 

Fig.S3a <- sum.coef %>%
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

Fig.S3b <- sum.coef %>%
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

plot_grid(Fig.S3a, Fig.S3b, labels = "AUTO")
ggsave("./figures/Fig.S3.png", height = 4, width = 9) 
ggsave("./figures/Fig.S3.pdf", height = 4, width = 9) 





