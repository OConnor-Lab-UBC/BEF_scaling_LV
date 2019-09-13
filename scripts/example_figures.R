pdf("./figures/landscape.pdf", height = 8, width = 8)
plot.igraph(graph, layout = layout.circle(graph), vertex.color = "forestgreen", vertex.label = NA, vertex.size = 7)
dev.off()


N.df <- data.frame(N = c(hold[[1]]), species = 1:100, patch = rep(1:100, each = 100), global_S = rep(1:100, each = 80*100))
dim(N.df)

N.df_sum <- N.df %>% 
  group_by(global_S, species) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  group_by(global_S) %>% 
  summarise(S = sum(N>0), N = sum(N))

ggplot(N.df_sum, aes(x = S, y = N))+
  geom_point()

env.df<-data.frame(env = c(phase.partnered(n = 80,gamma = 1, mu = 0.5, sigma = 0.25)$timeseries[,1],
                           phase.partnered(n = 80,gamma = 0, mu = 0.5, sigma = 0.25)$timeseries[,1],
                           phase.partnered(n = 80,gamma = 2, mu = 0.5, sigma = 0.25)$timeseries[,1]),
                   gamma = rep(c(1,0,2),each = 80),
                   patch = 1:80)

acf(phase.partnered(n = 80,gamma = -2, mu = 0.5, sigma = 0.25)$timeseries[,1])
acf(phase.partnered(n = 80,gamma = 0, mu = 0.5, sigma = 0.25)$timeseries[,1])
acf(phase.partnered(n = 80,gamma = 1, mu = 0.5, sigma = 0.25)$timeseries[,1])
acf(phase.partnered(n = 80,gamma = 2, mu = 0.5, sigma = 0.25)$timeseries[,1])

ggplot(env.df, aes(x = patch, y = env, color = gamma, group = gamma))+
  geom_line(size = 1)+
  theme_classic()+
  facet_grid(gamma~.)+
  scale_color_gradient(low = "grey", high = "red1", name = "gamma")+
  ylab("environment")
ggsave("./figures/env_space.pdf", height = 4, width = 5.5)

ggplot(env.df, aes(x = patch, y = env, color = gamma, group = gamma))+
  geom_line(size = 1)+
  theme_classic()+
  scale_color_gradient(low = "grey", high = "red1", name = "\u03B3")+
  ylab("environment")
ggsave("./figures/env_space_all.pdf", height = 4, width = 5)


hold<-meta_dyn_model_sp_pool(spatial_env = TRUE, env_gamma = 2, disp = 0)

env<-(decostand(c(hold[[2]][1,],-0.2,1.2),method = "range")[1:80])*100

comp_C <- data.frame(N = c(hold[[1]][,,100]), species = rep(1:100, each = 80), patch = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = patch, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F)+
  geom_line(data = data.frame(species = env, patch = 1:80), aes(size = NULL),color = "red", size = 1, alpha = 0.4)+
  ggtitle("gamma = 2")

hold2<-meta_dyn_model_sp_pool(spatial_env = TRUE, env_gamma = 1, disp = 0)
env<-(decostand(c(hold2[[2]][1,],-0.2,1.2),method = "range")[1:80])*100

comp_B <- data.frame(N = c(hold2[[1]][,,100]), species = rep(1:100, each = 80), patch = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = patch, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1),guide = F)+
  geom_line(data = data.frame(species = env, patch = 1:80), aes(size = NULL),color = "red", size = 1, alpha = 0.4)+
  ggtitle("gamma = 1")

hold3<-meta_dyn_model_sp_pool(spatial_env = TRUE, env_gamma = 0, disp = 0)
env<-(decostand(c(hold3[[2]][1,],-0.2,1.2),method = "range")[1:80])*100

comp_A <- data.frame(N = c(hold3[[1]][,,100]), species = rep(1:100, each = 80), patch = 1:80) %>%
  filter(N >0) %>% 
  ggplot(aes(x = patch, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1), guide = F)+
  geom_line(data = data.frame(species = env, patch = 1:80), aes(size = NULL),color = "red", size = 1, alpha = 0.4)+
  ggtitle("gamma = 0")

hold4<- meta_dyn_model_sp_pool_time(time_gamma = 0)
env<-(decostand(c(hold4[[2]],-0.2,1.2),method = "range")[501:1000])*100

comp_D <- data.frame(N = c(t(hold4[[1]][,501:1000,100])), species = rep(1:100, each = 500), time = 1:500) %>%
  filter(N >0) %>% 
  ggplot(aes(x = time, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1), guide = F)+
  geom_line(data = data.frame(species = env, time = 1:500), aes(size = NULL),color = "red", alpha = 0.4)+
  ggtitle("gamma = 0")

hold5<- meta_dyn_model_sp_pool_time(time_gamma = 1)
env<-(decostand(c(hold5[[2]],-0.2,1.2),method = "range")[501:1000])*100

comp_E <- data.frame(N = c(t(hold5[[1]][,501:1000,100])), species = rep(1:100, each = 500), time = 1:500) %>%
  filter(N >0) %>%  
  ggplot(aes(x = time, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1), guide = F)+
  geom_line(data = data.frame(species = env, time = 1:500), aes(size = NULL),color = "red", alpha = 0.4)+
  ggtitle("gamma = 1")

hold6<- meta_dyn_model_sp_pool_time(time_gamma = 2)
env<-(decostand(c(hold6[[2]],-0.2,1.2),method = "range")[501:1000])*100

comp_F <- data.frame(N = c(t(hold6[[1]][,501:1000,100])), species = rep(1:100, each = 500), time = 1:500) %>%
  filter(N >0) %>% 
  ggplot(aes(x = time, y = species, size = N))+
  geom_point()+
  ylab("species ID")+
  scale_size_continuous(range = c(0.1,1), guide = F)+
  geom_line(data = data.frame(species = env, time = 1:500), aes(size = NULL),color = "red", alpha = 0.4)+
  ggtitle("gamma = 2")


plot_grid(comp_A, comp_B, comp_C, comp_D,comp_E,comp_F,  nrow = 2, labels = c("a)", "b)", "c)", "d)", "e)","f)"))
ggsave("./figures/compositional_turnover_space.png", width = 10*1.2, height = 6*1.2)  
