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
b.df <- data.frame()
for(i in c(1,2,3,5,10,20,30,40,60,80)){
  for(g_sel in c(0,1,2)){
    hold.data = filter(hold.raw.data, rep == 2, gamma == g_sel, scale == i)
    b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
                                                  data = hold.data, 
                                                  start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2], 
                                   scale = i,
                                   gamma = g_sel))
  }}

A <- hold.raw.data %>% 
  filter(rep == 2, gamma == 0) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  geom_vline(data = filter(b.df, gamma == 0), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 0")+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

B <- hold.raw.data %>% 
  filter(rep == 2, gamma == 1) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  geom_vline(data = filter(b.df, gamma == 1), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 1")+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))

C <- hold.raw.data %>% 
  filter(rep == 2, gamma == 2) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "spatial\nscale", trans = "log10")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  geom_vline(data = filter(b.df, gamma == 2), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 2")+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylim(c(0,max(hold.raw.data$bmass)))+
  xlim(c(0,100))


load("./data/LV_time.RData")

b.df <- data.frame()
for(i in c(1,2,3,5,10,20,40,100,200,500)){
  for(g_sel in c(0,1,2)){
    hold.data = filter(hold.raw.data, rep == 2, gamma == g_sel, t_scale == i)
    b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
                                                  data = hold.data, 
                                                  start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2], 
                                   t_scale = i,
                                   gamma = g_sel))
  }}

D <- hold.raw.data %>% 
  filter(rep == 2, gamma == 0) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = t_scale, group = t_scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  #geom_smooth(method = "nls", 
  #            formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  geom_vline(data = filter(b.df, gamma == 0), aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 0")+
  ylim(c(0,2700))+
  xlim(c(0,100))

E <- hold.raw.data %>% 
  filter(rep == 2, gamma == 1) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = t_scale, group = t_scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  #geom_smooth(method = "nls", 
  #            formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  geom_vline(data = filter(b.df, gamma == 1), aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 1")+
  ylim(c(0,2700))+
  xlim(c(0,100))

FigF <- hold.raw.data %>% 
  filter(rep == 2, gamma == 2) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass, fill = t_scale, group = t_scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", trans = "log10")+
  ylab("cumulative biomass")+
  xlab("species richness")+
  #geom_smooth(method = "nls", 
  #            formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  geom_vline(data = filter(b.df, gamma == 2), aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 2")+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylim(c(0,2700))+
  xlim(c(0,100))

plot_grid(A,B,C,D,E,FigF, labels = c("a)", "b)", "c)", "d)", "e)", "f)"),nrow = 2)

ggsave("./figures/raw_BEF.png", height = 10, width = 16)


#averaged version####
load("./data/LV_data.RData")
b.df <- data.frame()
for(i in c(1,2,3,5,10,20,30,40,60,80)){
  for(g_sel in c(0,1,2)){
    hold.data = filter(hold.raw.data, rep == 2, gamma == g_sel, scale == i)
    b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
                                                  data = hold.data, 
                                                  start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2], 
                                   scale = i,
                                   gamma = g_sel))
  }}

A <- hold.raw.data %>% 
  filter(rep == 2, gamma == 0) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass/scale, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass/scale")+
  xlab("species richness")+
  geom_vline(data = filter(b.df, gamma == 0), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 0")+
  ylim(c(0,max(hold.raw.data$bmass)/max(hold.raw.data$scale)))+
  xlim(c(0,100))

B <- hold.raw.data %>% 
  filter(rep == 2, gamma == 1) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass/scale, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass/scale")+
  xlab("species richness")+
  geom_vline(data = filter(b.df, gamma == 1), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 1")+
  ylim(c(0,max(hold.raw.data$bmass)/max(hold.raw.data$scale)))+
  xlim(c(0,100))

C <- hold.raw.data %>% 
  filter(rep == 2, gamma == 2) %>%
  filter(scale %in% c(1,2,3,5,10,20,30,40,60,80)) %>% 
  ggplot(aes(x = SR, y = bmass/scale, fill = scale, group = scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "spatial\nscale", trans = "log10")+
  ylab("cumulative biomass/scale")+
  xlab("species richness")+
  geom_vline(data = filter(b.df, gamma == 2), aes(xintercept = b, color = factor(scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 2")+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylim(c(0,max(hold.raw.data$bmass)/max(hold.raw.data$scale)))+
  xlim(c(0,100))


load("./data/LV_time.RData")

b.df <- data.frame()
for(i in c(1,2,3,5,10,20,40,100,200,500)){
  for(g_sel in c(0,1,2)){
    hold.data = filter(hold.raw.data, rep == 2, gamma == g_sel, t_scale == i)
    b.df <- rbind(b.df, data.frame(b = coef(nlsLM(formula = bmass ~ (a * SR)/ (SR + b), 
                                                  data = hold.data, 
                                                  start = c(a = max(hold.data$bmass), b = max(hold.data$bmass)/2)))[2], 
                                   t_scale = i,
                                   gamma = g_sel))
  }}

D <- hold.raw.data %>% 
  filter(rep == 2, gamma == 0) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass/t_scale, fill = t_scale, group = t_scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass/scale")+
  xlab("species richness")+
  geom_vline(data = filter(b.df, gamma == 0), aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 0")+
  ylim(c(0,5.4))+
  xlim(c(0,100))

E <- hold.raw.data %>% 
  filter(rep == 2, gamma == 1) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass/t_scale, fill = t_scale, group = t_scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", guide = FALSE, trans = "log10")+
  ylab("cumulative biomass/scale")+
  xlab("species richness")+
  geom_vline(data = filter(b.df, gamma == 1), aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 1")+
  ylim(c(0,5.4))+
  xlim(c(0,100))

FigF <- hold.raw.data %>% 
  filter(rep == 2, gamma == 2) %>%
  filter(t_scale %in% c(1,2,3,5,10,20,40,100,200,500)) %>% 
  ggplot(aes(x = SR, y = bmass/t_scale, fill = t_scale, group = t_scale))+
  scale_fill_viridis_c(end = 0.8, option = "B", name = "temporal\nscale", trans = "log10")+
  ylab("cumulative biomass/scale")+
  xlab("species richness")+
  geom_vline(data = filter(b.df, gamma == 2), aes(xintercept = b, color = factor(t_scale), fill = NULL), linetype = 2)+
  geom_point(pch = 21, size = 2.5)+
  geom_smooth(method = "nls", 
              formula = y ~ a * x / (b + x), se = FALSE, aes(color = factor(t_scale)))+
  scale_color_viridis_d(end = 0.8, option = "B", guide = FALSE)+
  scale_linetype(guide = FALSE)+
  ggtitle("gamma = 2")+
  theme(legend.justification=c(1,0), legend.position=c(1,0))+
  ylim(c(0,5.4))+
  xlim(c(0,100))

plot_grid(A,B,C,D,E,FigF, labels = c("a)", "b)", "c)", "d)", "e)", "f)"),nrow = 2)

ggsave("./figures/raw_BEF_scaled.png", height = 10, width = 16)

