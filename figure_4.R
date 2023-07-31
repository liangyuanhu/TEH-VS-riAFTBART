library(tidyverse)
read_csv("result/figure_4.csv") %>% 
  mutate(method = factor(method, levels = c( "riAFT-BART","PEAMM", "FrailtyPenal", "FrailtyHL" ,"riCox") )) %>%  
  ggplot(aes(x = method, y = value))+
  geom_boxplot()+
  # scale_y_continuous(limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3))+
  labs(x = "", y = "Cross-validated concordance statistic")+
  theme_bw()

ggsave("plot/figure_4.jpg", width = 4, height = 3)