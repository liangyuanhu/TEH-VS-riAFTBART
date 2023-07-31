library(tidyverse)

read_csv("result/result_power_variable_nPH.csv") %>% 
  mutate(mod = ifelse(mod == "riGAPH", "PEAMM", mod)) %>% 
  mutate(mod = factor(mod, levels = c( "riAFT-BART","PEAMM", "FrailtyPenal", "FrailtyHL" ,"riAH","riCox") )) %>% 
  filter(mod != "riAH") %>% 
  mutate(ph = "nPH") %>% 
  mutate(mod = paste0(mod, "-nPH")) %>% 
  bind_rows(read_csv("result/result_power_variable_PH.csv") %>% 
              mutate(mod = ifelse(mod == "riGAPH", "PEAMM", mod)) %>% 
              mutate(mod = factor(mod, levels = c( "riAFT-BART","PEAMM", "FrailtyPenal", "FrailtyHL" ,"riAH","riCox") )) %>% 
              filter(mod != "riAH") %>% 
              mutate(ph = "PH") %>% 
              mutate(mod = paste0(mod, "-PH"))) %>% 
  mutate(mod = factor(mod, levels = c("riAFT-BART-PH","PEAMM-PH", "FrailtyPenal-PH", "FrailtyHL-PH" ,"riAH-PH","riCox-PH","riAFT-BART-nPH","PEAMM-nPH", "FrailtyPenal-nPH", "FrailtyHL-nPH" ,"riAH-nPH","riCox-nPH"))) %>% 
  # count(ph)
  ggplot(aes(x = time, y = power, shape = mod, fill = mod))+
  geom_point()+
  scale_shape_manual(name = "", values = c( 21:25, 21:25))+
  scale_fill_manual(name = "", values = c(rep("black",5), rep("white",5)))+
  theme_bw()+
  theme(legend.position = "top", panel.grid.minor = element_blank())+
  scale_x_continuous(breaks = 1:10, labels = c(expression(italic(x)[1]),
                                               expression(italic(x)[2]),
                                               expression(italic(x)[3]),
                                               expression(italic(x)[4]),
                                               expression(italic(x)[5]),
                                               expression(italic(x)[6]),
                                               expression(italic(x)[7]),
                                               expression(italic(x)[8]),
                                               expression(italic(x)[9]),
                                               expression(italic(x)[10])))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1))+
  guides(shape=guide_legend(nrow=2,byrow=T),
         fill = guide_legend(nrow=2,byrow=T))+
  labs(y = "Power",x= "", shape = "", color = "")

ggsave("plot/variable_power.jpg", width = 7, height = 4)