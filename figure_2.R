library(tidyverse)
library(facetscales)

vnames <-list(
  
  'Number of useful variables selected' = 'Number of useful \n variables selected',
  'Number of noise variables selected' = 'Number of noise \n variables selected'
)
vlabeller <- function(variable,value){
  return(vnames[value])
}


scales_y <- list(
  "Number of noise variables selected" = scale_y_continuous(breaks = seq(0,7,1), limits = c(0,7)),
  "Number of useful variables selected" = scale_y_continuous(breaks = seq(0,8,2),limits = c(0,8))
)

read_csv("result/figure_3.csv") %>% 
  mutate(category = factor(category, levels = c("Number of useful variables selected", "Number of noise variables selected"))) %>% 
  mutate(category = factor(category, levels = c("Number of useful variables selected", "Number of noise variables selected"))) %>% 
  mutate(mod = factor(mod, levels = c( "riAFT-BART","PEAMM", "FrailtyPenal", "FrailtyHL" ,"riAH","riCox") )) %>%
  mutate(ph = factor(ph, levels = c("PH", "nPH"))) %>% 
  # mutate(selected = ifelse(mod %in% c("FrailtyHL") &ph == "PH" & category == "Number of noise variables selected"& selected>2,5, selected)) %>% 
  ggplot(aes(x = mod, y = selected, fill = ph))+
  # facet_grid(category~., labeller = vlabeller)+
  scale_fill_manual(values = c("#C6CACE", "#6D767D"))+
  facet_grid_sc(rows = vars(category), scales = list(y = scales_y), labeller = vlabeller)+
  geom_boxplot()+
  labs(x = "", y ="", fill = "")+
  scale_y_continuous(breaks = c(0,2,4,6,8,10), limits = c(0,10))+
  # scale_fill_manual(
  #   values = c("#619CFF", "#00BA38", "#F8766D"))+
  theme_bw()+
  theme(legend.position = "top", legend.text = element_text(size = 14), strip.text = element_text(size = 14),axis.text  =  element_text(size = 14))+guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave("plot/figure_2.jpg", width = 7, height = 4.5)