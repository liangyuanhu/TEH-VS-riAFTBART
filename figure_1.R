library(tidyverse)
library(facetscales)

vnames_2 <-list(
  # 'ATE12' = bquote(ATE["1,2"]~"(Truth = 2.11 months)"),
  # 'ATE13' = bquote(ATE["1,3"]~"(Truth = 0.89 months)"),
  # 'ATE23' = bquote(ATE["2,3"]~"(Truth = -1.21 months)"),
  # 'ATE12PH' = bquote(ATE["1,2"]~"& PH"),
  # 'ATE13PH' = bquote(ATE["1,3"]~"& PH"),
  # 'ATE23PH' = bquote(ATE["2,3"]~"& PH"),
  # 'ATE12nPH' = bquote(ATE["1,2"]~"& nPH"),
  # 'ATE13nPH' = bquote(ATE["1,3"]~"& nPH"),
  # 'ATE23nPH' = bquote(ATE["2,3"]~"& nPH"),
  '40%' = "Censor rate = 40%",
  '10%' = "Censor rate = 10%",
  'PH'= "PH",
  'nPH' = "nPH",
  # 'ATE12' = bquote(CATE["1,2"]),
  # 'ATE13' = bquote(CATE["1,3"]),
  # 'ATE23' = bquote(CATE["2,3"])
  'ATE12' = "Trt 1 vs. 2",
  'ATE13' = "Trt 1 vs. 3",
  'ATE23' = "Trt 2 vs. 3"
)
vlabeller_2 <- function(variable,value){
  return(vnames_2[value])
}

(rmst_relative_bias_moderate <- read_csv("result/figure_1_relative_bias.csv") %>% 
    ggplot(aes(x = method, y = value, fill = hs))+
    geom_boxplot()+
    # geom_violin()+
    facet_grid(ph~group,labeller = vlabeller_2)+
    labs(x = "", y = "Relative bias", fill = "Heterogeneity setting")+
    # scale_fill_manual(values = c("#C6CACE", "#6D767D"))+
    scale_y_continuous(labels = scales::percent, limits = c(-0.2,0.5), breaks = c(-0.2,0,0.2,0.4))+
    geom_hline(yintercept = 0, linetype = 2)+
    theme_bw()+
    theme(legend.position = "top",
          axis.title = element_text(size = 12),
          panel.margin = unit(1, "lines"),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text = element_text(size = 12),
          legend.text =  element_text(size = 12),
          legend.title  =  element_text(size = 12)))
# ggsave("plot/rmst_relative_bias_moderate.jpg", width = 6, height = 5)


(rmst_rmse_moderate <- read_csv("result/figure_1_result_rmst.csv") %>% 
    ggplot(aes(x = method, y = value, fill = hs))+
    geom_boxplot()+
    # geom_violin()+
    facet_grid(ph~group,labeller = vlabeller_2)+
    labs(x = "", y = "RMSE", fill = "Heterogeneity setting")+
    scale_y_continuous( limits = c(0,1.3), breaks = c(0,0.4,0.8,1.2))+
    geom_hline(yintercept = 0, linetype = 2)+
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_text(size = 12),
          panel.margin = unit(1, "lines"),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12,        
                                     angle = 90, vjust = 0.5, hjust=1),
          strip.text = element_text(size = 12),
          legend.text =  element_text(size = 12),
          legend.title  =  element_text(size = 12)))
# ggsave("plot/rmst_rmse_moderate.jpg", width = 6, height = 5)

p_rmst <- cowplot::plot_grid(
  rmst_relative_bias_moderate,
  rmst_rmse_moderate,
  nrow = 2,
  rel_heights = c(0.45, 0.55),
  align = "v",
  labels = "AUTO"
)

ggsave(p_rmst, file = "plot/figure_1.jpg", width = 6, height = 7)