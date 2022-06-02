# FTICRRR: fticr results in R
# Kaizad F. Patel
# October 2020

################################################## #

## `functions_vankrevelens.R`
## this script will load functions for plotting Van Krevelen diagrams
## source this file in the `fticr_drake_plan.R` file, do not run the script here.

################################################## #
################################################## #


theme_kp <- function() {  # this for all the elements common across plots
  theme_bw() %+replace%
    theme(legend.position = "top",
          legend.key=element_blank(),
          #legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.key.size = unit(1.5, 'lines'),
          legend.background = element_rect(colour = NA),
          panel.border = element_rect(color="black",size=1.5, fill = NA),
          
          plot.title = element_text(hjust = 0.00, size = 14),
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, face = "bold", color = "black"),
          
          # formatting for facets
          panel.background = element_blank(),
          strip.background = element_rect(colour="white", fill="white"), #facet formatting
          panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
          panel.spacing.y = unit(1.5, "lines"), #facet spacing for x axis
          strip.text.x = element_text(size=12, face="bold"), #facet labels
          strip.text.y = element_text(size=12, face="bold", angle = 270) #facet labels
    )
}
gg_vankrev <- function(data,mapping){
  ggplot(data,mapping) +
    # plot points
    geom_point(size=0.5, alpha = 0.5) + # set size and transparency
    # axis labels
    ylab("H/C") +
    xlab("O/C") +
    # axis limits
    xlim(0,1.25) +
    ylim(0,2.5) +
    # add boundary lines for Van Krevelen regions
    geom_segment(x = 0.0, y = 1.5, xend = 1.2, yend = 1.5,color="black",linetype="longdash") +
    geom_segment(x = 0.0, y = 0.7, xend = 1.2, yend = 0.4,color="black",linetype="longdash") +
    geom_segment(x = 0.0, y = 1.06, xend = 1.2, yend = 0.51,color="black",linetype="longdash") +
    guides(colour = guide_legend(override.aes = list(alpha=1, size = 1)))
}


# van krevelen plots ------------------------------------------------------
plot_vankrevelen_domains = function(fticr_meta){
  
  gg_vk_domains = 
    fticr_meta %>%     
    mutate(Class = factor(Class, levels = c("aliphatic", "unsaturated/lignin", "aromatic", "condensed aromatic"))) %>% 
    filter(!is.na(Class)) %>% 
    gg_vankrev(aes(x = OC, y = HC, color = Class))+
    scale_color_manual(values = PNWColors::pnw_palette("Sunset2", 4))+
    theme_kp()+
    guides(color=guide_legend(nrow=2, override.aes = list(size = 1, alpha = 1)))+
    NULL

  gg_vk_domains_nosc = 
    gg_vankrev(fticr_meta, aes(x = OC, y = HC, color = as.numeric(NOSC)))+
    scale_color_gradientn(colors = PNWColors::pnw_palette("Bay"))+
    labs(color = "NOSC")+
    theme_kp()
  
  list(gg_vk_domains = gg_vk_domains,
       gg_vk_domains_nosc = gg_vk_domains_nosc)
}

plot_vankrevelens = function(fticr_data_trt, fticr_meta){
  
  fticr_hcoc = 
    fticr_data_trt %>% 
    left_join(dplyr::select(fticr_meta, formula, HC, OC), by = "formula")
  
  gg_vankrev_all = 
    fticr_hcoc %>% 
    gg_vankrev(aes(x = OC, y = HC, color = Incubation_Temp))+
    stat_ellipse(level = 0.90, show.legend = F)+
    scale_color_manual(values = pnw_palette("Bay", 3))+
    facet_wrap(~A_priori_CUE)+
    theme_kp()+
    NULL
  
  list(gg_vankrev_all = gg_vankrev_all)
}

plot_vk_unique = function(fticr_data_trt, fticr_meta){
  fticr_hcoc = 
    fticr_data_trt %>% 
    left_join(dplyr::select(fticr_meta, formula, HC, OC), by = "formula")
  
  peaks_unique = 
    fticr_hcoc %>% 
    group_by(A_priori_CUE, formula) %>% 
    dplyr::mutate(n = n()) %>% 
    filter(n == 1) %>% 
    ungroup()
  
  peaks_unique %>% 
    gg_vankrev(aes(x = OC, y = HC, color = Incubation_Temp))+
    stat_ellipse(level = 0.90, show.legend = F)+
    facet_wrap(~A_priori_CUE)+
    scale_color_manual(values = pnw_palette("Bay", 3))+
    labs(title = "unique peaks")+
    theme_kp()+
    NULL
  
    
}



################################################## #
################################################## #


# NOSC figures ------------------------------------------------------------


make_nosc_figures = function(fticr_data_trt, fticr_meta){
  fticr_data_nosc = 
    fticr_data_trt %>% 
    left_join(dplyr::select(fticr_meta, formula, NOSC)) %>% 
    mutate(length = factor(length, levels = c("timezero", "30d", "90d", "150d")))
  
  nosc_by_drying = 
    fticr_data_nosc %>% 
    ggplot(aes(x = NOSC, fill = drying, color = drying))+
    geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.5)+
    facet_grid(Site + depth ~ length + saturation)+
    theme_kp()+
    NULL
  
  nosc_by_saturation = 
    fticr_data_nosc %>% 
    ggplot(aes(x = NOSC, fill = saturation, color = saturation))+
    geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.5)+
    facet_grid(Site + depth ~ length + drying)+
    theme_kp()+
    NULL
  
  list(nosc_by_drying = nosc_by_drying,
       nosc_by_saturation = nosc_by_saturation)
  
  
## fticr_data_nosc %>% 
##   distinct(formula, NOSC, saturation, depth) %>% 
##   ggplot(aes(x = NOSC, fill = saturation, color = saturation))+
##   geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.5)+
##   facet_grid(depth ~ .)+
##   theme_kp()+
##   NULL
## 
## fticr_data_nosc %>% 
##   distinct(formula, NOSC, saturation, drying, depth) %>% 
##   ggplot(aes(x = NOSC, fill = drying, color = drying))+
##   geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.5)+
##   facet_grid(depth ~ saturation)+
##   theme_kp()+
##   NULL
  
}
