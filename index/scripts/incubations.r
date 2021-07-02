library(readxl)
library(here)
library(tidyverse)
library(ggplot2)
library(reshape2)

# import raw data


IC_raw_data <- read_excel(path.expand(here("index","data", "IC_raw_data.xlsx")))
incubationgraph <- read_excel(here("index","data","incubationgraph.xlsx"))

####== Transform the data ==####

# middle duplicates for IC data, and calculate umol/L
IC_inc_data <- 
  IC_raw_data %>% 
  filter(str_detect(sample,'inc')) %>% 
  group_by(sample) %>% 
  mutate(P = (mean(Phosphate)*1000)/94.9714, 
         "SO" = mean(Sulphate)*1000/96.06, 
         NO = mean(Nitrate)*1000/62.004, 
         "Br" = mean(Bromide)*1000/79.9, 
         "Fl" = mean(Fluoride)*1000/19,
         "Cl" = mean(Chloride)*1000/35.45) %>% 
  distinct(sample, .keep_all = TRUE) 


# Join IC and photometric data, and add treated and oxic status for all cores.
inc_data <- 
  full_join(
    incubationgraph,
    IC_inc_data,
    by = c("samplenr" = "sample"),
    copy = FALSE,
    suffix = c(".photometric", ".IC"),
    keep = FALSE
  ) %>% 
  mutate(
    HS         = 1.1*(((H2SuM*1.13759603802245)/0.026915542)/1000),
    PFe        = Preal/FeTot,
    FeIIFeIII  = FeII/FeIII, 
    FeIIFetot  = FeII/FeTot,
    oxic_state = case_when( str_detect(samplenr,"2\\.")  ~ "oxic",
                            str_detect(samplenr,"3\\.")  ~ "anoxic"),
    treated    = case_when( str_detect(samplenr,"A|C")  ~ "treated",
                            str_detect(samplenr,"B|D")  ~ "nontreated"),
    Location   = case_when( str_detect(samplenr,"A")  ~ "A",
                            str_detect(samplenr,"-B")  ~ "B",
                            str_detect(samplenr,"-\\C")  ~ "C",
                            str_detect(samplenr,"-\\D")  ~ "D")) 

# calculate benthic fluxes
flux_data <- inc_data %>% 
  group_by(Core.photometric) %>%  
  summarise(
    treated = treated,
    Location = Location,
    oxic_state = oxic_state,
    time = time, 
    dP = c(0, diff((0.001*Preal*volume_ml)/(pi*9))/diff(time)),
    dFe = c(0, diff((0.001*FeTot*volume_ml)/(pi*9))/diff(time)),
    dNH = c(0, diff((0.001*NH4*volume_ml)/(pi*9))/diff(time)),
    dSH = c(0, diff((0.001*H2SuM*volume_ml)/(pi*9))/diff(time)),
    dNO = c(0, diff((0.001*NO*volume_ml)/(pi*9))/diff(time)),
    dSO = c(0, diff((0.001*SO*volume_ml)/(pi*9))/diff(time))
  ) %>% ungroup()

# reshape the data by melting all variables, and grouping them by parameter.

molten_inc <- inc_data %>% select("Core.photometric", 
                                  "samplenr",
                                  "time", 
                                  "HS", 
                                  "FeTot", 
                                  "Preal", 
                                  "NH4", 
                                  "SO",
                                  "NO",
                                  "oxic_state", 
                                  "treated", 
                                  "Location") %>% 
  melt(id.vars = c("Core.photometric","samplenr","time", "oxic_state", "treated","Location"), 
       measured_.vars = c( "HS", 
                           "FeTot", 
                           "Preal", 
                           "NH4", 
                           "SO",
                           "NO"), na.rm = TRUE)

molten_flux <-  melt(flux_data, id.vars = c("Core.photometric","time","oxic_state","treated","Location"))

####== create plots for all parameters ==####

# Lable function
Lables = as_labeller( c(  "HS" = "HS^'-'", 
                          "FeII" = "Fe^'2+'",
                          "FeTot" = "Fe[tot]", 
                          "FeIII" = "Fe^'3+'",  
                          "Preal" = "PO[4]^'3-'",
                          "NH4" = "NH[4]^'+'",
                          "P"  = "PO[4]^'3-'",                      
                          "SO" = "SO[4]^'2-'", 
                          "NO" = "NO[3]^'-'",
                          "Br" = "Br^'-'",
                          "Fl" = "F^'-'", 
                          "Cl" = "Cl^'-'",
                          "oxic" = "'Oxic incubations'",
                          "anoxic" = "'Anoxic incubations'"), default =  label_parsed )

d <- list(molten_inc, molten_flux)
for (i in 1:2) {
  
  graph <- 
    ggplot(transform(d[[i]],
                     variable = factor(variable, levels = c("FeTot","Preal","HS","NH4","SO","NO"))), 
           mapping = aes(x = time , 
                         y = value, 
                         color = Location, 
                         shape = Location, 
                         by = Core.photometric)
    ) + 
    scale_color_manual( values=c("A" = "coral3", "B" = "deepskyblue4", "C" = "coral3", "D" = "deepskyblue4")) +
    scale_shape_manual( values=c("A" = 16, "B" = 1, "C" = 17,  "D" = 2)) +
    geom_point() + 
    geom_line() + 
    theme_gray() +
    theme(axis.title=element_text(size=20),
          plot.title = element_text(size=25, face="bold", hjust = 0.5), #hjust is position (left/right/ middle =0.5)
          axis.text=element_text(size=14,angle = 0, hjust = 0.5),
          strip.text.y   = element_text(size=14, angle = 0),
          strip.text.x   = element_text(size=14, angle = 0),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 21)) + 
    facet_grid(variable ~ oxic_state, scales = "free_y",switch = "y", labeller = Lables) + labs(title = "Benthic flux experiment", x= "time (days)", y = expression(paste("concentration in ",mu,"mol/L")))
  
  ggsave(paste("IC_incubations_",i,".png",sep=""), plot =graph, path = path.expand(here("index","figures")),
         
         width =30, height = 27,units = "cm",dpi = 600)
}

####== Benthic flux calculations ==####

# slope of NH4 line

initinc <- inc_data %>% filter(time > 15 ) %>% mutate(NH4mol = 0.001*NH4*`volume_ml` ) 
initplot <- ggplot( initinc,
                    mapping = aes(x = time , 
                                  y = NH4mol, 
                                  color = treated, 
                                  shape = Location, 
                                  by = Core.photometric)
) + 
  scale_color_manual( values=c("treated" = "red", "nontreated" = "blue")) +
  scale_shape_manual( values=c("A" = 16, "B" = 1, "C" = 17,  "D" = 2)) +
  geom_point() + 
  geom_line() + 
  stat_smooth(method = "lm") +
  facet_grid(~ oxic_state)


initslopes <- initinc %>%  group_by(Core.photometric) %>% 
  summarize(NH4slope = coef(lm(NH4mol ~ time))[2]
  ) %>% filter(str_detect(Core.photometric,"3\\.")) %>%  add_column(
    O2slope = c(-0.172133,
                -0.404987,
                -0.360635,
                -0.0793467,
                -0.370374,
                -0.174982,
                -0.227892,
                -0.167105
    )) %>% mutate(C_est = NH4slope*6.6, O2_umol_day = -1000*(O2slope*24)/31.9988, C_mg_day = -30*(O2slope*24)/31.9988) 


initreg <- ggplot(initslopes, 
                  mapping = aes( 
                    y=NH4slope, 
                    x=O2_umol_day),
                  scale_x_discrete(position = 'top')) +
  geom_point(
    stat = "identity", position = position_stack(reverse=TRUE),color="black") 

####== P to Fe ratio calculations ==####

inc_data %>%  
  group_by(Core.photometric) %>% 
  summarise(PtoFe = mean(PFe, na.rm = TRUE), 
            FeIItoFeIII = mean(FeIIFeIII, na.rm = TRUE), 
            FeIItoFetot = mean(FeIIFetot, na.rm = TRUE)) 

PFe_bp <- ggplot(inc_data, mapping = aes(x = Core.photometric, y = PFe)) + geom_boxplot()
FeIIFeIII_bp <- ggplot(inc_data, mapping = aes(x = Core.photometric, y = FeIIFetot)) + geom_boxplot()


PFe_inc_plot <-  ggplot(inc_data, 
                        mapping = aes(x = time , 
                                      y = PFe, 
                                      color = treated, 
                                      shape = Location, 
                                      by = Core.photometric)
) + 
  scale_color_manual( values=c("treated" = "red", "nontreated" = "blue")) +
  scale_shape_manual( values=c("A" = 16, "B" = 1, "C" = 17,  "D" = 2)) +
  geom_point() + 
  geom_line() + facet_grid(Location ~ oxic_state, scales = "free_y",switch = "y")

FeIIFetot_inc_plot <-  ggplot(inc_data, 
                              mapping = aes(x = time , 
                                            y = FeIIFetot, 
                                            color = treated, 
                                            shape = Location, 
                                            by = Core.photometric)
) + 
  scale_color_manual( values=c("treated" = "red", "nontreated" = "blue")) +
  scale_shape_manual( values=c("A" = 16, "B" = 1, "C" = 17,  "D" = 2)) +
  geom_point() + 
  geom_line() + facet_grid(Location ~ oxic_state, scales = "free_y")
show(PFe_inc_plot)
