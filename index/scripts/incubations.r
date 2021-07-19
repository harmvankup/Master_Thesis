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
  mutate(
    dP = (0.001*Preal*volume_ml)/(pi*9),
    dFe = (0.001*FeTot*volume_ml)/(pi*9),
    dNH =  (0.001*NH4*volume_ml)/(pi*9),
    dSH = (0.001*H2SuM*volume_ml)/(pi*9),
    dNO = (0.001*NO*volume_ml)/(pi*9),
    dSO = (0.001*SO*volume_ml)/(pi*9)
  ) 

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

molten_flux <- flux_data %>% select("Core.photometric", 
                                    "time", 
                                    "dSH", 
                                    "dFe", 
                                    "dP", 
                                    "dNH", 
                                    "dSO",
                                    "dNO",
                                    "oxic_state", 
                                    "treated", 
                                    "Location") %>%  
  melt( id.vars = c("Core.photometric","time","oxic_state","treated","Location"))

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
l <- list(c("FeTot","Preal","HS","NH4","SO","NO"), c("dFe",  "dP", "dSH",  "dNH", "dSO", "dNO")) 
for (i in 1:2) {
  
  graph <- 
    ggplot(transform(d[[i]],
                     variable = factor(variable, levels = l[[i]])), 
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
  
  ggsave(paste("IC_incubations_",i,".eps",sep=""), plot =graph, path = path.expand(here("index","figures")),
         
         width =30, height = 27,units = "cm",dpi = 600)
}

####== Benthic flux calculations ==####

# slope of NH4 line

inc_flux <- inc_data %>% mutate(NH4mol = (0.001*NH4*`volume_ml`)/(pi*0.09),
                                                      Pmol = (0.001*Preal*`volume_ml`)/(pi*0.09),
                                                      Femol = (0.001*FeTot*`volume_ml`)/(pi*0.09),
                                                      HSmol = (0.001*HS*`volume_ml`)/(pi*0.09),
                                                      NOmol = (0.001*NO*`volume_ml`)/(pi*0.09),
                                                      SOmol = (0.001*SO*`volume_ml`)/(pi*0.09)) 
fluxmin <- c(0,0,7)
fluxmax <- c(10, 25,60)
fluxan <- list(c("NH4", "Fe"), c("P"),c("HS","NH4","SO","Fe") )
fluxox <- list(c("NH4"), c("P"), c("SO","NO") )
fluxslopes <- tibble(row.names = c("A","A","B","B","C","C","D","D"))
for (i in 1:3) {
  d <- inc_flux[inc_flux$time < fluxmax[[i]] & inc_flux$time > fluxmin[[i]] ,] %>%  
    group_by(Core.photometric) %>% 
    summarize(NH4 = coef(lm(NH4mol ~ time))[2],
              P   = coef(lm(Pmol ~ time))[2],
              Fe = coef(lm(Femol ~ time))[2],
              HS = coef(lm(HSmol ~ time))[2],
              NO = coef(lm(NOmol ~ time))[2],
              SO = coef(lm(SOmol ~ time))[2]
    )
 
  an <- filter(d, str_detect(Core.photometric,"3\\.")) %>% select(fluxan[[i]]) 
  colnames(an) <- paste0(fluxan[[i]], c(rep(fluxmax[[i]], length(fluxan[[i]]))), c(rep("An",length(fluxan[[i]]))))
  ox <- filter(d, str_detect(Core.photometric,"2\\.")) %>% select(fluxox[[i]])   
  colnames(ox) <- paste0(fluxox[[i]], c(rep(fluxmax[[i]], length(fluxox[[i]]))), c(rep("Ox",length(fluxox[[i]]))))
  
  fluxslopes  <- bind_cols(fluxslopes, an, ox)
  
  
}

fluxslopes <- relocate(fluxslopes, row.names, Fe10An, Fe60An, P25An, NH410An, NH460An, HS60An, SO60An, P25Ox, NH410Ox, SO60Ox, NO60Ox)

#view(fluxslopes)

print(xtable(fluxslopes, type = "latex"), file = "slopes.tex")

O2slope = c(-0.172133,
            -0.404987,
            -0.360635,
            -0.0793467,
            -0.370374,
            -0.174982,
            -0.227892,
            -0.167105
)

 O2_umol_day = -1000*(O2slope*24)/31.9988 
 C_mg_day = -30*(O2slope*24)/31.9988


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
