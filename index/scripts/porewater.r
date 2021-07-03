
library(readxl)
library(here)
library(tidyverse)
library(ggplot2)
library(reshape2)

# import raw data


IC_raw_data <- read_excel(path.expand(here("index","data", "IC_raw_data.xlsx")))
pw_photometric <-  read_excel(path.expand(here("index","data", "pw_P_photometric.xlsx")))



####=== Statistics of duplicates in IC data. ===####

stat_IC <- IC_raw_data[duplicated(IC_raw_data$sample)|duplicated(IC_raw_data$sample, fromLast=TRUE),] %>%
  group_by(sample) %>% 
  summarize(P = mean(Phosphate),
            "SO" = mean(Sulphate), 
            NO = mean(Nitrate), 
            "Br" = mean(Bromide), 
            "Fl" = mean(Fluoride), 
            "Cl" = mean(Chloride),
            DP = abs(diff(Phosphate)),
            "DSO" = abs(diff(Sulphate)), 
            DNO = abs(diff(Nitrate)), 
            "DBr" = abs(diff(Bromide)), 
            "DFl" = abs(diff(Fluoride)), 
            "DCl" = abs(diff(Chloride))) 
hist_IC <-  stat_IC %>% select("DP":"DCl") %>% melt() %>% ggplot( aes(x = value)) + geom_histogram() + facet_grid(~variable)


IC_parameters <- c("P", "SO", "NO", "Br", "Fl", "Cl")
stat_IC_table <-  data.frame(parameter = IC_parameters)
IC_diff <- list()
for (i in 1:length(IC_parameters)) {
  IC_diff[[i]] <- paste("D", IC_parameters[[i]], sep = "")
  m = mean(pull(stat_IC, IC_parameters[[i]]))
  d = mean(pull(stat_IC, IC_diff[[i]]))
  var = mean(.5*(pull(stat_IC, IC_diff[[i]])^2))
  sd = sqrt( mean(.5*pull(stat_IC, IC_diff[[i]])^2))
  
  stat_IC_table$Mean[i] <- m
  stat_IC_table$Mean_diff[i] <-  d
  stat_IC_table$Variance[i] <- var
  stat_IC_table$Mean_sd[i] <-  sd
  
}


####== Transform the raw data. ==####


ICplotdata  <- IC_raw_data %>% 
  filter(str_detect(sample,"pw|sw")) %>%                  #select only porewater from all IC data
  group_by(sample) %>% 
  mutate(P = (mean(Phosphate)*1000)/94.9714,              #middle duplicates and transform from mg/L to umol/L
         "SO" = mean(Sulphate)*1000/96.06, 
         NO = mean(Nitrate)*1000/62.004, 
         "Br" = mean(Bromide)*1000/79.9,
         "Fl" = mean(Fluoride)*1000/19, 
         "Cl" = mean(Chloride)*1000/35.45 ) %>%
  distinct(sample, .keep_all = TRUE) %>% 
  mutate(treated =                                        #add columns with treated status and befor/after incubation
           case_when(
             str_detect(sample,"A|C")  ~ "treated",
             str_detect(sample,"B|D")  ~ "nontreated"),
         Incubation = 
           case_when(
             str_detect(sample, "1\\.") ~ "before",
             str_detect(sample, "3\\.")  ~ "after",
             T ~ "blank" )) %>% 
  left_join(   pw_photometric,                           # combine IC data with photometry data
               by = c("sample" = "Sample"),
               copy = FALSE,
               suffix = c(".IC", ".photometric"),
               keep = FALSE
             ) %>%  
  group_by(location) %>% arrange(cm_below_swi) 

####== Create pw profile plots ==####

# reshape data.
Parameters <-  c(  
  "SO", 
  "NO",
  "Br",
  "Fl", 
  "Cl")
molten_IC_core_data <- ICplotdata %>% select("location", "sample", "cm_below_swi", "SO" , "NO", "treated", "Incubation", "P_phot":"NH","SH") %>% 
  melt(id.vars = c("location", "sample", "cm_below_swi", "treated", "Incubation"), 
       na.rm = TRUE) 

# plot profiles by parameter

# lable functions
Lables = as_labeller( c( "P_phot"  = "PO[4]^'3-'",                      
                         "SO" = "SO[4]^'2-'", 
                         "NO" = "NO[3]^'-'",
                         "Br" = "Br^'-'",
                         "Fl" = "F^'-'", 
                         "Cl" = "Cl^'-'",
                         "Fetot" = "Fe[tot]",
                         "FeII" = "Fe^'2-'",
                         "NH" = "NH[4]^'+'",
                         "SH" = "SH^'-'",
                         "before" = "'Before'\\1\n'incubation'",
                         "after" = "'After'\\1\n'incubation'"), default =  label_parsed )

# scale of the y axis 
max<-max(na.omit(molten_IC_core_data$cm_below_swi))
limitsx <- ceiling(max+2)

# create plot
IC_profile_plots <- 
  ggplot(  transform(molten_IC_core_data, 
                     Incubation = factor(Incubation, 
                                         levels = c("before","after"), 
                                         labels = c("Before incubation","After incubation")),
                     variable =   factor(variable, 
                                       levels = c("Fetot","P_phot","SH","NH","SO","NO") 
                                       )
                     ), 
           mapping = aes(
                          y = value,
                          x = cm_below_swi ,
                          color = location,
                          shape = location)
        ) +
  
  scale_color_manual( values=c("A" = "coral3", "B" = "deepskyblue4", "C" = "coral3", "D" = "deepskyblue4")) +
  scale_shape_manual( values=c("A" = 16, "B" = 1, "C" = 17,  "D" = 2)) +
  geom_point(size = 2) + 
  geom_path() +
  #scale_y_continuous(name=NULL, labels=NULL, breaks=NULL, limits=c(0,limitsy[i]),
  #                    sec.axis = sec_axis(~ . *1,name=expression(paste("(",mu,"mol/l",")")))) +
  xlab("Depth (cm)") +
  scale_x_reverse(breaks=seq(0,limitsx[1],by=2), labels=seq(0,limitsx[1],by=2), limits=c(limitsx[1],0))+
  coord_flip() +
  facet_grid(Incubation ~ variable, scales = "free_x", labeller = labeller(.rows = label_wrap_gen(width = 10), .cols = Lables )) +
  labs(  y = expression(paste("concentration in ",mu,"mol/L")), 
         x = "Depth in cm",
         title = "porewater profiles" ) +
    theme_bw()+theme(axis.title=element_text(size=20),
                     plot.title = element_text(size=25, face="bold", hjust = 0.5), #hjust is position (left/right/ middle =0.5)
                     axis.text=element_text(size=14,angle = 0, hjust = 0.5),
                     strip.text.y   = element_text(size=14, angle = 0),
                     strip.text.x   = element_text(size=14, angle = 0),
                     # panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.title = element_text(size = 25),
                     legend.text = element_text(size = 21))

ggsave(paste("profiles.png",sep=""), plot =IC_profile_plots, path = path.expand(here("index","figures")),
       width =40, height = 25,units = "cm",dpi = 600)
