library(readxl)
library(here)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(gt)


#import sequential extraction data.
Seq_extr_data <- read_excel(here("index","data","Seq_extr_data.xlsx"))

####==== Define parameters and functions for calculations ====####


# define detection limits as given by ICP report
Detlim_Fe_MgCl_A <- 0.014
Detlim_Fe_MgCl_B <- 0.023
Detlim_Fe_HCl <- 0.021
Detlim_Fe_CDB <- 0.1071
Detlim_Fe_HNO <- 0.035
Detlim_P_MgCl_A <- 0.101
Detlim_P_MgCl_B <- 0.065
Detlim_P_HCl <- 0.043
Detlim_P_CDB <- 0.0976
Detlim_P_HNO <- 0.057

# density of extractants and molar weight of Fe and P
d_MgCl <- 1.0467
d_pyroP <- 1.0402
d_HCl <- 1.0152
d_CDB <- 1.0887
d_HNO <- 1.414
mw_Fe <- 55.845
mw_P <- 30.973762


# define functions to calculate umol/g sediment from raw ICP and colourimetric data.
ppm_to_umol_p_g <- function(mol_weight, 
                            Dens, 
                            total_weight, 
                            sample_weight, 
                            Dil = 10, 
                            detlim, 
                            ppm){
  
  conc <- case_when(ppm/Dil < detlim ~ 0, 
                    T ~ ((total_weight/Dens)*ppm)/(sample_weight*mol_weight) 
  )
  return(conc)
}

mM_to_umol_p_g<- function(  Dens, 
                            total_weight, 
                            sample_weight,
                            mM){
  
  conc <- case_when(mM < 0 ~ 0,
                    T ~ ((total_weight/Dens)*mM)/(sample_weight) )
  
  return(conc)
}

######====== Calculate abundancies of P and Fe in the fractions in umol/g ======######

####-------- Calculate umol/g for all fractions ---------####

# Transform raw data.
Seq_extr_trans <-  Seq_extr_data %>% 
  mutate(
    Location = case_when(
      str_detect(Station,"A")  ~ "A",
      str_detect(Station,"-B")  ~ "B",
      str_detect(Station,"-C")  ~ "C",
      str_detect(Station,"-D")  ~ "D",
      T ~ "blank"),
    Incubation = case_when(
      str_detect(Station, "1\\.") ~ "before",
      str_detect(Station, "3\\.")  ~ "after",
      T ~ "blank" ),
    Core = str_sub(Station, end = 7),
    Fe_MgCl_A_calc = ppm_to_umol_p_g(mol_weight = mw_Fe, 
                                     Dens = d_MgCl, 
                                     Dil = 5,
                                     total_weight = weight_MgCl_A, 
                                     sample_weight = sample_weight_A, 
                                     detlim = Detlim_Fe_MgCl_A,
                                     ppm = ICP_MgCl_A_Fe),
    Fe_MgCl_B_calc = ppm_to_umol_p_g(mol_weight = mw_Fe, 
                                     Dens = d_MgCl, 
                                     Dil = 5,
                                     total_weight = weight_MgCl_B, 
                                     sample_weight = sample_weight_B, 
                                     detlim = Detlim_Fe_MgCl_B,
                                     ppm = ICP_MgCl_B_Fe),
    Fe_pyroP_calc = ppm_to_umol_p_g(mol_weight = mw_Fe, 
                                    Dens = d_pyroP,
                                    Dil = 20,
                                    total_weight = weight_pyroP, 
                                    sample_weight = sample_weight_B, 
                                    detlim = Detlim_Fe_HNO,
                                    ppm = ICP_pyroP_Fe),
    Fe_HCl_calc = ppm_to_umol_p_g(mol_weight = mw_Fe, 
                                  Dens = d_HCl, 
                                  total_weight = weight_HCl, 
                                  sample_weight = sample_weight_A, 
                                  detlim = Detlim_Fe_HCl,
                                  ppm = ICP_HCl_Fe),
    Fetot_HCl_calc = mM_to_umol_p_g(Dens = d_HCl,
                                    total_weight = weight_HCl,
                                    sample_weight = sample_weight_A,
                                    mM = Col_HCl_Fetot),
    FeII_HCl_calc = mM_to_umol_p_g(Dens = d_HCl,
                                   total_weight = weight_HCl,
                                   sample_weight = sample_weight_A,
                                   mM = Col_HCl_FeII),
    FeIII_HCl_calc = mM_to_umol_p_g(Dens = d_HCl,
                                    total_weight = weight_HCl,
                                    sample_weight = sample_weight_A,
                                    mM = Col_HCl_FeIII),
    Fe_CDB_calc = ppm_to_umol_p_g(mol_weight = mw_Fe, 
                                  Dens = d_CDB, 
                                  total_weight = weight_CDB, 
                                  sample_weight = sample_weight_A, 
                                  detlim = Detlim_Fe_CDB,
                                  ppm = ICP_CDB_Fe),
    Fe_HNO_calc = ppm_to_umol_p_g(mol_weight = mw_Fe, 
                                  Dens = d_HNO, 
                                  total_weight = weight_HNO, 
                                  sample_weight = sample_weight_A, 
                                  detlim = Detlim_Fe_HNO,
                                  ppm = ICP_HNO_Fe*10),
    P_MgCl_A_calc = ppm_to_umol_p_g(mol_weight = mw_P, 
                                    Dens = d_MgCl, 
                                    Dil = 5,
                                    total_weight = weight_MgCl_A, 
                                    sample_weight = sample_weight_A, 
                                    detlim = Detlim_P_MgCl_A,
                                    ppm = ICP_MgCl_A_P),
    P_MgCl_B_calc = ppm_to_umol_p_g(mol_weight = mw_P, 
                                    Dens = d_MgCl, 
                                    Dil = 5,
                                    total_weight = weight_MgCl_B, 
                                    sample_weight = sample_weight_B, 
                                    detlim = Detlim_P_MgCl_B,
                                    ppm = ICP_MgCl_B_P),
    P_HCl_calc = ppm_to_umol_p_g(mol_weight = mw_P, 
                                 Dens = d_HCl, 
                                 total_weight = weight_HCl, 
                                 sample_weight = sample_weight_A, 
                                 detlim = Detlim_P_HCl,
                                 ppm = ICP_HCl_P),
    P_CDB_calc = ppm_to_umol_p_g(mol_weight = mw_P, 
                                 Dens = d_CDB, 
                                 total_weight = weight_CDB, 
                                 sample_weight = sample_weight_A, 
                                 detlim = Detlim_P_CDB,
                                 ppm = ICP_CDB_P),
    P_HNO_calc = ppm_to_umol_p_g(mol_weight = mw_P, 
                                 Dens = d_HNO, 
                                 total_weight = weight_HNO, 
                                 sample_weight = sample_weight_A, 
                                 detlim = Detlim_P_HNO,
                                 ppm = ICP_HNO_P*10)
  ) 

####-------------- transform data into data to be plotted -----------####

# Middle all duplicates
calcdata <- c("Fe_MgCl_A_calc",
                 "Fe_MgCl_B_calc",
                 "FeII_HCl_calc",
                 "FeIII_HCl_calc",
                 "Fe_CDB_calc",
                 "Fe_HNO_calc",
                 "Fe_pyroP_calc",
                 "Fe_HCl_calc",
                 "Fetot_HCl_calc",
                 "P_MgCl_A_calc",
                 "P_MgCl_B_calc",
                 "P_HCl_calc",
                 "P_CDB_calc",
                 "P_HNO_calc")
meannames <- str_sub(calcdata[[i]], end = -6)
molnames <- 
mol <- function(x,y){mean(x) * mean(y) }

pw_data <- filter(ICplotdata, str_detect(sample, "pw")) %>% mutate(Station = str_replace(sample, "pw-" , ""),
                                                                   FeIII = Fetot-FeII) 

Seq_extr_dr <-Seq_extr_trans %>% 
  filter(Incubation != "blank") %>% 
  group_by(Station) %>% 
  mutate_at(calcdata, list(~mean(.), ~ mol(., y = tot_weight))) %>% 
  rename_with(~ str_replace_all(., c("_calc_mean" = "", "_calc_mol" = "_umol"))) %>% 
  mutate(
         Fe_MgCl = (mean(Fe_MgCl_A_calc) + mean(Fe_MgCl_B_calc))/2,
         Fe_MgCl_umol = ((mean(Fe_MgCl_A_calc) + mean(Fe_MgCl_B_calc))/2)*mean(tot_weight),
         Fe_HClminPP = mean(Fetot_HCl_calc)-mean(Fe_pyroP_calc),
         Fe_HClminPP_umol = (mean(Fetot_HCl_calc)-mean(Fe_pyroP_calc))*mean(tot_weight),
         P_MgCl = (mean(P_MgCl_A_calc) + mean(P_MgCl_B_calc))/2,
         PFe_MgCl = case_when(((mean(Fe_MgCl_A_calc) + mean(Fe_MgCl_B_calc))/2) > 1 & 
                              ((mean(P_MgCl_A_calc) + mean(P_MgCl_B_calc))/2) ~ 
                                                                    ((mean(P_MgCl_A_calc) + mean(P_MgCl_B_calc))/2)/
                                                                    ((mean(Fe_MgCl_A_calc) + mean(Fe_MgCl_B_calc))/2),
                              T ~ 0),
         PFe_HCl = mean(P_HCl_calc)/mean(Fetot_HCl_calc),
         PFe_CDB = mean(P_CDB_calc)/mean(Fe_CDB_calc),
         PFe_HNO = mean(P_HNO_calc)/mean(Fe_HNO_calc),
         FeP_HCl = mean(Fetot_HCl_calc)/mean(P_HCl_calc),
         FeP_CDB = mean(Fe_CDB_calc)/mean(P_CDB_calc)) %>% 
  distinct(Station, .keep_all = TRUE) %>% 
  ungroup %>% 
  left_join(pw_data, 
            copy = FALSE, 
            suffix = c(".seq", ".pw"),
                        keep = FALSE)



# create molten dataframes for figures 
seq_extr_A_Fe_molten <- Seq_extr_dr %>% 
  select("depth","Location","Incubation", "FeII", "FeIII", Fe_MgCl_umol, Fe_pyroP_umol,  Fe_HClminPP_umol, Fe_CDB_umol, Fe_HNO_umol) %>% 
  melt(id.vars = c("depth","Location","Incubation", "FeII", "FeIII"), variable.name = "Fraction" , na.rm = TRUE) 

seq_extr_B_Fe_molten <- Seq_extr_dr %>% 
  select("depth","Location","Incubation", "FeII", "FeIII", Fe_MgCl, Fe_pyroP, Fe_HClminPP, Fe_CDB, Fe_HNO) %>% 
  melt(id.vars = c("depth","Location","Incubation", "FeII", "FeIII"), variable.name = "Fraction" , na.rm = TRUE) 

seq_extr_A_P_molten <- Seq_extr_dr %>% 
  select("depth","Location","Incubation", "FeII", "FeIII", P_MgCl, P_HCl, P_CDB, P_HNO) %>% 
  melt(id.vars = c("depth","Location","Incubation", "FeII", "FeIII"), variable.name = "Fraction" , na.rm = TRUE) 

seq_extr_untreated_Fe_molten <- Seq_extr_dr %>% filter(Location == "B" | Location == "D") %>% 
  select("depth","Location","Incubation", "FeII", "FeIII", Fe_MgCl, Fe_pyroP, Fe_HClminPP, Fe_CDB, Fe_HNO) %>% 
  melt(id.vars = c("depth","Location","Incubation", "FeII", "FeIII"), variable.name = "Fraction" , na.rm = TRUE) 

seq_extr_PFe_molten <- Seq_extr_dr %>% 
  select("depth","Location","Incubation",  PFe_HCl, PFe_CDB, PFe_HNO) %>% 
  melt(id.vars = c("depth","Location","Incubation"), variable.name = "Fraction" , na.rm = TRUE) 

seq_extr_FeP_molten <- Seq_extr_dr %>% 
  select("depth","Location","Incubation",  FeP_HCl, FeP_CDB) %>% 
  melt(id.vars = c("depth","Location","Incubation"), variable.name = "Fraction" , na.rm = TRUE) 

######====== Create figures from data ======######



# create sequential extraction figures
list_sq_extr_plots <-  list(seq_extr_A_Fe_molten, 
                            seq_extr_B_Fe_molten, 
                            seq_extr_A_P_molten,
                            seq_extr_untreated_Fe_molten)

legend_levels <- list(c("Fe_MgCl_umol", "Fe_pyroP_umol",  "Fe_HClminPP_umol", "Fe_CDB_umol", "Fe_HNO_umol"),
                      c("Fe_MgCl", "Fe_pyroP", "Fe_HClminPP", "Fe_CDB", "Fe_HNO"),
                      c("P_MgCl", "P_HCl", "P_CDB", "P_HNO"))
legend_labels <- list(c("MgCl; loosly adsorbed Fe",
                        "Pyrophosphate; organic bound Fe",
                        "HCl; Reactive Fe minerals",
                        "CBD; Crystalline Fe oxides",
                        "HNO[3]; Pyrite"),
                      c("MgCl;\n loosly adsorbed Fe",
                        "Pyrophosphate;\n organic bound Fe",
                        "HCl;\n Reactive Fe minerals",
                        "CBD;\n Crystalline Fe oxides",
                        "HNO[3];\n Pyrite"),
                      c("MgCl; loosly adsorbed Fe",
                        "HCl; Reactive Fe minerals",
                        "CBD; Crystalline Fe oxides",
                        "HNO[3]; Pyrite"),
                      c("MgCl;\n loosly adsorbed Fe",
                        "Pyrophosphate;\n organic bound Fe",
                        "HCl;\n Reactive Fe minerals",
                        "CBD;\n Crystalline Fe oxides",
                        "HNO[3];\n Pyrite"))

colors <- list(c("orange","coral4", "darkorange2","darkgoldenrod1","goldenrod"),
               c("orange","coral4", "darkorange2","darkgoldenrod1","goldenrod"),
               c("chartreuse", "darkolivegreen3", "palegreen2", "olivedrab1"),
               c("orange","coral4", "darkorange2","darkgoldenrod1","goldenrod"))

axtitles <-list(expression(paste("Extracted Fe (",mu,"mol",")")),
                expression(paste("Extracted Fe (",mu,"mol gDW"^-1,")")),
                expression(paste("Extracted P (",mu,"mol gDW"^-1,")")),
                expression(paste("Extracted Fe (",mu,"mol gDW"^-1,")")))


for (i in 1:length(list_sq_extr_plots)) {
  
  j <- ggplot(transform( list_sq_extr_plots[[i]], 
                         Incubation = factor(Incubation, 
                                             levels = c("before","after"), 
                                             labels = c("Before incubation","After incubation"))), 
              aes(fill=Fraction, 
                  y=value, 
                  x=depth),
              scale_x_discrete(position = 'top'))+
    geom_area(
      aes(fill = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE))+
    geom_point(
      aes(),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    geom_line(
      aes(),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    scale_fill_manual(values = colors[[i]], labels = legend_labels[[i]], name = "Extractant; Target pool") +
    #geom_line( mapping = 
    #  aes(x = depth, y = FeII*10,), color = "sienna2", size = 1.5) +
    #geom_line( mapping = 
    #            aes(x = depth, y = FeIII*10,), color = "tomato3", size = 1.5)+
    coord_flip()+
    scale_x_reverse() + 
    scale_y_continuous(breaks=NULL,labels=waiver(),name=NULL, sec.axis = sec_axis(~ . *1,name=axtitles[[i]]))+
    labs(y = axtitles[[i]],
         x = "Depth (cm)"
         ) +
    theme_bw()+
    #scale_fill_brewer(palette = "Reds",direction = -1)+
    theme(
          axis.title=element_text(size=20),plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=16,angle = 0, hjust = 0.5),
          #panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          strip.text.y   = element_text(size=16, angle = 0),
          strip.text.x   = element_text(size=17, angle = 0),
          legend.title = element_text( size = 20),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.3, "cm")) +
    facet_grid(Location ~ Incubation)
  
  ggsave(paste("seq_extr_Fe_",i,".eps",sep=""), plot =j, path = path.expand(here("index","figures")),
         
         width =27, height = 40,units = "cm",dpi = 600)
  
}


# Create figure show P/Fe ratio
Lables = as_labeller( c(  "PFe_MgCl" = "MgCl", 
                          "PFe_HCl" = "HCl",
                          "PFe_CDB" = "CBD", 
                          "PFe_HNO" = "HNO[3]",
                          "after" = "After",
                          "before" = "Before") , default =  label_parsed )
ratioplots <- list(seq_extr_PFe_molten, seq_extr_FeP_molten)
p <- list()
for (i in 1:2) {
  
p[[i]] <-  ggplot(transform(ratioplots[[i]], 
                                      Incubation = factor(Incubation, levels = c("before","after"))), 
                            mapping = aes(
  x = value,
  y = depth ,
  color = Location,
  shape = Location)
) +
  
  scale_color_manual( values=c("A" = "red", "B" = "blue", "C" = "red", "D" = "blue")) +
  scale_shape_manual( values=c("A" = 16, "B" = 1, "C" = 17,  "D" = 2)) +
  scale_y_reverse() +
  geom_point() + 
  geom_path() +
  facet_grid(Incubation ~ Fraction, scales = "free_y", space = "free_y", labeller = Lables) +
  
  labs(  x = expression(paste("P/Fe")), 
         y = "Depth in cm",
         title = "P/Fe ratio per pool") +
  theme_classic() +
  theme( plot.title = element_text(size=20),
    axis.text.y   = element_text(size=14),
    axis.text.x   = element_text(size=14),
    axis.title.y  = element_text(size=14),
    axis.title.x  = element_text(size=14),
    strip.text.y   = element_text(size=14, angle = 0),
    strip.text.x   = element_text(size=14, angle = 0),
    panel.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
}
 
 PFe_profile_plots <- p[[1]]
 show(p[[2]])
show(PFe_profile_plots)

ggsave("seq_extr_PFe.eps", plot = PFe_profile_plots, path = path.expand(here("index","figures")),
       
       width =25, height = 25,units = "cm",dpi = 600)

#########========== statistics and P/Fe ratio calculations =========#########


# Create a table with P/Fe ratio statistics
PFe_range <- seq_extr_PFe_molten %>% mutate(F_range = case_when(
  Incubation == "before" & Fraction == "PFe_HCl" & depth < 8 & (Location == "A"| Location =="C") ~ "HClshallowT",
  Incubation == "before" & Fraction == "PFe_HCl" & depth < 8 & (Location == "B"| Location =="D") ~ "HClshallowN",
  Incubation == "before" & Fraction == "PFe_HCl" & depth > 20 & (Location == "A"| Location =="C") ~ "HCldeepT",
  Incubation == "before" & Fraction == "PFe_HCl" & depth > 20 & (Location == "B"| Location =="D") ~ "HCldeepN",
  Incubation == "before" & Fraction == "PFe_CDB" & depth < 8 & (Location == "A"| Location =="C") ~ "CBDshallowT",
  Incubation == "before" & Fraction == "PFe_CDB" & depth < 8 & (Location == "B"| Location =="D") ~ "CBDshallowN",
  Incubation == "before" & Fraction == "PFe_CDB" & depth > 20 & (Location == "A"| Location =="C") ~ "CBDdeepT",
  Incubation == "before" & Fraction == "PFe_CDB" & depth > 20 & (Location == "B"| Location =="D") ~ "CBDdeepN",
  Incubation == "before" & Fraction == "PFe_HNO" & depth < 8 & (Location == "A"| Location =="C") ~ "HNOshallowT",
  Incubation == "before" & Fraction == "PFe_HNO" & depth < 8 & (Location == "B"| Location =="D") ~ "HNOshallowN",     
  Incubation == "before" & Fraction == "PFe_HNO" & depth > 20 & (Location == "A"| Location =="C") ~ "HNOdeepT",
  Incubation == "before" & Fraction == "PFe_HNO" & depth > 20 & (Location == "B"| Location =="D") ~ "HNOdeepN",
  T ~ "A")) 
group_by(PFe_range, F_range) %>% summarise(mean = mean(value), sd = sd(value)) 


PFe_Table <- matrix(nrow = 6, ncol = 5)
l <-  list("HClshallow","HCldeep","CBDshallow","CBDdeep","HNOshallow","HNOdeep")
for (i in 1:6) {
  v = 
    t.test(  pull(PFe_range[PFe_range$F_range == paste(l[i],"T", sep = ""), ], value),
             pull(PFe_range[PFe_range$F_range == paste(l[i],"N", sep = ""), ], value)
    )$p.value
  m1 = mean(pull(PFe_range[PFe_range$F_range == paste(l[i],"T", sep = ""), ], value))
  m2 = mean(pull(PFe_range[PFe_range$F_range == paste(l[i],"N", sep = ""), ], value))
  sd1 = sd(pull(PFe_range[PFe_range$F_range == paste(l[i],"T", sep = ""), ], value))
  sd2 = sd(pull(PFe_range[PFe_range$F_range == paste(l[i],"N", sep = ""), ], value))
  
  PFe_Table[i,5] <- v
  PFe_Table[i,2] <- sd1
  PFe_Table[i,4] <- sd2
  PFe_Table[i,1] <- m1
  PFe_Table[i,3] <- m2
} 
row_names <- c("mean p/Fe ","st. dev.","mean p/Fe","st. dev.", "p value")
col_names <- c("HClshallow","HCldeep","CBDshallow","CBDdeep","HNOshallow","HNOdeep")
PFe_matrix <- t(PFe_Table)
colnames(PFe_matrix) <- col_names
rownames(PFe_matrix) <-  row_names
PFe_tib <- as_tibble( PFe_matrix, rownames = "depth:")



PFe_tbl <- 
  gt(PFe_tib)  %>% tab_header(
    title = "P/Fe ratios statistical parameters",
    subtitle = "for selected fractions obtained with sequential extractions"
  ) %>% 
  fmt_number(
    columns = 2:7,
    decimals = 4) %>% 
  tab_row_group(
      label = "treated:",
    rows = 1:2
    ) %>% 
  tab_row_group(
    label = "nontreated:",
    rows = 3:4
  ) %>% 
  tab_row_group(
    label = "t test:",
    rows = 5
  ) %>% 
  row_group_order(c("treated:","nontreated:","t test:")) %>% 
  tab_spanner(
    label = "HCl",
    columns = c("HClshallow","HCldeep")
  ) %>%
  tab_spanner(
    label = "CBD",
    columns = c("CBDshallow","CBDdeep")
  ) %>%
  tab_spanner(
    label = "HNO",
    columns = c("HNOshallow","HNOdeep")
  ) %>% 
  tab_source_note(
    source_note = "shallow: top 8 cm; deep: below 20 cm"
  ) %>% 
  cols_label( "HClshallow" = "shallow",
              "HCldeep" = "deep",
              "CBDshallow" = "shallow",
              "CBDdeep" = "deep",
              "HNOshallow" = "shallow",
              "HNOdeep" = "deep") %>% 
  tab_options(
    row_group.background.color = "#FFEFDB"
  )

gtsave(PFe_tbl, "PFe_table.png",  path = path.expand(here("index","figures")))


# some statistical analyses


eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(
    as.expression(
      substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))
    )
  )
}

eqdepth <- Seq_extr_dr %>% filter(depth <= 10 & depth > 0)

loc  <-  c("A","B","C","D")
frac <- c("Fe_MgCl", "Fe_pyroP", "Fe_HCl_col", "Fe_CDB", "Fe_HNO")
tab <- matrix(nrow = 5, ncol = 5)
for (i in 1:4) {
  for (j in 1:5) {
    
    
    t = t.test(
      pull(eqdepth[eqdepth$Incubation == "before" & eqdepth$Location == loc[[i]],], frac[[j]]),
      pull(eqdepth[eqdepth$Incubation == "after" & eqdepth$Location == loc[[i]],], frac[[j]]),
      paired = TRUE)$p.value
    
    tab[j,i] <- t
  }
}
colnames(tab) <- loc
rownames(tab) <- frac


pyr_bp = ggplot(eqdepth,
                mapping = aes(y = Fe_HNO,
                              x = Core)) + geom_boxplot()
crys_bp = ggplot(eqdepth,
                 mapping = aes(y = Fe_CDB,
                               x = Core)) + geom_boxplot()
HCl_bp = ggplot(eqdepth,
                mapping = aes(y = Fe_HCl_col,
                              x = Core)) + geom_boxplot()
org_bp = ggplot(eqdepth,
                mapping = aes(y = PyroP,
                              x = Core)) + geom_boxplot()
salt_bp = ggplot(eqdepth,
                 mapping = aes(y = Fe_MgCl_A,
                               x = Core)) + geom_boxplot()


eqdepth %>% group_by(Core) %>% summarise(Location = first(Location), Incubation = first(Incubation), Pyrite = mean(Fe_HNO), sumPyrite = sum(Fe_HNO))

# test colorimetric Fe data against ICP data
instrument_seq_extr_Fe <- ggplot(Seq_extr_dr, 
                      mapping = aes( 
                        y=Fe_HCl_ICP, 
                        x=Fe_HCl_col),
                      scale_x_discrete(position = 'top')) +
  geom_point(
    stat = "identity", position = position_stack(reverse=TRUE),color="black") +
  stat_smooth(method = "lm")+
  geom_text(x = 25, y = 300, label = eq(Seq_extr_dr$Fe_HCl_col, Seq_extr_dr$Fe_HCl_ICP), parse = TRUE)
show(instrument_seq_extr_Fe)
