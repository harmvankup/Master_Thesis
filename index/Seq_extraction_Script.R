library(ggplot2)
library(data.table)
library(nlstools)
library(nlsMicrobio)
getwd()
require(graphics)
library(readxl)
library(reshape2)
library(dplyr)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(egg)
library(patchwork)
library(tibble)
#=================================================================
#Defining densities of the liquids
Dens_MgCl=1.069
Dens_HCl=1.015
Dens_Ac=1.0724
Dens_CD=1.0943
Dens_Bicab=1.05
#For total P calculation
P_MW=30.97
##=======Reading Excel Data into the R-Script=============########
##Reading tube weights into the file

my_data <- read_excel("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/SEDEX_CoreInc_ToRead.xlsx",
                      sheet = "KB", range = NULL, col_names = TRUE,col_types = NULL)
Data.df<-data.frame(my_data)
#====================================================================
## Reading Calibrations into the file
CalHCl<- read_excel("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/SEDEX_CoreInc_ToRead.xlsx",
                      sheet = "HCl_calc", range = NULL, col_names = TRUE,col_types = NULL)
CalMgCl<- read_excel("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/SEDEX_CoreInc_ToRead.xlsx",
                      sheet = "MgCl_calc", range = NULL, col_names = TRUE,col_types = NULL)
CalAc<- read_excel("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/SEDEX_CoreInc_ToRead.xlsx",
                      sheet = "Ac_calc", range = NULL,col_names = TRUE,col_types = NULL)
CalBicab<- read_excel("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/SEDEX_CoreInc_ToRead.xlsx",
                    sheet = "BiCab_calc", range = NULL,col_names = TRUE,col_types = NULL)

colnames(CalHCl) <- c("PO4","Abs1")
colnames(CalMgCl) <- c("PO4","Abs1")
colnames(CalAc) <- c("PO4","Abs1")
colnames(CalBicab) <- c("PO4","Abs1")

model <- lm(Abs1 ~ 0+ PO4, data = CalHCl)
SlopeHCl<-as.numeric(model$coefficients[1])
summary(model)
model <- lm(Abs1 ~ 0+PO4, data = CalAc)
SlopeAc<-as.numeric(model$coefficients[1])
summary(model)
model <- lm(Abs1 ~ 0+PO4, data = CalMgCl)
SlopeMgCl<-as.numeric(model$coefficients[1])
summary(model)
model <- lm(Abs1 ~ 0+PO4, data = CalBicab)
SlopeBicab<-as.numeric(model$coefficients[1])
summary(model)

#Polishing the raw data, leaving only those data needed to make the calculations using the calibration
keeps <- c("Station","Week","Depth","SampleWeight","MgCl1_Ex","MgCl1_AU","BiCab_Ex","BiCab_AU","MgCl2_Ex","MgCl2_AU"
           ,"CD_Ex","CD_P","MgCl3_Ex","MgCl3_AU","Ac_Ex","Ac_AU","MgCl4_Ex","MgCl4_AU","HCl1_Ex","HCl1_AU",
           "HCl2_Ex","HCl2_AU")
colnames(Data.df)
Data.df<-Data.df[keeps]

#Creating Data frame to store converted values
#MgCl Dataframe
ValMgCl.df<-data.frame(
  Sample=Data.df$Station,
  Depth=Data.df$Depth,
  week=Data.df$Week,
  MgCl1=as.numeric(seq(0,0,length.out=length(Data.df[,1]))),
  MgCl2=as.numeric(seq(0,0,length.out=length(Data.df[,1]))),
  MgCl3=as.numeric(seq(0,0,length.out=length(Data.df[,1]))),
  MgCl4=as.numeric(seq(0,0,length.out=length(Data.df[,1])))
)
#Bicab Dataframe
ValBicab.df<-data.frame(
  Sample=Data.df$Station,
  Depth=Data.df$Depth,
  week=Data.df$Week,
  Bicab=seq(0,0,length.out=length(Data.df[,1]))
  )
#HCl Dataframe
ValHCl.df<-data.frame(
  Sample=Data.df$Station,
  Depth=Data.df$Depth,
  week=Data.df$Week,
  HCl1=seq(0,0,length.out=length(Data.df[,1])),
  HCl2=seq(0,0,length.out=length(Data.df[,1]))
)
#Ac Dataframe
ValAc.df<-data.frame(
  Sample=Data.df$Station,
  Depth=Data.df$Depth,
  week=Data.df$Week,
  Ac=seq(0,0,length.out=length(Data.df[,1]))
)
#DB Dataframe
ValCD.df<-data.frame(
  Sample=Data.df$Station,
  Depth=Data.df$Depth,
  week=Data.df$Week,
  CD=Data.df$CD_P
)
#Converting AU data in a micromol/g value
DF_MgCl=3.125		
DF_HCl=	12.5
DF_Ac=6.25
DF_Bicab=12.5
DF_CD=10
Extracted_MgCl<- function(AU,Dens=Dens_MgCl,Slope=SlopeMgCl,V,SampWeight,DF=DF_MgCl) {
  (((AU/Slope)*DF)*((V/Dens)/1000))/SampWeight
}
Extracted_BiCab<- function(AU,Dens=Dens_Bicab,Slope=SlopeBicab,V,SampWeight,DF=DF_Bicab) {
  (((AU/Slope)*DF)*((V/Dens)/1000))/SampWeight
}
Extracted_HCl<- function(AU,Dens=Dens_HCl,Slope=SlopeHCl,V,SampWeight,DF=DF_HCl) {
  (((AU/Slope)*DF)*((V/Dens)/1000))/SampWeight
}
Extracted_Ac<- function(AU,Dens=Dens_Ac,Slope=SlopeAc,V,SampWeight,DF=DF_Ac) {
  (((AU/Slope)*DF)*((V/Dens)/1000))/SampWeight
}

#DB is special as it is not converted from AU but rather from ICP mg/L value
Extracted_CD<-function(P_Weight=P_MW,Dens=Dens_CD,V,SampWeight,DF=DF_CD,CD_conc) {
  ((V/Dens/1000)*((CD_conc*DF)/(P_Weight/1000)))/SampWeight
}

#Converting AU data in a micromol/g value
for (i in 1:length(ValMgCl.df[,1] )){
  ValMgCl.df[i,4]<-Extracted_MgCl(AU=Data.df$MgCl1_AU[i],V=Data.df$MgCl1_Ex[i],SampWeight=Data.df$SampleWeight[i])
  ValMgCl.df[i,5]<-Extracted_MgCl(AU=Data.df$MgCl2_AU[i],V=Data.df$MgCl2_Ex[i],SampWeight=Data.df$SampleWeight[i])
  ValMgCl.df[i,6]<-Extracted_MgCl(AU=Data.df$MgCl3_AU[i],V=Data.df$MgCl3_Ex[i],SampWeight=Data.df$SampleWeight[i])
  ValMgCl.df[i,7]<-Extracted_MgCl(AU=Data.df$MgCl4_AU[i],V=Data.df$MgCl4_Ex[i],SampWeight=Data.df$SampleWeight[i])
}
for (i in 1:length(ValHCl.df[,1])){
  ValHCl.df[i,4]<-Extracted_HCl(AU=Data.df$HCl1_AU[i],V=Data.df$HCl1_Ex[i],SampWeight=Data.df$SampleWeight[i])
  ValHCl.df[i,5]<-Extracted_HCl(AU=Data.df$HCl2_AU[i],V=Data.df$HCl2_Ex[i],SampWeight=Data.df$SampleWeight[i])
}
for (i in 1:length(ValAc.df[,1])){
  ValAc.df[i,4]<-Extracted_Ac(AU=Data.df$Ac_AU[i],V=Data.df$Ac_Ex[i],SampWeight=Data.df$SampleWeight[i])
}
for (i in 1:length(ValCD.df[,1])){
  ValCD.df[i,4]<-Extracted_CD(V=Data.df$CD_Ex[i],SampWeight=Data.df$SampleWeight[i],CD_conc=Data.df$CD_P[i])
}
for (i in 1:length(ValCD.df[,1])){
  ValBicab.df[i,4]<-Extracted_BiCab(AU=Data.df$BiCab_AU[i],V=Data.df$BiCab_Ex[i],SampWeight=Data.df$SampleWeight[i])
}


#Creating Data frame to store Final values
Extracted.df<-data.frame(
  Sample=Data.df$Station,
  Depth=Data.df$Depth,
  week=Data.df$Week,
  Exch_P=ValMgCl.df[,4],
  Bicab_P=as.numeric(seq(0,0,length.out=length(Data.df[,1]))),
  CD_P=as.numeric(seq(0,0,length.out=length(Data.df[,1]))),
  Ca_P=as.numeric(seq(0,0,length.out=length(Data.df[,1]))),
  Detr_P=ValHCl.df[,4],
  Org_P=ValHCl.df[,5]
)
for (i in 1:length(ValAc.df[,1])){
  Extracted.df$Bicab_P[i]<-ValBicab.df$Bicab[i]+ValMgCl.df$MgCl2[i]
  Extracted.df$CD_P[i]<-ValCD.df$CD[i]+ValMgCl.df$MgCl3[i]
  Extracted.df$Ca_P[i]<-ValAc.df$Ac[i]+ValMgCl.df$MgCl4[i]
}
maxy<-c(seq(0,0,length.out = length(Extracted.df[,1])))
for (i in 1:length(Extracted.df[,1])){
  maxy[i]<-c(as.numeric(Extracted.df[i,4]+Extracted.df[i,5]+Extracted.df[i,6]+Extracted.df[i,7]+Extracted.df[i,8]+Extracted.df[i,9]))
}
#Averaging the reference 1 and 2 for week 0
Ref_av.df<-data.frame(matrix(ncol=14,nrow=10))
colnames(Ref_av.df)<-c("Sample","Depth","Ex_av","Ex_sd","Bicab_av","Bicab_sd","CD_av","CD_sd",
                       "Ac_ac","Ac_sd","Detr_av","Detr_sd","Org_av","Org_sd")
for (i in 1:10){
  for(x in 4:ncol(Extracted.df)){
    Extracted.df[i,x]<-mean(c(Extracted.df[i,x]),Extracted.df[(i+10),x])
    Ref_av.df[i,(2*(x-3)+1)]<-mean(c(Extracted.df[i,x],Extracted.df[(i+10),x]))
    Ref_av.df[i,(2*(x-3)+2)]<-sd(c(Extracted.df[i,x],Extracted.df[(i+10),x]))
    }
}
for(i in 1:length(Extracted.df$Sample)){
  if (Extracted.df$Sample[i]=="Reference 1"){
    Extracted.df$Sample[i]<-"Reference"
  }
  if (Extracted.df$Sample[i]=="Reference 2"){
    Extracted.df$Sample[i]<-"Reference"
  }
}
Extracted.df<-Extracted.df[-c(11:20),]
#Changing the week numbers into characters:

for (i in 1:length(Extracted.df$Sample)){
  if (Extracted.df$week[i]=="0"){
    Extracted.df$week[i]<-"Week 0"
  }
  if (Extracted.df$week[i]=="1"){
    Extracted.df$week[i]<-"Week 1"
  }
  if (Extracted.df$week[i]=="2"){
    Extracted.df$week[i]<-"Week 2"
  }
  if (Extracted.df$week[i]=="3"){
    Extracted.df$week[i]<-"Week 3"
  }
}
#First we separate the treatments into separate lists. 
#Next we will plot the individual weekly results from these lists. 
SEDEX_List_Reference<-list()
SEDEX_List_Fesalt<-list()
SEDEX_List_Alsalt<-list()
SEDEX_List_Febyprod<-list()
Weeks<-unique(Extracted.df$week)
Stations<-unique(Extracted.df$Sample)
for (i in 1:length(Stations)){
  df<-Extracted.df[Extracted.df$Sample %in% Stations[i],]
  if (Stations [i]=="AlSalt"){
    for (x in 1:length(Weeks)){
      SEDEX_List_Alsalt[[x]]<-df[df$week %in% Weeks[x],]}
  }
  if (Stations [i]=="FeSalt"){
    for (x in 1:length(Weeks)){
      SEDEX_List_Fesalt[[x]]<-df[df$week %in% Weeks[x],]}
  }
  if (Stations [i]=="FeByprod"){
    for (x in 1:length(Weeks)){
      SEDEX_List_Febyprod[[x]]<-df[df$week %in% Weeks[x],]}
  }
  if (Stations [i]=="Reference"){
    for (x in 1:length(Weeks)){
      SEDEX_List_Reference[[x]]<-df[df$week %in% Weeks[x],]}
  }
}
SEDEX_Graph_List_Reference<-list()
SEDEX_Graph_List_Fesalt<-list()
SEDEX_Graph_List_Alsalt<-list()
SEDEX_Graph_List_Febyprod<-list()
#========================================================================================================================
#REF
for (i in 1:length(SEDEX_List_Reference)){
  df<-SEDEX_List_Reference[[i]]
  df<-subset(df, select = -c(Sample,week))
  row.names(df)<-df$Depth
  df<-reshape2::melt(df, id.vars = "Depth")
  colnames(df)<-c("Depth", "Fraction","Extracted")
  SEDEX_List_Reference [[i]] <- df # save your dataframes into the list

  SEDEX_Graph_List_Reference[[i]]<-ggplot(SEDEX_List_Reference [[i]], aes(fill=Fraction, y=Extracted, x=Depth),scale_x_discrete(position = 'top'))+
          geom_area(
            aes(color = Fraction),
            stat = "identity", position = position_stack(reverse=TRUE))+
          geom_point(
            aes(color = Fractionn),
            stat = "identity", position = position_stack(reverse=TRUE),color="black")+
          geom_line(
            aes(color = Fraction),
            stat = "identity", position = position_stack(reverse=TRUE),color="black")+
          coord_flip()+scale_x_reverse()+xlab("Depth (mm)")+ylim(0,max(maxy))+
          scale_y_continuous(breaks=NULL,labels=NULL,name=NULL,sec.axis = sec_axis(~ . *1,name=expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")"))))+
          ylab(expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")")))+theme_classic()+
          scale_color_grey()+ scale_fill_brewer(palette = "Greens",direction = -1)+
          theme(axis.text.y   = element_text(size=14),
                axis.text.x   = element_text(size=14),
                axis.title.y  = element_text(size=14),
                axis.title.x  = element_text(size=14),
                panel.background = element_blank(),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"),
                panel.border = element_rect(colour = "black", fill=NA, size=1)
          )+ggtitle(paste(Weeks[i]))
}

totplot<-SEDEX_Graph_List_Reference[[1]]+SEDEX_Graph_List_Reference[[2]]+
  SEDEX_Graph_List_Reference[[3]]+SEDEX_Graph_List_Reference[[4]]+
  plot_layout(guides = "collect")
setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
file_name = paste("Reference_Total_4Weekincubated.png", sep="")
png(filename=file_name, width = 960 , height = 960,
    units = "px")
print(totplot)
dev.off()

for (i in 1:length(SEDEX_Graph_List_Reference)){
setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
  file_name = paste("Seq_extract_Ref_CoreIncWeek_",i,".png", sep="")
  png(filename=file_name, width = 480, height = 480,
    units = "px")
  print(SEDEX_Graph_List_Reference[[i]])
  dev.off()
}
#===================================================================================================================
#FeSalt

for (i in 1:length(SEDEX_List_Fesalt)){
  df<-SEDEX_List_Fesalt[[i]]
  df<-subset(df, select = -c(Sample,week))
  row.names(df)<-df$Depth
  df<-reshape2::melt(df, id.vars = "Depth")
  colnames(df)<-c("Depth", "Fraction","Extracted")
  SEDEX_List_Fesalt [[i]] <- df # save your dataframes into the list
  
  SEDEX_Graph_List_Fesalt[[i]]<-ggplot(SEDEX_List_Fesalt [[i]], aes(fill=Fraction, y=Extracted, x=Depth),scale_x_discrete(position = 'top'))+
    geom_area(
      aes(color = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE))+
    geom_point(
      aes(color = Fractionn),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    geom_line(
      aes(color = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    coord_flip()+scale_x_reverse()+xlab("Depth (mm)")+ylim(0,max(maxy))+
    scale_y_continuous(breaks=NULL,labels=NULL,name=NULL,sec.axis = sec_axis(~ . *1,name=expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")"))))+
    ylab(expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")")))+theme_classic()+
    scale_color_grey()+ scale_fill_brewer(palette = "Greens",direction = -1)+
    theme(axis.title=element_text(size=28),plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=24,angle = 0, hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.title = element_text( size = 20),
          legend.text = element_text(size = 20))+
  ggtitle(paste(Weeks[i]))
}

totplot<-SEDEX_Graph_List_Fesalt[[2]]+SEDEX_Graph_List_Fesalt[[3]]+
  SEDEX_Graph_List_Fesalt[[4]]+
  plot_layout(guides = "collect")
setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
file_name = paste("Fesalt_Total_4Weekincubated.png", sep="")
png(filename=file_name, width = 960 , height = 480,
    units = "px")
print(totplot)
dev.off()

for (i in 2:length(SEDEX_Graph_List_Fesalt)){
  setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
  file_name = paste("Seq_extract_FeSalt_CoreIncWeek_",i,".png", sep="")
  png(filename=file_name, width = 480, height = 480,
      units = "px")
  print(SEDEX_Graph_List_Fesalt[[i]])
  dev.off()
}

#===========================================================================================================
#Alsalt

for (i in 2:length(SEDEX_List_Alsalt)){
  df<-SEDEX_List_Alsalt[[i]]
  df<-subset(df, select = -c(Sample,week))
  row.names(df)<-df$Depth
  df<-reshape2::melt(df, id.vars = "Depth")
  colnames(df)<-c("Depth", "Fraction","Extracted")
  SEDEX_List_Alsalt [[i]] <- df # save your dataframes into the list
  
  SEDEX_Graph_List_Alsalt[[i]]<-ggplot(SEDEX_List_Alsalt [[i]], aes(fill=Fraction, y=Extracted, x=Depth),scale_x_discrete(position = 'top'))+
    geom_area(
      aes(color = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE))+
    geom_point(
      aes(color = Fractionn),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    geom_line(
      aes(color = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    coord_flip()+scale_x_reverse()+xlab("Depth (mm)")+ylim(0,max(maxy))+
    scale_y_continuous(breaks=NULL,labels=NULL,name=NULL,sec.axis = sec_axis(~ . *1,name=expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")"))))+
    ylab(expression(paste("Extracted P (",mu,"mol*g DW"^-1,")")))+theme_classic()+
    scale_color_grey()+ scale_fill_brewer(palette = "Greens",direction = -1)+
    theme(axis.text.y   = element_text(size=14),
          axis.text.x   = element_text(size=14),
          axis.title.y  = element_text(size=14),
          axis.title.x  = element_text(size=14),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+ggtitle(paste(Weeks[i]))
}

totplot<-SEDEX_Graph_List_Alsalt[[2]]+SEDEX_Graph_List_Alsalt[[3]]+
  SEDEX_Graph_List_Alsalt[[4]]+
  plot_layout(guides = "collect")
setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
file_name = paste("Alsalt_Total_4Weekincubated.png", sep="")
png(filename=file_name, width = 960 , height = 480,
    units = "px")
print(totplot)
dev.off()
for (i in 2:length(SEDEX_List_Alsalt)){
  setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
  file_name = paste("Seq_extract_Alsalt_CoreIncWeek_",i,".png", sep="")
  png(filename=file_name, width = 480, height = 480,
      units = "px")
  print(SEDEX_Graph_List_Fesalt[[i]])
  dev.off()
}




#========================================================================================================
#Febyprod
for (i in 2:length(SEDEX_List_Febyprod)){
  df<-SEDEX_List_Febyprod[[i]]
  df<-subset(df, select = -c(Sample,week))
  row.names(df)<-df$Depth
  df<-reshape2::melt(df, id.vars = "Depth")
  colnames(df)<-c("Depth", "Fraction","Extracted")
  SEDEX_List_Febyprod [[i]] <- df # save your dataframes into the list
  
  SEDEX_Graph_List_Febyprod[[i]]<-ggplot(SEDEX_List_Febyprod [[i]], aes(fill=Fraction, y=Extracted, x=Depth),scale_x_discrete(position = 'top'))+
    geom_area(
      aes(color = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE))+
    geom_point(
      aes(color = Fractionn),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    geom_line(
      aes(color = Fraction),
      stat = "identity", position = position_stack(reverse=TRUE),color="black")+
    coord_flip()+scale_x_reverse()+xlab("Depth (mm)")+ylim(0,max(maxy))+
    scale_y_continuous(breaks=NULL,labels=NULL,name=NULL,sec.axis = sec_axis(~ . *1,name=expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")"))))+
    ylab(expression(paste("Extracted P (",mu,"mol'*g DW"^-1,")")))+theme_classic()+
    scale_color_grey()+ scale_fill_brewer(palette = "Greens",direction = -1)+
    theme(axis.text.y   = element_text(size=14),
          axis.text.x   = element_text(size=14),
          axis.title.y  = element_text(size=14),
          axis.title.x  = element_text(size=14),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    )+ggtitle(paste(Weeks[i]))
}

totplot<-SEDEX_Graph_List_Febyprod[[2]]+SEDEX_Graph_List_Febyprod[[3]]+
  SEDEX_Graph_List_Febyprod[[4]]+
  plot_layout(guides = "collect")
setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
file_name = paste("Febyprod_Total_4Weekincubated.png", sep="")
png(filename=file_name, width = 960 , height = 480,
    units = "px")
print(totplot)
dev.off()
for (i in 2:length(SEDEX_Graph_List_Febyprod)){
  setwd("F:/PhD_Bayreuth/Experimental/KleinerBrombachsee/October2020_3rdVisit/ComparativeSedimentIncubation/SequentialExtractions/P_SEDEX/Graphs")
  file_name = paste("Seq_extract_FeByprod_CoreIncWeek_",i,".png", sep="")
  png(filename=file_name, width = 480, height = 480,
      units = "px")
  print(SEDEX_Graph_List_Febyprod[[i]])
  dev.off()
}


