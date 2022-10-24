##### Figure 3 #####
#### Load packages ####
library("here")
library("tidyverse")
library("extrafont")
library("RColorBrewer")
library("scales")
library("data.table")
library("Amelia")
library("survival")
library("survminer")
library("patchwork")
library("brms")
library("tidybayes")
library("modelr")
library("gghalves") 
library("broom.mixed") 

#### Load data ####
weekly_monitoring<-read.table("Fig3_weekly_monitoring.txt",sep="\t",stringsAsFactors = F, 
                              header=T, check.names=F,as.is=T)
GSIS<-read.table("Fig3_GSIS.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)
Mice<-read.table("Fig3_animalTable.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)
blood<-read.table("Fig3_blood_percentages.txt",sep="\t",stringsAsFactors = F, 
                  header=T, check.names=F,as.is=T)
blood_cell_counts<-read_delim("Fig3_blood_cell_counts.txt", name_repair = "universal")
randomFed<-read.table("Fig3_randomFed.txt",sep="\t",stringsAsFactors = F, 
                      header=T, check.names=F,as.is=T)

GroupPalette<-c("#89adc5", #light grey-blue
                "#0158e2", #bright blue
                "#5de522", #bright green
                "#002c65",#dark blue
                "#ff52f6",#pink
                "#ff913f" #orange
) 

GroupShapes<-c(23,25,24,5,21,22)

#### Wrangle data ####
###### Mouse metadata #####
Mice$Termination<-as.POSIXct(strptime(Mice$Termination,"%Y/%m/%d",tz="America/Vancouver"))
Mice<-Mice%>%mutate(Group=case_when(Group=="0.5e6 CAR-T cells"~"CARlo",
                                    Group=="1e6 CAR-T cells + 2e6 PBMCs"~"CARmed_PBMClo",
                                    Group=="1e6 CAR-T cells + 9e6 PBMCs"~"CARmed_PBMChi",
                                    Group=="1e6 CAR-T cells"~"CARmed",
                                    Group=="3e6 CAR-T cells"~"CARhi",
                                    Group=="A2+ Human Islets"~"PBS",
                                    TRUE~Group))

Mice$InjDays <- sapply(1:nrow(Mice),function(i){
  if(Mice$Group[i]=="CARhi"){
    injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))} else if(Mice$Group[i]!="A2- Human Islets"){
      injDate<-as.POSIXct(strptime("2019/08/30","%Y/%m/%d",tz="America/Vancouver"))} else if (Mice$Group[i]=="PBS"){
        injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))
      } else {
        injDate<-as.POSIXct(strptime("2019/10/11","%Y/%m/%d",tz="America/Vancouver"))
      }
  d<-Mice$Termination[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

###### Weekly blood glucose and body weight measurements ######
weekly_monitoring$Group<-sapply(weekly_monitoring$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})
weekly_monitoring<-filter(weekly_monitoring,Group!="A2- Human Islets")
#we had another group of mice that received A2- human islets, but they developed xGvHD and their human C-peptide levels dropped
#we planned to repeat the experiment but then the COVID-19 pandemic hit

#these four mice received 3 million CAR-T cells at 18 weeks post-islet transplant
weekly_monitoring$Group[weekly_monitoring$AnimalID%in%c("CRISPR02M.13","CRISPR02M.27","CRISPR02M.31","CRISPR02M.35")&
                          weekly_monitoring$Date>"2019/12/11"]<-"CARhi"
weekly_monitoring$Date<-as.POSIXct(strptime(weekly_monitoring$Date,"%Y/%m/%d",tz="America/Vancouver"))
weekly_monitoring$Days <- sapply(1:nrow(weekly_monitoring),function(i){
  if(weekly_monitoring$Group[i]!="A2- Human Islets"){
    surgDate<-as.POSIXct(strptime("2019/08/10","%Y/%m/%d",tz="America/Vancouver"))} else{
      surgDate<-as.POSIXct(strptime("2019/09/21","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-weekly_monitoring$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
weekly_monitoring$Weeks<-round(weekly_monitoring$Days/7)
weekly_monitoring$Group<-factor(weekly_monitoring$Group,levels = c("PBS",
                                                                   "CARlo","CARmed",
                                                                   "CARhi",
                                                                   "CARmed_PBMClo",
                                                                   "CARmed_PBMChi"))

weekly_monitoring<-weekly_monitoring[order(weekly_monitoring$Date,
                                           weekly_monitoring$Group,weekly_monitoring$AnimalID),]
weekly_monitoring$AnimalID <- factor(weekly_monitoring$AnimalID, levels = as.character(unique(weekly_monitoring$AnimalID)))
weekly_monitoring$InjDays <- sapply(1:nrow(weekly_monitoring),function(i){
  if(weekly_monitoring$Group[i]=="CARhi"){
    injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))} else if(weekly_monitoring$Group[i]!="A2- Human Islets"){
      injDate<-as.POSIXct(strptime("2019/08/30","%Y/%m/%d",tz="America/Vancouver"))} else if (weekly_monitoring$Group[i]=="PBS"){
        injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))
      } else {
        injDate<-as.POSIXct(strptime("2019/10/11","%Y/%m/%d",tz="America/Vancouver"))
      }
  d<-weekly_monitoring$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})
weekly_monitoring$InjWeeks<-round(weekly_monitoring$InjDays/7)

names(GroupPalette)<-levels(weekly_monitoring$Group)
names(GroupShapes)<-levels(weekly_monitoring$Group)

weekly_monitoring$Factor<-interaction(weekly_monitoring$Group,weekly_monitoring$Date)

#No missing values
weeklyTab<-weekly_monitoring[weekly_monitoring$State=="Random",]%>%
  group_by(Factor)%>%
  summarise(BG.UQ=quantile(BG,0.75,na.rm=T), BG.LQ=quantile(BG,0.25,na.rm=T),
            BG.Median=median(BG,na.rm=T),
            BW.UQ=quantile(BW,0.75,na.rm=T), BW.LQ=quantile(BW,0.25,na.rm=T),
            BW.Median=median(BW,na.rm=T))

weeklyTab$Group<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$Group[weekly_monitoring$Factor==f]))
weeklyTab$Group<-factor(weeklyTab$Group,levels=levels(weekly_monitoring$Group))
weeklyTab$Days<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$Days[weekly_monitoring$Factor==f]))
weeklyTab$Weeks<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$Weeks[weekly_monitoring$Factor==f]))
weeklyTab$InjDays<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$InjDays[weekly_monitoring$Factor==f]))
weeklyTab$InjWeeks<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$InjWeeks[weekly_monitoring$Factor==f]))

weekly_monitoring$AnimalID<-as.character(weekly_monitoring$AnimalID)
weekly_monitoring$AnimalID[weekly_monitoring$Group=="CARhi"]<-paste0(weekly_monitoring$AnimalID[weekly_monitoring$Group=="CARhi"],"a")

dups<-weekly_monitoring%>%filter(AnimalID%in%c("CRISPR02M.13","CRISPR02M.27","CRISPR02M.31","CRISPR02M.35")&
                                   Days==104)
dups$Group<-"CARhi"
dups$Group<-factor(dups$Group,levels = c("PBS",
                                         "CARlo","CARmed",
                                         "CARhi",
                                         "CARmed_PBMClo",
                                         "CARmed_PBMChi"))
dups$InjDays<-round(as.numeric(difftime("2019/11/25","2019/12/11",units=c("days"))))
dups$InjWeeks<-round(dups$InjDays/7)
dups$AnimalID<-paste0(dups$AnimalID,"a")

weekly_monitoring<-bind_rows(weekly_monitoring,dups)

###### Graft versus Host Disease tracking ######
weekly_gvhd<-weekly_monitoring
weekly_gvhd$AnimalID<-as.character(weekly_gvhd$AnimalID)

###### Glucose stimulated insulin secretion ######
GSIS$Date<-as.POSIXct(strptime(GSIS$Date,"%Y/%m/%d",tz="America/Vancouver"))
GSIS$Group<-sapply(GSIS$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})
GSIS<-filter(GSIS,Group!="A2- Human Islets")
GSIS$Group[GSIS$AnimalID%in%c("CRISPR02M.13","CRISPR02M.27","CRISPR02M.31","CRISPR02M.35")&
             GSIS$Date>"2019/12/11"]<-"CARhi"

GSIS$Group<-factor(GSIS$Group,levels = levels(weekly_monitoring$Group))

GSIS$Days <- sapply(1:nrow(GSIS),function(i){
  if(GSIS$Group[i]!="A2- Human Islets"){
    surgDate<-as.POSIXct(strptime("2019/08/10","%Y/%m/%d",tz="America/Vancouver"))} else {
      surgDate<-as.POSIXct(strptime("2019/09/21","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-GSIS$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
GSIS$Weeks<-floor(GSIS$Days/7)
GSIS<-GSIS[order(GSIS$Weeks,GSIS$Group,GSIS$AnimalID),]
GSIS$AnimalID <- factor(GSIS$AnimalID, levels = as.character(unique(GSIS$AnimalID)))
GSIS$Factor <- paste(GSIS$Group,GSIS$Delivery,GSIS$Weeks,"weeks",GSIS$Time,"min",sep=".")
GSIS$Delivery<-factor(GSIS$Delivery)

GSIS$InjDays <- sapply(1:nrow(GSIS),function(i){
  if(GSIS$Group[i]=="CARhi"){
    injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))} else if(GSIS$Group[i]!="A2- Human Islets"){
      injDate<-as.POSIXct(strptime("2019/08/30","%Y/%m/%d",tz="America/Vancouver"))} else if (GSIS$Group[i]=="PBS"){
        injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))
      } else {
        injDate<-as.POSIXct(strptime("2019/10/11","%Y/%m/%d",tz="America/Vancouver"))
      }
  d<-GSIS$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

GSIS$InjWeeks<-round(GSIS$InjDays/7)
GSIS$InjWeeks[GSIS$Date=="2019-11-25 PST"&GSIS$Group=="PBS"]<-13
GSIS$InjWeeks[GSIS$Date=="2019-11-25 PST"&GSIS$Group=="A2- Human Islets"]<-7
GSIS$InjWeeks[GSIS$Date=="2019-09-30 PST"&GSIS$Group=="A2- Human Islets"]<--1

#C-peptide
GSIS$CP<-ifelse(GSIS$CP.End==GSIS$CP.Start,GSIS$CP.End,NA)
GSIS$BG<-ifelse(GSIS$BG.End==GSIS$BG.Start,GSIS$BG.End,NA)

#use priors for human C-peptide below limit of detection (22.5 pg/mL)
pr <- matrix(c(0,16,log(0.001),log(22.5),0.99,
               0,17,33.3,147,0.9999999),byrow=T, nrow=2,ncol = 5)

set.seed(12345)
GSIS.imp <- amelia(GSIS, m = 10, p2s=0, idvars = c("Date","Days","Delivery","InjDays","Weeks",
                                                   "CP.End","CP.Start",
                                                   "BG.End","BG.Start",
                                                   "Factor"),
                   ts="Time",cs="AnimalID",intercs = T,
                   noms = c("Group"),polytime=1, #linear effect of time
                   logs="CP",
                   priors = pr,
                   multicore=8,
                   empri = .01*nrow(GSIS), max.resample = 5000)
summary(GSIS.imp)
plot(GSIS.imp,which.vars = 16:17)

tabList<-lapply(GSIS.imp$imputations,function(imp){
  tab<-imp%>%
    group_by(Factor)%>% 
    summarise(BG.UQ=quantile(BG,0.75,na.rm=T), BG.LQ=quantile(BG,0.25,na.rm=T),
              BG.Median=median(BG,na.rm=T),
              CP.UQ=quantile(CP,0.75,na.rm=T), CP.LQ=quantile(CP,0.25,na.rm=T),
              CP.Median=median(CP,na.rm=T))})

GSIStab<-do.call(rbind, tabList)

GSIStab<-GSIStab%>%
  group_by(Factor)%>% 
  summarise(BG.UQ=mean(BG.UQ,na.rm=T), BG.LQ=mean(BG.LQ,na.rm=T),
            BG.Median=mean(BG.Median,na.rm=T),
            CP.UQ=mean(CP.UQ,na.rm=T), CP.LQ=mean(CP.LQ,na.rm=T),
            CP.Median=mean(CP.Median,na.rm=T))

GSIStab$Group<-sapply(GSIStab$Factor,function(f) unique(GSIS$Group[GSIS$Factor==f]))
GSIStab$Group<-factor(GSIStab$Group,levels=levels(GSIS$Group))
GSIStab$Time<-sapply(GSIStab$Factor,function(f) unique(GSIS$Time[GSIS$Factor==f]))
GSIStab$Weeks<-sapply(GSIStab$Factor,function(f) unique(GSIS$Weeks[GSIS$Factor==f]))
GSIStab$InjWeeks<-sapply(GSIStab$Factor,function(f) unique(GSIS$InjWeeks[GSIS$Factor==f]))

dups<-GSIS%>%filter(AnimalID%in%c("CRISPR02M.13","CRISPR02M.27","CRISPR02M.31","CRISPR02M.35")&
                      InjWeeks==12)
dups$Group<-"CARhi"
dups$Group<-factor(dups$Group,levels = c("PBS",
                                         "CARlo","CARmed",
                                         "CARhi",
                                         "CARmed_PBMClo",
                                         "CARmed_PBMChi"))
dups$InjDays<-round(as.numeric(difftime("2019/11/25","2019/12/11",units=c("days"))))
dups$InjWeeks<-round(dups$InjDays/7)
dups$Factor <- paste(dups$Group,dups$Weeks,"weeks",dups$Time,"min",sep=".")

GSIS<-rbind(GSIS,dups)

dups<-GSIStab[GSIStab$Group=="PBS"&GSIStab$InjWeeks==12,]
dups$Group<-"CARhi"
dups$Factor <- paste(dups$Group,dups$Weeks,"weeks",dups$Time,"min",sep=".")
dups$InjWeeks<--2

GSIStab<-rbind(GSIStab, dups)

###### Random fed human C-peptide ######
randomFed$Group<-sapply(randomFed$AnimalID,function(a){Mice$Group[Mice$AnimalID==a]})
randomFed<-filter(randomFed,Group!="A2- Human Islets")

randomFed$Group <- factor(randomFed$Group, levels = levels(weekly_monitoring$Group))
randomFed<-randomFed[order(randomFed$Group,randomFed$AnimalID),]

randomFed$TransplantDays <- sapply(1:nrow(randomFed),function(i){
  if(randomFed$Group[i]!="A2- Human Islets"){
    surgDate<-as.POSIXct(strptime("2019/08/10","%Y/%m/%d",tz="America/Vancouver"))} else{
      surgDate<-as.POSIXct(strptime("2019/09/21","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-randomFed$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})

randomFed$TransplantWeeks<-ceiling(randomFed$TransplantDays/7)

randomFed$InjDays <- sapply(1:nrow(randomFed),function(i){
  injDate<-case_when(grepl("A2-",randomFed$Group[i])~as.POSIXct(strptime("2019/10/11","%Y/%m/%d",tz="America/Vancouver")),
                     TRUE~as.POSIXct(strptime("2019/08/30","%Y/%m/%d",tz="America/Vancouver")))
  d<-randomFed$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})
randomFed$InjWeeks<-ceiling(randomFed$InjDays/7)

randomFed<-randomFed[order(randomFed$InjDays,randomFed$Group,randomFed$AnimalID),]
# randomFed$AnimalID <- factor(randomFed$AnimalID, levels = levels(weekly_monitoring$AnimalID))
randomFed$Factor<-paste(randomFed$Group,randomFed$InjDays,"days",sep=".")
randomFed$CP<-ifelse(randomFed$CP.End==randomFed$CP.Start,randomFed$CP.End,NA)

set.seed(12345)
randomFed.imp <- amelia(randomFed, m = 10, p2s=0, idvars = c("Date","TransplantDays","TransplantWeeks","InjDays",
                                                             "Factor","InjWeeks","AnimalID"),
                        noms = c("Group"),
                        empri = .05*nrow(randomFed), max.resample = 5000) #

plot(randomFed.imp,which.vars = 11)

tabList<-lapply(randomFed.imp$imputations,function(imp){
  tab<-imp%>%
    group_by(Factor)%>% 
    summarise(CP.UQ=quantile(CP,0.75,na.rm=T), CP.LQ=quantile(CP,0.25,na.rm=T),
              CP.Median=median(CP,na.rm=T))})

randomFedtab<-do.call(rbind, tabList)

randomFedtab<-randomFedtab%>%
  group_by(Factor)%>% 
  summarise(CP.UQ=mean(CP.UQ,na.rm=T), CP.LQ=mean(CP.LQ,na.rm=T),
            CP.Median=mean(CP.Median,na.rm=T))

randomFedtab$Group<-sapply(randomFedtab$Factor,function(f) unique(randomFed$Group[randomFed$Factor==f]))
randomFedtab$Group<-factor(randomFedtab$Group,levels=levels(randomFed$Group))
randomFedtab$TransplantWeeks<-sapply(randomFedtab$Factor,function(f) unique(randomFed$TransplantWeeks[randomFed$Factor==f]))
randomFedtab$InjWeeks<-sapply(randomFedtab$Factor,function(f) unique(randomFed$InjWeeks[randomFed$Factor==f]))

###### Percentages of immune cells in blood ######
blood$Group<-sapply(blood$AnimalID,function(a){Mice$Group[Mice$AnimalID==a]})
blood<-blood%>%filter(Group!="A2- Human Islets")
blood$Group <- factor(blood$Group, levels = levels(weekly_monitoring$Group))
blood<-blood[order(blood$Group,blood$AnimalID),]

blood$TransplantDays <- sapply(1:nrow(blood),function(i){
  if(blood$Group[i]!="A2- Human Islets"){
    surgDate<-as.POSIXct(strptime("2019/08/10","%Y/%m/%d",tz="America/Vancouver"))} else{
      surgDate<-as.POSIXct(strptime("2019/09/21","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-blood$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})

blood$InjDays <- sapply(1:nrow(blood),function(i){
  if(blood$Group[i]=="CARhi"){
    injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))} else if(blood$Group[i]!="A2- Human Islets"){
      injDate<-as.POSIXct(strptime("2019/08/30","%Y/%m/%d",tz="America/Vancouver"))} else if (blood$Group[i]=="PBS"){
        injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))
      } else {
        injDate<-as.POSIXct(strptime("2019/10/11","%Y/%m/%d",tz="America/Vancouver"))
      }
  d<-blood$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})
blood$InjWeeks<-round(blood$InjDays/7)

blood$CP<-vapply(1:nrow(blood),function(i){
  aID<-as.character(blood$AnimalID[i])
  d<-blood$TransplantDays[i]
  w=floor(d/7)
  CP<-ifelse(is_empty(GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Weeks==w&GSIS$Time==15]),
             NA,
             GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Weeks==w&GSIS$Time==15]/1000)
},double(1))

blood<-blood[order(blood$InjDays,blood$Group,blood$AnimalID),]
blood$Factor<-paste(blood$Group,blood$InjDays,"days",sep=".")

blood<-blood%>%filter(rowSums(blood[,3:10],na.rm=T)>0)

bloodTab<-blood%>%
  group_by(Factor)%>%
  summarise(percent_CD4posMycposCD45pos.UQ=quantile(percent_CD4posMycposCD45pos,0.75,na.rm=T), 
            percent_CD4posMycposCD45pos.LQ=quantile(percent_CD4posMycposCD45pos,0.25,na.rm=T),
            percent_CD4posMycposCD45pos.Median=median(percent_CD4posMycposCD45pos,na.rm=T),
            
            percent_CD8posMycposCD45pos.UQ=quantile(percent_CD8posMycposCD45pos,0.75,na.rm=T), 
            percent_CD8posMycposCD45pos.LQ=quantile(percent_CD8posMycposCD45pos,0.25,na.rm=T),
            percent_CD8posMycposCD45pos.Median=median(percent_CD8posMycposCD45pos,na.rm=T),
            
            percent_CD4pos_of_CD45posMycpos.UQ=quantile(percent_CD4pos_of_CD45posMycpos,0.75,na.rm=T), 
            percent_CD4pos_of_CD45posMycpos.LQ=quantile(percent_CD4pos_of_CD45posMycpos,0.25,na.rm=T),
            percent_CD4pos_of_CD45posMycpos.Median=median(percent_CD4pos_of_CD45posMycpos,na.rm=T),
            
            percent_CD8pos_of_CD45posMycpos.UQ=quantile(percent_CD8pos_of_CD45posMycpos,0.75,na.rm=T), 
            percent_CD8pos_of_CD45posMycpos.LQ=quantile(percent_CD8pos_of_CD45posMycpos,0.25,na.rm=T),
            percent_CD8pos_of_CD45posMycpos.Median=median(percent_CD8pos_of_CD45posMycpos,na.rm=T),
            
            percent_CD45pos.UQ=quantile(percent_CD45pos,0.75,na.rm=T), 
            percent_CD45pos.LQ=quantile(percent_CD45pos,0.25,na.rm=T),
            percent_CD45pos.Median=median(percent_CD45pos,na.rm=T),
            
            percent_Mycpos_hCD45pos.UQ=quantile(percent_Mycpos_hCD45pos,0.75,na.rm=T), 
            percent_Mycpos_hCD45pos.LQ=quantile(percent_Mycpos_hCD45pos,0.25,na.rm=T),
            percent_Mycpos_hCD45pos.Median=median(percent_Mycpos_hCD45pos,na.rm=T),
            
            percent_CD4posCD45pos.UQ=quantile(percent_CD4posCD45pos,0.75,na.rm=T), 
            percent_CD4posCD45pos.LQ=quantile(percent_CD4posCD45pos,0.25,na.rm=T),
            percent_CD4posCD45pos.Median=median(percent_CD4posCD45pos,na.rm=T),
            
            percent_CD8posCD45pos.UQ=quantile(percent_CD8posCD45pos,0.75,na.rm=T), 
            percent_CD8posCD45pos.LQ=quantile(percent_CD8posCD45pos,0.25,na.rm=T),
            percent_CD8posCD45pos.Median=median(percent_CD8posCD45pos,na.rm=T))

bloodTab$Group<-sapply(bloodTab$Factor,function(f) unique(blood$Group[blood$Factor==f]))
bloodTab$Group<-factor(bloodTab$Group,levels=levels(blood$Group))
bloodTab$InjDays<-sapply(bloodTab$Factor,function(f) unique(blood$InjDays[blood$Factor==f]))
bloodTab$TransplantDays<-sapply(bloodTab$Factor,function(f) unique(blood$TransplantDays[blood$Factor==f]))

###### Counts of immune cells in blood ######
blood_cell_counts$Group<-sapply(blood_cell_counts$AnimalID,function(a){Mice$Group[Mice$AnimalID==a]})
blood_cell_counts<-blood_cell_counts%>%filter(Group!="A2- Human Islets")
blood_cell_counts$Group <- factor(blood_cell_counts$Group, levels = levels(weekly_monitoring$Group)) #,"S5","S6"
blood_cell_counts<-blood_cell_counts[order(blood_cell_counts$Group,blood_cell_counts$AnimalID),]

blood_cell_counts$TransplantDays <- sapply(1:nrow(blood_cell_counts),function(i){
  if(blood_cell_counts$Group[i]!="A2- Human Islets"){
    surgDate<-as.POSIXct(strptime("2019/08/10","%Y/%m/%d",tz="America/Vancouver"))} else{
      surgDate<-as.POSIXct(strptime("2019/09/21","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-blood_cell_counts$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})

blood_cell_counts$InjDays <- sapply(1:nrow(blood_cell_counts),function(i){
  if(blood_cell_counts$Group[i]=="CARhi"){
    injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))} else if(blood_cell_counts$Group[i]!="A2- Human Islets"){
      injDate<-as.POSIXct(strptime("2019/08/30","%Y/%m/%d",tz="America/Vancouver"))} else if (blood_cell_counts$Group[i]=="PBS"){
        injDate<-as.POSIXct(strptime("2019/12/11","%Y/%m/%d",tz="America/Vancouver"))
      } else {
        injDate<-as.POSIXct(strptime("2019/10/11","%Y/%m/%d",tz="America/Vancouver"))
      }
  d<-blood_cell_counts$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})
blood_cell_counts$InjWeeks<-round(blood_cell_counts$InjDays/7)

blood_cell_counts$CP<-vapply(1:nrow(blood_cell_counts),function(i){
  aID<-as.character(blood_cell_counts$AnimalID[i])
  d<-blood_cell_counts$TransplantDays[i]
  w=floor(d/7)
  CP<-ifelse(is_empty(GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Weeks==w&GSIS$Time==15]),
             NA,
             GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Weeks==w&GSIS$Time==15]/1000)
},double(1))

blood_cell_counts<-blood_cell_counts[order(blood_cell_counts$InjDays,
                                           blood_cell_counts$Group,blood_cell_counts$AnimalID),]
blood_cell_counts$Factor<-paste(blood_cell_counts$Group,blood_cell_counts$InjDays,"days",sep=".")

blood_cell_counts<-blood_cell_counts%>%filter(rowSums(blood_cell_counts[,4:11],na.rm=T)>0)

bloodcellTab<-blood_cell_counts%>%
  group_by(Factor)%>%
  summarise(normal.count_CD4pos_of_CD45posMycpos.UQ=quantile(normal.count_CD4pos.of.CD45posMycpos,0.75,na.rm=T), 
            normal.count_CD4pos_of_CD45posMycpos.LQ=quantile(normal.count_CD4pos.of.CD45posMycpos,0.25,na.rm=T),
            normal.count_CD4pos_of_CD45posMycpos.Median=median(normal.count_CD4pos.of.CD45posMycpos,na.rm=T),
            normal.count_CD8pos_of_CD45posMycpos.UQ=quantile(normal.count_CD8pos.of.CD45posMycpos,0.75,na.rm=T), 
            normal.count_CD8pos_of_CD45posMycpos.LQ=quantile(normal.count_CD8pos.of.CD45posMycpos,0.25,na.rm=T),
            normal.count_CD8pos_of_CD45posMycpos.Median=median(normal.count_CD8pos.of.CD45posMycpos,na.rm=T),
            normal.count_CD45pos.UQ=quantile(normal.count_CD45pos,0.75,na.rm=T), 
            normal.count_CD45pos.LQ=quantile(normal.count_CD45pos,0.25,na.rm=T),
            normal.count_CD45pos.Median=median(normal.count_CD45pos,na.rm=T),
            normal.count_Mycpos_hCD45pos.UQ=quantile(normal.count_Mycpos_hCD45pos,0.75,na.rm=T), 
            normal.count_Mycpos_hCD45pos.LQ=quantile(normal.count_Mycpos_hCD45pos,0.25,na.rm=T),
            normal.count_Mycpos_hCD45pos.Median=median(normal.count_Mycpos_hCD45pos,na.rm=T))

bloodcellTab$Group<-sapply(bloodcellTab$Factor,function(f) unique(blood_cell_counts$Group[blood$Factor==f]))
bloodcellTab$Group<-factor(bloodcellTab$Group,levels=levels(blood_cell_counts$Group))
bloodcellTab$InjDays<-sapply(bloodcellTab$Factor,function(f) unique(blood_cell_counts$InjDays[blood$Factor==f]))
bloodcellTab$TransplantDays<-sapply(bloodcellTab$Factor,function(f) unique(blood_cell_counts$TransplantDays[blood$Factor==f]))

##### Statistical Analysis #####
#Note that no inferential statistical analyses were performed on Figure 3B, 3C or 3D -
#The body weight and blood glucose values are clearly not different, and
#Figure 3D is a subset of the same data as Figure 3E, so we did not analyze the same data
#in two separate ways

#Figure 3E
GSIS$AnimalID<-as.character(GSIS$AnimalID)
GSIS$AnimalID[GSIS$Group=="CARhi"]<-paste(GSIS$AnimalID[GSIS$Group=="CARhi"],"a",sep = "")

Total_CP<-GSIS%>%as_tibble()%>%
  group_by(AnimalID,InjDays)%>%summarise(TotalCP=sum(CP.End))

Total_CP$Rejected<-case_when(Total_CP$TotalCP<=(22.5*3/0.8)~1, #Total of limit of detection at 3 time points
                             TRUE~0)
Total_CP$Group<-sapply(Total_CP$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})

RF_CP<-randomFed%>%as_tibble()%>%
  group_by(AnimalID,InjDays)%>%summarise(TotalCP=sum(CP.End))

RF_CP$Rejected<-case_when(RF_CP$TotalCP<=(22.5/0.8)~1, #Total of limit of detection at 1 time point
                          TRUE~0)
RF_CP$Group<-sapply(RF_CP$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})

Total_CP2<-rbind(Total_CP,RF_CP)

sumsSurv<-tibble("AnimalID"=unique(Total_CP2$AnimalID),
                 "Time"=sapply(unique(Total_CP2$AnimalID),function(id) {
                   ifelse(sum(Total_CP2$Rejected[Total_CP2$AnimalID==id])>0,
                          min(Total_CP2$InjDays[Total_CP2$AnimalID==id&
                                                  Total_CP2$Rejected==1])/7, 
                          Mice$InjDays[Mice$AnimalID==id]/7)}),
                 "Rejected"=sapply(unique(Total_CP2$AnimalID),function(id) {
                   case_when(sum(Total_CP2$Rejected[Total_CP2$AnimalID==id])>0~1,
                             sum(Total_CP2$Rejected[Total_CP2$AnimalID==id])==0~0)}),
                 "Group"=unique(Total_CP2$AnimalID) %>% 
                   map_chr(function(id) unique(Total_CP2$Group[Total_CP2$AnimalID==id]))
)

pairwise_survdiff(Surv(time=Time, event=Rejected,type="right") ~ Group,
                  data = sumsSurv, p.adjust.method = "BH",
                  na.action=options()$na.action,
                  rho = 0)

#Figure 3F
gvhdSurv<-tibble("AnimalID"=unique(weekly_gvhd$AnimalID),
                 "Time"=sapply(unique(weekly_gvhd$AnimalID),function(id) {
                   ifelse(sum(weekly_gvhd$GvHD[weekly_gvhd$AnimalID==id])>0,
                          min(weekly_gvhd$InjDays[weekly_gvhd$AnimalID==id&
                                                    weekly_gvhd$GvHD>1])/7, 
                          max(Mice$InjDays[Mice$AnimalID==id]/7))}),
                 "GvHD"=sapply(unique(weekly_gvhd$AnimalID),function(id) {
                   ifelse(sum(weekly_gvhd$GvHD[weekly_gvhd$AnimalID==id])>0,
                          1, 0)}),
                 "Group"=unique(weekly_gvhd$AnimalID) %>% 
                   map_chr(function(id) as.character(unique(weekly_gvhd$Group[weekly_gvhd$AnimalID==id])))
)

pairwise_survdiff(Surv(time=Time, event=GvHD,type="right") ~ Group,
                  data = gvhdSurv, p.adjust.method = "BH",
                  na.action=options()$na.action,
                  rho = 0)

#Figure 3G
#other models were tested and assessed for fit using pp_check() and loo_compare()
CD45mod_slope <- brm(percent_CD45pos~InjStage*Group+(InjStage|AnimalID),
                     data=blood%>%
                       mutate(InjStage=factor(InjWeeks),
                              percent_CD45pos=percent_CD45pos+0.001)%>%
                       filter(Group!="PBS",InjDays>0,!is.na(percent_CD45pos)),
                     prior=c(set_prior("normal(-5,5)", class = "Intercept"),
                             set_prior("normal(-5,5)", class = "b")),
                     iter=10000,
                     warmup=1000,
                     family="lognormal",
                     control = list(adapt_delta=0.975,
                                    max_treedepth=12),
                     backend = "cmdstanr",
                     threads = threading(parallel::detectCores()),
                     silent=0,
                     save_pars = save_pars(all = TRUE))

pp_check(CD45mod_slope)+scale_x_log10()

cd45_res_tab<-rbind(hypothesis(CD45mod_slope,"exp(Intercept+InjStage2+GroupCARmed_PBMChi+InjStage2:GroupCARmed_PBMChi)-
                               exp(Intercept+InjStage2+GroupCARmed+InjStage2:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3+GroupCARmed+InjStage3:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage2+GroupCARmed_PBMChi+InjStage2:GroupCARmed_PBMChi)-
                               exp(Intercept+InjStage2+GroupCARhi+InjStage2:GroupCARhi)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3+GroupCARhi+InjStage3:GroupCARhi)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage2+GroupCARmed_PBMChi+InjStage2:GroupCARmed_PBMChi)-
                               exp(Intercept+InjStage2+GroupCARmed_PBMClo+InjStage2:GroupCARmed_PBMClo)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                               exp(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage2+GroupCARmed_PBMChi+InjStage2:GroupCARmed_PBMChi)-
                               exp(Intercept+InjStage2)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARhi+InjStage3:GroupCARhi)-
                               exp(Intercept+InjStage3+GroupCARmed+InjStage3:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage4+GroupCARhi+InjStage4:GroupCARhi)-
                               exp(Intercept+InjStage4+GroupCARmed+InjStage4:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage6+GroupCARhi+InjStage6:GroupCARhi)-
                               exp(Intercept+InjStage6+GroupCARmed+InjStage6:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage8+GroupCARhi+InjStage8:GroupCARhi)-
                               exp(Intercept+InjStage9+GroupCARmed+InjStage9:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage13+GroupCARhi+InjStage13:GroupCARhi)-
                               exp(Intercept+InjStage13+GroupCARmed+InjStage13:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARhi+InjStage3:GroupCARhi)-
                               exp(Intercept+InjStage3)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage4+GroupCARhi+InjStage4:GroupCARhi)-
                               exp(Intercept+InjStage4)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage6+GroupCARhi+InjStage6:GroupCARhi)-
                               exp(Intercept+InjStage6)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage8+GroupCARhi+InjStage8:GroupCARhi)-
                               exp(Intercept+InjStage9)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage12+GroupCARhi+InjStage12:GroupCARhi)-
                               exp(Intercept+InjStage12)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage13+GroupCARhi+InjStage13:GroupCARhi)-
                               exp(Intercept+InjStage13)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage8+GroupCARhi+InjStage8:GroupCARhi)-
                               exp(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)<0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage3+GroupCARmed+InjStage3:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage4+GroupCARmed+InjStage4:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage5+GroupCARmed_PBMClo+InjStage5:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage5+GroupCARmed+InjStage5:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage6+GroupCARmed+InjStage6:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage7+GroupCARmed_PBMClo+InjStage7:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage7+GroupCARmed+InjStage7:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage9+GroupCARmed+InjStage9:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage12+GroupCARmed_PBMClo+InjStage12:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage12+GroupCARmed+InjStage12:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage13+GroupCARmed+InjStage13:GroupCARmed)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage3)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage4)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage5+GroupCARmed_PBMClo+InjStage5:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage5)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage6)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage7+GroupCARmed_PBMClo+InjStage7:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage7)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage9)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage12+GroupCARmed_PBMClo+InjStage12:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage12)>0")$hypothesis,
                    hypothesis(CD45mod_slope,"exp(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage13)>0")$hypothesis)

#Figure 3H
#other models were tested and assessed for fit using pp_check() and loo_compare()
MycCD45mod_noRE <- brm(percent_Mycpos_hCD45pos~InjStage*Group,
                       data=blood%>%
                         mutate(InjStage=factor(InjWeeks),
                                percent_Mycpos_hCD45pos=percent_Mycpos_hCD45pos)%>%
                         filter(Group%in%c("CARhi","CARmed_PBMChi",
                                           "CARmed_PBMClo"),
                                InjDays>6,!is.na(percent_Mycpos_hCD45pos)),
                       prior=c(set_prior("normal(100,15)", class = "Intercept"),
                               set_prior("normal(0,50)", class = "b")),
                       iter=10000,
                       warmup=1500,
                       family="skew_normal",
                       control = list(adapt_delta=0.95),
                       backend = "cmdstanr",
                       threads = threading(parallel::detectCores()),
                       silent=0,
                       save_pars = save_pars(all = TRUE))

pp_check(MycCD45mod_noRE)

MycCD45_res_tab<-rbind(hypothesis(MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept)<0")$hypothesis,
                       hypothesis(MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept+GroupCARmed_PBMClo)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMClo)-
                               (Intercept)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                      (Intercept+InjStage3)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                      (Intercept+InjStage4)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                      (Intercept+InjStage6)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                      (Intercept+InjStage8)<0")$hypothesis,
                      hypothesis(MycCD45mod_noRE,"(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)-
                      (Intercept+InjStage13)<0")$hypothesis)
#Figure 3I
CD4MycCD45mod_noRE <- brm(percent_CD4pos_of_CD45posMycpos~InjStage*Group,
                          data=blood%>%
                            mutate(InjStage=factor(InjWeeks),
                                   percent_CD4pos_of_CD45posMycpos=percent_CD4pos_of_CD45posMycpos)%>%
                            filter(Group%in%c("CARhi","CARmed_PBMChi",
                                              "CARmed_PBMClo"),
                                   InjDays>6,!is.na(percent_CD4pos_of_CD45posMycpos)),
                          prior=c(set_prior("normal(50,25)", class = "Intercept"),
                                  set_prior("normal(0,50)", class = "b")),
                          iter=10000,
                          warmup=1500,
                          family="skew_normal",
                          control = list(adapt_delta=0.95),
                          backend = "cmdstanr",
                          threads = threading(parallel::detectCores()),
                          silent=0,
                          save_pars = save_pars(all = TRUE))

CD4MycCD45_res_tab<-rbind(hypothesis(CD4MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept)<0")$hypothesis,
                          hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3)<0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept+GroupCARmed_PBMClo)>0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)>0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMClo)-
                               (Intercept)<0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                      (Intercept+InjStage3)<0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                      (Intercept+InjStage4)<0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                      (Intercept+InjStage6)<0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                      (Intercept+InjStage8)<0")$hypothesis,
                      hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)-
                      (Intercept+InjStage13)<0")$hypothesis)

#we also want to know if the percentages are different from the injected percentages, ~50%
CD4MycCD45_res_tab2<-rbind(hypothesis(CD4MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage5+GroupCARmed_PBMClo+InjStage5:GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage7+GroupCARmed_PBMClo+InjStage7:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage12+GroupCARmed_PBMClo+InjStage12:GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)<50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage3)>50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage4)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage6)>50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage8)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage11)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage13)=50")$hypothesis,
                           hypothesis(CD4MycCD45mod_noRE,"(Intercept+InjStage14)=50")$hypothesis
)

#Figure 3J
CD8MycCD45mod_noRE <- brm(percent_CD8pos_of_CD45posMycpos~InjStage*Group,
                          data=blood%>%
                            mutate(InjStage=factor(InjWeeks),
                                   percent_CD8pos_of_CD45posMycpos=percent_CD8pos_of_CD45posMycpos)%>%
                            filter(Group%in%c("CARhi","CARmed_PBMChi",
                                              "CARmed_PBMClo"),
                                   InjDays>6,!is.na(percent_CD8pos_of_CD45posMycpos)),
                          prior=c(set_prior("normal(50,25)", class = "Intercept"),
                                  set_prior("normal(0,50)", class = "b")),
                          iter=10000,
                          warmup=1500,
                          family="skew_normal",
                          control = list(adapt_delta=0.95),
                          backend = "cmdstanr",
                          threads = threading(parallel::detectCores()),
                          silent=0,
                          save_pars = save_pars(all = TRUE))

pp_check(CD8MycCD45mod_noRE)

CD8MycCD45_res_tab<-rbind(hypothesis(CD8MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept)>0")$hypothesis,
                          hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3)>0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept+GroupCARmed_PBMClo)<0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)<0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMClo)-
                               (Intercept)>0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                      (Intercept+InjStage3)>0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                      (Intercept+InjStage4)>0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                      (Intercept+InjStage6)>0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                      (Intercept+InjStage8)>0")$hypothesis,
                      hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)-
                      (Intercept+InjStage13)>0")$hypothesis)

CD8MycCD45_res_tab2<-rbind(hypothesis(CD8MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMChi)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+GroupCARmed_PBMClo)>50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage5+GroupCARmed_PBMClo+InjStage5:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)>50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage7+GroupCARmed_PBMClo+InjStage7:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)>50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage12+GroupCARmed_PBMClo+InjStage12:GroupCARmed_PBMClo)>50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage13+GroupCARmed_PBMClo+InjStage13:GroupCARmed_PBMClo)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage3)<50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage4)<50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage6)<50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage8)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage11)=50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage13)<50")$hypothesis,
                           hypothesis(CD8MycCD45mod_noRE,"(Intercept+InjStage14)=50")$hypothesis
)

#### Figures ####
# Figure 3A was made on Biorender

###### Figure 3B #####
Fig3B<-ggplot(weekly_monitoring[!is.na(weekly_monitoring$BW)&weekly_monitoring$State=="Random",],
               aes(x=Days/7,y=BW,group=Group))+
  #geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=2,alpha=0.7)+
  #facet_grid(Species~.,scales="free_x",space="free_x")+
  geom_ribbon(aes(x=Days/7,y=BW.Median,ymin = BW.LQ,ymax = BW.UQ, fill = Group),
              data=weeklyTab[!is.na(weeklyTab$BW.Median),],alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+ 
  geom_line(aes(y=BW.Median,colour = Group,group=Group), 
            data=weeklyTab[!is.na(weeklyTab$BW.Median),],size=0.5,alpha=0.8)+
  geom_point(aes(y=BW.Median,colour = Group,group=Group,shape=Group,fill=Group), #
             data=weeklyTab[!is.na(weeklyTab$BW.Median),],size=1,alpha=0.8)+
  geom_vline(xintercept=as.numeric(difftime("2019/08/30","2019/08/10",units="weeks")),
             colour=gray(0.20),linetype=2)+
  # geom_label(x=as.numeric(difftime("2019/08/30","2019/08/10",units="weeks")),
  #            y=38,label="Injection",family="Arial",fill="white",
  #            size=4.90,hjust=0.5,vjust=-0.2,colour="black")+
  geom_vline(xintercept=as.numeric(difftime("2019/12/11","2019/08/10",units="weeks")),
             colour=gray(0.20),linetype=2)+
  # geom_label(x=as.numeric(difftime("2019/12/11","2019/08/10",units="weeks")),
  #            y=38,label="Injection",family="Arial",fill="white",
  #            size=4.90,hjust=0.5,vjust=-0.2,colour="black")+
  labs(x="Weeks Post-Islet Transplant",y="Body\nWeight (g)")+
  #scale_x_datetime(labels = date_format("%m-%d-%Y"),breaks=unique(weekly_monitoring$Date))+
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_color_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  scale_y_continuous(breaks=pretty_breaks(n=5))+
  scale_x_continuous(breaks=pretty_breaks(n=6))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = c(0.235,0.235),
        legend.key.size = unit(0.35,"cm"))+
  theme(strip.text = element_text(family="Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3B

###### Figure 3C ######
Fig3C<-ggplot(weekly_monitoring[!is.na(weekly_monitoring$BG)&weekly_monitoring$State=="Random",],
               aes(x=Days/7,y=BG,group=Group))+
  geom_ribbon(aes(y=BG.Median,ymin = BG.LQ,ymax = BG.UQ,fill=Group), 
              data=weeklyTab[!is.na(weeklyTab$BG.Median),],alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  geom_line(aes(y=BG.Median,colour = Group,group=Group), #
            data=weeklyTab[!is.na(weeklyTab$BG.Median),],size=0.5,alpha=0.8)+
  geom_point(aes(y=BG.Median,colour = Group,group=Group,shape=Group,fill=Group), 
             data=weeklyTab[!is.na(weeklyTab$BG.Median),],size=1,alpha=0.8)+
  scale_x_continuous(breaks=pretty_breaks(n=6))+
  geom_vline(xintercept=as.numeric(difftime("2019/08/30","2019/08/10",units="weeks")),
             colour=gray(0.20),linetype=2)+
  geom_vline(xintercept=as.numeric(difftime("2019/12/11","2019/08/10",units="weeks")),
             colour=gray(0.20),linetype=2)+
  labs(x="Weeks Post-Islet Transplant",y="Blood\nGlucose (mM)")+
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_color_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = c(0.235,0.235),
        legend.key.size = unit(0.35,"cm"))+
  theme(strip.text = element_text(family="Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3C

###### Figure 3D ######
breaks_fun <- function(x) {
  if (max(x) > 3) {
    c(-1,1,3,5,7,9,11)
  } else {
    c(-1,1,3)
  }
}

Fig3D<-ggplot(GSIS[GSIS$Time%in%c(15),],
                     aes(x=InjWeeks,y=CP.End/3020.29,colour=Group,fill=Group))+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  facet_grid(.~Group,scales = "free_x",space="free_x",drop=T)+
  geom_line(aes(x=InjWeeks,y=CP.Median/3020.29,group=Group),
            data=GSIStab[GSIStab$Time%in%c(15),],
            size=0.5,alpha=1) +
  geom_point(aes(x=InjWeeks,y=CP.Median/3020.29,group=Group,shape=Group),
             data=GSIStab[GSIStab$Time%in%c(15),],
             size=1,alpha=0.8) +
  geom_linerange(aes(ymin=CP.Start/3020.29,ymax=CP.End/3020.29,group=AnimalID),
                 linetype=3,size=0.3,alpha=0.8)+
  geom_ribbon(aes(x=InjWeeks,y=CP.Median/1000,ymin = CP.LQ/3020.29, ymax = CP.UQ/3020.29,
                  group=Group),colour=NA,data=GSIStab[GSIStab$Time%in%c(15),],
              size=0.5,alpha=0.2)+
  scale_fill_manual(name="Group", values=GroupPalette,
                    labels=c("PBS",
                             "CAR-lo",
                             "CAR-med",
                             "CAR-hi",
                             "CAR-med PBMC-lo",
                             "CAR-med PBMC-hi"))+
  scale_color_manual(name="Group", values=GroupPalette,
                     labels=c("PBS",
                              "CAR-lo",
                              "CAR-med",
                              "CAR-hi",
                              "CAR-med PBMC-lo",
                              "CAR-med PBMC-hi"))+
  scale_shape_manual(name="Group", values=GroupShapes,
                     labels=c("PBS",
                              "CAR-lo",
                              "CAR-med",
                              "CAR-hi",
                              "CAR-med PBMC-lo",
                              "CAR-med PBMC-hi"))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y="Stimulated\nhuC-peptide (nM)")+
  geom_hline(yintercept = 4.5*5/3020.29,linetype=3,alpha=0.7)+
  scale_x_continuous(limits=c(-3,NA),breaks = breaks_fun)+
  scale_y_continuous(breaks=pretty_breaks(3))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "bottom")+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3D

###### Figure 3E ######
survLine<-c("41","3212","11","1","4121","6121")
names(survLine)<-paste0("Group=",c("PBS",
                                   "CARlo","CARmed",
                                   "CARhi",
                                   "CARmed_PBMClo",
                                   "CARmed_PBMChi"))

SurvPal<-GroupPalette
names(SurvPal)<-paste0("Group=",c("PBS",
                                  "CARlo","CARmed",
                                  "CARhi",
                                  "CARmed_PBMClo",
                                  "CARmed_PBMChi"))

fit <- survfit(Surv(time=Time, event=Rejected,type="right") ~ Group,
               data = sumsSurv)

Fig3E<-ggsurvplot(fit, data = sumsSurv, pval = F,
                       break.time.by = 2,
                       risk.table = TRUE,
                       risk.table.col = "strata",
                       risk.table.height = 0.5,
                       risk.table.y.text=T,
                       palette = SurvPal,
                       linetype=c("Group"),
                       legend="none",
                       ylab="Fraction\nw/o\nRejection",
                       xlab="Weeks Post-\nA2-CAR T Cell Injection")
risk.table.ylab <- ""
Fig3E$plot<-Fig3E$plot  + scale_x_continuous(breaks = seq(-2,14,2))+
  scale_linetype_discrete(breaks=survLine)
Fig3E

###### Figure 3F #####
gvhd_fit <- survfit(Surv(Time, GvHD,type="right") ~ Group,
                    data = gvhdSurv)

Fig3F<-ggsurvplot(gvhd_fit, data = gvhdSurv, pval = F,
                     break.time.by = 1,
                     xlim=c(1,14),
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     risk.table.height = 0.5,
                     risk.table.y.text=T,
                     palette = SurvPal,
                     linetype=c("strata"),
                     legend="none",
                     ylab="Fraction\nw/o GvHD",
                     xlab="Weeks Post-\nA2-CAR T Cell Injection")
risk.table.ylab <- ""
Fig3F$plot<-Fig3F$plot + scale_x_continuous(breaks = seq(-2,14,2))+
  scale_linetype_discrete(breaks=survLine)
Fig3F$table <- Fig3F$table + labs(y = risk.table.ylab)
Fig3F

###### Figure 3G ######
my_labeller = as_labeller(
  c(`PBS` = "PBS",
    `CARlo`="`A2-CAR`^lo",
    `CARmed` = "`A2-CAR`^med", 
    `CARhi` = "`A2-CAR`^hi",
    `CARmed_PBMClo` = "`A2-CAR`^med*PBMC^lo",
    `CARmed_PBMChi` = "`A2-CAR`^med*PBMC^hi"), 
  default = label_parsed)
ylab=expression(atop("Percent huCD45"^{"+"}," Cells in Blood (%)"))
Fig3G<-ggplot(blood%>%
                 filter(Group!="PBS",InjDays>0,!is.na(percent_CD45pos)),
               aes(x=InjDays/7,y=percent_CD45pos,colour=Group,group=Group,fill=Group))+
  geom_ribbon(aes(x=InjDays/7,y=percent_CD45pos.Median,
                  ymin = percent_CD45pos.LQ,ymax = percent_CD45pos.UQ,fill=Group),
              data=bloodTab%>%filter(Group!="PBS",InjDays>0,
                                     !is.na(percent_CD45pos.Median)),
              colour=NA,alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  geom_point(aes(x=InjDays/7,y=percent_CD45pos.Median,
                 colour = Group,group=Group,shape=Group), #
             data=bloodTab%>%filter(Group!="PBS",InjDays>0,
                                    !is.na(percent_CD45pos.Median)),
             size=1,alpha=0.8)+
  geom_line(aes(x=InjDays/7,y=percent_CD45pos.Median,colour = Group,group=Group), #
            data=bloodTab%>%filter(Group!="PBS",InjDays>0,
                                   !is.na(percent_CD45pos.Median)),
            size=0.5,alpha=0.8)+
  scale_x_continuous(limits=c(0,14),
                     breaks=pretty_breaks(n=8))+
  labs(x="",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[-1],
                      labels=c("CAR-lo",
                               "CAR-med",
                               "CAR-hi",
                               "CAR-med PBMC-lo",
                               "CAR-med PBMC-hi"))+
  scale_shape_manual(name="Group", values=GroupShapes[-1],
                     labels=c("CAR-lo",
                              "CAR-med",
                              "CAR-hi",
                              "CAR-med PBMC-lo",
                              "CAR-med PBMC-hi"))+
  scale_fill_manual(name="Group", values=GroupPalette[-1],
                    labels=c("CAR-lo",
                             "CAR-med",
                             "CAR-hi",
                             "CAR-med PBMC-lo",
                             "CAR-med PBMC-hi"))+
  geom_hline(yintercept = mean(blood$percent_CD45pos[blood$Group=="PBS"|
                                                       blood$InjDays<0],na.rm=T),colour="black",linetype=3)+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position ="bottom")+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3G

###### Figure 3H ######
my_labeller = as_labeller(
  c(`PBS` = "PBS",
    `CARlo`="`A2-CAR`^lo",
    `CARmed` = "`A2-CAR`^med", 
    `CARhi` = "`A2-CAR`^hi",
    `CARmed_PBMClo` = "`A2-CAR`^med*PBMC^lo",
    `CARmed_PBMChi` = "`A2-CAR`^med*PBMC^hi"), 
  default = label_parsed)
ylab=expression(paste("Percent Myc"^{"+"}," of huCD45"^{"+"}," in Blood (%)"))
Fig3H<-ggplot(blood%>%
                    filter(Group%in%c("CARhi","CARmed_PBMChi",
                                      "CARmed_PBMClo"),
                           InjDays>6,!is.na(percent_Mycpos_hCD45pos)),
                  aes(x=InjDays/7,y=percent_Mycpos_hCD45pos,colour=Group,group=Group,fill=Group))+
  geom_ribbon(aes(x=InjDays/7,y=percent_Mycpos_hCD45pos.Median,
                  ymin = percent_Mycpos_hCD45pos.LQ,ymax = percent_Mycpos_hCD45pos.UQ,fill=Group),
              data=bloodTab%>%
                filter(Group%in%c("CARhi","CARmed_PBMChi",
                                  "CARmed_PBMClo"),InjDays>6,
                       !is.na(percent_Mycpos_hCD45pos.Median)),
              colour=NA,alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  geom_point(aes(x=InjDays/7,y=percent_Mycpos_hCD45pos.Median,
                 colour = Group,group=Group,shape=Group), #
             data=bloodTab%>%
               filter(Group%in%c("CARhi","CARmed_PBMChi",
                                 "CARmed_PBMClo"),InjDays>6,
                      !is.na(percent_Mycpos_hCD45pos.Median)),
             size=1,alpha=0.8)+
  geom_line(aes(x=InjDays/7,y=percent_Mycpos_hCD45pos.Median,colour = Group,group=Group),
            data=bloodTab%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),InjDays>6,
                     !is.na(percent_Mycpos_hCD45pos.Median)),
            size=0.5,alpha=0.8)+
  scale_x_continuous(limits=c(0,14),
                     breaks=pretty_breaks(n=8))+
  labs(x="",y=ylab)+
scale_colour_manual(name="Group", values=GroupPalette[-c(1,2,3)],
                    labels=c("CAR-hi",
                             "CAR-med PBMC-lo",
                             "CAR-med PBMC-hi"))+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1,2,3)],
                     labels=c("CAR-hi",
                              "CAR-med PBMC-lo",
                              "CAR-med PBMC-hi"))+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1,2,3)],
                    labels=c("CAR-hi",
                             "CAR-med PBMC-lo",
                             "CAR-med PBMC-hi"))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position ="bottom")+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3H

###### Figure 3I ######
my_labeller = as_labeller(
  c(`PBS` = "PBS",
    `CARlo`="`A2-CAR`^lo",
    `CARmed` = "`A2-CAR`^med", 
    `CARhi` = "`A2-CAR`^hi",
    `CARmed_PBMClo` = "`A2-CAR`^med*PBMC^lo",
    `CARmed_PBMChi` = "`A2-CAR`^med*PBMC^hi"), 
  default = label_parsed)
ylab=expression(atop("Percent CD4"^{"+"},"of A2-CAR T Cells\nin Blood (%)"))
Fig3I<-ggplot(blood%>%
                       filter(Group%in%c("CARhi","CARmed_PBMChi",
                                         "CARmed_PBMClo"),
                              InjDays>6,!is.na(percent_CD4pos_of_CD45posMycpos)),
                     aes(x=InjDays/7,y=percent_CD4pos_of_CD45posMycpos,colour=Group,group=Group,fill=Group))+
  geom_ribbon(aes(x=InjDays/7,y=percent_CD4pos_of_CD45posMycpos.Median,
                  ymin = percent_CD4pos_of_CD45posMycpos.LQ,ymax = percent_CD4pos_of_CD45posMycpos.UQ,fill=Group),
              data=bloodTab%>%filter(Group%in%c("CARhi","CARmed_PBMChi",
                                                "CARmed_PBMClo"),
                                     InjDays>6,
                                     !is.na(percent_CD4posMycposCD45pos.Median)),
              colour=NA,alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  geom_line(aes(x=InjDays/7,y=percent_CD4pos_of_CD45posMycpos.Median,colour = Group,group=Group), #
            data=bloodTab%>%filter(Group%in%c("CARhi","CARmed_PBMChi",
                                              "CARmed_PBMClo"),
                                   InjDays>6,
                                   !is.na(percent_CD4posMycposCD45pos.Median)),
            size=0.5,alpha=0.8)+
  geom_point(aes(x=InjDays/7,y=percent_CD4pos_of_CD45posMycpos.Median,
                 colour = Group,group=Group,shape=Group), #
             data=bloodTab%>%filter(Group%in%c("CARhi","CARmed_PBMChi",
                                               "CARmed_PBMClo"),
                                    InjDays>6,
                                    !is.na(percent_CD4posMycposCD45pos.Median)),
             size=1,alpha=0.8)+
  scale_x_continuous(limits=c(0,15),
                     breaks=pretty_breaks(n=6))+
  labs(x="Weeks Post-\nA2-CAR T Cell Injection",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1,2,3)],
                      labels=c("CAR-hi",
                               "CAR-med PBMC-lo",
                               "CAR-med PBMC-hi"))+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1,2,3)],
                     labels=c("CAR-hi",
                              "CAR-med PBMC-lo",
                              "CAR-med PBMC-hi"))+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1,2,3)],
                    labels=c("CAR-hi",
                             "CAR-med PBMC-lo",
                             "CAR-med PBMC-hi"))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position ="bottom")+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3I

###### Figure 3J ######
my_labeller = as_labeller(
  c(`PBS` = "PBS",
    `CARlo`="`A2-CAR`^lo",
    `CARmed` = "`A2-CAR`^med", 
    `CARhi` = "`A2-CAR`^hi",
    `CARmed_PBMClo` = "`A2-CAR`^med*PBMC^lo",
    `CARmed_PBMChi` = "`A2-CAR`^med*PBMC^hi"), 
  default = label_parsed)
ylab=expression(atop("Percent CD8"^{"+"},"of A2-CAR T Cells\nin Blood (%)"))
Fig3J<-ggplot(blood%>%
                       filter(Group%in%c("CARhi","CARmed_PBMChi",
                                         "CARmed_PBMClo"),
                              InjDays>6,
                              !is.na(percent_CD8pos_of_CD45posMycpos)),
                     aes(x=InjDays/7,y=percent_CD8pos_of_CD45posMycpos,colour=Group,group=Group,fill=Group))+
  geom_ribbon(aes(x=InjDays/7,y=percent_CD8pos_of_CD45posMycpos.Median,
                  ymin = percent_CD8pos_of_CD45posMycpos.LQ,ymax = percent_CD8pos_of_CD45posMycpos.UQ,fill=Group),
              data=bloodTab%>%filter(Group%in%c("CARhi","CARmed_PBMChi",
                                                "CARmed_PBMClo"),
                                     InjDays>6,
                                     !is.na(percent_CD8pos_of_CD45posMycpos.Median)),
              colour=NA,alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  geom_line(aes(x=InjDays/7,y=percent_CD8pos_of_CD45posMycpos.Median,colour = Group,group=Group), #
            data=bloodTab%>%filter(Group%in%c("CARhi","CARmed_PBMChi",
                                              "CARmed_PBMClo"),
                                   InjDays>6,
                                   !is.na(percent_CD8pos_of_CD45posMycpos.Median)),
            size=0.5,alpha=0.8)+
  geom_point(aes(x=InjDays/7,y=percent_CD8pos_of_CD45posMycpos.Median,
                 colour = Group,group=Group,shape=Group),
             data=bloodTab%>%filter(Group%in%c("CARhi","CARmed_PBMChi",
                                               "CARmed_PBMClo"),
                                    InjDays>6,
                                    !is.na(percent_CD8pos_of_CD45posMycpos.Median)),
             size=1,alpha=0.8)+
  scale_x_continuous(limits=c(0,15),
                     breaks=pretty_breaks(n=6))+
  labs(x="Weeks Post-\nA2-CAR T Cell Injection",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1,2,3)],
                      labels=c("CAR-hi",
                               "CAR-med PBMC-lo",
                               "CAR-med PBMC-hi"))+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1,2,3)],
                     labels=c("CAR-hi",
                              "CAR-med PBMC-lo",
                              "CAR-med PBMC-hi"))+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1,2,3)],
                    labels=c("CAR-hi",
                             "CAR-med PBMC-lo",
                             "CAR-med PBMC-hi"))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+ 
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position ="bottom")+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig3J

#### Use patchwork to assemble panels ####
#use Fig3B as a placeholder for the schematic; edit in Inkscape
(Figure3<-(Fig3B)/(Fig3B+ theme(legend.position = "None") | Fig3C+plot_layout(guides = 'collect')) /
   (Fig3D+ theme(legend.position = "None"))/
   (Fig3E$plot+theme(axis.text.x = element_text(family = "Arial",color="black",size=12),
                          axis.text.y = element_text(family = "Arial",color="black",size=12),
                          axis.title.x = element_text(family = "Arial",color="black",size=12),
                          axis.title.y = element_text(family = "Arial",color="black",size=12)) | 
      Fig3F$plot+theme(axis.text.x = element_text(family = "Arial",color="black",size=12),
                          axis.text.y = element_text(family = "Arial",color="black",size=12),
                          axis.title.x = element_text(family = "Arial",color="black",size=12),
                          axis.title.y = element_text(family = "Arial",color="black",size=12)))/
   (Fig3G+ theme(legend.position = "None",
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) | Fig3H+ theme(legend.position = "None",
                                                                     axis.text.x = element_blank(),
                                                                     axis.ticks.x = element_blank())) /
   (Fig3I+ theme(legend.position = "None") | Fig3J+ theme(legend.position = "None")) + 
   plot_annotation(tag_levels = 'A') &
   theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('Figure3.svg', Figure3, width=8.5, height=9, units = "in", dpi=500)
ggsave('Figure3_tall.svg', Figure3, width=8.5, height=12, units = "in", dpi=500)
#use the taller version panels in the shorter version - patchwork squishes them a little