##### Supplemental Figures S4, S5, S6 and S7 #####
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
#No statistical analyses were performed on Figures S4 or S5 (C-peptide data is analyzed with survival analysis)
#Statistical analyses for Figure S6 is shown on Figure 3
#Figure S7A
CD45mod_slope_counts <- brm(normal.count_CD45pos~InjStage*Group+(InjStage|AnimalID),
                            data=blood_cell_counts%>%
                              mutate(InjStage=factor(InjWeeks),
                                     normal.count_CD45pos=normal.count_CD45pos+0.001)%>%
                              filter(Group!="PBS",InjDays>0,!is.na(normal.count_CD45pos)),
                            prior=c(set_prior("normal(-5,5)", class = "Intercept"),
                                    set_prior("normal(-5,5)", class = "b")),
                            iter=10000,
                            warmup=1000,
                            family="lognormal",
                            control = list(adapt_delta=0.975,
                                           max_treedepth=12),
                            #chains=1,
                            backend = "cmdstanr",
                            threads = threading(parallel::detectCores()),
                            silent=0,
                            save_pars = save_pars(all = TRUE))

cd45_res_tab_counts<-rbind(hypothesis(CD45mod_slope_counts,"exp(Intercept+GroupCARmed_PBMChi)-
                               exp(Intercept+GroupCARmed)>0")$hypothesis,
                           hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3+GroupCARmed+InjStage3:GroupCARmed)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+GroupCARmed_PBMChi)-
                               exp(Intercept+GroupCARmed_PBMClo)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                               exp(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+GroupCARmed_PBMChi)-
                               exp(Intercept)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage3+GroupCARmed+InjStage3:GroupCARmed)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage4+GroupCARmed+InjStage4:GroupCARmed)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage5+GroupCARmed_PBMClo+InjStage5:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage5+GroupCARmed+InjStage5:GroupCARmed)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage6+GroupCARmed+InjStage6:GroupCARmed)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage9+GroupCARmed+InjStage9:GroupCARmed)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage3+GroupCARmed_PBMClo+InjStage3:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage3)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage4+GroupCARmed_PBMClo+InjStage4:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage4)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage5+GroupCARmed_PBMClo+InjStage5:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage5)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage6+GroupCARmed_PBMClo+InjStage6:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage6)>0")$hypothesis,
                      hypothesis(CD45mod_slope_counts,"exp(Intercept+InjStage9+GroupCARmed_PBMClo+InjStage9:GroupCARmed_PBMClo)-
                               exp(Intercept+InjStage9)>0")$hypothesis)

#Figure S7B
MycCD45mod_counts_ZIP <- brm(normal.count_Mycpos_hCD45pos~InjStage*Group+(InjStage|AnimalID),
                             data=blood_cell_counts%>%
                               mutate(InjStage=factor(InjWeeks),
                                      normal.count_Mycpos_hCD45pos=round(normal.count_Mycpos_hCD45pos))%>%
                               filter(Group%in%c("CARhi","CARmed_PBMChi",
                                                 "CARmed_PBMClo"),
                                      InjDays>6,!is.na(normal.count_Mycpos_hCD45pos)),
                             prior=c(set_prior("normal(0,5)", class = "Intercept"),
                                     set_prior("normal(-5,5)", class = "b")),
                             iter=8000,
                             warmup=1000,
                             family="zero_inflated_poisson",
                             control = list(adapt_delta=0.99),
                             #chains=1,
                             backend = "cmdstanr",
                             threads = threading(parallel::detectCores()),
                             silent=1,
                             save_pars = save_pars(all = TRUE))

#super small number of divergent counts (2 of 28,000), but ESS and Rhats are perfect, so ignore
percentMycCD45_mod_counts<-blood_cell_counts%>%
  filter(Group%in%c("CARhi","CARmed_PBMChi",
                    "CARmed_PBMClo"),
         InjDays>6,!is.na(normal.count_Mycpos_hCD45pos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=factor(Group,levels=c("CARmed_PBMChi",
                                     "CARmed_PBMClo")))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(MycCD45mod_counts_ZIP,value="normal.count_Mycpos_hCD45pos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))

percentMycCD45_mod_counts%>%filter(Group=="CARmed_PBMChi")

MycCD45_res_tab_counts<-rbind(hypothesis(MycCD45mod_counts_ZIP,"exp(Intercept+GroupCARmed_PBMChi)-
                               exp(Intercept)>0")$hypothesis,
                              hypothesis(MycCD45mod_counts_ZIP,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3)>0")$hypothesis)

#Figure S7C
CD4MycCD45mod_counts_ZIP <- brm(normal.count_CD4pos.of.CD45posMycpos~InjStage*Group+(InjStage|AnimalID),
                                data=blood_cell_counts%>%
                                  mutate(InjStage=factor(InjWeeks),
                                         normal.count_CD4pos.of.CD45posMycpos=round(normal.count_CD4pos.of.CD45posMycpos,digits=0))%>%
                                  filter(Group%in%c("CARhi","CARmed_PBMChi",
                                                    "CARmed_PBMClo"),
                                         InjDays>6,!is.na(normal.count_CD4pos.of.CD45posMycpos)),
                                prior=c(set_prior("normal(5,10)", class = "Intercept"),
                                        set_prior("normal(-5,10)", class = "b")),
                                # prior=prior_CD4MycCD45_mix,
                                iter=8000,
                                warmup=1000,
                                family="zero_inflated_poisson",
                                control = list(adapt_delta=0.98,
                                               max_treedepth=12),
                                #chains=1,
                                backend = "cmdstanr",
                                threads = threading(parallel::detectCores()),
                                silent=0,
                                save_pars = save_pars(all = TRUE))
#super small number of divergent counts (35 of 28,000), but ESS and Rhats are perfect, so ignore

CD4MycCD45_res_tab_counts<-rbind(hypothesis(CD4MycCD45mod_counts_ZIP,"exp(Intercept+GroupCARmed_PBMChi)-
                               exp(Intercept)>0")$hypothesis,
                               hypothesis(CD4MycCD45mod_counts_ZIP,"exp(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      exp(Intercept+InjStage3)>0")$hypothesis)

#Fig S7D
#works up until here
CD8MycCD45mod_counts_ZIP <- brm(normal.count_CD8pos.of.CD45posMycpos~InjStage*Group+(InjStage|AnimalID),
                                data=blood_cell_counts%>%
                                  mutate(InjStage=factor(InjWeeks),
                                         normal.count_CD8pos.of.CD45posMycpos=round(normal.count_CD8pos.of.CD45posMycpos,digits=0))%>%
                                  filter(Group%in%c("CARhi","CARmed_PBMChi",
                                                    "CARmed_PBMClo"),
                                         InjDays>6,!is.na(normal.count_CD8pos.of.CD45posMycpos)),
                                prior=c(set_prior("normal(5,10)", class = "Intercept"),
                                        set_prior("normal(-5,10)", class = "b")),
                                # prior=prior_CD4MycCD45_mix,
                                iter=8000,
                                warmup=1000,
                                family="zero_inflated_poisson",
                                control = list(adapt_delta=0.98,
                                               max_treedepth=12),
                                #chains=1,
                                backend = "cmdstanr",
                                threads = threading(parallel::detectCores()),
                                silent=0,
                                save_pars = save_pars(all = TRUE))

CD8MycCD45_res_tab_counts<-rbind(hypothesis(CD8MycCD45mod_counts_ZIP,"(Intercept+GroupCARmed_PBMChi)-
                               (Intercept)>0")$hypothesis,
                               hypothesis(CD8MycCD45mod_counts_ZIP,"(Intercept+InjStage3+GroupCARmed_PBMChi+InjStage3:GroupCARmed_PBMChi)-
                      (Intercept+InjStage3)>0")$hypothesis)

#### Figures ####
###### Figure S4A #####
breaks<-c(0,15,60)
limits<-c(-3,65)
my_labeller = as_labeller(
  c(`-2`="`-2`",
    `-1`="`-1`",
    `1`="`1`",
    `3`="3",
    `5`="5",
    `7`="7",
    `9`="9",
    `11`="11",
    `12`="12",
    `PBS` = "PBS",
    `CARlo`="`A2-CAR`^lo",
    `CARmed` = "`A2-CAR`^med", 
    `CARhi` = "`A2-CAR`^hi",
    `CARmed_PBMClo` = "`A2-CAR`^med*PBMC^lo",
    `CARmed_PBMChi` = "`A2-CAR`^med*PBMC^hi"), 
  default = label_parsed)

FigS4A<-ggplot(GSIS[GSIS$Delivery==f&GSIS$Time%in%breaks,],
              aes(x=Time,y=BG.Start,colour=Group,group=Group,fill=Group))+
  facet_grid(Group~InjWeeks,scales = "free_y",
             labeller=my_labeller)+
  geom_point(aes(x=Time,y=BG.Median,shape=Group),
             data=GSIStab[GSIStab$Time%in%breaks&GSIStab$Delivery==f,],
             size=1,alpha=1,position=pd) +
  geom_line(aes(x=Time,y=BG.Median),data=GSIStab[GSIStab$Time%in%breaks&GSIStab$Delivery==f,],
            size=0.5,alpha=1,position=pd) +
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  scale_fill_manual(values=GroupPalette)+
  scale_color_manual(values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  geom_line(aes(colour=Group,group=AnimalID),position=pd,linetype=2,alpha=0.5)+
  labs(x="Time (minutes)",y="Blood Glucose (mM)")+
  geom_hline(data=data.frame(yint = 33.3,Species="Mouse"),aes(yintercept=yint),colour="#666666",linetype=3)+
  scale_y_continuous(limits=c(0, 35),breaks=c(0,10,20,30))+
  scale_x_continuous(breaks = breaks)+
  geom_ribbon(aes(x=Time,y=BG.Median,
                  ymin = BG.LQ, ymax = BG.UQ,fill=Group),colour=NA,
              data=GSIStab[GSIStab$Time%in%breaks&GSIStab$Weeks!=8&GSIStab$Delivery==f,],
              size=0.5,alpha=0.2)+
  theme(legend.position = c(0.5,0.10))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "bottom")+ 
  guides(fill=guide_legend(ncol=3),
         colour=guide_legend(ncol=3),
         shape=guide_legend(ncol=3))+
  theme(strip.text.x = element_text(family="Arial",color="black",size=12),
        strip.text.y = element_text(family="Arial",color="black",size=6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Weeks Post-A2-CAR T Cell Injection")+
  theme(plot.title = element_text(family = "Arial",color="black", size=11, hjust=0.5))
FigS4A

###### Figure S4B #####
breaks<-c(0,15,60)
limits<-c(-3,65)
my_labeller = as_labeller(
  c(`-2`="`-2`",
    `-1`="`-1`",
    `1`="`1`",
    `3`="3",
    `5`="5",
    `7`="7",
    `9`="9",
    `11`="11",
    `12`="12",
    `PBS` = "PBS",
    `CARlo`="`A2-CAR`^lo",
    `CARmed` = "`A2-CAR`^med", 
    `CARhi` = "`A2-CAR`^hi",
    `CARmed_PBMClo` = "`A2-CAR`^med*PBMC^lo",
    `CARmed_PBMChi` = "`A2-CAR`^med*PBMC^hi"), 
  default = label_parsed)

FigS4B<-ggplot(GSIS[GSIS$Time%in%breaks,],
                     aes(x=Time,y=CP.End/3020.29,colour=Group,group=Group,fill=Group))+
  facet_grid(Group~InjWeeks,scales = "free",labeller = my_labeller)+
  geom_point(aes(x=Time,y=CP.Median/3020.29,shape=Group),data=GSIStab,
             size=1,alpha=1,position=pd) +
  geom_line(aes(x=Time,y=CP.Median/3020.29),data=GSIStab,
            size=0.5,alpha=1,position=pd) +
  geom_linerange(aes(ymin=CP.Start/3020.29,ymax=CP.End/3020.29,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  geom_ribbon(aes(x=Time,y=CP.Median/3020.29,ymin = CP.LQ/3020.29, ymax = CP.UQ/3020.29,
                  fill=Group),colour=NA,data=GSIStab[GSIStab$Time%in%breaks,],
              size=0.5,alpha=0.2,position=pd)+
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_color_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  geom_line(aes(colour=Group,group=AnimalID),position=pd,size=0.3,linetype=2,alpha=0.5)+
  labs(x="Time (minutes)",y="Human C-peptide (nM)")+
  geom_hline(yintercept = 4.5*5/3020.29,linetype=3,alpha=0.7)+
  scale_x_continuous(limits=limits,breaks = breaks)+
  scale_y_continuous(breaks=pretty_breaks(3))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "bottom")+
  guides(fill=guide_legend(ncol=3),
         colour=guide_legend(ncol=3),
         shape=guide_legend(ncol=3))+
  theme(strip.text.x = element_text(family="Arial",color="black",size=12),
        strip.text.y = element_text(family="Arial",color="black",size=6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Weeks Post-A2-CAR T Cell Injection")+
  theme(plot.title = element_text(family = "Arial",color="black", size=12, hjust=0.5))
FigS4B

###### Figure S5 #####
FigS5<-ggplot(randomFed,
                    aes(x=InjWeeks,y=CP.End/3020.29,colour=Group,group=Group,fill=Group))+
  geom_point(aes(shape=Group),size=1,alpha=0.5)+
  geom_linerange(aes(ymin=CP.Start/3020.29,ymax=CP.End/3020.29,group=AnimalID),
                 linetype=3,size=0.3,alpha=0.8)+
  scale_fill_manual(name="Group", values=GroupPalette,
                    labels=c("PBS",
                             bquote("A2-CAR"^lo),
                             bquote("A2-CAR"^med),
                             bquote("A2-CAR"^hi),
                             bquote("A2-CAR"^med*" PBMC"^lo),
                             bquote("A2-CAR"^med*" PBMC"^hi)))+
  scale_color_manual(name="Group", values=GroupPalette,
                     labels=c("PBS",
                              bquote("A2-CAR"^lo),
                              bquote("A2-CAR"^med),
                              bquote("A2-CAR"^hi),
                              bquote("A2-CAR"^med*" PBMC"^lo),
                              bquote("A2-CAR"^med*" PBMC"^hi)))+
  scale_shape_manual(name="Group", values=GroupShapes,
                     labels=c("PBS",
                              bquote("A2-CAR"^lo),
                              bquote("A2-CAR"^med),
                              bquote("A2-CAR"^hi),
                              bquote("A2-CAR"^med*" PBMC"^lo),
                              bquote("A2-CAR"^med*" PBMC"^hi)))+
  geom_line(aes(colour=Group,group=AnimalID),size=0.8,linetype=1,alpha=0.8)+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y="Human C-peptide (nM)")+
  geom_hline(yintercept = 4.5*5/3020.29,linetype=3,alpha=0.7)+
  scale_x_continuous(breaks=unique(randomFed$InjWeeks))+
  theme(plot.title = element_text(family = "Arial", color="black",  size=32)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "bottom")+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(strip.text = element_text(family="Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS5

###### Figure S6A #####
percentCD45_mod<-blood%>%
  filter(Group!="PBS",InjDays>0,!is.na(percent_CD45pos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=as.character(Group))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="1"~1,
                            InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(CD45mod_slope,value="percent_CD45pos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

percentCD45_mod%>%filter(Group=="CARmed_PBMChi")

ylab=expression(atop("Percent huCD45"^{"+"}," Cells in Blood (%)"))
FigS6A<-ggplot(percentCD45_mod%>%
                     filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                   aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = percent_CD45pos),alpha=0.6) +
  geom_line(aes(y=percent_CD45pos,
                group=AnimalID,colour=Group),
            data=blood%>%
              filter(Group!="PBS",InjDays>0,
                     !is.na(percent_CD45pos)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(limits=c(0,14),
                     breaks=pretty_breaks(n=8))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[-1],
                      labels=c(bquote("A2-CAR"^lo),
                               bquote("A2-CAR"^med),
                               bquote("A2-CAR"^hi),
                               bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  geom_hline(yintercept = mean(blood$percent_CD45pos[blood$Group=="PBS"|
                                                       blood$InjDays<0],na.rm=T),colour="black",linetype=3)+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS6A

###### Figure S6B #####
percentMycCD45_mod<-blood%>%
  filter(Group%in%c("CARhi","CARmed_PBMChi",
                    "CARmed_PBMClo"),
         InjDays>6,!is.na(percent_Mycpos_hCD45pos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=factor(Group,levels=c("CARhi","CARmed_PBMChi",
                                     "CARmed_PBMClo")))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(MycCD45mod_noRE,value="percent_Mycpos_hCD45pos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

percentMycCD45_mod%>%filter(Group=="CARmed_PBMChi")

ylab=expression(paste("Percent Myc"^{"+"}," of huCD45"^{"+"}," in Blood (%)"))
FigS6B<-ggplot(percentMycCD45_mod%>%
                        filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                      aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = percent_Mycpos_hCD45pos),alpha=0.6) +
  geom_line(aes(y=percent_Mycpos_hCD45pos,
                group=AnimalID,colour=Group),
            data=blood%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),
                     InjDays>6,!is.na(percent_Mycpos_hCD45pos))%>%
              mutate(InjStage=factor(InjWeeks)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(limits=c(0,14),
                     breaks=pretty_breaks(n=8))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[-1],
                      labels=c(bquote("A2-CAR"^lo),
                               bquote("A2-CAR"^med),
                               bquote("A2-CAR"^hi),
                               bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS6B

###### Figure S6C #####
percentCD4MycCD45_mod<-blood%>%
  filter(Group%in%c("CARhi","CARmed_PBMChi",
                    "CARmed_PBMClo"),
         InjDays>6,!is.na(percent_CD4pos_of_CD45posMycpos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=factor(Group,levels=c("CARhi","CARmed_PBMChi",
                                     "CARmed_PBMClo")))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(CD4MycCD45mod_noRE,value="percent_CD4pos_of_CD45posMycpos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

percentCD4MycCD45_mod%>%filter(Group=="CARmed_PBMChi")

ylab=expression(atop("Percent CD4"^{"+"},"of A2-CAR T Cells in Blood (%)"))
FigS6C<-ggplot(percentCD4MycCD45_mod%>%
                           filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                         aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = percent_CD4pos_of_CD45posMycpos),alpha=0.6) +
  geom_line(aes(y=percent_CD4pos_of_CD45posMycpos,
                group=AnimalID,colour=Group),
            data=blood%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),
                     InjDays>6,!is.na(percent_CD4pos_of_CD45posMycpos))%>%
              mutate(InjStage=factor(InjWeeks)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(limits=c(0,14),
                     breaks=pretty_breaks(n=8))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[-1],
                      labels=c(bquote("A2-CAR"^lo),
                               bquote("A2-CAR"^med),
                               bquote("A2-CAR"^hi),
                               bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS6C

###### Figure S6D #####
percentCD8MycCD45_mod<-blood%>%
  filter(Group%in%c("CARhi","CARmed_PBMChi",
                    "CARmed_PBMClo"),
         InjDays>6,!is.na(percent_CD8pos_of_CD45posMycpos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=factor(Group,levels=c("CARhi","CARmed_PBMChi",
                                     "CARmed_PBMClo")))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(CD8MycCD45mod_noRE,value="percent_CD8pos_of_CD45posMycpos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

percentCD8MycCD45_mod%>%filter(Group=="CARmed_PBMChi")

ylab=expression(atop("Percent CD8"^{"+"},"of A2-CAR T Cells in Blood (%)"))
FigS6D<-ggplot(percentCD8MycCD45_mod%>%
                           filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                         aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = percent_CD8pos_of_CD45posMycpos),alpha=0.6) +
  geom_line(aes(y=percent_CD8pos_of_CD45posMycpos,
                group=AnimalID,colour=Group),
            data=blood%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),
                     InjDays>6,!is.na(percent_CD8pos_of_CD45posMycpos))%>%
              mutate(InjStage=factor(InjWeeks)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(limits=c(0,14),
                     breaks=pretty_breaks(n=8))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[-1],
                      labels=c(bquote("A2-CAR"^lo),
                               bquote("A2-CAR"^med),
                               bquote("A2-CAR"^hi),
                               bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS6D

###### Figure S7A #####
CD45_mod_counts<-blood_cell_counts%>%
  filter(Group!="PBS",InjDays>0,!is.na(normal.count_CD45pos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=as.character(Group))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="1"~1,
                            InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(CD45mod_slope_counts,value="normal.count_CD45pos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

CD45_mod_counts%>%filter(Group=="CARmed_PBMChi")

ylab=expression(atop("Number of huCD45"^{"+"}," Cells in Blood"))
FigS7A<-ggplot(CD45_mod_counts%>%
                            filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                          aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = normal.count_CD45pos),alpha=0.6) +
  geom_line(aes(y=normal.count_CD45pos,
                group=AnimalID,colour=Group),
            data=blood_cell_counts%>%
              filter(Group!="PBS",InjDays>0,
                     !is.na(normal.count_CD45pos)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(breaks=unique(blood_cell_counts$InjWeeks))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[-1],
                      labels=c(bquote("A2-CAR"^lo),
                               bquote("A2-CAR"^med),
                               bquote("A2-CAR"^hi),
                               bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  geom_hline(yintercept = mean(blood_cell_counts$normal.count_CD45pos[blood_cell_counts$Group=="PBS"|
                                                                        blood_cell_counts$InjDays<0],na.rm=T),colour="black",linetype=3)+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS7A

###### Figure S7B #####
ylab=expression(paste("Number of Myc"^{"+"}," huCD45"^{"+"}," Cells in Blood"))
FigS7B<-ggplot(percentMycCD45_mod_counts%>%
                               filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                             aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = normal.count_Mycpos_hCD45pos),alpha=0.6) +
  geom_line(aes(y=normal.count_Mycpos_hCD45pos,
                group=AnimalID,colour=Group),
            data=blood_cell_counts%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),
                     InjDays>6,!is.na(normal.count_Mycpos_hCD45pos))%>%
              mutate(InjStage=factor(InjWeeks)),
            linetype=2,size=0.7,alpha=0.6)+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_x_continuous(breaks=unique(blood_cell_counts$InjWeeks))+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[c(4:6)],
                      labels=c(bquote("A2-CAR"^hi),
                               bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS7B

###### Figure S7C #####
percentCD4MycCD45_mod_counts<-blood_cell_counts%>%
  filter(Group%in%c("CARhi","CARmed_PBMChi",
                    "CARmed_PBMClo"),
         InjDays>6,!is.na(normal.count_CD4pos.of.CD45posMycpos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=factor(Group,levels=c("CARmed_PBMChi",
                                     "CARmed_PBMClo")))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(CD4MycCD45mod_counts_ZIP,value="normal.count_CD4pos.of.CD45posMycpos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

percentCD4MycCD45_mod_counts%>%filter(Group=="CARmed_PBMClo")

ylab=expression(atop("Number of CD4"^{"+"}," A2-CAR T Cells in Blood"))
FigS7C<-ggplot(percentCD4MycCD45_mod_counts%>%
                                  filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                                aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = normal.count_CD4pos.of.CD45posMycpos),alpha=0.6) +
  geom_line(aes(y=normal.count_CD4pos.of.CD45posMycpos,
                group=AnimalID,colour=Group),
            data=blood_cell_counts%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),
                     InjDays>6,!is.na(normal.count_CD4pos.of.CD45posMycpos))%>%
              mutate(InjStage=factor(InjWeeks)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(breaks=unique(blood_cell_counts$InjWeeks))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[c(5,6)],
                      labels=c(bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS7C

###### Figure S7D #####
percentCD8MycCD45_mod_counts<-blood_cell_counts%>%
  filter(Group%in%c("CARhi","CARmed_PBMChi",
                    "CARmed_PBMClo"),
         InjDays>6,!is.na(normal.count_CD8pos.of.CD45posMycpos))%>%
  mutate(InjStage=factor(InjWeeks),
         Group=factor(Group,levels=c("CARmed_PBMChi",
                                     "CARmed_PBMClo")))%>%
  data_grid(
    InjStage,
    Group
  )%>%group_by(Group)%>%
  mutate(InjWeeks=case_when(InjStage=="2"~2,
                            InjStage=="3"~3,
                            InjStage=="4"~4,
                            InjStage=="5"~5,
                            InjStage=="6"~6,
                            InjStage=="7"~7,
                            InjStage=="8"~8,
                            InjStage=="9"~9,
                            InjStage=="11"~11,
                            InjStage=="12"~12,
                            InjStage=="13"~13,
                            InjStage=="14"~14))%>%
  add_epred_rvars(CD8MycCD45mod_counts_ZIP,value="normal.count_CD8pos.of.CD45posMycpos",
                  re_formula = NA)%>%
  filter(!(InjWeeks>3&Group=="CARmed_PBMChi"))%>%
  filter((Group!="CARhi"&InjWeeks%in%c(2,3,4,5,6,7,9,12,13))|
           (Group=="CARhi"&InjWeeks%in%c(2,3,4,6,8,11,13,14)))

percentCD8MycCD45_mod_counts%>%filter(Group=="CARmed_PBMChi")

ylab=expression(atop("Number of CD8"^{"+"}," A2-CAR T Cells in Blood"))
FigS7D<-ggplot(percentCD8MycCD45_mod_counts%>%
                                  filter(!(InjWeeks>3&Group=="CARmed_PBMChi")),
                                aes(x=InjWeeks,group=Group,colour=Group))+
  stat_lineribbon(aes(dist = normal.count_CD8pos.of.CD45posMycpos),alpha=0.6) +
  geom_line(aes(y=normal.count_CD8pos.of.CD45posMycpos,
                group=AnimalID,colour=Group),
            data=blood_cell_counts%>%
              filter(Group%in%c("CARhi","CARmed_PBMChi",
                                "CARmed_PBMClo"),
                     InjDays>6,!is.na(normal.count_CD8pos.of.CD45posMycpos))%>%
              mutate(InjStage=factor(InjWeeks)),
            linetype=2,size=0.7,alpha=0.6)+
  scale_x_continuous(breaks=unique(blood_cell_counts$InjWeeks))+
  labs(x="Weeks Post-A2-CAR T Cell Injection",y=ylab)+
  scale_fill_brewer(palette = "Greys") +
  scale_colour_manual(name="Group", values=GroupPalette[c(5,6)],
                      labels=c(bquote("A2-CAR"^med*" PBMC"^lo),
                               bquote("A2-CAR"^med*" PBMC"^hi)))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+ #,hjust=1,angle = 45
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_text(family="Arial",color="black",size=12),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS7D

#### Use patchwork to assemble panels ####
(FigureS4<-(FigS4A+ theme(legend.position = "None")) /
  (FigS4B+ theme(legend.position = "None"))+ 
  #plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('FigureS4.svg', FigureS4, width=6, height=12, units = "in", dpi=600)

(FigureS5<-FigS5+
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('FigureS5.svg', FigureS5, width=4.5, height=4.5, units = "in", dpi=600)

(FigureS6<-(FigS6A | FigS6B)/
    (FigS6C | FigS6D)+ 
    plot_layout(guides = 'collect')+
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('FigureS6.svg', FigureS6, width=9, height=4.5, units = "in", dpi=600)

(FigureS7<-(FigS7A | FigS7B)/
    (FigS7C | FigS7D)+ 
    plot_layout(guides = 'collect')+
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('../CRISPR02/Figures/20220824_FigureS7.svg', FigureS7, width=9, height=4.5, units = "in", dpi=600)