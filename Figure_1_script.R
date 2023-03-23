##### Figure 1 #####
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

##### Load data #####
weekly_monitoring<-read.table("./Figure1_data/Fig1_weekly_monitoring.txt",sep="\t",stringsAsFactors = F, 
                              header=T, check.names=F,as.is=T)
GSIS<-read.table("./Figure1_data/Fig1_GSIS.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)
#human C-peptide values are in pg/mL
#Blood glucose values are in mM
Mice<-read.table("./Figure1_data/Fig1_animalTable.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)

blood<-read.table("./Figure1_data/Fig1_blood_percentages.txt",sep="\t",stringsAsFactors = F, 
                  header=T, check.names=F,as.is=T)
blood_cell_counts<-read_delim("Fig1_blood_cell_counts.txt", name_repair = "universal")

GroupPalette<-c("#000000","#E69F00","#56B4E9","#009E73",
                "#F0E442","#0072B2","#D55E00","#CC79A7")[c(1,6,4,5,7)]

GroupShapes<-c(25,24,21,22,23)

#### Wrangle data ####
###### Mouse metadata ######
Mice$Termination<-as.POSIXct(strptime(Mice$Termination,"%Y/%m/%d",tz="America/Vancouver"))
Mice$InjDays <- sapply(1:nrow(Mice),function(i){
  if(Mice$Group[i]=="2.5e6 CAR-T: 19 days"){
    injDate<-as.POSIXct(strptime("2019/04/16","%Y/%m/%d",tz="America/Vancouver"))} else {
      injDate<-as.POSIXct(strptime("2019/04/02","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-Mice$Termination[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

Mice$Group<-str_replace(Mice$Group,"2.5e6 CAR-T: 5 days","CAR-D5-hi")%>%
  str_replace("2.5e6 CAR-T: 19 days","CAR-D19-hi")%>%
  str_replace("1e6 CAR-T: 5 days","CAR-D5-lo")%>%
  str_replace("0e6 CAR-T: 5 days","PBMC")%>%
  str_replace("Control","PBS")

###### Weekly blood glucose and body weight measurements ######
weekly_monitoring$Group<-sapply(weekly_monitoring$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})
weekly_monitoring$Date<-as.POSIXct(strptime(weekly_monitoring$Date,"%Y/%m/%d",tz="America/Vancouver"))
weekly_monitoring$Days <- sapply(1:nrow(weekly_monitoring),function(i){
  surgDate<-as.POSIXct(strptime("2019/03/28","%Y/%m/%d",tz="America/Vancouver"))
  d<-weekly_monitoring$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
weekly_monitoring$Weeks<-round(weekly_monitoring$Days/7)

weekly_monitoring$InjDays <- sapply(1:nrow(weekly_monitoring),function(i){
  if(weekly_monitoring$Group[i]=="CAR-D19-hi"){
    injDate<-as.POSIXct(strptime("2019/04/16","%Y/%m/%d",tz="America/Vancouver"))} else {
      injDate<-as.POSIXct(strptime("2019/04/02","%Y/%m/%d",tz="America/Vancouver"))
    }
  d<-weekly_monitoring$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

weekly_monitoring$Group<-factor(weekly_monitoring$Group,levels = c("PBS","PBMC","CAR-D5-lo",
                                                                   "CAR-D5-hi",
                                                                   "CAR-D19-hi"))
weekly_monitoring<-weekly_monitoring[order(weekly_monitoring$Weeks,
                                           weekly_monitoring$Group,weekly_monitoring$AnimalID),]
weekly_monitoring$AnimalID <- factor(weekly_monitoring$AnimalID, levels = as.character(unique(weekly_monitoring$AnimalID)))

names(GroupPalette)<-levels(weekly_monitoring$Group)
names(GroupShapes)<-levels(weekly_monitoring$Group)

weekly_monitoring$Factor<-interaction(weekly_monitoring$Group,weekly_monitoring$Date)

#No missing values
weeklyTab<-weekly_monitoring%>%
  group_by(Factor)%>%
  summarise(BG.UQ=quantile(BG,0.75,na.rm=T), BG.LQ=quantile(BG,0.25,na.rm=T),
            BG.Median=median(BG,na.rm=T),
            BW.UQ=quantile(BW,0.75,na.rm=T), BW.LQ=quantile(BW,0.25,na.rm=T),
            BW.Median=median(BW,na.rm=T))

weeklyTab$Group<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$Group[weekly_monitoring$Factor==f]))
weeklyTab$Group<-factor(weeklyTab$Group,levels=levels(weekly_monitoring$Group))
weeklyTab$Days<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$Days[weekly_monitoring$Factor==f]))
weeklyTab$InjDays<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$InjDays[weekly_monitoring$Factor==f]))
weeklyTab$Weeks<-sapply(weeklyTab$Factor,function(f) unique(weekly_monitoring$Weeks[weekly_monitoring$Factor==f]))

###### Glucose stimulated insulin secretion ######
GSIS$Date<-as.POSIXct(strptime(GSIS$Date,"%Y/%m/%d",tz="America/Vancouver"))
GSIS$Group<-sapply(GSIS$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})
GSIS$Group<-factor(GSIS$Group,levels = levels(weekly_monitoring$Group))

GSIS$Days <- sapply(1:nrow(GSIS),function(i){
  surgDate<-as.POSIXct(strptime("2019/03/28","%Y/%m/%d",tz="America/Vancouver"))
  d<-GSIS$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})
GSIS$Weeks<-round(GSIS$Days/7)
GSIS<-GSIS[order(GSIS$Weeks,GSIS$Group,GSIS$AnimalID),]
GSIS$AnimalID <- factor(GSIS$AnimalID, levels = as.character(unique(GSIS$AnimalID)))
GSIS$Factor <- paste(GSIS$Group,GSIS$Days,"days",GSIS$Time,"min",sep=".")

GSIS$InjDays <- sapply(1:nrow(GSIS),function(i){
  injDate<-case_when(grepl("19",GSIS$Group[i])~as.POSIXct(strptime("2019/04/16","%Y/%m/%d",tz="America/Vancouver")),
                     TRUE~as.POSIXct(strptime("2019/04/02","%Y/%m/%d",tz="America/Vancouver")))
  d<-GSIS$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

#C-peptide
GSIS$CP<-ifelse(GSIS$CP.End==GSIS$CP.Start,GSIS$CP.End,NA)
GSIS$BG<-ifelse(GSIS$BG.End==GSIS$BG.Start,GSIS$BG.End,NA)

#missing human C-peptide values are below the limit of detection - 45 pg/mL
pr <- matrix(c(0,15,log(0.01),log(45),0.99),byrow=T, nrow=1,ncol = 5)

set.seed(12345)
GSIS.imp <- amelia(GSIS, m = 20, p2s=0, idvars = c("Date","Days","Delivery","InjDays","Weeks",
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
plot(GSIS.imp,which.vars = 15)

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
GSIStab$Days<-sapply(GSIStab$Factor,function(f) unique(GSIS$Days[GSIS$Factor==f]))
GSIStab$InjDays<-sapply(GSIStab$Factor,function(f) unique(GSIS$InjDays[GSIS$Factor==f]))
GSIStab$Weeks<-sapply(GSIStab$Factor,function(f) unique(GSIS$Weeks[GSIS$Factor==f]))

###### Percentages of immune cells in blood ######
blood$Group<-sapply(blood$AnimalID,function(a){Mice$Group[Mice$AnimalID==a]})
blood$Group <- factor(blood$Group, levels = levels(weekly_monitoring$Group))
blood<-blood[order(blood$Group,blood$AnimalID),]

blood$TransplantDays <- sapply(1:nrow(blood),function(i){
  surgDate<-as.POSIXct(strptime("2019/03/28","%Y/%m/%d",tz="America/Vancouver"))
  d<-blood$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})

blood$InjDays <- sapply(1:nrow(blood),function(i){
  injDate<-case_when(grepl("19",blood$Group[i])~as.POSIXct(strptime("2019/04/16","%Y/%m/%d",tz="America/Vancouver")),
                     TRUE~as.POSIXct(strptime("2019/04/02","%Y/%m/%d",tz="America/Vancouver")))
  d<-blood$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

blood$CP<-vapply(1:nrow(blood),function(i){
  aID<-as.character(blood$AnimalID[i])
  d<-blood$TransplantDays[i]
  CP<-ifelse(is_empty(GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Days==d&GSIS$Time==30]),
             NA,
             GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Days==d&GSIS$Time==30]/1000)
},double(1))

blood<-blood[order(blood$InjDays,blood$Group,blood$AnimalID),]
blood$AnimalID <- factor(blood$AnimalID, levels = levels(weekly_monitoring$AnimalID))
blood$Factor<-paste(blood$Group,blood$InjDays,"days",sep=".")

bloodTab<-blood%>%
  group_by(Factor)%>%
  summarise(percent_CD4posMycposCD45pos.UQ=quantile(percent_CD4posMycposCD45pos,0.75,na.rm=T), 
            percent_CD4posMycposCD45pos.LQ=quantile(percent_CD4posMycposCD45pos,0.25,na.rm=T),
            percent_CD4posMycposCD45pos.Median=median(percent_CD4posMycposCD45pos,na.rm=T),
            
            percent_CD8posMycposCD45pos.UQ=quantile(percent_CD8posMycposCD45pos,0.75,na.rm=T), 
            percent_CD8posMycposCD45pos.LQ=quantile(percent_CD8posMycposCD45pos,0.25,na.rm=T),
            percent_CD8posMycposCD45pos.Median=median(percent_CD8posMycposCD45pos,na.rm=T),
            
            percent_CD45pos.min=quantile(percent_CD45pos,0.025,na.rm=T),
            percent_CD45pos.max=quantile(percent_CD45pos,0.95,na.rm=T),
            percent_CD45pos.UQ=quantile(percent_CD45pos,0.75,na.rm=T), 
            percent_CD45pos.LQ=quantile(percent_CD45pos,0.25,na.rm=T),
            percent_CD45pos.Median=median(percent_CD45pos,na.rm=T),
            
            percent_Mycpos_hCD45pos.min=quantile(percent_Mycpos_hCD45pos,0.025,na.rm=T), 
            percent_Mycpos_hCD45pos.max=quantile(percent_Mycpos_hCD45pos,0.95,na.rm=T),
            percent_Mycpos_hCD45pos.UQ=quantile(percent_Mycpos_hCD45pos,0.75,na.rm=T), 
            percent_Mycpos_hCD45pos.LQ=quantile(percent_Mycpos_hCD45pos,0.25,na.rm=T),
            percent_Mycpos_hCD45pos.Median=median(percent_Mycpos_hCD45pos,na.rm=T)
  )


bloodTab$Group<-sapply(bloodTab$Factor,function(f) unique(blood$Group[blood$Factor==f]))
bloodTab$Group<-factor(bloodTab$Group,levels=levels(blood$Group))
bloodTab$InjDays<-sapply(bloodTab$Factor,function(f) unique(blood$InjDays[blood$Factor==f]))
bloodTab$TransplantDays<-sapply(bloodTab$Factor,function(f) unique(blood$TransplantDays[blood$Factor==f]))

###### Counts of immune cells in blood ######
blood_cell_counts$Group<-sapply(blood_cell_counts$AnimalID,function(a){Mice$Group[Mice$AnimalID==a]})
blood_cell_counts$Group <- factor(blood_cell_counts$Group, levels = levels(weekly_monitoring$Group)) #,"S5","S6"
blood_cell_counts<-blood_cell_counts[order(blood_cell_counts$Group,blood_cell_counts$AnimalID),]

blood_cell_counts$TransplantDays <- sapply(1:nrow(blood_cell_counts),function(i){
  surgDate<-as.POSIXct(strptime("2019/03/28","%Y/%m/%d",tz="America/Vancouver"))
  d<-blood_cell_counts$Date[i]
  t<-round(as.numeric(difftime(d,surgDate,units="days")))
})

blood_cell_counts$InjDays <- sapply(1:nrow(blood_cell_counts),function(i){
  injDate<-case_when(grepl("19",blood_cell_counts$Group[i])~as.POSIXct(strptime("2019/04/16","%Y/%m/%d",tz="America/Vancouver")),
                     TRUE~as.POSIXct(strptime("2019/04/02","%Y/%m/%d",tz="America/Vancouver")))
  d<-blood_cell_counts$Date[i]
  t<-round(as.numeric(difftime(d,injDate,units="days")))
})

blood_cell_counts$CP<-vapply(1:nrow(blood_cell_counts),function(i){
  aID<-as.character(blood_cell_counts$AnimalID[i])
  d<-blood$TransplantDays[i]
  CP<-ifelse(is_empty(GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Days==d&GSIS$Time==30]),
             NA,
             GSIS$CP.Start[GSIS$AnimalID==aID&GSIS$Days==d&GSIS$Time==30]/1000)
},double(1))

blood_cell_counts<-blood_cell_counts[order(blood_cell_counts$InjDays,
                                           blood_cell_counts$Group,blood_cell_counts$AnimalID),]
blood_cell_counts$AnimalID <- factor(blood_cell_counts$AnimalID, levels = levels(weekly_monitoring$AnimalID))
blood_cell_counts$Factor<-paste(blood_cell_counts$Group,blood_cell_counts$InjDays,"days",sep=".")

bloodTab_counts<-blood_cell_counts%>%
  group_by(Factor)%>%
  summarise(normal.count_CD4posMycposCD45pos.UQ=quantile(normal.count_CD4pos.of.CD45posMycpos,0.75,na.rm=T), 
            normal.count_CD4posMycposCD45pos.LQ=quantile(normal.count_CD4pos.of.CD45posMycpos,0.25,na.rm=T),
            normal.count_CD4posMycposCD45pos.Median=median(normal.count_CD4pos.of.CD45posMycpos,na.rm=T),
            normal.count_CD8posMycposCD45pos.UQ=quantile(normal.count_CD4pos.of.CD45posMycpos,0.75,na.rm=T), 
            normal.count_CD8posMycposCD45pos.LQ=quantile(normal.count_CD4pos.of.CD45posMycpos,0.25,na.rm=T),
            normal.count_CD8posMycposCD45pos.Median=median(normal.count_CD4pos.of.CD45posMycpos,na.rm=T),
            normal.count_CD45pos.min=quantile(normal.count_CD45pos,0.025,na.rm=T),
            normal.count_CD45pos.max=quantile(normal.count_CD45pos,0.95,na.rm=T),
            normal.count_CD45pos.UQ=quantile(normal.count_CD45pos,0.75,na.rm=T), 
            normal.count_CD45pos.LQ=quantile(normal.count_CD45pos,0.25,na.rm=T),
            normal.count_CD45pos.Median=median(normal.count_CD45pos,na.rm=T),
            normal.count_Mycpos_hCD45pos.min=quantile(normal.count_Mycpos_hCD45pos,0.025,na.rm=T), 
            normal.count_Mycpos_hCD45pos.max=quantile(normal.count_Mycpos_hCD45pos,0.95,na.rm=T),
            normal.count_Mycpos_hCD45pos.UQ=quantile(normal.count_Mycpos_hCD45pos,0.75,na.rm=T), 
            normal.count_Mycpos_hCD45pos.LQ=quantile(normal.count_Mycpos_hCD45pos,0.25,na.rm=T),
            normal.count_Mycpos_hCD45pos.Median=median(normal.count_Mycpos_hCD45pos,na.rm=T)
  )


bloodTab_counts$Group<-sapply(bloodTab_counts$Factor,function(f) unique(blood_cell_counts$Group[blood_cell_counts$Factor==f]))
bloodTab_counts$Group<-factor(bloodTab_counts$Group,levels=levels(blood_cell_counts$Group))
bloodTab_counts$InjDays<-sapply(bloodTab_counts$Factor,function(f) unique(blood_cell_counts$InjDays[blood_cell_counts$Factor==f]))
bloodTab_counts$TransplantDays<-sapply(bloodTab_counts$Factor,function(f) unique(blood_cell_counts$TransplantDays[blood_cell_counts$Factor==f]))

##### Statistical Analysis #####
#Note that no inferential statistical analyses were performed on Figure 1B, 1C or 1D -
#The body weight and blood glucose values are clearly not different, and
#Figure 1D is a subset of the same data as Figure 1E, so we did not analyze the same data
#in two separate ways

GSIS$AnimalID<-as.character(GSIS$AnimalID)

Total_CP<-GSIS%>%as_tibble()%>%
  group_by(AnimalID,InjDays)%>%summarise(TotalCP=sum(CP.End))

Total_CP$Group<-sapply(Total_CP$AnimalID,function(a){
  Mice$Group[Mice$AnimalID==a]})

Total_CP<-Total_CP[Total_CP$Group!="PBS"&Total_CP$InjDays>1,]

Total_CP$Rejected<-case_when(Total_CP$TotalCP<=(45*2/0.8)&Total_CP$Group!="CAR-D19-hi"~1, #Total of limit of detection at 2 time points
                             Total_CP$TotalCP<=(45*3/0.8)&Total_CP$Group=="CAR-D19-hi"~1, #Total of limit of detection at 3 time points
                             TRUE~0)

sumsSurv<-tibble("AnimalID"=unique(Total_CP$AnimalID),
                 "Time"=sapply(unique(Total_CP$AnimalID),function(id) {
                   fifelse(sum(Total_CP$Rejected[Total_CP$AnimalID==id])>0,
                           min(Total_CP$InjDays[Total_CP$AnimalID==id&
                                                  Total_CP$Rejected==1]), 
                           13)}),
                 "Rejected"=sapply(unique(Total_CP$AnimalID),function(id) {
                   case_when(sum(Total_CP$Rejected[Total_CP$AnimalID==id])>0~1,
                             sum(Total_CP$Rejected[Total_CP$AnimalID==id])==0~0)}),
                 "Group"=unique(Total_CP$AnimalID) %>% 
                   map_chr(function(id) unique(Total_CP$Group[Total_CP$AnimalID==id])))
sumsSurv$Time[sumsSurv$Group=="CAR-D19-hi"]<-9

sumsSurv$Group<-factor(sumsSurv$Group,levels = c("PBMC","CAR-D5-lo",
                                                 "CAR-D5-hi",
                                                 "CAR-D19-hi"))

#Figure 1E
pairwise_survdiff(Surv(time=Time, event=Rejected,type="right") ~ Group,
                  data = sumsSurv, p.adjust.method = "BH",
                  na.action=options()$na.action,
                  rho = 0)

#Figure 1F-I
#drop rows where no cells were counted - technical errors
blood<-blood%>%filter(rowSums(blood[,3:6],na.rm=T)>0)%>%
  mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                            InjDays==16~"End"))

blood2<-blood

blood2%>%dplyr::filter(AnimalID=="CRISPR01M.02")
blood%>%dplyr::filter(AnimalID=="CRISPR01M.02")

mod1<-lm(Value~Parameter*Group*InjStage,
         data=blood%>%
           filter(InjDays>3)%>%
           filter(Group!="PBS")%>%
           mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                     InjDays==16~"End"))%>%
           pivot_longer(-c(AnimalID,Date,Group,
                           TransplantDays,Factor,CP,InjDays,InjStage),
                        names_to="Parameter",
                        values_to = "Value"))
mod1
an<-aov(mod1)
dat<-TukeyHSD(an, ordered = TRUE)
dat<-as.data.frame(dat$`Parameter:Group:InjStage`)
dat<-dat%>%filter((str_count(row.names(dat),"Early")==2|
                     str_count(row.names(dat),"End")==2)&(
                       str_count(row.names(dat),"percent_CD8posMycposCD45pos")==2|
                         str_count(row.names(dat),"percent_CD4posMycposCD45pos")==2|
                         str_count(row.names(dat),"percent_CD45pos")==2|
                         str_count(row.names(dat),"percent_Mycpos_hCD45pos")==2))

dat%>%filter(`p adj`<0.05)

g="PBMC"
blood_sub<-blood%>%dplyr::filter(Group==g)
blood_sub_early<-blood_sub%>%
  dplyr::filter(InjStage=="Early")
blood_sub_end<-blood_sub%>%
  dplyr::filter(InjStage=="End")
res_early_CD4<-t.test(blood_sub_early$percent_CD4posMycposCD45pos, 
                      mu = 50, alternative = "two.sided")
res_early_CD8<-t.test(blood_sub_early$percent_CD8posMycposCD45pos, 
                      mu = 50, alternative = "two.sided")
res_end_CD4<-t.test(blood_sub_end$percent_CD4posMycposCD45pos, 
                    mu = 50, alternative = "two.sided")
res_end_CD8<-t.test(blood_sub_end$percent_CD8posMycposCD45pos, 
                    mu = 50, alternative = "two.sided")

t_res<-tibble(pvalues=c(
  res_early_CD4$p.value,res_early_CD8$p.value,
  res_end_CD4$p.value,res_end_CD8$p.value
),
Estimate=c(res_early_CD4$estimate,res_early_CD8$estimate,
           res_end_CD4$estimate,res_end_CD8$estimate),
CI=c(paste(res_early_CD4$conf.int[1],res_early_CD4$conf.int[2],sep=","),
     paste(res_early_CD8$conf.int[1],res_early_CD8$conf.int[2],sep=","),
     paste(res_end_CD4$conf.int[1],res_end_CD4$conf.int[2],sep=","),
     paste(res_end_CD8$conf.int[1],res_end_CD8$conf.int[2],sep=",")),
Group=g,
Cell_Type=c("CD4","CD8","CD4","CD8"),
InjStage=c("Early","Early","End","End"))

for(g in c("CAR-D5-lo","CAR-D5-hi","CAR-D19-hi")){
  blood_sub<-blood%>%dplyr::filter(Group==g)
  blood_sub_early<-blood_sub%>%
    dplyr::filter(InjStage=="Early")
  blood_sub_end<-blood_sub%>%
    dplyr::filter(InjStage=="End")
  res_early_CD4<-t.test(blood_sub_early$percent_CD4posMycposCD45pos, 
                        mu = 50, alternative = "two.sided")
  res_early_CD8<-t.test(blood_sub_early$percent_CD8posMycposCD45pos, 
                        mu = 50, alternative = "two.sided")
  res_end_CD4<-t.test(blood_sub_end$percent_CD4posMycposCD45pos, 
                      mu = 50, alternative = "two.sided")
  res_end_CD8<-t.test(blood_sub_end$percent_CD8posMycposCD45pos, 
                      mu = 50, alternative = "two.sided")
  
  t_res_temp<-tibble(pvalues=c(
    res_early_CD4$p.value,res_early_CD8$p.value,
    res_end_CD4$p.value,res_end_CD8$p.value
  ),
  Estimate=c(res_early_CD4$estimate,res_early_CD8$estimate,
             res_end_CD4$estimate,res_end_CD8$estimate),
  CI=c(paste(res_early_CD4$conf.int[1],res_early_CD4$conf.int[2],sep=","),
       paste(res_early_CD8$conf.int[1],res_early_CD8$conf.int[2],sep=","),
       paste(res_end_CD4$conf.int[1],res_end_CD4$conf.int[2],sep=","),
       paste(res_end_CD8$conf.int[1],res_end_CD8$conf.int[2],sep=",")),
  Group=g,
  Cell_Type=c("CD4","CD8","CD4","CD8"),
  InjStage=c("Early","Early","End","End")
  )
  t_res<-bind_rows(t_res,t_res_temp)
}
t_res$p_adj<-p.adjust(t_res$pvalues, method = "BH")
t_res

#### Figures ####
# Figure 1A was made on Biorender
###### Figure 1B ######
Fig1B<-ggplot(weekly_monitoring[!is.na(weekly_monitoring$BW),]%>%
                 filter(Days>-7),
               aes(x=Days,y=BW,group=Group))+ 
  geom_ribbon(aes(y=BW.Median,ymin = BW.LQ,ymax = BW.UQ, fill = Group),
              data=weeklyTab[!is.na(weeklyTab$BW.Median),]%>%
                filter(Days>-7),alpha=0.2)+ 
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,linewidth=0.3,alpha=0.6)+ 
  geom_line(aes(y=BW.Median,colour = Group,group=Group), 
            data=weeklyTab[!is.na(weeklyTab$BW.Median),]%>%
              filter(Days>-7),size=0.5,alpha=0.8)+
  geom_point(aes(y=BW.Median,colour = Group,group=Group,shape=Group,fill=Group), #
             data=weeklyTab[!is.na(weeklyTab$BW.Median),]%>%
               filter(Days>-7),size=1,alpha=0.8)+
  geom_vline(xintercept=as.numeric(difftime("2019/04/16","2019/03/27",units="days")),
             colour=gray(0.20),linetype=2)+
  geom_vline(xintercept=as.numeric(difftime("2019/04/02","2019/03/27",units="days")),
             colour=gray(0.20),linetype=2)+
  labs(x="Days Post-Islet Transplantation",y="Body\nWeight (g)")+
  scale_colour_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_y_continuous(breaks=c(26,30,34))+
  scale_x_continuous(breaks=pretty_breaks(n=5))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = c(0.25,0.22))+
  theme(strip.text = element_text(family="Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1B

###### Figure 1C ######
Fig1C<-ggplot(weekly_monitoring[!is.na(weekly_monitoring$BG),]%>%
                 filter(Days>-7),
               aes(x=Days,y=BG,group=Group))+
  geom_ribbon(aes(y=BG.Median,ymin = BG.LQ,ymax = BG.UQ,fill=Group), 
              data=weeklyTab[!is.na(weeklyTab$BG.Median),]%>%
                filter(Days>-7),alpha=0.2)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,linewidth=0.3,alpha=0.6)+ 
  geom_line(aes(y=BG.Median,colour = Group,group=Group), #
            data=weeklyTab[!is.na(weeklyTab$BG.Median),]%>%
              filter(Days>-7),linewidth=0.5,alpha=0.8)+
  geom_point(aes(y=BG.Median,colour = Group,group=Group,shape=Group,fill=Group), 
             data=weeklyTab[!is.na(weeklyTab$BG.Median),]%>%
               filter(Days>-7),size=1,alpha=0.8)+
  scale_x_continuous(breaks=pretty_breaks(n=5))+
  geom_vline(xintercept=as.numeric(difftime("2019/04/16","2019/03/27",units="days")),
             colour=gray(0.20),linetype=2)+
  geom_vline(xintercept=as.numeric(difftime("2019/04/02","2019/03/27",units="days")),
             colour=gray(0.20),linetype=2)+
  labs(x="Days Post-Islet Transplantation",y="Blood\nGlucose (mM)")+
  scale_colour_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  scale_fill_manual(name="Group", values=GroupPalette)+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = c(0.35,0.22))+
  theme(strip.text = element_text(family="Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1C

###### Figure 1D ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
#divide pg/mL values by 3020.29 to convert to nM

Fig1D<-ggplot(GSIS%>%filter(Time==30)%>%filter(!Group%in%c("PBS")),
                     aes(x=InjDays,y=CP.End/3020.29,colour=Group,fill=Group))+
  geom_point(aes(x=InjDays,y=CP.Median/3020.29,group=Group,shape=Group),
             data=GSIStab%>%filter(Time==30)%>%
               filter(!Group%in%c("PBS")),size=1,alpha=0.8)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  facet_grid(.~Group,scales = "free_x",space="free_x",drop=T,
             labeller = my_labeller)+
  geom_line(aes(x=InjDays,y=CP.Median/3020.29,group=Group),
            data=GSIStab%>%filter(Time==30)%>%
              filter(!Group%in%c("PBS")),
            size=0.5,alpha=1) +
  geom_linerange(aes(ymin=CP.Start/3020.29,ymax=CP.End/3020.29,group=AnimalID),
                 linetype=3,size=0.3,alpha=0.8)+
  geom_ribbon(aes(x=InjDays,y=CP.Median/1000,ymin = CP.LQ/3020.29, ymax = CP.UQ/3020.29,
                  group=Group),colour=NA,data=GSIStab%>%filter(Time==30)%>%
                filter(!Group%in%c("PBS")),
              size=0.5,alpha=0.2)+
  scale_fill_manual(breaks=names(GroupPalette)[-c(1)],
                    name="Group", values=GroupPalette)+
  scale_color_manual(breaks=names(GroupPalette)[-c(1)],
                     name="Group", values=GroupPalette)+
  scale_shape_manual(breaks=names(GroupPalette)[-c(1)],
                     name="Group", values=GroupShapes)+
  labs(x="Days Post-A2-CAR T Cell Injection",y="Human C-peptide (nM)")+
  geom_hline(yintercept = 4.5*10/3020.29,linetype=3,alpha=0.7)+
  scale_y_continuous(breaks=pretty_breaks(4))+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = c(0.6,0.85))+
  guides(fill=guide_legend(ncol=2),
         colour=guide_legend(ncol=2),
         shape=guide_legend(ncol=2))+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank() ,
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1D

###### Figure 1E ######
fit <- survfit(Surv(time=Time, event=Rejected,type="right") ~ Group,
               data = sumsSurv)

SurvPal<-GroupPalette[2:5]
names(SurvPal)<-paste0("Group=",levels(sumsSurv$Group))

survLine<-c(1:4)
names(survLine)<-paste0("Group=",levels(sumsSurv$Group))

Fig1E<-ggsurvplot(fit, 
                       data = sumsSurv, pval = F,
                       break.time.by = 2,
                       risk.table = TRUE,
                       risk.table.col = "strata",
                       risk.table.height = 0.5,
                       risk.table.y.text=T,
                       palette = SurvPal,
                       linetype = survLine,
                       legend="none",
                       ylab="Fraction Without Rejection",
                       xlab="Days Post-A2-CAR T Cell Injection")
risk.table.ylab <- ""
Fig1E$plot<-Fig1E$plot + scale_x_continuous(breaks = seq(-2,16,2))
Fig1E$table <- Fig1E$table + labs(y = risk.table.ylab)
Fig1E

#Ideally the survival table should have been included in the manuscript

###### Figure 1F ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(paste("Percent huCD45"^{"+"}," Cells","in Blood (%)"))
Fig1F<-ggplot(blood%>%
                 filter(Group!="PBS",
                        InjDays>3,!is.na(percent_CD45pos))%>%
                 mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                           InjDays==16~"End")),
               aes(x=InjStage,y=percent_CD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(aes(x=InjStage),
               alpha = 1/4)+
  scale_y_continuous(limits = c(0,100))+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  labs(x="",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1)])+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1)])+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1)])+
  geom_hline(yintercept = mean(blood$percent_CD45pos[blood$Group=="PBS"|
                                                       blood$InjDays<0]),colour="black",linetype=3)+
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
  theme(strip.text = element_blank(),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1F

###### Figure 1G ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(paste("Percent Myc"^{"+"}," of huCD45"^{"+"}," in Blood (%)"))
Fig1G<-ggplot(blood%>%
                    filter(Group!="PBS",
                           InjDays>3,!is.na(percent_Mycpos_hCD45pos))%>%
                    mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                              InjDays==16~"End")),
                  aes(x=InjStage,y=percent_Mycpos_hCD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(alpha = 1/4)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  labs(x="",y=ylab)+
  scale_y_continuous(limits = c(0,100))+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1)])+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1)])+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1)])+
  geom_hline(yintercept = mean(blood$percent_CD45pos[blood$Group=="PBS"|
                                                       blood$InjDays<0]),colour="black",linetype=3)+
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
  theme(strip.text = element_blank(),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1G

###### Figure 1H ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(atop("Percent CD4"^{"+"},"of A2-CAR T Cells in Blood (%)"))
Fig1H<-ggplot(blood%>%
                       filter(Group!="PBS",Group!="PBMC",
                              InjDays>3,!is.na(percent_CD4posMycposCD45pos))%>%
                       mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                                 InjDays==16~"End")),
                     aes(x=InjStage,y=percent_CD4posMycposCD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(alpha = 1/4)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  labs(x="",y=ylab)+
  scale_y_continuous(limits=c(0,100),
                     breaks=pretty_breaks(n=6))+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1:2)])+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1:2)])+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1:2)])+
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
  theme(strip.text = element_blank(),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1H

###### Figure 1I ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(atop("Percent CD8"^{"+"},"of A2-CAR T Cells in Blood (%)"))
Fig1I<-ggplot(blood%>%
                       filter(Group!="PBS",Group!="PBMC",
                              InjDays>3,!is.na( percent_CD8posMycposCD45pos))%>%
                       mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                                 InjDays==16~"End")),
                     aes(x=InjStage,y=percent_CD8posMycposCD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(aes(x=InjStage),alpha = 1/4,outlier.size = 0)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  scale_y_continuous(limits=c(0,100),
                     breaks=pretty_breaks(n=6))+
  labs(x="",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1:2)])+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1:2)])+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1:2)])+
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
  theme(strip.text = element_blank(),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig1I

#### Use patchwork to assemble panels ####
#use Fig1B as a placeholder for the schematic; edit in Inkscape
(Figure1<-((Fig1B+ theme(legend.position = "None")) | plot_spacer() |((Fig1B+ theme(legend.position = "None"))/  
                                                                         (Fig1C+theme(legend.position = "None")))) /
    ((Fig1D)|
       (Fig1E$plot+theme(axis.text.x = element_text(family = "Arial",color="black",size=12),
                              axis.text.y = element_text(family = "Arial",color="black",size=12),
                              axis.title.x = element_text(family = "Arial",color="black",size=12),
                              axis.title.y = element_text(family = "Arial",color="black",size=12)))) / 
    ((Fig1F+ theme(legend.position = "None",
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank())) | (Fig1G+ theme(legend.position = "None",
                                                                         axis.text.x = element_blank(),
                                                                         axis.ticks.x = element_blank()))) /
    ((Fig1H+ theme(legend.position = "None")) | (Fig1I+ theme(legend.position = "None")))  + 
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('Figure1.svg', Figure1, width=7.5, height=9, units = "in", dpi=600)
