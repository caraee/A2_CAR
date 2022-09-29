#### Supplemental Figures S2 and S3 ####
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
weekly_monitoring<-read.table("Fig1_weekly_monitoring.txt",sep="\t",stringsAsFactors = F, 
                              header=T, check.names=F,as.is=T)
GSIS<-read.table("Fig1_GSIS.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)
#human C-peptide values are in pg/mL
#Blood glucose values are in mM
Mice<-read.table("Fig1_animalTable.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)

blood<-read.table("Fig1_blood_percentages.txt",sep="\t",stringsAsFactors = F, 
                  header=T, check.names=F,as.is=T)
blood_cell_counts<-read_delim("Fig1_blood_cell_counts.txt", name_repair = "universal")

GroupPalette<-c("#332288","#117733","#44AA99","#88CCEE",
                "#DDCC77","#CC6677","#AA4499","#882255")[c(4,2,8,1,6)]
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
blood$Group <- factor(blood$Group, levels = levels(weekly_monitoring$Group)) #,"S5","S6"
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

#### Statistical Analyses ####
#No inferential statistical analyses were performed on Figure S2

library("ggpubr")
ggqqplot(blood_cell_counts$normal.count_CD4posMycposCD45pos[blood_cell_counts$Group!="PBMC"&
                                                              blood_cell_counts$Group!="PBS"&blood_cell_counts$InjDays>3],
         ggtheme = theme_minimal())


blood_cell_counts<-blood_cell_counts%>%filter(rowSums(blood_cell_counts[,4:11],na.rm=T)>0)%>%
  mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                            InjDays==16~"End"))

mod1<-lm(Value~Parameter*Group*InjStage,
         data=blood_cell_counts%>%
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
                       str_count(row.names(dat),"normal.count_CD8posMycposCD45pos")==2|
                         str_count(row.names(dat),"normal.count_CD4posMycposCD45pos")==2|
                         str_count(row.names(dat),"normal.count_CD45pos")==2|
                         str_count(row.names(dat),"normal.count_Mycpos_hCD45pos")==2))

dat%>%filter(`p adj`<0.05)
dat%>%filter(str_count(row.names(dat),"normal.count_Mycpos_hCD45pos")>1)

#### Supplemental Figure S2 ####
my_labeller = as_labeller(
  c(`-15`="`-15`",
    `-1`="`-1`",
    `9`="9",
    `13`="13",
    `23`="23",
    PBS="PBS",
    PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
pd <- position_dodge(width=2)
breaks<-c(0,15,30)
limits<-c(-3,35)

###### Figure S2A ######
FigS2A<-ggplot(GSIS[GSIS$Time%in%breaks,],
              aes(x=Time,y=BG.Start,colour=Group,group=Group,fill=Group))+
  facet_grid(Group~InjDays,scales = "free_y",
             labeller=my_labeller)+
  geom_point(aes(x=Time,y=BG.Median,shape=Group),
             data=GSIStab[GSIStab$Time%in%breaks,],
             size=1,alpha=1,position=pd) +
  geom_line(aes(x=Time,y=BG.Median),data=GSIStab[GSIStab$Time%in%breaks,],
            size=0.5,alpha=1,position=pd) +
  geom_linerange(aes(ymin=BG.Start,ymax=BG.End,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  scale_fill_manual(values=GroupPalette)+
  scale_color_manual(values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  geom_line(aes(colour=Group,group=AnimalID),position=pd,linetype=2,alpha=0.5)+
  labs(x="Time (minutes)",y="Blood Glucose (mM)")+
  geom_hline(data=data.frame(yint = 33.3),aes(yintercept=yint),colour="#666666",linetype=3)+
  scale_y_continuous(limits=c(0, 35),breaks=c(0,10,20,30))+
  scale_x_continuous(breaks = breaks)+
  geom_ribbon(aes(x=Time,y=BG.Median,
                  ymin = BG.LQ, ymax = BG.UQ,fill=Group),colour=NA,
              data=GSIStab[GSIStab$Time%in%breaks&GSIStab$Weeks!=8,],
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
        strip.text.y = element_text(family="Arial",color="black",size=12))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle("Days Post-A2-CAR T Cell Injection")+
  theme(plot.title = element_text(family = "Arial",color="black", size=11, hjust=0.5))
FigS2A

###### Figure S2B ######
my_labeller = as_labeller(
  c(`-15`="`-15`",
    `-1`="`-1`",
    `9`="9",
    `13`="13",
    PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)

FigS2B<-ggplot(GSIS[GSIS$Time%in%breaks&GSIS$Group!="PBS",]%>%
                       filter(Days>-7),
                     aes(x=Time,y=CP.End/3020.29,colour=Group,group=Group,fill=Group))+
  facet_grid(Group~InjDays,scales = "fixed",
             labeller = my_labeller)+
  geom_point(aes(shape=Group),position=pd,size=1.75,alpha=0.5)+
  geom_line(aes(x=Time,y=CP.Median/3020.29),data=GSIStab[GSIStab$Group!="PBS",],
            size=0.5,alpha=1,position=pd) +
  geom_linerange(aes(ymin=CP.Start/3020.29,ymax=CP.End/3020.29,group=AnimalID),
                 linetype=3,position=pd,
                 size=0.3,alpha=0.8)+
  geom_ribbon(aes(x=Time,y=CP.Median/3020.29,ymin = CP.LQ/3020.29, ymax = CP.UQ/3020.29,
                  fill=Group),colour=NA,data=GSIStab[GSIStab$Time%in%breaks&
                                                       GSIStab$Group!="PBS",],
              size=0.5,alpha=0.2,position=pd)+
  ggtitle("Days Post-A2-CAR T Cell Injection")+
  theme(plot.title = element_text(family = "Arial",color="black", size=11, hjust=0.5))+
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_color_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=GroupShapes)+
  geom_line(aes(colour=Group,group=AnimalID),position=pd,size=0.3,linetype=2,alpha=0.5)+
  labs(x="Time (minutes)",y="Human C-peptide (nM)")+
  geom_hline(yintercept = 4.5*10/3020.29,linetype=3,alpha=0.7)+
  scale_x_continuous(limits=limits,breaks = breaks)+ 
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "None")+ 
  guides(fill=guide_legend(ncol=3),
         colour=guide_legend(ncol=3),
         shape=guide_legend(ncol=3))+
  theme(strip.text.x = element_text(family="Arial",color="black",size=12),
        strip.text.y = element_text(family="Arial",color="black",size=12),)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS2B

#### Supplementary Figure S3 ####
###### Figure S3A ####
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(paste("Number of huCD45"^{"+"}," Cells","in Blood"))
FigS3A<-ggplot(blood_cell_counts%>%
                        filter(Group!="PBS",
                               InjDays>3,!is.na(normal.count_CD45pos))%>%
                        mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                                  InjDays==16~"End")),
                      aes(x=InjStage,y=normal.count_CD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(aes(x=InjStage),
               alpha = 1/4)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  labs(x="",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[2:4])+
  scale_shape_manual(name="Group", values=GroupShapes[2:4])+
  scale_fill_manual(name="Group", values=GroupPalette[2:4])+
  geom_hline(yintercept = mean(blood_cell_counts$normal.count_CD45pos[blood_cell_counts$Group=="PBS"|
                                                                        blood_cell_counts$InjDays<0]),colour="black",linetype=3)+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks.y = element_line(colour="black"),
        axis.ticks.x = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank())+
  theme(strip.text = element_blank(),
        strip.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
FigS3A

###### Figure S3B ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(paste("Number of Myc"^{"+"}," huCD45"^{"+"}," in Blood"))
FigS3B<-ggplot(blood_cell_counts%>%
                           filter(Group!="PBS",
                                  InjDays>3,!is.na(normal.count_Mycpos_hCD45pos))%>%
                           mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                                     InjDays==16~"End")),
                         aes(x=InjStage,y=normal.count_Mycpos_hCD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(alpha = 1/4)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  labs(x="",y=ylab)+
  scale_colour_manual(name="Group", values=GroupPalette[-c(1)])+
  scale_shape_manual(name="Group", values=GroupShapes[-c(1)])+
  scale_fill_manual(name="Group", values=GroupPalette[-c(1)])+
  geom_hline(yintercept = mean(blood_cell_counts$normal.count_CD45pos[blood_cell_counts$Group=="PBS"|
                                                                        blood_cell_counts$InjDays<0]),colour="black",linetype=3)+
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
FigS3B

###### Figure S3C ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(atop("Number of CD4"^{"+"}," A2-CAR T Cells in Blood"))
FigS3C<-ggplot(blood_cell_counts%>%
                              filter(Group!="PBS",Group!="PBMC",
                                     InjDays>3,!is.na(normal.count_CD4posMycposCD45pos))%>%
                              mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                                        InjDays==16~"End")),
                            aes(x=InjStage,y=normal.count_CD4posMycposCD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(alpha = 1/4)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
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
FigS3C

###### Figure S3D ######
my_labeller = as_labeller(
  c(PBMC = "PBMC", 
    `CAR-D5-lo` = "`A2-CAR`^lo*`-D5`",
    `CAR-D5-hi` = "`A2-CAR`^hi*`-D5`",
    `CAR-D19-hi` = "`A2-CAR`^hi*`-D19`"), 
  default = label_parsed)
ylab=expression(atop("Number of CD8"^{"+"}," A2-CAR T Cells in Blood"))
FigS3D<-ggplot(blood_cell_counts%>%
                              filter(Group!="PBS",Group!="PBMC",
                                     InjDays>3,!is.na( normal.count_CD8posMycposCD45pos))%>%
                              mutate(InjStage=case_when(InjDays%in%c(9,13)~"Early",
                                                        InjDays==16~"End")),
                            aes(x=InjStage,y=normal.count_CD8posMycposCD45pos,colour=Group,fill=Group))+
  geom_point(aes(shape=Group,color=Group,fill=Group),na.rm=T,size=1,alpha=0.8,
             position=position_dodge2(width=0.1,preserve="total"))+
  geom_boxplot(aes(x=InjStage),alpha = 1/4,outlier.size = 0)+
  facet_grid(~Group,
             labeller = my_labeller)+
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
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
FigS3D

#### Assemble figures in patchwork ####
(FigureS2<-(FigS2A+ theme(legend.position = "None")) /
  (FigS2B+ theme(legend.position = "None"))+ 
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('FigureS2.svg', FigureS2, width=6, height=10.5, units = "in", dpi=600)

(FigureS3<-(FigS3A | FigS3B+ theme(legend.position = "None"))/
    (FigS3C+ theme(legend.position = "None") | FigS3D+ theme(legend.position = "None"))+ 
    plot_layout(guides = 'collect')+
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(family = "Arial",color="black",size=20)))

ggsave('FigureS3.svg', FigureS3, width=6, height=12, units = "in", dpi=600)