library("here")
library("tidyverse")
library("extrafont")
library("RColorBrewer")
library("scales")
library("grid")
library("gtable")
library("emmeans")
library("lme4")
#library("lmerTest")
library("data.table")
library("Amelia")


GSIS<-read.table("GSIS.txt",sep="\t",stringsAsFactors = F, 
                 header=T, check.names=F,as.is=T)

GroupPalette<-c("#c41b38",
                "#0142a9")

GSIS$Factor <- paste(GSIS$Group,GSIS$InjDay,"Days",sep=".")

GSIStab<-GSIS%>%
  group_by(Factor)%>%
  summarise(CP.UQ=quantile(CP,0.75,na.rm=T), CP.LQ=quantile(CP,0.25,na.rm=T),
            CP.Median=median(CP,na.rm=T))
GSIStab$InjDay<-sapply(GSIStab$Factor,function(f) unique(GSIS$InjDay[GSIS$Factor==f]))
GSIStab$Group<-sapply(GSIStab$Factor,function(f) unique(GSIS$Group[GSIS$Factor==f]))

my_labeller = as_labeller(
  c(`HLA-A2-` = "`HLA-A2`^—", 
    `HLA-A2+` = "`HLA-A2`^+"), 
  default = label_parsed)

group_labs<-c(`HLA-A2-`=expression(paste("HLA-A2"^neg)),
              `HLA-A2+`=expression(paste("HLA-A2"^pos)))

#group_labs<-c(`HLA-A2-`=expression(paste("HLA-A2"^"—")),
#              `HLA-A2+`=expression(paste("HLA-A2"^"+")))

Fig5B<-ggplot(GSIS,aes(x=InjDay,
                              y=CP/3020.29,colour=Group,fill=Group))+
  geom_point(aes(x=InjDay,y=CP.Median/3020.29,group=Group,shape=Group),
             data=GSIStab,
             size=1,alpha=0.8) +
  geom_line(aes(group=AnimalID,colour=Group),linetype=2,size=0.3,alpha=0.6)+
  geom_line(aes(x=InjDay,y=CP.Median/3020.29,group=Group),
            data=GSIStab,
            size=0.5,alpha=1) +
  geom_ribbon(aes(x=InjDay,y=CP.Median/1000,ymin = CP.LQ/3020.29, ymax = CP.UQ/3020.29,
                  group=Group),colour=NA,data=GSIStab,
              size=0.5,alpha=0.2)+
  scale_fill_manual(name="Group", values=GroupPalette,
                    labels=group_labs)+
  scale_color_manual(name="Group", values=GroupPalette,
                     labels=group_labs)+
  scale_shape_manual(name="Group", values=c(19,17),
                     labels=group_labs)+
  scale_y_continuous(limits=c(0,0.2),breaks=c(0,0.05,0.1,0.15,0.2))+
  labs(x="Days Post-\nA2-CAR T Cell Injection",y="Human C-peptide (nM)")+
  #facet_grid(~Weeks,scales="free_x",space="free_x")+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = c(0.25,0.22))+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank() ,
        strip.background = element_blank(),
        legend.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig5B

ggsave("Fig5B.svg",width = 3, height = 2.25, units = "in")
#ggsave("CP_symbs.png",path="./Figures/",width = 7.5, height = 7.5, units = "cm")

write_csv(GSIS,file="GSIS.csv")


mod<-lmer(CP~Group*InjDay+(1|AnimalID),
          data=GSIS%>%mutate(InjDay=factor(InjDay,levels=unique(InjDay))))
emm = emmeans(mod, ~ Group*InjDay)
pairs(emm)
pwpp(emm)
pwpm(emm)
