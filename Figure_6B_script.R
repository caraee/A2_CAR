library("tidyverse")
library("pzfx")

pzfx_tables("/home/cara/RProjects/MM_ACE02/MM_CP_hCD45_dat.pzfx")
df_CP <- read_pzfx("/home/cara/RProjects/MM_ACE02/MM_CP_hCD45_dat.pzfx", table="C-peptide in ACE-1 study (pMol/L)")%>%
  rename("A2-CAR"="Pherripheral A2-CAR",
         "PBS"="No CAR T cell")
head(df_CP)
df_CP<-df_CP%>%pivot_longer(cols=everything(),
                      names_to = "Group",
                      values_to="CP")%>%
  filter(Group%in%c("PBS","A2-CAR"))%>%
  mutate(Group=factor(Group,levels=c("PBS","A2-CAR")))

GroupPalette<-c("#D81B60",
                         "#1E88E5")
                         
Fig6B<-ggplot(df_CP,aes(x=Group,
                       y=CP,colour=Group,fill=Group))+
  geom_boxplot(alpha=0.8)+
  geom_point(aes(group=Group,shape=Group),colour="black",
             size=1,alpha=0.8,position = position_dodge2(width=0.1)) +
  geom_hline(yintercept = 2.2,colour="black",linetype=3)+
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_color_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=c(21,23))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=6))+
  labs(x="",y="Human C-peptide\n(pM)")+
  #facet_grid(~Weeks,scales="free_x",space="free_x")+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "None")+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank() ,
        strip.background = element_blank(),
        legend.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig6B
ggsave("Fig6B.png",path="./Figures/",width = 2, height = 2.5, units = "in")

df_hCD45 <- read_pzfx("/home/cara/RProjects/MM_ACE02/MM_CP_hCD45_dat.pzfx", table="Engraftment oif hCD45")%>%
  rename("Days"="ROWTITLE")
head(df_hCD45)
df_hCD45<-df_hCD45%>%pivot_longer(cols=2:10,
                      names_to = "Mouse",
                      values_to="hCD45")%>%
  mutate(Group=str_split_fixed(Mouse," [(]",2)[,1])

df_hCD45$Group<-df_hCD45$Group%>%str_replace("Peripheral A2-CAR","A2-CAR")%>%
  str_replace("No cell recipients","PBS")
df_hCD45<-df_hCD45%>%
  filter(Group%in%c("PBS","A2-CAR"))%>%
  mutate(Group=factor(Group,levels=c("PBS","A2-CAR")))

CD45tab<-df_hCD45%>%
  group_by(Days,Group)%>%
  summarise(CD45.UQ=quantile(hCD45,0.75,na.rm=T), CD45.LQ=quantile(hCD45,0.25,na.rm=T),
            CD45.Median=median(hCD45,na.rm=T))
GSIStab$InjDay<-sapply(GSIStab$Factor,function(f) unique(GSIS$InjDay[GSIS$Factor==f]))
GSIStab$Group<-sapply(GSIStab$Factor,function(f) unique(GSIS$Group[GSIS$Factor==f]))

Fig6C<-ggplot(df_hCD45,aes(x=Days,
                     y=hCD45,colour=Group,fill=Group))+
  geom_line(aes(group=Mouse,colour=Group),linetype=2,linewidth=0.3,alpha=0.6)+
  geom_line(aes(y=CD45.Median,group=Group),
            data=CD45tab,
            size=0.5,alpha=1) +
  geom_ribbon(aes(y=CD45.Median,ymin = CD45.LQ, ymax = CD45.UQ,
                  group=Group),colour=NA,data=CD45tab,
              size=0.5,alpha=0.2)+
  geom_point(aes(group=Group,shape=Group),colour="black",
             size=2,alpha=0.8,position = position_dodge2(width=0.1)) +
  scale_fill_manual(name="Group", values=GroupPalette)+
  scale_color_manual(name="Group", values=GroupPalette)+
  scale_shape_manual(name="Group", values=c(21,23))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=4))+
  labs(x="Days Post-A2-CAR\nT Cell Injection",y="hCD45/µL of Blood")+
  #facet_grid(~Weeks,scales="free_x",space="free_x")+
  theme(axis.title = element_text(family = "Arial", color="black", size=12),
        axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial",color="black",size=12))+
  theme(axis.text.y = element_text(family = "Arial",color="black",size=12))+
  theme(legend.text = element_text(family = "Arial",color="black",size=12), 
        legend.title = element_blank(),
        legend.position = "None")+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank() ,
        strip.background = element_blank(),
        legend.background = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
Fig6C
ggsave("Fig6C.png",path="./Figures/",width = 2.25, height = 2.75, units = "in")
