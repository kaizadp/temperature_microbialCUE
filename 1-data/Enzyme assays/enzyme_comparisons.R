

rm(list=ls())

library(ggplot2)



###################
######PCOA Analysis
###################
###################
setwd("~/OneDrive - PNNL/Documents/temperature experiment/Enzyme assays")

enzyme = read.table("enzyme_R_sheet.txt", sep="\t", header=TRUE,row.names=1)


enzyme$CUE = factor(enzyme$CUE, levels=c("High","Mid","Low","Blank"))

ggplot(enzyme, aes(x=CUE,y=uM.mass))+
  geom_boxplot()+
  facet_wrap(~Time)+
  theme_bw()

ggplot(enzyme, aes(x=Temperature,y=uM.mass))+
  geom_boxplot()+
  facet_wrap(~Time)+
  theme_bw()


T2_enzyme = enzyme[ which(enzyme$Time=="T2"),]
T1_enzyme = enzyme[ which(enzyme$Time=="T1"),]

T2_aov = aov(uM.mass~CUE, data=T2_enzyme)
summary(T2_aov)
TukeyHSD(T2_aov)

T1_aov = aov(uM.mass~CUE, data=T1_enzyme)
summary(T1_aov)
TukeyHSD(T1_aov)

ggplot(enzyme, aes(x=Temperature,y=uM.mass))+
  geom_boxplot()+
  facet_wrap(~Time)+
  theme_bw()

T2_aov = aov(uM.mass~Temperature, data=T2_enzyme)
summary(T2_aov)
TukeyHSD(T2_aov)

T1_aov = aov(uM.mass~Temperature, data=T1_enzyme)
summary(T1_aov)
TukeyHSD(T1_aov)


ggplot(T2_enzyme, aes(x=CUE, y=uM.mass))+
  geom_boxplot()+
  facet_wrap(~Temperature)+
  theme_bw()

T2_aov = aov(uM.mass~CUE+Temperature+Temperature*CUE, data=T2_enzyme)
summary(T2_aov)
TukeyHSD(T2_aov)

ggplot(T2_enzyme, aes(x=Temperature, y=uM.mass))+
  geom_boxplot()+
  facet_wrap(~CUE)+
  theme_bw()
