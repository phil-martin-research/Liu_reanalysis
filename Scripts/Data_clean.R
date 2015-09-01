##################################
#script to tidy up Liu et al data#
###to remove sites without age####
##################################

#read in data
library(gdata)
Liu<-read.csv("Data/Liu_et_al.csv",header=T)

#subset to give only columns of interest
Liu_sub<-Liu[-c(2:4,6,18:22)]
Liu_sub<-Liu_sub[-c(898:903),]
colnames(Liu_sub)<-c("ID","Site","Lat","Long","Mean_T","Mean_precip","AGB","L_AGB","T_AGB","Age","A_L_Ratio","A_T_Ratio","Ref")
Liu_sub<-subset(Liu_sub,!is.na(AGB))
#set levels for sites
levels(Liu_sub$Site)<-seq(1:449)

write.csv(Liu_sub,"Data/Liu_Aged.csv")
