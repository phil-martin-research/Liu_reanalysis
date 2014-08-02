##################################
#script to tidy up Liu et al data#
###to remove sites without age####
##################################

#read in data
<<<<<<< HEAD
library(gdata)
Liu <- read.xls("Data/geb12113-sup-0001-ts1.xlsx",sheet=1)
Liu<-read.csv("Data/Liu_et_al.csv",header=T)
=======
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Publications/Liu_et_al/Liu_reanalysis")
Liu<-read.csv("Liu_et_al.csv")
>>>>>>> a8cf70aa1cf6548b832224134b17a7e862e1a2cd

#subset to give only columns of interest
head(Liu)
as.factor(Liu$Site)
factor(Liu$Site)
Liu_sub<-Liu[-c(2:4,6,18:22)]
Liu_sub<-Liu_sub[-c(898:903),]
colnames(Liu_sub)<-c("ID","Site","Lat","Long","Mean_T","Mean_precip","AGB","L_AGB","T_AGB","Age","A_L_Ratio","A_T_Ratio","Ref")
Liu_sub<-subset(Liu_sub,AGB>0)
head(Liu_sub)
levels(Liu_sub$Site)<-seq(1:449)
Liu_sub$Age_sq<-Liu_sub$Age^2
Liu_sub2<-subset(Liu_sub,Age>0)

<<<<<<< HEAD
write.csv(Liu_sub2,"Liu_Aged.csv")
=======
write.csv(Liu_sub2,"Liu_Aged.csv")
>>>>>>> a8cf70aa1cf6548b832224134b17a7e862e1a2cd
