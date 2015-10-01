Studies<-read.csv("Data/Studies_table.csv",)

library(plyr)

head(Studies)
Studies_sum<-ddply(Studies,.(Reference.1),summarise,M_Temp=mean(Temp,na.rm = T),M_Precip=mean(Precip,),m_Age=mean(Age,na.rm = T))
Studies2<-subset(Studies_sum,!is.na(m_Age))
unique(Studies)
write.csv(Studies2,"Results/Studies.csv",row.names=F)
