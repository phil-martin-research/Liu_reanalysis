#script to look at results from Lio et al 2013 paper in GEB

library(ggplot2)
library(nlme)
library(GGally)
library(ape)
library(cshapes)
library(MuMIn)
library(mgcv)

#read in data
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Publications/Liu_et_al")
Liu<-read.csv("Liu_et_al.csv")

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


#first we need to account for spatial autocorrelation
#and random effects from differences between studies
#fit a dummy random variable
Liu_sub$dummy<-rep(1,893)
#now build in spatial autocorrelation
head(Liu_sub)
Liu_sub$Lat_J<-Liu_sub$Lat+(rnorm(length(Liu_sub$Lat),0,0.01))
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu_sub)
corMatrix(cs1Exp)[1:10, 1:4]


#start models - these models are similar to those of liu et al
#considering everything independantly to each other
#but using our random variable structure

#first we fit a dummy random variable
Liu_sub2$dummy<-rep(1,572)
#now build in spatial autocorrelation
head(Liu_sub2)
Liu_sub2$Lat_J<-Liu_sub2$Lat+(rnorm(length(Liu_sub2$Lat),0,0.01))

cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu_sub2)
corMatrix(cs1Exp)[1:10, 1:4]

null.model<-lme(log(AGB)~1,data=Liu_sub2,random=~1|dummy)
null.model2<-lme(log(AGB)~1,data=Liu_sub2,random=~1|Ref)
null.model3<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long))
AICc(null.model,null.model2,null.model3)
#fitting a model that accounts for between study differences is better
#and we NEED to account for spatial autocorrelation for our results to
#be statistically valid

#now let's fit the global models of Liu et al for: 
#precipitation
Precip_model<-lme(log(AGB)~Mean_precip,data=Liu_sub2,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Temperature
Temp_model<-lme(log(AGB)~Mean_T,data=Liu_sub2,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Age
Age_model<-lme(log(AGB)~Age,data=Liu_sub2,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#all variables interacting
All_model<-lme(log(AGB)~Age*Mean_precip*Mean_T*Age_sq,data=Liu_sub2,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))
#additive model
Add_model<-lme(log(AGB)~Age+Mean_precip+Mean_T+Age_sq,data=Liu_sub2,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))


AICc(Precip_model,Temp_model,Age_model,All_model,Add_model)

r.squaredGLMM(Add_model)

Liu_sub2$Pred<-exp(predict(Add_model))

#now a global model to test all possibilities

ggplot(Liu_sub2,aes(Long,Lat,colour=Pred))+geom_point(alpha=0.8)+scale_colour_gradient(low="grey",high="dark green")

global.model<-update(null.model3,~Mean_precip*Mean_T*Age+Mean_precip*Mean_T*Age_sq)

plot(global.model)
summary(global.model)


#now we can dredge the model so that the value of each variable in predicting biomass
#can be assessed rather than using them in isolation
MS1<-dredge(global.model,evaluate=T,rank=AICc,trace=T,subset=dc(Age,Age_sq),REML=F)
poss_mod<-get.models(MS1,subset=delta<7)
modsumm<- model.sel(poss_mod, rank = "AICc")
modsumm2<-subset(modsumm,modsumm$delta<7)
modsumm2
averaged<-model.avg(modsumm)

top_model<-lme(log(AGB)~Age+Mean_precip+Mean_T+Age_sq+Mean_precip*Mean_T,data=Liu_sub2,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))
r.squaredGLMM(top_model)
AICc(top_model,Precip_model,Temp_model,Age_model)




Age<-seq(100,1200,1)
pred<-3.351+(0.003874*Age)+(-2.927e-06*Age^2)+(0.04635*4.4)+(0.0006104*984)+(6698.1*-2.613e-05)
precip<-seq(213,4000,1)
temp<-seq(-18,28,.1)
pred2<-3.351+(0.0006104*precip)+(162*0.003874)+(162*-2.927e-06)+(0.04635*4.4)+((precip*4.4)*-2.613e-05)
pred3<-3.351+(0.0006104*984)+(162*0.003874)+(162*-2.927e-06)+(0.04635*temp)+((984*temp)*-2.613e-05)


plot(Liu_sub2$Age,Liu_sub2$AGB)
lines(Age,exp(pred))

plot(Liu_sub2$Mean_precip,Liu_sub2$AGB)
lines(precip,exp(pred2))

plot(Liu_sub2$Mean_T,Liu_sub2$AGB)
lines(temp,exp(pred3))


Liu_sub2$pred<-(exp(predict(averaged)))[1:572]

head(Liu_sub2)

ggplot(data=Liu_sub2,aes(x=Mean_precip,y=Mean_T,size=AGB,colour=Age))+geom_point(shape=15,alpha=0.9)

Pred<-(exp(predict(averaged)))[1:572]
ggplot(data=Liu_sub2,aes(x=Age,y=AGB))+geom_point()+facet_wrap(~Ref)



#plot a map of model residuals
world <- cshp(date=as.Date("2008-1-1"))
world.points <- fortify(world, region='COWCODE')
p <- ggplot(world.points, aes(long,lat,group=group)) + geom_polygon(fill="grey")
p+geom_point(data=Liu_sub2,aes(Long,Lat,colour=Pred,group=NULL))


ggplot(data=Liu_sub2,aes(x=abs(Long)))+geom_histogram()

ggplot(data=Liu_sub2,aes(x=Mean_T,y=Mean_precip,size=AGB))+geom_point(alpha=0.5)
ggplot(data=Liu_sub2,aes(x=Mean_T,y=Age,size=AGB))+geom_point(alpha=0.5)
ggplot(data=Liu_sub2,aes(x=Mean_precip,y=Age,size=AGB))+geom_point(alpha=0.5)

#set a categorical temp bin
Liu_sub2$Temp_Cat<-NA
Temp_cat<-seq(from=-20,to=30,by=5)
for (i in 1:length(Temp_cat)) { 
Liu_sub2$Temp_Cat<-ifelse(test=Liu_sub2$Mean_T>=Temp_cat[i]&Liu_sub2$Mean_T<=Temp_cat[i+1],yes=Temp_cat[i],no=Liu_sub2$Temp_Cat)
}

#set a categorical precip bin
Liu_sub2$Precip_cat<-NA
Precip_cat<-seq(from=0,to=3000,by=500)
for (i in 1:length(Precip_cat)) { 
Liu_sub2$Precip_cat<-ifelse(test=Liu_sub2$Mean_precip>=Precip_cat[i]&Liu_sub2$Mean_precip<=Precip_cat[i+1],yes=Precip_cat[i+1],no=Liu_sub2$Precip_cat)
}

ggplot(data=Liu_sub2,aes(x=Age,y=AGB))+geom_point(shape=1,size=4)+facet_grid(Precip_cat~Temp_Cat,scales="free")



