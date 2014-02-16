# this is a script to reanalyse the 
# data from Lio et al 2013 paper in GEB

library(ggplot2)
library(nlme)
library(GGally)
library(ape)
library(cshapes)
library(MuMIn)
library(mgcv)

#read in data
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Publications/Liu_et_al/Liu_reanalysis/Data")
Liu<-read.csv("Liu_Aged.csv")
head(Liu)

#start models - these models are similar to those of liu et al
#considering everything independantly to each other
#but using our random variable structure to account for
#spatial autocorrelation and any systematic differences amongst studies

#first we fit a dummy random variable
Liu$dummy<-rep(1,572)

#now build in spatial autocorrelation

#change coordinates slightly since some sites 
#have exactly the same coordinates
Liu$Lat_J<-Liu$Lat+(rnorm(length(Liu$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu)
corMatrix(cs1Exp)[1:10, 1:4]


#first run some null models to check our random
#variable structure is appropriate
null.model<-lme(log(AGB)~1,data=Liu,random=~1|dummy)
null.model2<-lme(log(AGB)~1,data=Liu,random=~1|Ref)
null.model3<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long))
AICc(null.model,null.model2,null.model3)

#fitting a model that accounts for between study differences is better
#and we NEED to account for spatial autocorrelation for our results to
#be statistically valid

#now let's fit the models that Liu et al used for: 
#precipitation - a model with a squared term for mean_precip
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Temperature
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2)+I(Mean_T^2),data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Age
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#now a global model for use in model averaging
#that contains all varibles that are needed
#I haven't included the squared and cubed terms for temp
#and precipitation becuase I don't think they make biological
#sense
All_model<-lme(log(AGB)~Age*Mean_precip*Mean_T*Age_sq,data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

plot(All_model)

#now we can dredge the model so that the value of each variable in predicting biomass
#can be assessed rather than using them in isolation
MS1<-dredge(global.model,evaluate=T,rank=AICc,trace=T,subset=dc(Age,Age_sq),REML=F)
poss_mod<-get.models(MS1,subset=delta<7)
modsumm<- model.sel(poss_mod, rank = "AICc",fit=T)
modsumm2<-subset(modsumm,modsumm$delta<7)
modsumm2
averaged<-model.avg(modsumm)

top_model<-lme(log(AGB)~Age+Mean_precip+Mean_T+Age_sq+Mean_precip*Mean_T,data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))
r.squaredGLMM(top_model)

AICc_res<-AICc(top_model,Precip_model,Temp_model,Age_model)

AICc_table<-data.frame(Model=row.names(AICc_res),AICc=AICc_res$AICc)
AICc_table$delta<-AICc_table$AICc-AICc_table$AICc[1]

#

Age<-seq(min(Liu$Age),max(Liu$Age),.1)
pred<-3.351+(0.003874*Age)+(-2.927e-06*Age^2)+(0.04635*4.4)+(0.0006104*984)+(6698.1*-2.613e-05)
Age_pred<-data.frame(Age,pred)
precip<-seq(213,4000,1)
pred2<-3.351+(0.0006104*precip)+(162*0.003874)+(162*-2.927e-06)+(0.04635*4.4)+((precip*4.4)*-2.613e-05)
Precip_pred<-data.frame(precip,pred2)

temp<-seq(-18,28,.1)

pred3<-3.351+(0.0006104*984)+(162*0.003874)+(162*-2.927e-06)+(0.04635*temp)+((984*temp)*-2.613e-05)
Temp_pred<-data.frame(temp,pred3)

a<-ggplot(Liu,aes(Age,AGB))+geom_point(shape=1,alpha=0.8)
a+geom_line(data=Age_pred,aes(Age,exp(pred)),size=1)

a<-ggplot(Liu,aes(Mean_precip,AGB))+geom_point(shape=1,alpha=0.8)
a+geom_line(data=Precip_pred,aes(precip,exp(pred2)),size=1)

a<-ggplot(Liu,aes(Mean_T,AGB))+geom_point(shape=1,alpha=0.8)
a+geom_line(data=Temp_pred,aes(temp,exp(pred3)),size=1)

PT_int<-expand.grid(Temp=seq(min(Liu$Mean_T),max(Liu$Mean_T),1),Precip=seq(min(Liu$Mean_precip),max(Liu$Mean_precip),10))
PT_int$int<-PT_int$Temp*PT_int$Precip
PT_int$pred<-3.351+(0.0006104*PT_int$Precip)+(162*0.003874)+(162*-2.927e-06)+(0.04635*PT_int$Temp)+((PT_int$Temp*PT_int$Precip)*-2.613e-05)
summary(Liu$Mean_precip*Liu$Mean_T)
PT_int2<-subset(PT_int,!(int>111900.0))
PT_int2<-subset(PT_int,!(Precip>3800))

hist(PT_int$int)

a<-ggplot(PT_int2,aes(Precip,Temp,fill=exp(pred)))+geom_raster()
a<-ggplot(Liu,aes(Mean_precip,Mean_T,fill=AGB))+geom_point()
a
a+geom_point(data=Liu,aes(Mean_precip,Mean_T,fill=NULL))


a+scale_fill_gradient(low="grey",high="dark green")
lines(Age,exp(pred))

plot()



plot(Liu$Mean_precip,Liu$AGB)
lines(precip,exp(pred2))

plot(Liu$Mean_T,Liu$AGB)
lines(temp,exp(pred3))


Liu$pred<-(exp(predict(averaged)))[1:572]

head(Liu)

ggplot(data=Liu,aes(x=Mean_precip,y=Mean_T,size=AGB,colour=Age))+geom_point(shape=15,alpha=0.9)

Pred<-(exp(predict(averaged)))[1:572]
ggplot(data=Liu,aes(x=Age,y=AGB))+geom_point()+facet_wrap(~Ref)



#plot a map of model residuals
world <- cshp(date=as.Date("2008-1-1"))
world.points <- fortify(world, region='COWCODE')
p <- ggplot(world.points, aes(long,lat,group=group)) + geom_polygon(fill="grey")
p+geom_point(data=Liu,aes(Long,Lat,colour=Pred,group=NULL))


ggplot(data=Liu,aes(x=abs(Long)))+geom_histogram()

ggplot(data=Liu,aes(x=Mean_T,y=Mean_precip,size=AGB))+geom_point(alpha=0.5)
ggplot(data=Liu,aes(x=Mean_T,y=Age,size=AGB))+geom_point(alpha=0.5)
ggplot(data=Liu,aes(x=Mean_precip,y=Age,size=AGB))+geom_point(alpha=0.5)

#set a categorical temp bin
Liu$Temp_Cat<-NA
Temp_cat<-seq(from=-20,to=30,by=5)
for (i in 1:length(Temp_cat)) { 
Liu$Temp_Cat<-ifelse(test=Liu$Mean_T>=Temp_cat[i]&Liu$Mean_T<=Temp_cat[i+1],yes=Temp_cat[i],no=Liu$Temp_Cat)
}

#set a categorical precip bin
Liu$Precip_cat<-NA
Precip_cat<-seq(from=0,to=3000,by=500)
for (i in 1:length(Precip_cat)) { 
Liu$Precip_cat<-ifelse(test=Liu$Mean_precip>=Precip_cat[i]&Liu$Mean_precip<=Precip_cat[i+1],yes=Precip_cat[i+1],no=Liu$Precip_cat)
}

ggplot(data=Liu,aes(x=Age,y=AGB))+geom_point(shape=1,size=4)+facet_grid(Precip_cat~Temp_Cat,scales="free")



