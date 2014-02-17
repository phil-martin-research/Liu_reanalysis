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

hist(Liu$Age)
hist(Liu$Mean_T)
hist(Liu$Mean_precip)

#looking at these histograms I wouldn't be happy to make predictions
#about biomass in forests of Age >500 years, Temp <-5 or >20, or precip >3000
#so I will subset the data to remove these
Liu_subset<-subset(Liu,Age<400&Mean_T>-5&Mean_T<20&Mean_precip<3000)
hist(Liu_subset$Age)
hist(Liu_subset$Mean_T)
hist(Liu_subset$Mean_precip)


#start models - these models are similar to those of liu et al
#considering everything independantly to each other
#but using our random variable structure to account for
#spatial autocorrelation and any systematic differences amongst studies

#first we fit a dummy random variable
Liu_subset$dummy<-rep(1,501)

#now build in spatial autocorrelation

#change coordinates slightly since some sites 
#have exactly the same coordinates
Liu_subset$Lat_J<-Liu$Lat+(rnorm(length(Liu$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu_subset)
corMatrix(cs1Exp)[1:10, 1:4]

#first run some null models to check our random
#variable structure is appropriate
null.model<-lme(log(AGB)~1,data=Liu_subset,random=~1|dummy)
null.model2<-lme(log(AGB)~1,data=Liu_subset,random=~1|Ref)
null.model3<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long))
AICc(null.model,null.model2,null.model3)

#fitting a model that accounts for between study differences is better
#and we NEED to account for spatial autocorrelation for our results to
#be statistically valid

#now let's fit the models that Liu et al used for: 
#precipitation - a model with a squared term for mean_precip
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu_subset,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Temperature
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2)+I(Mean_T^3),data=Liu_subset,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Age
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu_subset,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#now a global model for use in model averaging
#that contains all varibles that are needed
#I haven't included the squared and cubed terms for temp
#and precipitation becuase I don't think they make biological
#sense
All_model<-lme(log(AGB)~Age*Mean_precip*Mean_T*Age_sq,data=Liu_subset,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

plot(All_model)

#now we can dredge the model so that the value of each variable in predicting biomass
#can be assessed rather than using them in isolation
MS1<-dredge(All_model,evaluate=T,rank=AICc,trace=T,subset=dc(Age,Age_sq),REML=F)
poss_mod<-get.models(MS1,subset=delta<7)
modsumm<- model.sel(poss_mod, rank = "AICc",fit=T)
modsumm2<-subset(modsumm,modsumm$delta<7)
modsumm2
averaged<-model.avg(modsumm,fit=T)

#run the top model from the model averaging to get
#r squared statistic
top_model<-lme(log(AGB)~Age+Mean_precip+Mean_T+Age_sq+Mean_precip*Mean_T,data=Liu_subset,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))


#create table of different model AICc, AICc delta and marginal R squared
AICc_res<-AICc(top_model,Precip_model,Temp_model,Age_model)
AICc_table<-data.frame(Model=row.names(AICc_res),AICc=AICc_res$AICc)
AICc_table$delta<-AICc_table$AICc-AICc_table$AICc[1]
AICc_table$R_squared<-c(r.squaredGLMM(top_model)[1],r.squaredGLMM(Precip_model)[1],r.squaredGLMM(Temp_model)[1],r.squaredGLMM(Age_model)[1])
AICc_table

#create dataframe for predictions so paramaters can be plotted
#first - Age
nseq <- function(x, len = length(x)) seq(min(x, na.rm = TRUE),
    max(x, na.rm=TRUE), length = len)
head(Liu_sub2)
newdata_Age<- as.data.frame(lapply(lapply(Liu_subset[c(5,6,10,14)], mean), rep, 300))
newdata_Age$Age<- nseq(Liu_sub2$Age, nrow(newdata_Age))
newdata_Age$Age_sq<- newdata_Age$Age^2
newdata_Age$pred<-predict(top_model,newdata_Age,level=0)
#add predictions from model which is most similar to that of Liu et al
newdata_Age$Liu_pred<-predict(Age_model,newdata_Age,level=0)
#create design matrix
Designmat <- model.matrix(eval(eval(top_model$call$fixed)[-2]), newdata_Age[-ncol(newdata_Age)])

#compute standard error for predictions
predvar <- diag(Designmat %*% top_model$varFix %*% t(Designmat))
newdata_Age$SE <- sqrt(predvar) 
newdata_Age$SE2 <- sqrt(predvar+top_model$sigma^2)

#now plot the predictions for Age
theme_set(theme_bw(base_size=20))
Age_plot<-ggplot(data=newdata_Age,aes(x=Age,y=exp(pred)))+geom_line()
Age_plot2<-Age_plot+geom_line(data=newdata_Age,aes(y=exp(Liu_pred)),colour="red")
Age_plot3<-Age_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Age_plot4<-Age_plot3+ylab(expression(paste("Aboveground biomass (Mg ",ha^-1,")")))+xlab("Age (Years)")

#now make predictions for Temp - with a range of different precipitation levels
#to show interaction between the two
plot(Liu_sub2$Mean_T,Liu_sub2$Mean_precip)
#predictions for 1000, 2000 and 3000mm seem sensible
head(Liu_sub2)
newdata_Temp1<- as.data.frame(lapply(lapply(Liu_sub2[c(5,6,10,14)], mean), rep, 300))
newdata_Temp1$Mean_T<- nseq(Liu_sub2$Mean_T, nrow(newdata_Temp1))
newdata_Temp1$Mean_precip<-1000
newdata_Temp2<- as.data.frame(lapply(lapply(Liu_sub2[c(5,6,10,14)], mean), rep, 300))
newdata_Temp2$Mean_T<- nseq(Liu_sub2$Mean_T, nrow(newdata_Temp2))
newdata_Temp2$Mean_precip<-2000
newdata_Temp3<- as.data.frame(lapply(lapply(Liu_sub2[c(5,6,10,14)], mean), rep, 300))
newdata_Temp3$Mean_T<- nseq(Liu_sub2$Mean_T, nrow(newdata_Temp3))
newdata_Temp3$Mean_precip<-3000
newdata_Temp<-rbind(newdata_Temp1,newdata_Temp2,newdata_Temp3)

#now use predict
newdata_Temp$pred<-predict(top_model,newdata_Temp,level=0)
#add predictions from model which is most similar to that of Liu et al
newdata_Temp$Liu_pred<-predict(Temp_model,newdata_Temp,level=0)

#now plot the predictions for Temp
theme_set(theme_bw(base_size=20))
Temp_plot<-ggplot(data=newdata_Temp,aes(x=Mean_T,y=exp(pred),group=as.factor(Mean_precip),colour=as.factor(Mean_precip)))+geom_line()
Temp_plot2<-Temp_plot+geom_line(data=newdata_Temp,aes(y=exp(Liu_pred)),colour="red")
Temp_plot3<-Temp_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Temp_plot4<-Temp_plot3+ylab(expression(paste("Aboveground biomass (Mg ",ha^-1,")")))+xlab("Mean annual temperature")
Temp_plot4

#now predictions for changes in precipitation
newdata_P1<- as.data.frame(lapply(lapply(Liu_sub2[c(5,6,10,14)], mean), rep, 300))
newdata_P1$Mean_precip<- nseq(Liu_sub2$Mean_precip, nrow(newdata_P1))
newdata_P1$Mean_T<-0
newdata_P2<- as.data.frame(lapply(lapply(Liu_sub2[c(5,6,10,14)], mean), rep, 300))
newdata_P2$Mean_precip<- nseq(Liu_sub2$Mean_precip, nrow(newdata_P2))
newdata_P2$Mean_T<-10
newdata_P3<- as.data.frame(lapply(lapply(Liu_sub2[c(5,6,10,14)], mean), rep, 300))
newdata_P3$Mean_precip<- nseq(Liu_sub2$Mean_precip, nrow(newdata_P3))
newdata_P3$Mean_T<-20
newdata_P<-rbind(newdata_P1,newdata_P2,newdata_P3)

summary(newdata_P)


#now use predict
newdata_P$pred<-predict(top_model,newdata_P,level=0)
plot(newdata_P$Mean_precip,exp(newdata_P$pred))

#add predictions from model which is most similar to that of Liu et al
newdata_P$Liu_pred<-predict(Precip_model,newdata_P,level=0)

#now plot
theme_set(theme_bw(base_size=20))
Temp_plot<-ggplot(data=newdata_P,aes(x=Mean_precip,y=exp(pred),group=as.factor(Mean_T),colour=as.factor(Mean_T)))+geom_line()
Temp_plot
Temp_plot2<-Temp_plot+geom_line(data=newdata_Temp,aes(y=exp(Liu_pred)),colour="red")
Temp_plot3<-Temp_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Temp_plot4<-Temp_plot3+ylab(expression(paste("Aboveground biomass (Mg ",ha^-1,")")))+xlab("Mean annual temperature")
Temp_plot4