# this is a script to reanalyse the 
# data from Lio et al 2013 paper in GEB

library(ggplot2)
library(nlme)
library(GGally)
library(ape)
library(cshapes)
library(MuMIn)
library(mgcv)
library(geoR)


#read in data
setwd("C:/Users/Phil/Dropbox/Work/PhD/Publications, Reports and Responsibilities/Publications/Liu_et_al/Liu_reanalysis/Data")
Liu<-read.csv("Liu_aged.csv")
head(Liu)

#subset data to only include sites with ages
Liu_sub<-subset(Liu,Age>0)

Liu_sub$Biome<-ifelse(abs(Liu_sub$Lat)<23.5,"Tropical",NA)
Liu_sub$Biome<-ifelse(abs(Liu_sub$Lat)>23.5&abs(Liu_sub$Lat)<50,"Temperate",Liu_sub$Biome)
Liu_sub$Biome<-ifelse(abs(Liu_sub$Lat)>50,"Boreal",Liu_sub$Biome)

ggplot(data=Liu_sub,aes(y=Age,x=Biome,colour=Biome))+geom_boxplot(size=1)+ coord_trans(y = "log10")+scale_y_continuous(breaks=c(100,200,300,400,500,600,700,800,1200))

head(Liu)

hist(Liu$Age)
hist(Liu$Mean_T)
hist(Liu$Mean_precip)

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


#produce variogram to examine spatial autocorrelation
head(Liu)
dists<-dist(cbind(Liu[,5],Liu[,17]))
summary(dists)
breaks<-seq(0, 325,by = 5)
v1 <- variog(coords = cbind(Liu[,5],Liu[,17]), data = Liu[,8], breaks = breaks)


plot(v1, type="b")

v1.summary <- cbind(v1$v, v1$n)
colnames(v1.summary) <- c("lag", "semi-variance", "# of pairs")

v1.summary
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
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2)+I(Mean_T^3),data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#Age
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

#now a global model for use in model averaging
#that contains all varibles that are needed
#I haven't included the squared and cubed terms for temp
#and precipitation becuase I don't think they make biological
#sense
All_model<-lme(log(AGB)~Age*Mean_precip*Mean_T+Age_sq,data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))

plot(All_model)

#now we can dredge the model so that the value of each variable in predicting biomass
#can be assessed rather than using them in isolation
MS1<-dredge(All_model,evaluate=T,rank=AICc,trace=T,subset=dc(Age,Age_sq),REML=F)
poss_mod<-get.models(MS1,subset=delta<7)
modsumm<- model.sel(poss_mod, rank = "AICc",fit=T)
modsumm2<-subset(modsumm,modsumm$delta<7)
modsumm2
averaged<-model.avg(modsumm2,fit=T)

#run the top model from the model averaging to get
#r squared statistic
top_model<-lme(log(AGB)~Age+Mean_precip*Mean_T+Age_sq,data=Liu,random=~1|Ref,correlation = corExp(1, form = ~ Lat_J + Long))


#create table of different model AICc, AICc delta and marginal R squared
AICc_res<-AICc(top_model,Precip_model,Temp_model,Age_model)
AICc_table<-data.frame(Model=row.names(AICc_res),AICc=AICc_res$AICc)
AICc_table$delta<-AICc_table$AICc-AICc_table$AICc[1]
AICc_table$R_squared<-c(r.squaredGLMM(top_model)[1],r.squaredGLMM(Precip_model)[1],r.squaredGLMM(Temp_model)[1],r.squaredGLMM(Age_model)[1])
AICc_table

setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Publications/Liu_et_al/Liu_reanalysis/Results")
write.csv(AICc_table,"Model_comp.csv",row.names=F)


#create dataframe for predictions so paramaters can be plotted
#first - Age
nseq <- function(x, len = length(x)) seq(min(x, na.rm = TRUE),
    max(x, na.rm=TRUE), length = len)
str(Liu_subset)
newdata_Age1<- as.data.frame(lapply(lapply(Liu[c(6,7,11,15)], mean), rep, 300))
newdata_Age1$Age<- nseq(Liu$Age, nrow(newdata_Age1))
newdata_Age1$Age_sq<- newdata_Age$Age^2
newdata_Age1$Mean_T<-0
newdata_Age1$pred<-predict(top_model,newdata_Age,level=0)
newdata_Age2<- as.data.frame(lapply(lapply(Liu[c(6,7,11,15)], mean), rep, 300))
newdata_Age2$Age<- nseq(Liu$Age, nrow(newdata_Age1))
newdata_Age2$Age_sq<- newdata_Age$Age^2
newdata_Age2$Mean_T<-5
newdata_Age2$pred<-predict(top_model,newdata_Age,level=0)
newdata_Age3<- as.data.frame(lapply(lapply(Liu[c(6,7,11,15)], mean), rep, 300))
newdata_Age3$Age<- nseq(Liu$Age, nrow(newdata_Age1))
newdata_Age3$Age_sq<- newdata_Age$Age^2
newdata_Age3$Mean_T<-10
newdata_Age3$pred<-predict(top_model,newdata_Age,level=0)

newdata_Age<-rbind(newdata_Age1,newdata_Age2,newdata_Age3)

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
setwd("C:/Users/Phil/Documents/My Dropbox/Work/PhD/Publications, Reports and Responsibilities/Publications/Liu_et_al/Liu_reanalysis/Figures")
theme_set(theme_bw(base_size=12))
Age_plot<-ggplot(data=newdata_Age,aes(x=Age,y=exp(pred),group=as.factor(Mean_T),colour=as.factor(Mean_T)))+geom_line(size=2)
Age_plot2<-Age_plot+geom_line(data=newdata_Age,aes(y=exp(Liu_pred),group=NULL),colour="black",lty=2,size=2)
Age_plot3<-Age_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Age_plot4<-Age_plot3+ylab(expression(paste("Aboveground \nbiomass (Mg ",ha^-1,")")))+xlab("Age (Years)")
Age_plot5<-Age_plot4+scale_colour_brewer(expression("Mean annual\ntemperature("*degree*C*")"),palette="Set1")
Age_plot5
ggsave(filename="Age_temp.png",width=9,height=4,units="in",dpi=400)

#now make predictions for Temp
plot(Liu_subset$Age,Liu_subset$Mean_T)
newdata_Temp1<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Temp1$Mean_T<- nseq(Liu_subset$Mean_T, nrow(newdata_Temp1))
newdata_Temp1$Age<-100
newdata_Temp1$Age_sq<-newdata_Temp1$Age^2
newdata_Temp2<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Temp2$Mean_T<- seq(from=-5,to=15,length.out=nrow(newdata_Temp3))
newdata_Temp2$Age<-150
newdata_Temp2$Age_sq<-newdata_Temp2$Age^2
newdata_Temp3<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Temp3$Mean_T<- seq(from=-5,to=5, length.out=nrow(newdata_Temp3))
newdata_Temp3$Age<-200
newdata_Temp3$Age_sq<-newdata_Temp3$Age^2
newdata_Temp<-rbind(newdata_Temp1,newdata_Temp2,newdata_Temp3)


#now use predict
newdata_Temp$pred<-predict(top_model,newdata_Temp,level=0)
#add predictions from model which is most similar to that of Liu et al
newdata_Temp$Liu_pred<-predict(Temp_model,newdata_Temp,level=0)

summary(newdata_Temp)

plot(newdata_Temp$Mean_T,exp(newdata_Temp$pred))

#now plot the predictions for Temp
theme_set(theme_bw(base_size=10))
Temp_plot<-ggplot(data=newdata_Temp,aes(x=Mean_T,y=exp(pred),group=as.factor(Age),colour=as.factor(Age)))+geom_line(size=2)
Temp_plot2<-Temp_plot+geom_line(data=newdata_Temp,aes(y=exp(Liu_pred)),colour="black",lty=2,size=2)
Temp_plot3<-Temp_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Temp_plot4<-Temp_plot3+ylab(expression(paste("Aboveground \nbiomass (Mg ",ha^-1,")")))+xlab("Mean annual temperature")
Temp_plot4+scale_colour_brewer("Age \n(Years)",palette="Set1")
ggsave(filename="Temp_age.png",width=9,height=4,units="in",dpi=400)

#now predictions for changes in precipitation
newdata_P<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_P$Mean_precip<- nseq(Liu_subset$Mean_precip, nrow(newdata_P1))

#now use predict
newdata_P$pred<-predict(top_model,newdata_P,level=0)

#add predictions from model which is most similar to that of Liu et al
newdata_P$Liu_pred<-predict(Precip_model,newdata_P,level=0)

#now plot
theme_set(theme_bw(base_size=10))
Temp_plot<-ggplot(data=newdata_P,aes(x=Mean_precip,y=exp(pred)))+geom_line(size=2)
Temp_plot
Temp_plot2<-Temp_plot+geom_line(data=newdata_P,aes(y=exp(Liu_pred)),colour="black",lty=2,size=2)
Temp_plot2
Temp_plot3<-Temp_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Temp_plot4<-Temp_plot3+ylab(expression(paste("Aboveground\n biomass (Mg ",ha^-1,")")))+xlab("Mean annual precipitation (mm)")
Temp_plot4
ggsave(filename="Precipitation.png",width=9,height=4,units="in",dpi=400)
