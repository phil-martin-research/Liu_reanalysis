#### Preparation and standard Package Loading ####
rm(list=ls())
library(ggplot2)
library(nlme)
library(MuMIn)
library(mgcv)
library(AICcmodavg)
library(scales)
library(lattice)
library(ncf)
library(plyr)
library(raster)
library(sp)
library(maptools)

#### Data loading and format ####
# Load in the supplementary data from LIU et al. DOI: 10.1111/geb.12113
Liu <- read.csv("Data/Liu_aged.csv")

# Null-models to test if ether random variables or Spatial Autocorrelation 
# are appropriate on the Above-Ground-Biomass

null.model<-lme(log(AGB)~1,data=Liu,random=~1|dummy,method="ML") # Without Random Structure - equivalent to original OLS 
null.model2<-lme(log(AGB)~1,data=Liu,random=~1|Ref,method="ML") # With Random Structure
null.model3<-lme(log(AGB)~1,data=Liu,random=~1|Ref/Site,method="ML") # Hierarchical nested random Structure 
null.model4<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long),method="ML") # The same as above but including the Spatial Autocorrelation matrix
null.model5<- update(null.model3, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

# Null model comparison
(nm_out <- aictab(list(null.model,null.model2,null.model3,null.model4,null.model5),sort=T,second.ord=F,
                  modnames=c("Null - w/o Random","Null - w. Random","Null - w. nested Random","Null - w. Random + SAC","Null - w. nested Random + SAC")) )
# Models that account for differences between studies and / or spatial autocorrelation
# outperform models that have no such structure (such as Liu et al. original OLS)

#### Analysis - Model with Interactions ####
# First we build a global model for use in model averaging that contains all varibles that are needed
# Squared and cubed terms temperature and precipitation are not included due to missing biological sense

Liu$Mean_T2<-(Liu$Mean_T-mean(Liu$Mean_T))/sd(Liu$Mean_T)
Liu$Precip2<-(Liu$Mean_precip-mean(Liu$Mean_precip))/sd(Liu$Mean_precip)
Liu$Age2<-(Liu$Age-mean(Liu$Age))/sd(Liu$Age)
Liu$Age2sq<-Liu$Age2^2

mymodel<-lme(log(AGB)~Mean_T2*Precip2+Age2*Precip2+Age*Mean_T2,
             data=Liu,
             random=~1|Ref/Block,
             correlation = corExp(1, form = ~ Lat_J + Long),
             method="ML")

plot(mymodel)


# Check for heteroskedasticity
plot(mymodel, which=2) # Somewhat greater spread at higher AGB

qplot(Liu$Age,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_T2,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_precip,resid(mymodel))+geom_smooth()
qqnorm(mymodel,abline = c(0, 1))
# Heteroskedasticity is present in some cases, 
# likely due to small sample sizes with high variances in extreme regions (few samples in tropics)

# Now we dredge the model so that the value of each variable in predicting biomass
# can be assessed rather than using them in isolation
# Use second-order Information Criterion and keep Age as explanatory variable
MS1 <- dredge(mymodel,evaluate=T,rank=AICc,trace=T,extra=c("R^2","adjR^2"))
poss_mod <- get.models(MS1,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
averaged <- model.avg(modsumm2,fit=T,subset=delta<7)
summary(averaged)

#### ReAnalysis - Model predictions ####
# Check the value distribution per rounded data level
Liu$Agebin <- ( round_any(Liu$Age,100) )
Liu$Precipbin <-( round_any(Liu$Mean_precip,1000) )
Liu$Tempbin <- ( round_any(Liu$Mean_T,10) )

# Finally run the best performing model including interactions from the model averaging process
# Use Liu et al. original models in comparison which consider those factors in isolation,
# but control for random structure and spatial autocorrelation
top_model<-lme(averaged$formula,data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),weights= varFunc(~I(Mean_T2)),method="ML")
Precip_model<-lme(AGB~Mean_precip+I(Mean_precip^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Temp_model<-lme(AGB~Mean_T+I(Mean_T^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Age_model<-lme(AGB~Age+I(Age^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(out <- aictab(list(top_model,Precip_model,Temp_model,Age_model),sort=T,second.ord=T,
               modnames=c("Model including Interactions","Precip. only","Temp. only","Age only")) )
write.csv(as.data.frame(out),"Results/FinalModelComparison.csv")

# Our Best performing model includes several interactions between all used Predictors
# outperforms all other of Liu et al.s Models even if accounted for Spatial Autocorrelation 
# and random Structure. Thus indicates that Interactions are indeed important to determine Above Ground Biomass

#subset dataset to give only sites where precipitation bin is 1000-3000


