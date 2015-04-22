#### Preparation and standard Package Loading ####
rm(list=ls())
library(gdata)
library(ggplot2)
library(nlme)
library(MuMIn)
library(mgcv)
library(AICcmodavg)
library(scales)
library(lattice)
library(ncf)
library(plyr)



#### Data loading formation ####
# Load in the supplementary data from LIU et al. DOI: 10.1111/geb.12113
Liu <- read.csv("Data/geb12113-sup-0001-ts1.csv")

# Format and cleanup Data for the reanalysis
Liu_sub<-Liu[-c(2:4,6,18:22)] # Get Columns of Interest
colnames(Liu_sub)<-c("ID","Site","Lat","Long","Mean_T","Mean_precip","AGB","L_AGB","T_AGB","Age","A_L_Ratio","A_T_Ratio","Ref")
Liu_sub<-subset(Liu_sub,AGB>0&Age>0) # Only Site with given Age and AGB
levels(Liu_sub$Site)<-seq(1:449) # Correct Site Level
Liu_sub$Age_sq<-Liu_sub$Age^2
Liu <- Liu_sub
rm(Liu_sub)

# Investigate the Structure
hist(Liu$Age)
hist(Liu$Mean_T)
hist(Liu$Mean_precip)

# Histograms indicate a general paucity at Age >500 years, Temp <-5 or >20, or precip >3000

# We fit a dummy variable to be used for models without random factor
Liu$dummy <- rep(1,nrow(Liu))

#before analyses some data exploration
#first look at it spatially
world_map <- map_data("world")#Get world map info
p <- ggplot() + coord_fixed()#Create a base plot
base_world <- p + geom_polygon(data=world_map,aes(x=long,y=lat,group=group))#Add map to base plot
base_world+geom_point(data=Liu,aes(x=Long,y=Lat,colour=Ref),alpha=0.2)+facet_wrap(~Ref)
#references may focus on particular areas of the globe but the big ones (Luo, 1996; Ma, 2012; Luyassaert et al 2007)
#come from quite a spread of different locations so I'm not super sure we need to include random effects
ncf.cor <- correlog(Liu$Long, Liu$Lat, Liu$AGB,increment=500, resamp=500,latlon = T)
qplot(x=ncf.cor$mean.of.class,y=ncf.cor$correlation)+geom_line()
#not sure what's going on with the very high Moran's I values but this looks
#like a classic case of sites closer to each other being more similar in AGB

ggplot(Liu,aes(x=Age,y=AGB))+geom_point()+geom_smooth()#very high AGB of Australian forest at ~500yrs
ggplot(Liu,aes(x=Mean_T,y=AGB))+geom_point()+geom_smooth()#high AGB again at ~10 degrees C
ggplot(Liu,aes(x=Mean_precip,y=AGB))+geom_point()+geom_smooth()#positive trend but litte data >2000mm

#### ReAnalysis - Nullmodels####

# We start modeling and build models similar to those of liu et al (OLS) mentioned in the Appendix
# We considering everything independently to each other but using our random variable structure to account for
# spatial autocorrelation and any systematic differences amongst studies

# First build a spatial correlation Matrix
# We change coordinates slightly since some sites have exactly the same coordinates
Liu$Lat_J <- Liu$Lat+(rnorm(length(Liu$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu)
corMatrix(cs1Exp)[1:10, 1:4] # Looks good

# Then we run some null-models to test if ether random variables or Spatial Autocorrelation 
# are appropriate on the Above-Ground-Biomass

null.model<-lme(log(AGB)~1,data=Liu,random=~1|dummy,method="ML") # Without Random Strucutre
null.model2<-lme(log(AGB)~1,data=Liu,random=~1|Ref,method="ML") # With Random Structure
null.model3<-lme(log(AGB)~1,data=Liu,random=~1|Ref/Site,method="ML") # Hierarchical nested random Structure 
null.model4<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long),method="ML") # The same as above but including the Spatial Autocorrelation matrix
null.model5<- update(null.model3, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(nm_out <- aictab(list(null.model,null.model2,null.model3,null.model4,null.model5),sort=T,second.ord=F,
                  modnames=c("Null - w/o Random","Null - w. Random","Null - w. nested Random","Null - w. Random + SAC","Null - w. nested Random + SAC")) )



# Fitting a model that accounts for between study differences is better
# and we NEED to account for spatial autocorrelation for our results to
# be statistically valid. Models that account for Spatial Autocorrelation and random Structure perform best

#### ReAnalysis - Liu original Models ####
# The Following part is to show that models with included random factor and taking SAC
# into account the models outperform models without

# Precipitation - a model with a squared term for mean_precip
Precip_model<-lme(AGB~Mean_precip+I(Mean_precip^2),data=Liu,random=~1|dummy,method="ML")
Precip_model_sac<-lme(AGB~Mean_precip+I(Mean_precip^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

# Temperature
Temp_model<-lme(AGB~Mean_T+I(Mean_T^2),data=Liu,random=~1|dummy,method="ML")
Temp_model_sac<-lme(AGB~Mean_T+I(Mean_T^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

# Age
Age_model<-lme(AGB~Age+I(Age^2),data=Liu,random=~1|dummy,method="ML")
Age_model_sac<-lme(AGB~Age+I(Age^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(liuold_out <- aictab(list(Precip_model,Precip_model_sac,Temp_model,Temp_model_sac,Age_model,Age_model_sac),sort=T,second.ord=F,
                      modnames=c("AGB - Precip","AGB - Precip + SAC","AGB ~ Temp","AGB - Temp + SAC","AGB - Age","AGB - Age + SAC")) )

mean(liuold_out$Delta_AIC[3:6])

#look at spatial autocorrelation from Liu models

ncf.cor_precip<- correlog(Liu$Long, Liu$Lat, resid(Precip_model),increment=500, resamp=500,latlon = T)
ncf.cor_precip_ourmodel<- correlog(Liu$Long, Liu$Lat, resid(Precip_model_sac),increment=500, resamp=500,latlon = T)

qplot(x=ncf.cor_precip$mean.of.class,y=ncf.cor_precip$correlation)+geom_line()
qplot(x=ncf.cor_precip_ourmodel$mean.of.class,y=ncf.cor_precip_ourmodel$correlation)+geom_line()

ncf.cor_temp<- correlog(Liu$Long, Liu$Lat, resid(Temp_model),increment=500, resamp=500,latlon = T)
ncf.cor_temp_ourmodel<- correlog(Liu$Long, Liu$Lat, resid(Temp_model_sac),increment=500, resamp=500,latlon = T)
qplot(x=ncf.cor_temp$mean.of.class,y=ncf.cor_temp$correlation)+geom_line()
qplot(x=ncf.cor_temp_ourmodel$mean.of.class,y=ncf.cor_temp_ourmodel$correlation)+geom_line()

ncf.cor_age<- correlog(Liu$Long, Liu$Lat, resid(Age_model),increment=500, resamp=500,latlon = T)
ncf.cor_age_ourmodel<- correlog(Liu$Long, Liu$Lat, resid(Age_model_sac),increment=500, resamp=500,latlon = T)
qplot(x=ncf.cor_age$mean.of.class,y=ncf.cor_age$correlation)+geom_line()
qplot(x=ncf.cor_age_ourmodel$mean.of.class,y=ncf.cor_age_ourmodel$correlation)+geom_line()


# All models with SAC and random structure perform better than Liu et al. original models.
# Furthermore due to colinearity between Predictors we need to consider interactions between
# predictors

#### ReAnalysis - Model with Interactions ####
# First we build a global model for use in model averaging that contains all varibles that are needed
# Squared and cubed terms temperature and precipitation are not included due to missing biological sense

vf1<-varFunc(~Mean_T)
Liu$Mean_T2<-Liu$Mean_T+17
Liu$logAge<-log(Liu$Age)
Liu$Age_sq<-Liu$Age^2

mymodel<-lme(AGB~Age*Mean_precip+Age*Mean_T2+Mean_T2*Mean_precip+Age_sq+logAge*Mean_precip+logAge*Mean_T2,data=Liu,random=~1|Ref/Site,weights=varFixed(~Mean_T2),correlation = corExp(1, form = ~ Lat_J + Long),method="ML")


qplot(Liu$Age,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_T2,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_precip,resid(mymodel))+geom_smooth()
qplot(x=fitted(mymodel),y=resid(mymodel))+geom_smooth()
qqnorm(mymodel,abline = c(0, 1))
#generally all of this seems to be ok, residuals look fine


# Now we dredge the model so that the value of each variable in predicting biomass
# can be assessed rather than using them in isolation
# Use second-order Information Criterion and keep Age as explanatory variable
MS1 <- dredge(mymodel,evaluate=T,rank=AICc,trace=T,REML=F,subset = !(Age && logAge)&&dc(Age,Age_sq))
poss_mod <- get.models(MS1,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
modsumm2
averaged <- model.avg(modsumm2,fit=T,subset=delta<7)

#since the model without log terms comes out best, rerun the model averaging
#routine without this the log term
mymodel2<-lme(log(AGB)~Age*Mean_precip+Age*Mean_T2+Mean_T2*Mean_precip,data=Liu,random=~1|Ref/Site,weights=varFixed(~Mean_T2),correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
plot(mymodel2)

MS1 <- dredge(mymodel2,evaluate=T,rank=AICc,trace=T,REML=F)
poss_mod <- get.models(MS1,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
modsumm2
averaged <- model.avg(modsumm2,fit=T,subset=delta<7)

averaged$av

# Finally run the best performing model including interactions from the model averaging process
# Use Liu et al. original models in comparison which consider those factors in isolation,
# but control for random structure and spatial autocorrelation
top_model<-lme(averaged$formula,data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")



summary(Liu$Mean_precip)

library(plyr)

Liu$Agebin<-round_any(Liu$Age,100)
Liu$Precipbin<-round_any(Liu$Mean_precip,1000)
Liu$Tempbin<-round_any(Liu$Mean_T,10)

ddply(Liu,.(Agebin),summarize,minp=min(Mean_precip),maxp=max(Mean_precip),no=length(Mean_precip))
ddply(Liu,.(Precipbin),summarize,minp=min(Age),max=max(Age),no=length(Age))
ddply(Liu,.(Tempbin),summarize,minp=min(Age),maxp=max(Age),no=length(Age))


#figure for interaction between age, precipitation and AGB
AP_pred<-data.frame(rbind(data.frame(Age=seq(80,795,1),Mean_precip=1000,Mean_T2=mean(Liu$Mean_T2)),
               data.frame(Age=seq(87,1200,1),Mean_precip=2000,Mean_T2=mean(Liu$Mean_T2)),
               data.frame(Age=seq(85,750,1),Mean_precip=3000,Mean_T2=mean(Liu$Mean_T2))))
AP_pred$Age_sq<-AP_pred$Age^2
AP_pred$Pred<-predict(averaged,AP_pred,level=0,se.fit=T)$fit
AP_pred$UCI<-AP_pred$Pred+(predict(averaged,AP_pred,level=0,se.fit=T)$se.fit*2)
AP_pred$LCI<-AP_pred$Pred-(predict(averaged,AP_pred,level=0,se.fit=T)$se.fit*2)

#subset dataset to give only sites where precipitation bin is 1000-3000
Liu_precip<-subset(Liu,Precipbin>0&Precipbin<=3000)
Liu_precip$Mean_precip<-Liu_precip$Precipbin

#plot predictions
theme_set(theme_bw(base_size=12))
Age_precip1<-ggplot(AP_pred,aes(Age,exp(Pred),ymax=exp(UCI),ymin=exp(LCI),group=as.factor(Mean_precip),fill=as.factor(Mean_precip)))+geom_line()+geom_ribbon(alpha=0.2)
Age_precip2<-Age_precip1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")+facet_wrap(~Mean_precip)
Age_precip3<-Age_precip2+scale_fill_brewer(palette = "Set1")+geom_rug(data=Liu_precip,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL))+geom_point(data=Liu_precip,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL),shape=1,alpha=0.2)
Age_precip3+ylab(expression(paste("Aboveground biomass (Mg ",ha^-1,")",sep="")))+xlab("Estimated forest age")
ggsave("Figures/Age_Precip.png",height=5,width=8,dpi=800,units="in")

#now age and temperature
Liu_Temp<-subset(Liu,Tempbin>=0&Tempbin<=20)
Liu_Temp$Mean_T<-Liu_Temp$Tempbin

AT_pred<-data.frame(rbind(data.frame(Age=seq(80,290,1),Mean_precip=mean(Liu$Mean_precip),Mean_T=0),
                          data.frame(Age=seq(80,750,1),Mean_precip=mean(Liu$Mean_precip),Mean_T=10),
                          data.frame(Age=seq(80,200,1),Mean_precip=mean(Liu$Mean_precip),Mean_T=20)))
AT_pred$Age_sq<-AT_pred$Age^2
AT_pred$Mean_T2<-AT_pred$Mean_T+17
AT_pred$Pred<-predict(averaged,AT_pred,level=0,se.fit=T)$fit
AT_pred$UCI<-AT_pred$Pred+(predict(averaged,AT_pred,level=0,se.fit=T)$se.fit*2)
AT_pred$LCI<-AT_pred$Pred-(predict(averaged,AT_pred,level=0,se.fit=T)$se.fit*2)

Temp_Age1<-ggplot(AT_pred,aes(Age,exp(Pred),ymax=exp(UCI),ymin=exp(LCI),group=as.factor(Mean_T),fill=as.factor(Mean_T)))+geom_line()+geom_ribbon(alpha=0.5)
Temp_Age2<-Temp_Age1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")+facet_wrap(~Mean_T)
Temp_Age3<-Temp_Age2+scale_fill_brewer(palette = "Set1")+geom_rug(data=Liu_Temp,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL))+geom_point(data=Liu_Temp,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL),shape=1,alpha=0.2)
Temp_Age3+ylab(expression(paste("Aboveground biomass (Mg ",ha^-1,")",sep="")))+xlab("Estimated forest age")
ggsave("Figures/Age_Temp.png",height=3,width=8,dpi=800,units="in")

#now temperature and precipitation
Liu_Temp<-subset(Liu,Tempbin>=10&Tempbin<=20)
Liu_Temp$Mean_T<-Liu_Temp$Tempbin

ddply(Liu,.(Tempbin),summarize,minp=min(Mean_precip),maxp=max(Mean_precip),no=length(Mean_precip))

AP_pred<-data.frame(rbind(data.frame(Age=mean(Liu$Age),Mean_precip=seq(352,1937,1),Mean_T=0),
                          data.frame(Age=mean(Liu$Age),Mean_precip=seq(355,3669),Mean_T=10),
                          data.frame(Age=mean(Liu$Age),Mean_precip=seq(707,2500),Mean_T=20)))
AP_pred$Age_sq<-AP_pred$Age^2
AP_pred$Mean_T2<-AP_pred$Mean_T+17
AP_pred$Pred<-predict(averaged,AP_pred,level=0,se.fit=T)$fit
AP_pred$UCI<-AP_pred$Pred+(predict(averaged,AP_pred,level=0,se.fit=T)$se.fit*2)
AP_pred$LCI<-AP_pred$Pred-(predict(averaged,AP_pred,level=0,se.fit=T)$se.fit*2)

#now plot this
theme_set(theme_bw(base_size=12))
Temp_precip1<-ggplot(AP_pred,aes(Mean_precip,exp(Pred),ymax=exp(UCI),ymin=exp(LCI),group=as.factor(Mean_T),fill=as.factor(Mean_T)))+geom_line()+geom_ribbon(alpha=0.5)
Temp_precip2<-Temp_precip1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")+facet_wrap(~Mean_T)
Temp_precip3<-Temp_precip2+scale_fill_brewer(palette = "Set1")+geom_rug(data=Liu_Temp,aes(x=Mean_precip,y=AGB,ymax=NULL,ymin=NULL,fill=NULL))+geom_point(data=Liu_Temp,aes(x=Mean_precip,y=AGB,ymax=NULL,ymin=NULL,fill=NULL),shape=1,alpha=0.2)+xlim(0,4000)
Temp_precip3+ylab(expression(paste("Aboveground biomass (Mg ",ha^-1,")",sep="")))+xlab("Mean annual precipitation (mm)")
ggsave("Figures/Temp_Precip.png",height=5,width=8,dpi=800,units="in")


ggplot()

require(maps)


head(Liu_subset)

base_world+geom_point(data=Liu_subset,aes(x=Long,y=Lat,size=sqrt(Resid^2),colour=Resid))

(out <- aictab(list(top_model,Precip_model,Temp_model,Age_model),sort=T,second.ord=T,
               modnames=c("Model including Interactions","Precip. only","Temp. only","Age only")) )

out$adj.r.squared <- unlist(lapply(list(top_model,Temp_model,Age_model,Precip_model),FUN=function(x)attr(r.squaredLR(x,null=null.model5),"adj.r.squared")) )
(out <- as.data.frame(out))

# Our Best performing model includes several interactions between all used Predictors
# outperforms all other of Liu et al.s Models even if accounted for Spatial Autocorrelation 
# and random Structure.
# Thus indicates that Interactions are indeed important to determine Above Ground Biomass
