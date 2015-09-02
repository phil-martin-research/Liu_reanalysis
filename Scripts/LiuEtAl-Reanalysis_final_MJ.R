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
temp<-raster("Data/bio1.bil")
precip<-raster("Data/bio12.bil")
forest<-raster("Data/Forest.tif")


#do spatial analysis to look at representativeness
#of forests used  in our study
Liu_coords<-SpatialPoints(cbind(Liu[,5],Liu[,4]))
Liu_clim<-data.frame(precip=extract(precip,Liu_coords),temp=extract(temp/10,Liu_coords),data="Our data")

#create grid with 0.5 degree resolution
Grid<-expand.grid(x=seq(-180,180,by = 0.5),y=seq(-90,90,by=0.5))
coordinates(Grid)<-c("x", "y")
gridded(Grid) <- TRUE
Forest_climate<-data.frame(precip=extract(precip,Grid),temp=extract(temp,Grid),forest=extract(forest,Grid))
Forest_climate2<-subset(Forest_climate,forest==1)
all_data<-data.frame(precip=as.numeric(Forest_climate2$precip),
                     temp=as.numeric(Forest_climate2$temp/10),
                     data="Global data")

#stick these data from sites and all forests together
Climate<-rbind(Liu_clim,all_data)
Climate<-subset(Climate,!is.na(precip)&!is.na(temp))
Climate$Climate_precip_bin<-as.numeric(cut(Climate$precip,
    breaks=(seq(min(Climate$precip),max(Climate$precip),by=200)),
    labels=seq(min(Climate$precip),max(Climate$precip),by=200)[-1]))
Climate$Climate_temp_bin<-as.numeric(cut(Climate$temp,
                        breaks=(seq(min(Climate$temp),max(Climate$temp),by=1)),
                        labels=seq(min(Climate$temp),max(Climate$temp),by=1)[-1]))


ddply(Climate,.(data),summarise,
      Temp_perc=(Climate_temp_bin/length(Climate_temp_bin))*100,
      Clim_perc=(Climate_precip_bin/length(Climate_precip_bin))*100)


#find convex hull for these data
find_hull <- function(Climate) Climate[chull(Climate$precip,Climate$temp), ]
hulls <- ddply(Climate, "data", find_hull)

plot <- ggplot(data = Climate, aes(x = precip, y = temp, colour=data)) +
   geom_point(alpha=0.5)+geom_polygon(data=hulls,fill=NA)+
  labs(x = "Precipitation", y = "Temperature")+facet_wrap(~data)
plot




bb <- extent(-180, 180, -90, 90)
precip2 <- setExtent(precip, bb, keepres=TRUE)
forest2 <- setExtent(forest, bb, keepres=TRUE)

plot(precip2)
plot(forest)

Precip_mask<-mask(precip2,forest2)

#look at biases in the data that may influence results

#age
theme_set(theme_bw(base_size=12))
Geom_hist<-ggplot(Liu,aes(x=Age))+geom_histogram()+ylab("number of sites")+xlab("Estimate age (Years)")+
theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA))+geom_hline(y=0)


#location
theme_set(theme_bw(base_size=12))
world_map <- map_data("world")#Get world map info
p <- ggplot() + coord_fixed()#Create a base plot
base_world <- p + geom_polygon(data=world_map,aes(x=long,y=lat,group=group),fill="light grey")#Add map to base plot
Location<-base_world + geom_point(data=Liu,aes(x=Long,y=Lat,size=Age),colour="black",alpha=0.2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())




# ---------------------------------------------------- #
# Spatial autocorrelation of Liu et al. original data

ncf.cor <- correlog(Liu$Long, Liu$Lat, Liu$Age,increment=50, resamp=50,latlon = T)
qplot(x=ncf.cor$mean.of.class,y=ncf.cor$correlation)+geom_smooth()

# So there is some autocorrelation present in smaller and larger geographic scales
# in the Age dataset

#### ReAnalysis - Nullmodels ####

# We fit a dummy variable to be used for models without random factor
Liu$dummy <- rep(1,nrow(Liu))

# First build a spatial correlation Matrix
# We change coordinates slightly since some sites have exactly the same coordinates
Liu$Lat_J <- Liu$Lat+(rnorm(length(Liu$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu)
corMatrix(cs1Exp)[1:10, 1:4] # Looks good

# Then we run some null-models to test if ether random variables or Spatial Autocorrelation 
# are appropriate on the Above-Ground-Biomass

null.model<-lme(log(AGB)~1,data=Liu,random=~1|dummy,method="ML") # Without Random Structure - equivalent to original OLS 
null.model2<-lme(log(AGB)~1,data=Liu,random=~1|Ref,method="ML") # With Random Structure
null.model3<-lme(log(AGB)~1,data=Liu,random=~1|Ref/Site,method="ML") # Hierarchical nested random Structure 
null.model4<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long),method="ML") # The same as above but including the Spatial Autocorrelation matrix
null.model5<- update(null.model3, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

# Null model comparison
(nm_out <- aictab(list(null.model,null.model2,null.model3,null.model4,null.model5),sort=T,second.ord=F,
                  modnames=c("Null - w/o Random","Null - w. Random","Null - w. nested Random","Null - w. Random + SAC","Null - w. nested Random + SAC")) )

write.csv(as.data.frame(nm_out),"Results/NullModel.Comparison.csv")
# Models that account for differences between studies and / or spatial autocorrelation
# outperform models that have no such structure (such as Liu et al. original OLS)

#### ReAnalysis - Liu original Models ####

# We start modeling and build models similar to those of liu et al (OLS) mentioned in the Appendix
# We considering everything independently to each other but using our random variable structure to account for
# spatial autocorrelation and any systematic differences amongst studies

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
write.csv(as.data.frame(liuold_out),"Results/LIUmodels_SAC.csv")

mean(liuold_out$Delta_AIC[3:6])

# Look at spatial autocorrelation in the model residuals for both models
sar.df <- data.frame()
# Precip
ncf.cor_precip<- correlog(Liu$Long, Liu$Lat, resid(Precip_model),increment=50, resamp=50,latlon = T)
ncf.cor_precip_ourmodel<- correlog(Liu$Long, Liu$Lat, resid(Precip_model_sac),increment=50, resamp=50,latlon = T)
sar.df <- rbind(sar.df,data.frame(mean.of.class=ncf.cor_precip$mean.of.class,correlation=ncf.cor_precip$correlation,type="Precipitation model (with SAR)"))
sar.df <- rbind(sar.df,data.frame(mean.of.class=ncf.cor_precip_ourmodel$mean.of.class,correlation=ncf.cor_precip_ourmodel$correlation,type="Precipitation model (corrected)") )

# Temperature
ncf.cor_temp<- correlog(Liu$Long, Liu$Lat, resid(Temp_model),increment=50, resamp=50,latlon = T)
ncf.cor_temp_ourmodel<- correlog(Liu$Long, Liu$Lat, resid(Temp_model_sac),increment=50, resamp=50,latlon = T)
sar.df <- rbind(sar.df,data.frame(mean.of.class=ncf.cor_temp$mean.of.class,correlation=ncf.cor_temp$correlation,type="Temperature model (with SAR)"))
sar.df <- rbind(sar.df,data.frame(mean.of.class=ncf.cor_temp_ourmodel$mean.of.class,correlation=ncf.cor_temp_ourmodel$correlation,type="Temperature model (corrected)") )

# Age
ncf.cor_age<- correlog(Liu$Long, Liu$Lat, resid(Age_model),increment=50, resamp=50,latlon = T)
ncf.cor_age_ourmodel<- correlog(Liu$Long, Liu$Lat, resid(Age_model_sac),increment=50, resamp=50,latlon = T)
sar.df <- rbind(sar.df,data.frame(mean.of.class=ncf.cor_age$mean.of.class,correlation=ncf.cor_age$correlation,type="Age model (with SAR)"))
sar.df <- rbind(sar.df,data.frame(mean.of.class=ncf.cor_age_ourmodel$mean.of.class,correlation=ncf.cor_age_ourmodel$correlation,type="Age model (corrected)") )

# Plotting
g <- ggplot(sar.df,aes(x=mean.of.class,y=correlation))
g <- g + geom_line() + facet_wrap(~type,nrow = 3,as.table = T)
g <- g + labs(x="Distance class",y="Correlation",title="Correcting for spatial autocorrelation in models")
ggsave("Figures/SpatialAutocorrelationOfOriginalModelResiduals.png",plot=g,height=fig_h,width=fig_w,dpi=fig_dpi,units=fig_units,scale=fig_scale)

# All models with SAC and random structure perform better than Liu et al. original models
# and reduce the spatial autocorrelation especially at larger scales!
# Furthermore due to colinearity between Predictors we need to consider interactions between
# predictors

#### ReAnalysis - Model with Interactions ####
# First we build a global model for use in model averaging that contains all varibles that are needed
# Squared and cubed terms temperature and precipitation are not included due to missing biological sense

Liu$Mean_T2<-Liu$Mean_T+17 
Liu$logAge<-log(Liu$Age)
Liu$Age_sq<-Liu$Age^2

mymodel<-lme(log(AGB)~Age*Mean_precip+Age*Mean_T2+Mean_T2*Mean_precip+Age_sq+logAge*Mean_precip+logAge*Mean_T2,
             data=Liu,
             weights= varFunc(~I(Mean_T2)), # To approach homoscedasticity
             random=~1|Ref/Site,
             correlation = corExp(1, form = ~ Lat_J + Long),
             method="ML")

# Check for heteroskedasticity
plot(mymodel, which=2) # Somewhat greater spread at higher AGB
plot(ranef(mymodel)) # Random effects seem okay

qplot(Liu$Age,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_T2,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_precip,resid(mymodel))+geom_smooth()
qqnorm(mymodel,abline = c(0, 1))
# Heteroskedasticity is present in some cases, 
# likely due to small sample sizes with high variances in extreme regions (few samples in tropics)

# Now we dredge the model so that the value of each variable in predicting biomass
# can be assessed rather than using them in isolation
# Use second-order Information Criterion and keep Age as explanatory variable
MS1 <- dredge(mymodel,evaluate=T,rank=AICc,trace=T,subset = !(Age&&logAge) && dc(Age,Age_sq) ,extra=c("R^2","adjR^2"))
poss_mod <- get.models(MS1,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
modsumm2
averaged <- model.avg(modsumm2,fit=T,subset=delta<7)

# since the model without log terms comes out best, rerun the model averaging
# routine without this the log term
mymodel2<-lme(log(AGB)~Age*Mean_precip+Age*Mean_T2+Mean_T2*Mean_precip,
              data=Liu,
              weights= varFunc(~I(Mean_T2)), # To approach homoscedasticity
              random=~1|Ref/Site,
              correlation = corExp(1, form = ~ Lat_J + Long),
              method="ML")
plot(mymodel2) 
anova(mymodel,mymodel2)
MS2 <- dredge(mymodel2,evaluate=T,rank=AICc,trace=T,REML=F)
poss_mod <- get.models(MS2,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
modsumm2
averaged <- model.avg(modsumm2,fit=T,subset=delta<7)

averaged$formula

#### ReAnalysis - Model predictions ####
# Check the value distribution per rounded data level
Liu$Agebin <- ( round_any(Liu$Age,100) )
Liu$Precipbin <-( round_any(Liu$Mean_precip,1000) )
Liu$Tempbin <- ( round_any(Liu$Mean_T,10) )

# Finally run the best performing model including interactions from the model averaging process
# Use Liu et al. original models in comparison which consider those factors in isolation,
# but control for random structure and spatial autocorrelation
top_model<-lme(averaged$formula,data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),weights= varFunc(~I(Mean_T2)),method="ML")
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(out <- aictab(list(top_model,Precip_model,Temp_model,Age_model),sort=T,second.ord=T,
               modnames=c("Model including Interactions","Precip. only","Temp. only","Age only")) )
write.csv(as.data.frame(out),"Results/FinalModelComparison.csv")

# Our Best performing model includes several interactions between all used Predictors
# outperforms all other of Liu et al.s Models even if accounted for Spatial Autocorrelation 
# and random Structure. Thus indicates that Interactions are indeed important to determine Above Ground Biomass

#subset dataset to give only sites where precipitation bin is 1000-3000
Liu_precip<-subset(Liu,Precipbin>0&Precipbin<=3000)
Liu_precip$Mean_precip<-Liu_precip$Precipbin

ddply(Liu,.(Precipbin),summarize,minp=min(Age),max=max(Age),no=length(Age))

# Figure for interaction between age, precipitation and AGB
AP_pred<-data.frame(rbind(data.frame(Age=seq(80,795,1),Mean_precip=1000,Mean_T2=mean(Liu$Mean_T2)),
               data.frame(Age=seq(80,1200,1),Mean_precip=2000,Mean_T2=mean(Liu$Mean_T2)),
               data.frame(Age=seq(80,750,1),Mean_precip=3000,Mean_T2=mean(Liu$Mean_T2))))
AP_pred$Age_sq<-AP_pred$Age^2
AP_pred$Pred<-predict(top_model,AP_pred,level=0,se.fit=T,backtransform=T)$fit 
AP_pred$UCI<-AP_pred$Pred+(predict(top_model,AP_pred,level=0,se.fit=T)$se.fit*2)
AP_pred$LCI<-AP_pred$Pred-(predict(top_model,AP_pred,level=0,se.fit=T)$se.fit*2)

# Plot predictions
theme_set(theme_bw(base_size=12))
Age_precip1 <- ggplot(AP_pred,aes(Age,exp(Pred),ymax=exp(UCI),ymin=exp(LCI),group=as.factor(Mean_precip),fill=as.factor(Mean_precip)))+geom_line()+geom_ribbon(alpha=0.2)
Age_precip2 <- Age_precip1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")+facet_wrap(~Mean_precip)
Age_precip3 <- Age_precip2+scale_fill_brewer(palette = "Set1")+geom_rug(data=Liu_precip,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL))+geom_point(data=Liu_precip,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL),shape=1,alpha=0.5)
Age_precip3 <- Age_precip3+labs(y=expression(paste("Aboveground biomass (Mg ",ha^-1,")",sep="")),
                                x="Estimated forest age")
ggsave("Figures/Age_Precip.png",plot = Age_precip3,height=fig_h,width=fig_w,dpi=fig_dpi,units=fig_units,scale=fig_scale)

# Now age and temperature
Liu_Temp <- subset(Liu,Tempbin>=0&Tempbin<=20)
Liu_Temp$Mean_T <- Liu_Temp$Tempbin

ddply(Liu,.(Tempbin),summarize,minp=min(Age),maxp=max(Age),no=length(Age))

AT_pred <- data.frame(rbind(data.frame(Age=seq(80,1200,1),Mean_precip=mean(Liu$Mean_precip),Mean_T=0),
                          data.frame(Age=seq(80,1000,1),Mean_precip=mean(Liu$Mean_precip),Mean_T=10),
                          data.frame(Age=seq(80,200,1),Mean_precip=mean(Liu$Mean_precip),Mean_T=20)))
AT_pred$Age_sq<-AT_pred$Age^2
AT_pred$Mean_T2<-AT_pred$Mean_T+17
AT_pred$Pred<-predict(top_model,AT_pred,level=0,se.fit=T)$fit
AT_pred$UCI<-AT_pred$Pred+(predict(top_model,AT_pred,level=0,se.fit=T)$se.fit*2)
AT_pred$LCI<-AT_pred$Pred-(predict(top_model,AT_pred,level=0,se.fit=T)$se.fit*2)

Temp_Age1 <- ggplot(AT_pred,aes(Age,exp(Pred),ymax=exp(UCI),ymin=exp(LCI),group=as.factor(Mean_T),fill=as.factor(Mean_T)))+geom_line()+geom_ribbon(alpha=0.5)
Temp_Age2 <- Temp_Age1+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")+facet_wrap(~Mean_T)
Temp_Age3 <- Temp_Age2+scale_fill_brewer(palette = "Set1")+geom_rug(data=Liu_Temp,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL))+geom_point(data=Liu_Temp,aes(x=Age,y=AGB,ymax=NULL,ymin=NULL,fill=NULL),shape=1,alpha=0.2)
Temp_Age3 <- Temp_Age3 + labs(y = expression(paste("Aboveground biomass (Mg ",ha^-1,")",sep="")),
                              x = "Estimated forest age")
ggsave("Figures/Age_Temp.png",plot = Temp_Age3,height=fig_h,width=fig_w,dpi=fig_dpi,units=fig_units,scale=fig_scale)

#now temperature and precipitation
Liu_Temp<-subset(Liu,Tempbin>=0&Tempbin<=20)
Liu_Temp$Mean_T <- Liu_Temp$Tempbin

ddply(Liu,.(Tempbin),summarize,minp=min(Mean_precip),maxp=max(Mean_precip),no=length(Mean_precip))

AP_pred<-data.frame(rbind(data.frame(Age=mean(Liu$Age),Mean_precip=seq(0,3000,1),Mean_T=0),
                          data.frame(Age=mean(Liu$Age),Mean_precip=seq(0,3700,1),Mean_T=10),
                          data.frame(Age=mean(Liu$Age),Mean_precip=seq(0,5800,1),Mean_T=20)))

AP_pred$Age_sq<-AP_pred$Age^2
AP_pred$Mean_T2<-AP_pred$Mean_T+17
AP_pred$Pred<-predict(top_model,AP_pred,level=0,se.fit=T)$fit
AP_pred$UCI<-AP_pred$Pred+(predict(top_model,AP_pred,level=0,se.fit=T)$se.fit*2)
AP_pred$LCI<-AP_pred$Pred-(predict(top_model,AP_pred,level=0,se.fit=T)$se.fit*2)

#now plot this
theme_set(theme_bw(base_size=12))
Temp_precip1 <- ggplot(AP_pred,aes(Mean_precip,exp(Pred),ymax=exp(UCI),ymin=exp(LCI),group=as.factor(Mean_T),fill=as.factor(Mean_T)))+geom_line()+geom_ribbon(alpha=0.5)
Temp_precip2 <- Temp_precip1 + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_rect(size=1.5,colour="black",fill=NA))+ theme(legend.position="none")+facet_wrap(~Mean_T)
Temp_precip3 <- Temp_precip2 + scale_fill_brewer(palette = "Set1")+geom_rug(data=Liu_Temp,aes(x=Mean_precip,y=AGB,ymax=NULL,ymin=NULL,fill=NULL))+geom_point(data=Liu_Temp,aes(x=Mean_precip,y=AGB,ymax=NULL,ymin=NULL,fill=NULL),shape=1,alpha=0.2)
Temp_precip3 <- Temp_precip3 + labs(y=expression(paste("Aboveground biomass (Mg ",ha^-1,")",sep="")),
                                   x= "Mean annual precipitation (mm)")
ggsave("Figures/Temp_Precip.png",plot=Temp_precip3,height=fig_h,width=fig_w,dpi=fig_dpi,units=fig_units,scale=fig_scale)

# Spatial look at the residuals
#r <- residuals(top_model)
#base_world+geom_point(data=Liu,aes(x=Long,y=Lat,size=sqrt(r^2)),color="blue")
