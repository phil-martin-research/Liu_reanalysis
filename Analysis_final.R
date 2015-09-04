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

#### Output Figure extent parameters ####
fig_h = 5
fig_w = 8
fig_dpi = 400
fig_units = "in"
fig_scale = 1.2

#### Data loading and format ####
# Load in the supplementary data from LIU et al. DOI: 10.1111/geb.12113
# and save a comma seperated file
Liu <- read.csv("Data/geb12113-sup-0001-ts1.csv")

# Format and cleanup Data for the reanalysis
Liu_sub<-Liu[-c(2:4,6,18:22)] # Get Columns of Interest
colnames(Liu_sub)<-c("ID","Site","Lat","Long","Mean_T","Mean_precip","AGB","L_AGB","T_AGB","Age","A_L_Ratio","A_T_Ratio","Ref")
Liu_sub<-subset(Liu_sub,AGB>0&Age>0) # Only Site with given Age and AGB
levels(Liu_sub$Site)<-seq(1:449) # Correct Site Level
Liu_sub$Age_sq<-Liu_sub$Age^2 # squared age
Liu <- Liu_sub
rm(Liu_sub)

#### Explanatory analysis ####
# Investigate the Structure
hist(Liu$Age)
hist(Liu$Mean_T)
hist(Liu$Mean_precip)

ggplot(Liu,aes(x=Age,y=AGB))+geom_point()+geom_smooth() # very high AGB of Australian forest at ~500yrs
ggplot(Liu,aes(x=Mean_T,y=AGB))+geom_point()+geom_smooth()# high AGB again at ~10 degrees C
ggplot(Liu,aes(x=Mean_precip,y=AGB))+geom_point()+geom_smooth()# positive trend but litte data >2000mm

# Histograms indicate a general paucity at Age >500 years, Temp <-5 or >20, or precip >3000

# Look at it spatially
world_map <- map_data("world")#Get world map info
p <- ggplot() + coord_fixed()#Create a base plot
base_world <- p + geom_polygon(data=world_map,aes(x=long,y=lat,group=group))#Add map to base plot
base_world + geom_point(data=Liu,aes(x=Long,y=Lat,colour=Ref),alpha=0.5)+facet_wrap(~Ref,5)
ggsave("Figures/DataDistribution.png",height=fig_h,width=fig_w,dpi=fig_dpi,units=fig_units,scale=fig_scale)

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

#### ReAnalysis without Interactions ####

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

mean(liuold_out$Delta_AIC[3:6]) # average delta AIC between the models

#### Spatial Autocorrelation - Correlog ####
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

# Model diagnostics
# Check for heteroskedasticity
plot(mymodel, which=2) # Somewhat greater spread at higher AGB
plot(ranef(mymodel)) # Random effects seem okay

qplot(Liu$Age,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_T2,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_precip,resid(mymodel))+geom_smooth()
qqnorm(mymodel,abline = c(0, 1))
# The variance function decreased Heteroskedasticity tremendously, 

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

MS2 <- dredge(mymodel2,evaluate=T,rank=AICc,trace=T,REML=F)
poss_mod <- get.models(MS2,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
modsumm2
averaged <- model.avg(modsumm2,fit=T,subset=delta<7)

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
