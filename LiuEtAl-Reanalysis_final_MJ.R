#### Preperation and standard Package Loading ####
rm(list=ls())
library(gdata)
library(ggplot2)
library(nlme)
library(GGally)
library(ape)
library(MuMIn)
library(mgcv)
library(influence.ME)
library(AICcmodavg)
library(scales)

#### Data loading formation ####
# Load in the supplementary data from LIU et al. DOI: 10.1111/geb.12113
si <- "http://onlinelibrary.wiley.com/store/10.1111/geb.12113/asset/supinfo/geb12113-sup-0001-ts1.xlsx?v=1&s=5429cca13645c52cec45abf5e842af1e69d9d740"
Liu <- read.xls(si,sheet=1)

# Format and cleanup Data for the reanalysis
Liu_sub<-Liu[-c(2:4,6,18:22)] # Get Columns of Interest
#Liu_sub<-Liu_sub[-c(898:903),]
colnames(Liu_sub)<-c("ID","Site","Lat","Long","Mean_T","Mean_precip","AGB","L_AGB","T_AGB","Age","A_L_Ratio","A_T_Ratio","Ref")
Liu_sub<-subset(Liu_sub,AGB>0&Age>0) # Only Site with given Age and AGB
levels(Liu_sub$Site)<-seq(1:449) # Correct Site Level
Liu_sub$Age_sq<-Liu_sub$Age^2
Liu <- Liu_sub
rm(si,Liu_sub)

#### GGMap of study points ####
library(maps)
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
d <- SpatialPointsDataFrame(cbind(Liu$Long,Liu$Lat),data=data.frame(Liu$Age,Liu$AGB))
plot(d,color="red",cex=d$Age,pch=21,col="red",add=T)

library(ggmap)
br <- c(100,300,500,900,1200)
theme_set(theme_classic(base_size=12,base_family = "sans"))
mapWorld <- borders("world", colour="gray50",fill="gray50") # create a layer of borders
map <- ggplot() + mapWorld
map <- map + geom_path() + scale_y_continuous(breaks=(-2:2) * 30) +
  scale_x_continuous(breaks=(-4:4) * 45)
map <- map + coord_map("vandergrinten",xlim=c(-180,180),ylim=c(-60,70))
# Add rectangluar shape for belt of cancer 23.5 North and South
map <- map + geom_rect(aes(xmin=-180,xmax=180,ymin=-23.5,ymax=23.5),fill="orange",alpha=.3)
map <- map + ylab("Latitude (°)") + xlab("Longitude (°)")
# Add Coordinates of Study Points
coord <- subset(Liu,select=c(Long,Lat,AGB,Age))
map <- map + geom_point(data=coord,aes(x=Long, y=Lat,size=Age),alpha = .5,color="darkblue")  + scale_size_continuous(breaks=c(100,300,500,900,1200),range=c(1,6))# + scale_colour_gradient(low="blue",high="darkred",guide="colourbar") 
map <- map + scale_size_continuous(guide_legend(title = "Age (in years)"))
map
ggsave(filename="Figures/LIU_StudyMap-Size.png",plot=map,width=9,height=4,units="in",dpi=400)

# How many sites are in the tropics
Liu_trop <- subset(Liu,(Lat>-23.5&Lat<23.5))
nrow(Liu_trop)
max(Liu_trop$Age)

# Investigate the Structure
hist(Liu$Age)
hist(Liu$Mean_T)
hist(Liu$Mean_precip)

# Histograms indicate paucity at Age >500 years, Temp <-5 or >20, or precip >3000
# Thus we will sub the data for further analysis
Liu_subset<-subset(Liu,Age<500&Mean_T>-5&Mean_T<20&Mean_precip<3000)

# We fit a dummy variable to be used for models without random factor
Liu_subset$dummy<-rep(1,nrow(Liu_subset))


#### ReAnalysis - Nullmodels####

# We start modeling and build models imilar to those of liu et al (OLS) mentioned in the Appendix
# We considering everything independantly to each other but using our random variable structure to account for
# spatial autocorrelation and any systematic differences amongst studies

# First build a spatial correlation Matrix
# We change coordinates slightly since some sites have exactly the same coordinates
Liu_subset$Lat_J<-Liu_subset$Lat+(rnorm(length(Liu_subset$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu_subset)
corMatrix(cs1Exp)[1:10, 1:4] # Looks good

# Then we run some null-models to test if ether random variables or Spatial Autocorrelation 
# are appropriate on the Above-Ground-Biomass

null.model<-lme(log(AGB)~1,data=Liu_subset,random=~1|dummy,method="ML") # Without Random Strucutre
null.model2<-lme(log(AGB)~1,data=Liu_subset,random=~1|Ref,method="ML") # With Random Structure
null.model3<-lme(log(AGB)~1,data=Liu_subset,random=~1|Ref/Site,method="ML") # Hierarchical nested random Structure 
null.model4<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long),method="ML") # The same as above but including the Spatial Autocorrelation matrix
null.model5<- update(null.model3, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(nm_out <- aictab(list(null.model,null.model2,null.model3,null.model4,null.model5),sort=T,second.ord=F,
                  modnames=c("Null - w/o Random","Null - w. Random","Null - w. nested Random","Null - w. Random + SAC","Null - w. nested Random + SAC")) )

# Fitting a model that accounts for between study differences is better
# and we NEED to account for spatial autocorrelation for our results to
# be statistically valid. Models that account for Spatial Autocorrelation and random Structure perform best

#### ReAnalysis - Liu original Models ####
#The Following part is to show that models with included random factor and taking SAC
# into account the models outperform models without

#precipitation - a model with a squared term for mean_precip
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu_subset,random=~1|dummy,method="ML")
Precip_model_sac<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

#Temperature
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2),data=Liu_subset,random=~1|dummy,method="ML")
Temp_model_sac<-lme(log(AGB)~Mean_T+I(Mean_T^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

#Age
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu_subset,random=~1|dummy,method="ML")
Age_model_sac<-lme(log(AGB)~Age+I(Age^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(liuold_out <- aictab(list(Precip_model,Precip_model_sac,Temp_model,Temp_model_sac,Age_model,Age_model_sac),sort=T,second.ord=F,
                      modnames=c("AGB - Precip","AGB - Precip + SAC","AGB ~ Temp","AGB - Temp + SAC","AGB - Age","AGB - Age + SAC")) )

# All models with SAC and random structure perform better than Liu et al. original models.
# Furthermore due to colinearity between Predictors we need to consider interactions between
# predictors

#### ReAnalysis - Model with Interactions ####
# First we build a global model for use in model averaging that contains all varibles that are needed
# Squared and cubed terms temperature and precipitation are not included due to missing biological sense

All_model<-lme(log(AGB)~Age*Mean_precip*Mean_T*Age_sq,data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

# Now we dredge the model so that the value of each variable in predicting biomass
# can be assessed rather than using them in isolation
# Use second-order Information Criterion and keep Age as explanatory variable
MS1 <- dredge(All_model,evaluate=T,rank=AICc,trace=T,subset=dc(Age,Age_sq))
poss_mod <- get.models(MS1,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
modsumm2 <- subset(modsumm,modsumm$delta<7)
modsumm2
averaged <- model.avg(modsumm2,cumsum(weight)<=.95,fit=T)

# Finally run the best performing model including interactions from the model averaging process
# Use Liu et al. original models in comparison which consider those factors in isolation,
# but control for random structure and spatial autocorrelation
top_model<-lme(averaged$formula,data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(out <- aictab(list(top_model,Precip_model,Temp_model,Age_model),sort=T,second.ord=T,
               modnames=c("Model including Interactions","Precip. only","Temp. only","Age only")) )

# Our Best performing model includes several interactions between all used Predictors
# outperforms all other of Liu et al.s Models even if accounted for Spatial Autocorrelation 
# and random Structure.
# Thus indicates that Interactions are indeed important to determine Above Ground Biomass
