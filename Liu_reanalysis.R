# this is a script to reanalyse the 
# data from Lio et al 2013 paper in GEB
rm(list=ls())
library(ggplot2)
library(nlme)
library(GGally)
library(ape)
library(MuMIn)
library(mgcv)
library(influence.ME)
library(AICcmodavg)
library(scales)
#read in data
Liu<-read.csv("Liu_Aged.csv")

hist(Liu$Age)
hist(Liu$Mean_T)
hist(Liu$Mean_precip)
hist(log(Liu$AGB))

#### GGMap of study points ####
library(maps)
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
d <- SpatialPointsDataFrame(cbind(Liu$Long,Liu$Lat),data=data.frame(Liu$Age,Liu$AGB))
plot(d,color="red",cex=d$Age,pch=21,col="red",add=T)

library(ggmap)
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
map
ggsave(filename="Figures/LIU_StudyMap-Size.png",plot=map,width=9,height=4,units="in",dpi=400)

# How many sites are in the tropics
Liu_trop <- subset(Liu,(Lat>-23.5&Lat<23.5))
nrow(Liu_trop)
max(Liu_trop$Age)

#looking at these histograms I wouldn't be happy to make predictions
#about biomass in forests of Age >500 years, Temp <-5 or >20, or precip >3000
#so I will subset the data to remove these
Liu_subset<-subset(Liu,Age<400&Mean_T>-5&Mean_T<20&Mean_precip<3000)
hist(Liu_subset$Age)
hist(Liu_subset$Mean_T)
hist(Liu_subset$Mean_precip)
hist(Liu_subset$AGB)

# Correlation between Predictors
cor.test(Liu_subset$Mean_precip,Liu_subset$Mean_T,method="pearson",exact=T,conf.level=T)
cor.test(Liu$Mean_precip,Liu$Mean_T,method="pearson",exact=T,conf.level=T)

require(psych)
pairs.panels(Liu[,c(4:8,11)]) # Correlations between Predictors


# Although statistically reasonable this really limits our direct comparison with Liu et al.s 
# Models, which is why I will use the whole dataset in the following analysis and investigate 
# the influence of those poorly sampled regions, respectively 

#start models - these models are similar to those of liu et al
#considering everything independantly to each other
#but using our random variable structure to account for
#spatial autocorrelation and any systematic differences amongst studies

#first we fit a dummy random variable
Liu_subset$dummy<-rep(1,nrow(Liu_subset))

#now build in spatial autocorrelation

#change coordinates slightly since some sites 
#have exactly the same coordinates
Liu_subset$Lat_J<-Liu_subset$Lat+(rnorm(length(Liu_subset$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu_subset)
corMatrix(cs1Exp)[1:10, 1:4]

#first run some null models to check our random
#variable structure is appropriate
null.model<-lme(log(AGB)~1,data=Liu_subset,random=~1|dummy,method="ML")
null.model2<-lme(log(AGB)~1,data=Liu_subset,random=~1|Ref,method="ML")
null.model3<-lme(log(AGB)~1,data=Liu_subset,random=~1|Ref/Site,method="ML")
null.model4<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
null.model5<- update(null.model3, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
null.model6<- update(null.model, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(nm_out <- aictab(list(null.model,null.model2,null.model3,null.model4,null.model5,null.model6),sort=T,second.ord=F,
       modnames=c("Null - w/o Random","Null - w. Random","Null - w. nested Random","Null - w. Random + SAC","Null - w. nested Random + SAC","Null - w/o Random + SAC")) )

aictab(list(null.model,null.model2),modnames=c("null","Null+SAC"),second.order=F)
write.csv(nm_out,"Results/NullModel.Comparison.csv",row.names=F)

#fitting a model that accounts for between study differences is better
#and we NEED to account for spatial autocorrelation for our results to
#be statistically valid

#### LIUs Model ####
#Following part is to show that even with included random factor and taking SAC
#into account the models perform bad if interaction effects are not considered

#now let's fit the models that Liu et al used for: 
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
mean(c(liuold_out$Delta_AIC[6] - liuold_out$Delta_AIC[1],
liuold_out$Delta_AIC[5] - liuold_out$Delta_AIC[2], 
liuold_out$Delta_AIC[4] - liuold_out$Delta_AIC[3])
)
liuold_out$adj.r.squared =NA
liuold_out$adj.r.squared[c(4,5,6)]<- unlist(lapply(list(Precip_model,Age_model,Temp_model),
                                          FUN=function(x)attr(r.squaredLR(x,null=null.model3),"adj.r.squared")) )
# For LR-ratio comparison with SAC nullmodel
liuold_out$adj.r.squared[c(1,2,3)]<- unlist(lapply(list(Precip_model_sac,Temp_model_sac,Age_model_sac),
                                                   FUN=function(x)attr(r.squaredLR(x,null=null.model5),"adj.r.squared")) )

write.csv(liuold_out,"Results/LIUmodels_SAC.csv",row.names=F)

# Alternatively build a correlogramm with Moran's I using your corrected Lat_j and Long
# as input. Also we probably want to check if the residuals of Liu et al. original models
# without your correlation matrix are autocorrelated.
#require(pgirmess)
#corr_p <- correlog(cbind(Liu_subset$Lat_J,Liu_subset$Long),resid(Precip_model_sac));plot(corr_p);abline(h=0)
#corr_t <- correlog(cbind(Liu_subset$Lat_J,Liu_subset$Long),resid(Temp_model_sac));plot(corr_t);abline(h=0)
#corr_a <- correlog(cbind(Liu_subset$Lat_J,Liu_subset$Long),resid(Age_model_sac));plot(corr_a);abline(h=0)
# Not all of them look too serious at all, but SAC is existant. You can improve the Liu's model
# relative parsimony by incorporating /site as nested effect as well. Thus the low number of different sites
# and values (especially in the tropics) seems to even prevent more serious SAC.

#### Automated Model selection process ####
# now a global model for use in model averaging that contains all varibles that are needed
#I haven't included the squared and cubed terms for temp and precipitation becuase I don't think they make biological sense
All_model<-lme(log(AGB)~Age*Mean_precip*Mean_T*Age_sq,data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
#now we can dredge the model so that the value of each variable in predicting biomass
#can be assessed rather than using them in isolation

MS1<-dredge(All_model,evaluate=T,rank=AICc,trace=T,subset=dc(Age,Age_sq))
poss_mod<-get.models(MS1,subset=delta<7)
modsumm<- model.sel(poss_mod, rank = "AICc",fit=T)
modsumm2<-subset(modsumm,modsumm$delta<7)
modsumm2
averaged<-model.avg(modsumm2,cumsum(weight)<=.95,fit=T)

#run the top model from the model averaging to get
top_model<-lme(averaged$formula,data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Precip_model<-lme(log(AGB)~Mean_precip+I(Mean_precip^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Temp_model<-lme(log(AGB)~Mean_T+I(Mean_T^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")
Age_model<-lme(log(AGB)~Age+I(Age^2),data=Liu_subset,random=~1|Ref/Site,correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

(out <- aictab(list(top_model,Precip_model,Temp_model,Age_model),sort=T,second.ord=T,
              modnames=c("Model including Interactions","Precip. only","Temp. only","Age only")) )
# Our Model performs best.
# Include results of LR-ratio test against null model
out$adj.r.squared <- unlist(lapply(list(top_model,Precip_model,Temp_model,Age_model),FUN=function(x)attr(r.squaredLR(x,null=null.model5),"adj.r.squared")) )
write.csv(out,"Results/FinalModelComparison.csv",row.names=F)

#### Figures ####
Liu_subset$pred <- NA # For the Rug
#create dataframe for predictions so paramaters can be plotted
#first - Age
nseq <- function(x, len = length(x)) seq(min(x, na.rm = TRUE),
    max(x, na.rm=TRUE), length = len)
str(Liu_subset)
newdata_Age1<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Age1$Age<- nseq(Liu_subset$Age, nrow(newdata_Age1))
newdata_Age1$Age_sq<- newdata_Age1$Age^2
newdata_Age1$Mean_T<-0
newdata_Age1$pred<-predict(top_model,newdata_Age1,level=0)
newdata_Age2<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Age2$Age<- nseq(Liu_subset$Age, nrow(newdata_Age2))
newdata_Age2$Age_sq<- newdata_Age2$Age^2
newdata_Age2$Mean_T<-5
newdata_Age2$pred<-predict(top_model,newdata_Age2,level=0)
newdata_Age3<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Age3$Age<- nseq(Liu_subset$Age, nrow(newdata_Age3))
newdata_Age3$Age_sq<- newdata_Age3$Age^2
newdata_Age3$Mean_T<-10
newdata_Age3$pred<-predict(top_model,newdata_Age3,level=0)

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

names(newdata_Age)
#now plot the predictions for Age
theme_set(theme_classic(base_size=12,base_family = "sans"))
Age_plot <- ggplot(data=newdata_Age,aes(x=Age,y=exp(pred),group=as.factor(Mean_T),colour=as.factor(Mean_T)))
Age_plot <- Age_plot+geom_line(data=newdata_Age,aes(y=exp(Liu_pred),group=NULL),colour="black",lty=2,size=2,lineend="round") + geom_rug(data = Liu_subset,side="b",color="black",show_guide=F)
Age_plot2 <- Age_plot + geom_line(size=2,lineend="round")
#Age_plot3 <- Age_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Age_plot4 <- Age_plot2+ylab(expression(paste("Aboveground \nbiomass (Mg ",ha^-1,")")))+xlab("Age (Years)") 
Age_plot5 <- Age_plot4+scale_colour_brewer(expression("Mean annual\ntemperature("*degree*C*")"),palette="Set1")
Age_plot5
ggsave(filename="Figures/Age_temp.png",plot=Age_plot5,width=9,height=4,units="in",dpi=400)

#now make predictions for Temp
plot(Liu_subset$Age,Liu_subset$Mean_T)
newdata_Temp1<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Temp1$Mean_T<- nseq(Liu_subset$Mean_T, nrow(newdata_Temp1))
newdata_Temp1$Age<-100
newdata_Temp1$Age_sq<-newdata_Temp1$Age^2
newdata_Temp2<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Temp2$Mean_T<- nseq(Liu_subset$Mean_T, nrow(newdata_Temp2))
newdata_Temp2$Age<-150
newdata_Temp2$Age_sq<-newdata_Temp2$Age^2
newdata_Temp3<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Temp3$Mean_T<- nseq(Liu_subset$Mean_T, nrow(newdata_Temp3))
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
theme_set(theme_classic(base_size=12,base_family = "sans"))
Temp_plot<-ggplot(data=newdata_Temp,aes(x=Mean_T,y=exp(pred),group=as.factor(Age),colour=as.factor(Age)))
Temp_plot2<-Temp_plot+geom_line(data=newdata_Temp,aes(y=exp(Liu_pred)),colour="black",lty=2,size=2) +geom_line(size=2,guide="legend") + geom_rug(data = Liu_subset,side="b",color="black",show_guide=F)
#Temp_plot3<-Temp_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Temp_plot4 <- Temp_plot2+ylab(expression(paste("Aboveground \nbiomass (Mg ",ha^-1,")")))+xlab("Mean annual temperature (in °C)")
Temp_plot5 <- Temp_plot4 + scale_colour_brewer("Age \n(Years)",palette="Set1")
Temp_plot5 <- Temp_plot5 + ggtitle("Effect of Temperature")
Temp_plot5   
ggsave(filename="Figures/Temp_age.png",plot=Temp_plot5,width=9,height=4,units="in",dpi=400)

#now predictions for changes in precipitation
newdata_P<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_P$Mean_precip<- nseq(Liu_subset$Mean_precip, nrow(newdata_P))

#now use predict
newdata_P$pred<-predict(top_model,newdata_P,level=0)

#add predictions from model which is most similar to that of Liu et al
newdata_P$Liu_pred<-predict(Precip_model,newdata_P,level=0)

#now plot
theme_set(theme_classic(base_size=12,base_family = "sans"))
Temp_plot<-ggplot(data=newdata_P,aes(x=Mean_precip,y=exp(pred)))+geom_line(size=2)
Temp_plot2<-Temp_plot+geom_line(data=newdata_P,aes(y=exp(Liu_pred)),colour="black",lty=2,size=2) + geom_rug(data = Liu_subset,side="b",color="black",show_guide=F)
Temp_plot4<-Temp_plot2+ylab(expression(paste("Aboveground\n biomass (Mg ",ha^-1,")")))+xlab("Mean annual precipitation (mm)")
Temp_plot4
ggsave(filename="Figures/Precipitation.png",plot=Temp_plot4,width=9,height=4,units="in",dpi=400)


# Do the comparison for a age-precip design matrix
plot(Liu_subset$Age,Liu_subset$Mean_precip)
newdata_Precip1<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Precip1$Mean_precip<- nseq(Liu_subset$Mean_precip, nrow(newdata_Precip1))
newdata_Precip1$Age<-100
newdata_Precip1$Age_sq<-newdata_Precip1$Age^2
newdata_Precip2<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Precip2$Mean_precip<- nseq(Liu_subset$Mean_precip, nrow(newdata_Precip2))
newdata_Precip2$Age<-150
newdata_Precip2$Age_sq<-newdata_Precip2$Age^2
newdata_Precip3<- as.data.frame(lapply(lapply(Liu_subset[c(6,7,11,15)], mean), rep, 300))
newdata_Precip3$Mean_precip<- nseq(Liu_subset$Mean_precip, nrow(newdata_Precip3))
newdata_Precip3$Age<-200
newdata_Precip3$Age_sq<-newdata_Precip3$Age^2
newdata_Precip<-rbind(newdata_Precip1,newdata_Precip2,newdata_Precip3)

#now use predict
newdata_Precip$pred<-predict(top_model,newdata_Precip,level=0)
#add predictions from model which is most similar to that of Liu et al
newdata_Precip$Liu_pred<-predict(Precip_model,newdata_Precip,level=0)

summary(newdata_Precip)
plot(newdata_Precip$Mean_precip,exp(newdata_Precip$pred))

#now plot the predictions for Precipitation
theme_set(theme_classic(base_size=12,base_family = "sans"))
Precip_plot<-ggplot(data=newdata_Precip,aes(x=Mean_precip,y=exp(pred),group=as.factor(Age),colour=as.factor(Age)))
Precip_plot2<-Precip_plot+geom_line(data=newdata_Precip,aes(y=exp(Liu_pred)),colour="black",lty=2,size=2) +geom_line(size=2)+ geom_rug(data = Liu_subset,side="b",color="black",show_guide=F)
#Temp_plot3<-Temp_plot2+theme(panel.grid.major = element_line(colour =NA),panel.grid.minor = element_line(colour =NA),panel.border = element_rect(size=1.5,colour="black",fill=NA))
Precip_plot4 <- Precip_plot2+xlab("Mean annual precipitation (mm)")+ylab("")#+ylab(expression(paste("Aboveground \nbiomass (Mg ",ha^-1,")")))
Precip_plot5 <- Precip_plot4 + scale_colour_brewer("Age \n(Years)",palette="Set1") + theme(legend.position="bottom")
Precip_plot5 <- Precip_plot5 + ggtitle("Effect of Precipitation")
Precip_plot5  
ggsave(filename="Figures/Precip_age.png",plot=Precip_plot5,width=9,height=4,units="in",dpi=400)

# On one page with combined Age legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Precip_plot5)


library(gridExtra)
png("Figures/CombinedTemperaturePrecip.png",width=8,height=6,units="in",res=400)
grid.arrange(arrangeGrob(Temp_plot5 + theme(legend.position="none"),
                               Precip_plot5 + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))
dev.off()

