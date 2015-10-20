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
# Load in data taken from LIU et al. DOI: 10.1111/geb.12113
Liu <- read.csv("Data/Liu_aged.csv")

# Null-models to test if ether random variables or Spatial Autocorrelation 
# are appropriate on the Above-Ground-Biomass
# We fit a dummy variable to be used for models without random factor
Liu$dummy <- rep(1,nrow(Liu))

# First build a spatial correlation Matrix
# We change coordinates slightly since some sites have exactly the same coordinates
Liu$Lat_J <- Liu$Lat+(rnorm(length(Liu$Lat),0,0.00001)) 
cs1Exp <- corExp(1, form = ~ Lat_J + Long)
cs1Exp <- Initialize(cs1Exp, Liu)

null.model<-lme(AGB~1,data=Liu,random=~1|dummy,method="ML") # Without Random Structure - equivalent to original OLS 
null.model2<-lme(AGB~1,data=Liu,random=~1|Ref,method="ML") # With Random Structure
null.model3<-lme(AGB~1,data=Liu,random=~1|Ref/Site,method="ML") # Hierarchical nested random Structure 
null.model4<- update(null.model2, correlation = corExp(1, form = ~ Lat_J + Long),method="ML") # The same as above but including the Spatial Autocorrelation matrix
null.model5<- update(null.model3, correlation = corExp(1, form = ~ Lat_J + Long),method="ML")

# Null model comparison
(nm_out <- aictab(list(null.model,null.model2,null.model3,null.model4,null.model5),sort=T,second.ord=F,
                  modnames=c("Null - w/o Random","Null - w. Random","Null - w. nested Random","Null - w. Random + SAC","Null - w. nested Random + SAC")) )
# Models that account for differences between studies and / or spatial autocorrelation
# outperform models that have no such structure (such as Liu et al. original OLS)

#### Analysis - Model with Interactions ####
# First we build a global model for use in model averaging that contains all varibles that are needed

Liu$Mean_T2<-(Liu$Mean_T-mean(Liu$Mean_T))/sd(Liu$Mean_T)
Liu$Precip2<-(Liu$Mean_precip-mean(Liu$Mean_precip))/sd(Liu$Mean_precip)
Liu$Age2<-(Liu$Age-mean(Liu$Age))/sd(Liu$Age)
Liu$Age2sq<-Liu$Age2^2
Liu$Age3<-(log(Liu$Age)-log(mean(Liu$Age)))/sd(log(Liu$Age))

mymodel<-lme(log(AGB)~Age3*Mean_T2+Age3*Precip2+Precip2*Mean_T2,
             data=Liu,
             random=~1|Ref/Site,
             correlation = corExp(1, form = ~ Lat_J + Long),
             method="ML")

plot(mymodel)

qplot(Liu$Age2,resid(mymodel))+geom_smooth()
qplot(Liu$Mean_T2,resid(mymodel))+geom_smooth()
qplot(Liu$Precip2,resid(mymodel))+geom_smooth()
qqnorm(mymodel,abline = c(0, 1))

# Heteroskedasticity is present in some cases, 
# likely due to small sample sizes with high variances in extreme regions (few samples in tropics)

# Now we dredge the model so that the value of each variable in predicting biomass
# can be assessed rather than using them in isolation
# Use second-order Information Criterion and keep Age as explanatory variable
MS1 <- dredge(mymodel,evaluate=T,rank=AICc,trace=T)
Model_comp<-data.frame(model.sel(MS1, rank = "AICc",fit=T)) # Rank and select the best models
Model_comp<-round(Model_comp[,-1],2)
colnames(Model_comp)[1:6]<-c("Age","Temperature","Precipitation","Age*Temperature","Temperature*Precipitation","Age*Precipitation")
write.csv(Model_comp,"Tables/Model_comparison.csv",row.names=F)

poss_mod <- get.models(MS1,subset=delta<7)
modsumm <- model.sel(poss_mod, rank = "AICc",fit=T) # Rank and select the best models
Mod.avg<-model.avg(modsumm)
summary(Mod.avg)
r.squaredGLMM(mymodel)


#now create predictions
Liu$Mean_T_bins<-cut(Liu$Mean_T2,(quantile(Liu$Mean_T2)),labels=c("-16.8 ~*C - 0.2 ~*C", "0.2 ~*C - 3.8 ~*C", "3.8 ~*C - 8.0 ~*C", "8.0 ~*C - 26.2 ~*C"),include.lowest = T)
Liu$Mean_P_bins<-cut(Liu$Precip2,(quantile(Liu$Precip2)),labels=c("250-583mm", "583-796mm", "796-1079mm", "1079-5800mm"),include.lowest = T)


new.data<-expand.grid(Mean_T2=c(mean(c(-3.35,-0.66)),mean(c(-0.66,-0.09)),mean(c(-0.09,0.57)),mean(c(0.57,3.455))),
            Precip2=c(mean(c(-1.11,-0.61)),mean(c(-0.61,-0.28)),mean(c(-0.28,0.14)),mean(c(0.14,7.27))),
            Age3=seq(-0.7,8.79,by=0.1)
            )
new.data$AGB<-exp(predict(mymodel,newdata = new.data,level=0))
new.data$Mean_T_bins<-cut(new.data$Mean_T2,as.numeric(quantile(new.data$Mean_T2)),labels=c("-16.8 - 0.2", "0.2 - 3.8", "3.8 - 8.0", "8.0 - 26.2"),include.lowest = T)
new.data$Mean_P_bins<-cut(new.data$Precip2,as.numeric(quantile(new.data$Precip2)),labels=c("250 - 583", "583 - 796", "796 - 1079", "1079 - 5800"),include.lowest = T)


new.data2<-data.frame(Mean_T2=rep(c(mean(c(-3.35,-0.66)),mean(c(-0.66,-0.09)),mean(c(-0.09,0.57)),mean(c(0.57,3.455))),4),
                      Precip2=rep(c(mean(c(-1.11,-0.61)),mean(c(-0.61,-0.28)),mean(c(-0.28,0.14)),mean(c(0.14,7.27))),each = 4))
new.data2$Mean_T_bins<-cut(new.data2$Mean_T2,as.numeric(quantile(new.data2$Mean_T2)),labels=c("-16.8 - 0.2", "0.2 - 3.8", "3.8 - 8.0", "8.0 - 26.2"),include.lowest = T)
new.data2$Mean_P_bins<-cut(new.data2$Precip2,as.numeric(quantile(new.data2$Precip2)),labels=c("250-583mm", "583-796mm", "796-1079 mm", "1079-5800 mm"),include.lowest = T)


new.data3<-ddply(Liu,.(Mean_P_bins,Mean_T_bins),summarise,min_age=min(Age3),max_age=max(Age3),Mean_T2=mean(Mean_T2),Precip2=mean(Precip2))
new.data4<-merge(new.data2,new.data3,by=c("Mean_P_bins","Mean_T_bins"))
head(new.data3)

new.data5<-NULL
for (i in 1:nrow(new.data3)){
  data.sub<-data.frame(Mean_P_bins=new.data3$Mean_P_bins[i],
             Mean_T_bins=new.data3$Mean_T_bins[i],
             Mean_T2=new.data3$Mean_T2[i],
             Precip2=new.data3$Precip2[i],
             Age3=seq(new.data3$min_age[i],new.data3$max_age[i],0.01))
  new.data5<-rbind(new.data5,data.sub)
}

new.data5$AGB<-exp(predict(Mod.avg,newdata=new.data5,level=0))

(log(Liu$Age)-log(mean(Liu$Age)))/sd(log(Liu$Age))

theme_set(theme_bw(base_size=12))
Plot1<-ggplot(new.data5,aes(x=exp((Age3*sd(log(Liu$Age)))+log(mean(Liu$Age))),y=AGB))+geom_line(size=1)+facet_grid(Mean_T_bins~Mean_P_bins)
Plot1+geom_point(data=Liu,shape=1,alpha=0.2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA),
        legend.position="none")+xlab("Age (years)")+
  ylab(expression(paste("Aboveground biomass (Mg ", ha^-1,")",sep="")))
ggsave("Figures/AGB_change.pdf",width = 8,height=6,units="in",dpi=400)

