# Packages
library(maps)
library(sp)
library(ggmap)

#### Output Figure extent parameters ####
fig_h = 5
fig_w = 8
fig_dpi = 400
fig_units = "in"
fig_scale = 1.2

#### GGMap of study points ####
map("world", fill=TRUE, col="white", bg="lightblue", ylim=c(-60, 90), mar=c(0,0,0,0))
d <- SpatialPointsDataFrame(cbind(Liu$Long,Liu$Lat),data=data.frame(Liu$Age,Liu$AGB))
plot(d,color="red",cex=d$Age,pch=21,col="red",add=T)


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
ggsave(filename="Figures/LIU_StudyMap-Size.png",plot=map,height=fig_h,width=fig_w,dpi=fig_dpi,units=fig_units,scale=fig_scale)

# How many sites are in the tropics
Liu_trop <- subset(Liu,(Lat>-23.5&Lat<23.5))
nrow(Liu_trop)
max(Liu_trop$Age)


#### Figures Predictions ####
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
