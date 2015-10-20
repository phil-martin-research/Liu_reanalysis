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
library(cowplot)
library(gridExtra)

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
                     data="Global forests")

#stick these data from sites and all forests together
Climate<-rbind(Liu_clim,all_data)
Climate<-subset(Climate,!is.na(precip)&!is.na(temp))
Climate$Climate_precip_bin<-as.numeric(as.character(cut(Climate$precip,
    breaks=(seq(min(Climate$precip),max(Climate$precip),by=200)),
    labels=seq(min(Climate$precip),max(Climate$precip),by=200)[-1])))
Climate$Climate_temp_bin<-as.numeric(as.character(cut(Climate$temp,
                        breaks=(seq(min(Climate$temp),max(Climate$temp),by=1)),
                        labels=seq(min(Climate$temp),max(Climate$temp),by=1)[-1])))
Climate$Climate_temp_bin2<-as.factor(Climate$Climate_temp_bin)

Climate_sum1<-ddply(Climate,.(data),mutate,
      Total=length(Climate_temp_bin))
head(Climate_sum1)

Climate_sum2<-ddply(Climate_sum1,.(data,Climate_temp_bin,Climate_precip_bin),
                    summarise,mean_Total=mean(Total),
                    Perc=(length(Climate_temp_bin)/(mean(Total))*100))
head(Climate_sum2)

sum(unique(Climate_sum2$Perc))




#look at biases in the data that may influence results

#climate space
theme_set(theme_bw(base_size=12))
C_plot1<-ggplot(Climate_sum2,aes(x=Climate_precip_bin,y=Climate_temp_bin,fill=Perc*2))+geom_raster()+facet_wrap(~data)
C_plot2<-C_plot1+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA),legend.position="none")
C_plot3<-C_plot2+ylab("Mean annual temperature")+xlab("Total annual precipitation (mm)")+scale_fill_gradient("Percentage",low="light grey",high="black")
C_plot4<-C_plot3+geom_text(data=data.frame(x=400,y=30,text="(b)",data="Our data"),aes(x=x,y=y,label=text,fill=NULL))


#age
theme_set(theme_bw(base_size=12))
Geom_hist<-ggplot(Liu,aes(x=Age))+geom_histogram()+ylab("number of sites")+xlab("Estimate age (Years)")+
theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA))+
        coord_cartesian(xlim = c(-10,1050),ylim=c(0,200))+
        annotate(geom = "text",x = 30,y=190,label="(c)")


#location
forest_coarse <- aggregate(forest,fact=20)
forest_points<-rasterToPoints(forest_coarse)
df <- data.frame(forest_points)
colnames(df)<-c("Long","Lat","Forest")



theme_set(theme_bw(base_size=12))
world_map <- map_data("world")#Get world map info

p <- ggplot() + coord_equal(xlim = c(-170,180),ylim=c(-60,80))#Create a base plot

base_world <- p + geom_polygon(data=world_map,aes(x=long,y=lat,group=group),fill="light grey")#Add map to base plot
base_world2<-base_world+geom_raster(data=df,aes(y=Lat,x=Long,fill=Forest,alpha=Forest))+scale_fill_gradient(low="white",high="light green")
Location<-base_world2 + geom_point(data=Liu,aes(x=Long,y=Lat),colour="black",alpha=0.2)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.5,colour="black",fill=NA),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none")+annotate(geom = "text",x = -160,y=70,label="(a)")
Location

pdf("Figures/Bias.pdf",width = 8,height = 8,units = "in",res =400)
grid.arrange(Location,arrangeGrob(C_plot4,Geom_hist,ncol=2),ncol=1)
dev.off()
