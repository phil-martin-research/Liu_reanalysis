# Include BAAD data
#http://www.esajournals.org/doi/abs/10.1890/14-1889.1
library(dplyr)
library(ggplot2)
fd <- tempfile()
download.file("http://www.esapubs.org/archive/ecol/E096/128/baad_data.zip",destfile = fd)
b <- unzip(fd)
data <- read.csv("baad_data/baad_data.csv",header=T) # the data
meta <- read.csv("baad_data/baad_dictionary.csv",header=T) # the column explanations
# Filter those out without an age
data_a <- data %>% filter(!is.na(age))

summary(data_a$age)
qplot(data_a$age,geom="density")

# On a map
world_map <- map_data("world")#Get world map info
p <- ggplot() + coord_fixed()#Create a base plot
p <- p + geom_polygon(data=world_map,aes(x=long,y=lat,group=group),fill="white")#Add map to base plot
p <- p + geom_rect(aes(xmin=-180,xmax=180,ymin=-23.5,ymax=23.5),fill="orange",alpha=.3)
p <- p + ylab("Latitude (°)") + xlab("Longitude (°)")
p + geom_point(data=data_a,aes(x=longitude,y=latitude,size=age),alpha=0.5)


# After Brown and Lugo
# Aboveground biomass density (t/ha) = VOB * WD * BEF 
# AGB = r.st ...
