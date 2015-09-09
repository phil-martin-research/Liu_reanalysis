library(readxl)
library(dplyr)
library(taxize)
library(ggplot2)
library(stringr)
# Load data
data <- read_excel("Data/geb12113-sup-0001-ts1.xlsx",1)

# Format names
names(data) <- make.names(names(data))

# Select cols of interest and aggregate
d <- data %>% select(Number,Biome,Region,Above.ground.biomass.carbon.density..BCDa...Mg.C.ha.1,Stand.age,FullBibEntry) %>% 
  group_by(Biome,Region,FullBibEntry) %>% 
  summarise(AGB_avg = mean(Above.ground.biomass.carbon.density..BCDa...Mg.C.ha.1,na.rm=T),
            Age_avg = mean(Stand.age,na.rm=T)) %>% 
  # Reorder and rename
  select(Biome,Region,AGB_avg,Age_avg,FullBibEntry) %>% 
  rename("FullBibEntry" = "Reference")

write.csv2(d,"Results/StudyList.csv",row.names=T)
