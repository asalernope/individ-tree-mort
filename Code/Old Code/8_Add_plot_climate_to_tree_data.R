###########################################################
#  Add derived climate variables to tree dataset
#  Dec 2024
###########################################################

rm(list=ls())

computer <- "~"

# Directores
if(computer == "~"){
  Dir.Base <- "~/Desktop/Prague/mort-likelihood"
}else{
  Dir.Base <- "insert path"
}

Dir.Data <- file.path(Dir.Base, "Data")
Dir.Code <- file.path(Dir.Base, "Code")
Dir.Clim <- file.path(Dir.Data,"Climate")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim), function(x) if(!dir.exists(x)) dir.create(x))


# load climate
clim_df <- load(file=file.path(Dir.Clim,"6_plot_mean_census_climate.Rdata"))
clim_df <- get(clim_df)
names(clim_df)

# get tree plot data
# load(file=file.path(Dir.Data,"7_Trees_with_mortality.Rdata"))

 load("~/Desktop/Prague/mort-likelihood/Data/8_Trees_with_mortality.Rdata")
targets<-get(targets)


# combine trees and climate data for period of available climate (1958-2023)
trees <- merge(targets,clim_df, by=c("plotid"), all=FALSE)
str(trees)
summary(trees)
# nrow(trees[trees$shading==0,])
# nrow(trees[trees$shading==1,])

### save
save(trees,file=file.path(Dir.Data,"8_Trees_with_derived_climate_variables.Rdata"))


#####################################
# exploratory graphs of annual climate
Dir.Fig <- file.path(Dir.Base, "Figures")
if(!dir.exists(Dir.Fig)) {dir.create(Dir.Fig)}

ns <- length(levels(trees$stand))
colors <- rainbow(ns)
clab <- 1.3

# anomalies FAGSYL
quartz(height=6,width=6)
par(mar=c(5,5,3,1), oma=c(1,1,0,1))

plot(trees$temp_k_ann_diff[trees$SPCD %in% "FAGSYL"],trees$precip_ann_diff[trees$SPCD %in% "FAGSYL"],
	xlab="",ylab="",las=1,main="Fagus climate anomalies",col=colors[trees$stand])
mtext("Annual temperature diff (K)",side=1,line=3,outer=F,cex=clab)
mtext("Annual precipitation diff (mm)",side=2,line=4,outer=F,cex=clab)

setwd(Dir.Fig)
quartz.save("FAGSYL_climate_anomalies.pdf", type="pdf")


# anomalies all species
quartz(height=6,width=6)
par(mar=c(5,5,3,1), oma=c(1,1,0,1))

plot(trees$temp_k_ann_diff,trees$precip_ann_diff,col=colors[trees$stand],
	xlab="",ylab="",las=1,main="Anomalies all species")
mtext("Annual temperature diff (K)",side=1,line=3,outer=F,cex=clab)
mtext("Annual precipitation diff (mm)",side=2,line=4,outer=F,cex=clab)

quartz.save("All_species_climate_anomalies.pdf", type="pdf")


# census interval all species
quartz(height=6,width=6)
par(mar=c(5,5,3,1), oma=c(1,1,0,1))

plot(trees$mean_ann_temps_int,trees$precip_mean_ann_int,col=colors[trees$stand],
	xlab="",ylab="",las=1,main="All species census climate")
mtext("Annual census temperature (K)",side=1,line=3,outer=F,cex=clab)
mtext("Annual census precipitation (mm)",side=2,line=4,outer=F,cex=clab)

quartz.save(file=file.path(Dir.Fig,"All_species_census_climate.pdf"), type="pdf")

library(tidyverse)
library(ggplot2)

clim_df %>%
ggplot() +
  geom_point(aes(temp_k_summer_diff,WB_mm_mean_summer_diff))

