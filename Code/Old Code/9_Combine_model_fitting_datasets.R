###############################################################
#  Combine tree, climate, N deposition and disturbance datasets
#  Add variable for time since max and most recent disturbances
#  12.02.25
###############################################################

rm(list=ls())

#  set user computer
computer <- "~"

#  directores
if(computer == "~"){
  Dir.Base <- "~/Desktop/Prague/mort-likelihood"
}else{
  Dir.Base <- "insert path"
}
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Code <- file.path(Dir.Base, "Code")
Dir.Clim <- file.path(Dir.Data,"Climate")
Dir.Dist <- file.path(Dir.Data,"Disturbance")

Dir_list <- c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim, Dir.Dist)

for(index in 1:length(Dir_list))   # using loop instead of sapply
{  Dir_sel <- Dir_list[index]   # select directories one by one
   if(!dir.exists(Dir_sel)) { dir.create(Dir_sel) }
}




### get tree plot data
load("~/Desktop/Prague/mort-likelihood/Data/8_Trees_with_mortality.Rdata")
targets$plotid <- as.factor(targets$plotid)
str(targets)  # n=198 plotids




### load climate
clim_df <- load(file=file.path(Dir.Clim,"6_plot_mean_census_climate.Rdata"))
clim_df <- get(clim_df)
str(clim_df)  # n=198 plotids

# ### get N deposition
# load(file=file.path(Dir.Data, "7_plot_mean_census_nitrogen_deposition_v2.Rdata"))
# str(plot_Ndep)   # n=198 plotids


### load disturbance history variables (revised version)
#  max severity and timing, recent severity and timing, number of peaks, and mean decadal disturbance
Dist <- read.csv(file=file.path(Dir.Data,"Disturbance.csv"))
Dist$plotid <- as.factor(Dist$plotid)
str(Dist)   # n=207 plotids


####### check for missing data
dist_plots <- levels(Dist$plotid)
target_plots <- levels(targets$plotid)

setdiff(dist_plots,target_plots)
#  [1] "CRO_RIS_009_1" "CRO_RIS_009_2" "CRO_RIS_010_1" "CRO_RIS_010_2" "CRO_RIS_011_1"
#  [6] "CRO_RIS_011_2" "CRO_RIS_015_1" "CRO_RIS_015_2" "CRO_RIS_016_1" "CRO_RIS_016_2"
# [11] "CRO_RIS_017_1" "CRO_RIS_017_2" "CRO_RIS_022_2"

setdiff(target_plots,dist_plots)
#[1] "BOS_PER_202_2" "BUL_STE_487_2" "BUL_STE_488_1" "CRO_RIE_092_2"


### merge datasets
Dist <- Dist[,!names(Dist) %in% "stand"]

trees = Reduce(function(x,y) {merge(x,y, by=c("plotid"),all.x=T)}, 
               list(targets, clim_df, Dist) )


# add variables for time since disturbance
trees$time_since_dist_max <- trees$date_census_1 - trees$my_max_year
trees$time_since_dist_recent <- trees$date_census_1 - trees$my_rec_year

# save
save(trees, file=file.path(Dir.Data,"9_Trees_with_derived_climate_variables_v2.Rdata"))





### get tree plot data
load(file=file.path(Dir.Data,"8_Trees_with_mortality.Rdata"))
targets$plotid <- as.factor(targets$plotid)
str(targets)

### load climate
clim_df <- load(file=file.path(Dir.Clim,"6_plot_mean_census_climate.Rdata"))
clim_df <- get(clim_df)
str(clim_df)

### get N deposition
load(file=file.path(Dir.Data, "7_plot_mean_census_nitrogen_deposition.Rdata"))
str(plot_Ndep)


### load disturbance history variables
#  max severity and timing, recent severity and timing, number of peaks, and mean decadal disturbance
Dist <- read_csv("~/Desktop/Prague/mort-likelihood/Data/disturbance.csv")
Dist$plotid <- as.factor(Dist$plotid)
head(Dist)
str(Dist)  # 166 plotids ?
summary(Dist)

# Dist <- Dist %>%
#   filter(!stand %in% "Risnjak")

### merge datasets
Dist <- Dist[,!names(Dist) %in% "stand"]

trees = Reduce(function(x,y) {merge(x,y, by=c("plotid"),all.x=T)}, 
               list(targets, clim_df, Dist) )


# add variables for time since disturbance
trees$time_since_dist_max <- trees$date_census_1 - trees$my_max_year
trees$time_since_dist_recent <- trees$date_census_1 - trees$my_rec_year

str(trees)
summary(trees)


####### check for missing data
dist_plots <- levels(Dist$plotid)
target_plots <- levels(targets$plotid)

setdiff(dist_plots,target_plots)
#[1] "ROM_IZV_377_1" "ROM_IZV_402_1" "ROM_IZV_402_2" "ROM_IZV_403_2" "ROM_IZV_404_1"
#[6] "ROM_IZV_428_1" "ROM_IZV_428_2" "ROM_IZV_429_1" "ROM_IZV_429_2"
setdiff(target_plots,dist_plots)
# character(0)

# save
save(trees, file=file.path(Dir.Data,"9_Full_model_fitting_dataset.Rdata"))


# 
# #####################################
# # exploratory graphs of annual climate
# Dir.Fig <- file.path(Dir.Base, "Figures")
# if(!dir.exists(Dir.Fig)) {dir.create(Dir.Fig)}
# setwd(Dir.Fig)
# 
# 
# trees$stand <- as.factor(trees$stand)
# ns <- length(levels(trees$stand))
# colors <- rainbow(ns)
# clab <- 1.3
# 
# 
# #### anomalies all species
# quartz(height=6,width=6)
# par(mar=c(5,5,3,1), oma=c(1,1,0,1))
# 
# plot(trees$temp_k_summer_diff, trees$WD_mm_mean_summer_diff,col=colors[trees$stand],
# 	xlab="",ylab="",las=1,main="Anomalies all species")
# mtext("Summer temperature anomaly (K)",side=1,line=3,outer=F,cex=clab)
# mtext("Summer water deficit anomaly (mm)",side=2,line=4,outer=F,cex=clab)
# 
# quartz.save("All_species_temp_WD_summer_anomalies.pdf", type="pdf")
# 
# 
# ##### census interval all species
# quartz(height=6,width=6)
# par(mar=c(5,5,3,1), oma=c(1,1,0,1))
# 
# plot(trees$temp_k_mean_summer_int,trees$WB_mm_mean_seas_int,col=colors[trees$stand],
# 	xlab="",ylab="",las=1,main="All species census climate")
# mtext("Annual census temperature (K)",side=1,line=3,outer=F,cex=clab)
# mtext("Seasonal water supply (mm)",side=2,line=4,outer=F,cex=clab)
# 
# quartz.save(file=file.path(Dir.Fig,"All_species_census_temp_vs_water_deficit.pdf"), type="pdf")
# 
# 
# 











