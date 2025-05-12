#####################################################################################################
#  Mortality - disturbance project Balkans
#  Calculate census interval mean climate for each plot
#  calculate climate anomalies (interval - antecedent)
#  Dec 2024
#  12.02.25 Revised tree dataset; 
#           Added extreme climate values for census interval
#               based on Canham & Murphy 2017
#           Fixed error in for loop
#           Revised climate variable names
#  April 2025: -use expanded tree dataset; 
#              -census interval varies in some plots
#              -add climate terms based on edits in step 4 (extreme monthly values for growing season
#####################################################################################################

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
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim), function(x) if(!dir.exists(x)) dir.create(x))


# read tree dataset
# name <- "A2_Balkan_trees_no_cores_code"
# load(file=file.path(Dir.Data, paste0(name,".Rdata")))

# targets <- load("~/Desktop/Prague/mort-likelihood/Data/full_trees.Rdata")
# targets<-get(targets)
# head(targets)
# summary(targets)
# 
# # rename census date columns
# colnames(targets)[names(targets) %in% "date"] <- "date_census_1"
# colnames(targets)[names(targets) %in% "date_re"] <- "date_census_2"
# 
# 
# ##### combine census dates
# survey_dates <- subset(targets,select=c(plotid,date_census_1,date_census_2))
# survey_dates <- unique(survey_dates)
# 
# 
# # remove NA
# survey_dates <- survey_dates[!is.na(survey_dates$date_census_2),]
# survey_dates <- survey_dates[order(survey_dates$plotid),]
# which(duplicated(survey_dates$plotid))
# nrow(survey_dates)
# 
# 
# write.csv(survey_dates,file=file.path(Dir.Data,"survey_dates.csv"),row.names=F)
# 
# # census interval summary
# summary(survey_dates$date_census_2 - survey_dates$date_census_1)  # consistent 5 year intervals
# #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# #     5       5       5       5       5       5 
# 
# 
# ################
# # load annual-resolution climate data (plot scale)
# clim_df <- load(file=file.path(Dir.Clim,"5_plot_monthly_water_budget_variables.Rdata"))
# clim_df <- get(clim_df)
# head(clim_df)
# 
# 
# clim_array <- array(data=list(NULL),dim=nrow(survey_dates))
# 
# # set length of antecedent climate interval
# period_length <- 30
# 
# # loop through plot and calculate quantiles
# pb <- txtProgressBar(min=0, max=length(survey_dates$plotid), style=3)  # progress bar
# tind=0
# 
# # loop through plots and compute mean climate for census and prior intervals
# #i <- 1
# for(i in 1:length(survey_dates$plotid) )   
# {
# 	tind=tind+1
# 	setTxtProgressBar(pb, tind)
# 	
# 	p <- survey_dates$plotid[i]
# 	
# 	# get census interval
# 	start_yr <- survey_dates$date_census_1[which(survey_dates$plotid %in% p)]
# 	end_yr <- survey_dates$date_census_2[which(survey_dates$plotid %in% p)]
# 	
# 	# subset climate data to census interval
# 	dat.plot.int <- subset(clim_df, year>=start_yr &
# 								year<=(end_yr-5) &
# 								plotid %in% p )
# 	
# 	### calculate mean climate for census interval
# 	temp_k_mean_ann_int <- mean(dat.plot.int$temp_ann_K,na.rm=T)   
# 	precip_mean_ann_int <- mean(dat.plot.int$precip_ann_mm,na.rm=T)
# 	wind_mean_ann_int <- mean(dat.plot.int$ws_ann,na.rm=T)
# 	WD_mm_mean_ann_int <- mean(dat.plot.int$WD_ann_mm,na.rm=T) # water deficit (PET-AET)
# 	WB_mm_mean_seas_int <- mean(dat.plot.int$water_balance_seasonal_mm,na.rm=T) # effective growing season water supply
# 	
# 	temp_k_mean_summer_int <- mean(dat.plot.int$temp_seasonal_K,na.rm=T)   
# 	WD_mm_mean_summer_int <- mean(dat.plot.int$WD_seasonal_mm,na.rm=T) # water deficit (PET-AET)
# 	SMI_mm_mean_summer_int <- mean(dat.plot.int$SMI_seasonal,na.rm=T) # effective growing season water supply
# 
# 	
# 	
# 	# extreme census climate values  (See Canham & Murphy 2017)
# 	temp_k_max_summer_int <- max(dat.plot.int$temp_seasonal_K,na.rm=T)   
# 	WD_mm_max_summer_int <- max(dat.plot.int$WD_seasonal_mm,na.rm=T) # water deficit (PET-AET)
# 	WB_mm_min_seas_int <- min(dat.plot.int$water_balance_seasonal_mm,na.rm=T) # effective growing season water supply
# 	
# 	
# 	### calculate mean climate for prior interval (antecedent)
# 	# end_yr <- start_yr - period_length
# 
# 	dat.plot.ant <- subset(clim_df, year < 2000 &
# 							year>= 1970 &
# 							plotid %in% p )
# 	
# # 	mean_ann_temps_ant <- quantile(dat.plot.ant$temp_ann_K,probs=c(0.5),na.rm=T)
# #   mean_ann_precips_ant <- quantile(dat.plot.ant$precip_ann_mm,probs=c(0.5),na.rm=T)
# # 	mean_ann_WD_ann_ant <- quantile(dat.plot.ant$WD_ann_mm,probs=c(0.5),na.rm=T)
# # 	mean_ann_WB_ant <- quantile(dat.plot.ant$water_balance_seasonal_mm,probs=c(0.5),na.rm=T)
# 	temp_k_mean_summer_ant <- mean(dat.plot.ant$temp_seasonal_K,na.rm=T)
# 	temp_k_mean_ann_ant <- mean(dat.plot.ant$temp_ann_K,na.rm=T)
# 	WD_mm_mean_summer_ann_ant <- mean(dat.plot.ant$WD_seasonal_mm,na.rm=T)
# #	WD_mm_mean_summer_ann_ant <- mean(dat.plot.ant$WD_seasonal_mm,na.rm=T)
# 	WB_mm_mean_seas_ant <- mean(dat.plot.ant$water_balance_seasonal_mm,na.rm=T)
# 
# 	
# 	### calculate climate anomalies
# 	temp_k_summer_diff <- temp_k_mean_summer_int - temp_k_mean_summer_ant 
# 	temp_k_ann_diff <- temp_k_mean_ann_int - temp_k_mean_ann_ant 
# 	WD_mm_mean_summer_diff <- WD_mm_mean_summer_int - WD_mm_mean_summer_ann_ant 
# 	WB_mm_mean_summer_diff <- WB_mm_mean_seas_int - WB_mm_mean_seas_ant 
# 
# 	# compile
# 	clim_array[[i]] <- data.frame("plotid"=p, temp_k_mean_ann_int,WD_mm_mean_ann_int,WB_mm_mean_seas_int,
# 	temp_k_mean_summer_int, WD_mm_mean_summer_int, SMI_mm_mean_summer_int, 
# 	temp_k_max_summer_int, WD_mm_max_summer_int, WB_mm_min_seas_int,
# 	temp_k_mean_summer_ant, WD_mm_mean_summer_ann_ant, WB_mm_mean_seas_ant,
# 	temp_k_summer_diff, WD_mm_mean_summer_diff, WB_mm_mean_summer_diff,
# 	temp_k_ann_diff, precip_mean_ann_int, wind_mean_ann_int)
# 	
# 	rm(start_yr,end_yr)
# 	
# }
# 
# warnings() 
# # convert to dataframe
# plot_clim <- do.call(rbind,clim_array)
# plot_clim$plotid <- as.factor(plot_clim$plotid)
# str(plot_clim)
# summary(plot_clim)
# 
# # save census interval climate
# save(plot_clim,file=file.path(Dir.Clim,"6_plot_mean_census_climate.Rdata"))
# 	
# 


# read tree dataset

targets <- load("~/Desktop/Prague/mort-likelihood/Data/full_trees.Rdata")
targets<-get(targets)
head(targets)
# name <- "A2_Balkan_trees_no_cores_code_v2"
# load(file=file.path(Dir.Data, paste0(name,".Rdata")))

# rename census date columns
colnames(targets)[names(targets) %in% "date"] <- "date_census_1"
colnames(targets)[names(targets) %in% "date_re"] <- "date_census_2"


##### combine census dates
survey_dates <- subset(targets,select=c(plotid,date_census_1,date_census_2))
survey_dates <- unique(survey_dates)

# remove NA
survey_dates <- survey_dates[!is.na(survey_dates$date_census_2),]
survey_dates <- survey_dates[order(survey_dates$plotid),]
which(duplicated(survey_dates$plotid))
summary(survey_dates)
write.csv(survey_dates,file=file.path(Dir.Data,"survey_dates.csv"),row.names=F)

# census interval summary
summary(survey_dates$date_census_2 - survey_dates$date_census_1)  # max = 7 years


################
# load annual-resolution climate data (plot scale)
clim_df <- load(file=file.path(Dir.Clim,"5_plot_monthly_water_budget_variables.Rdata"))
clim_df <- get(clim_df)

clim_array <- array(data=list(NULL),dim=nrow(survey_dates))

# set length of antecedent climate interval
period_length <- 30

# loop through plot and calculate quantiles
pb <- txtProgressBar(min=0, max=length(survey_dates$plotid), style=3)  # progress bar
tind=0

# loop through plots and compute mean climate for census and prior intervals
i <- 1
for(i in 1:length(survey_dates$plotid) )   
{
  tind=tind+1
  setTxtProgressBar(pb, tind)
  
  p <- survey_dates$plotid[i]
  
  # get census interval
  start_yr <- survey_dates$date_census_1[which(survey_dates$plotid %in% p)]
  end_yr <- survey_dates$date_census_2[which(survey_dates$plotid %in% p)]
  
  # subset climate data to census interval
  # dat.plot.int <- subset(clim_df, year>=start_yr &
  #                          year<=end_yr &
  #                          plotid %in% p )
  
  	# subset climate data to census interval
  	dat.plot.int <- subset(clim_df, year>=start_yr &
  								year<=(end_yr-5) &
  								plotid %in% p )
  
  ### calculate mean climate for census interval
  temp_k_mean_ann_int <- mean(dat.plot.int$temp_K_ann_ave,na.rm=T)   
  temp_k_mean_seas_int <- mean(dat.plot.int$temp_K_seasonal_ave,na.rm=T)  
  WD_mm_mean_ann_int <- mean(dat.plot.int$WD_mm_ann,na.rm=T) # water deficit (PET-AET)
  WD_mm_mean_seas_int <- mean(dat.plot.int$WD_mm_seasonal,na.rm=T) # water deficit (PET-AET)
  WB_mm_mean_seas_int <- mean(dat.plot.int$water_balance_seasonal_mm,na.rm=T) # effective growing season water supply
  SMI_mm_mean_seas_int <- mean(dat.plot.int$SMI_seasonal_ave,na.rm=T) # effective growing season water supply
  ppt_mm_ann_int<- mean(dat.plot.int$ppt_mm_ann_sum,na.rm=T)
  ppt_mm_seas_int<- mean(dat.plot.int$ppt_mm_seasonal_sum,na.rm=T)
  terra_WD_mm_mean_ann_int <- mean(dat.plot.int$terra_WD_mm_ann,na.rm=T) # water deficit (PET-AET)
  terra_WD_mm_mean_seas_int <- mean(dat.plot.int$terra_WD_mm_seasonal,na.rm=T) # water deficit (PET-AET)
  
  ### extreme census climate values  (See Canham & Murphy 2017)
  temp_k_max_seas_month_int <- max(dat.plot.int$temp_K_seasonal_month_max,na.rm=T)   # max monthly_level temp during growing season
  ppt_mm_min_seas_month_int <- min(dat.plot.int$ppt_mm_seasonal_month_min,na.rm=T)  
  SMI_min_seas_month_int <- min(dat.plot.int$SMI_seasonal_month_min,na.rm=T)  
  WD_mm_max_seas_int <- max(dat.plot.int$WD_mm_seasonal,na.rm=T)
  WB_mm_min_seas_int <- min(dat.plot.int$water_balance_seasonal_mm,na.rm=T)
  
  
  ### calculate mean climate for prior interval (antecedent)
  # start_yr <- end_yr
  # end_yr <- end_yr - period_length
  
  # dat.plot.ant <- subset(clim_df, year<start_yr &
  #                          year>=end_yr &
  #                          plotid %in% p )
  dat.plot.ant <- subset(clim_df, year < 2000 &
                           					year>= 1970 &
                            					plotid %in% p)
                           
  
  temp_k_mean_seas_ant <- mean(dat.plot.ant$temp_K_seasonal_ave,na.rm=T)
  temp_k_mean_ann_ant <- mean(dat.plot.ant$temp_K_ann_ave,na.rm=T)
  WD_mm_mean_seas_ant <- mean(dat.plot.ant$WD_mm_seasonal,na.rm=T)
  WD_mm_mean_ann_ant <- mean(dat.plot.ant$WD_mm_ann,na.rm=T)
  WB_mm_mean_seas_ant <- mean(dat.plot.ant$water_balance_seasonal_mm,na.rm=T)
  ppt_mm_ann_ant<- mean(dat.plot.ant$ppt_mm_ann_sum,na.rm=T)
  ppt_mm_seas_ant<- mean(dat.plot.ant$ppt_mm_seasonal_sum,na.rm=T)
  terra_WD_mm_mean_ann_ant <- mean(dat.plot.ant$terra_WD_mm_ann,na.rm=T) # water deficit (PET-AET)
  terra_WD_mm_mean_seas_ant <- mean(dat.plot.ant$terra_WD_mm_seasonal,na.rm=T) # water deficit (PET-AET)
  
  
  ### calculate climate anomalies
  # temp_k_summer_diff <- temp_k_mean_summer_int - temp_k_mean_summer_ant 
  # WD_mm_mean_summer_diff <- WD_mm_mean_summer_int - WD_mm_mean_summer_ann_ant 
  # WB_mm_mean_summer_diff <- WB_mm_mean_seas_int - WB_mm_mean_seas_ant 
  
  # compile
  clim_array[[i]] <- data.frame("plotid"=p,  temp_k_mean_ann_int,temp_k_mean_seas_int,WD_mm_mean_seas_int,WD_mm_mean_ann_int,WB_mm_mean_seas_int,SMI_mm_mean_seas_int,ppt_mm_ann_int,ppt_mm_seas_int,temp_k_max_seas_month_int,ppt_mm_min_seas_month_int,SMI_min_seas_month_int,WD_mm_max_seas_int,WB_mm_min_seas_int, temp_k_mean_seas_ant,temp_k_mean_ann_ant,WD_mm_mean_seas_ant,WD_mm_mean_ann_ant,WB_mm_mean_seas_ant,ppt_mm_ann_ant, ppt_mm_seas_ant,terra_WD_mm_mean_ann_ant,terra_WD_mm_mean_seas_ant, terra_WD_mm_mean_ann_int, terra_WD_mm_mean_seas_int)
  
  rm(start_yr,end_yr)
}

# convert to dataframe
plot_clim <- do.call(rbind,clim_array)

summary(plot_clim)

# save census interval climate
save(plot_clim,file=file.path(Dir.Clim,"6_plot_mean_census_climate.Rdata"))





	

	
	















