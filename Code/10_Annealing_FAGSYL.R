###########################################################
#  Disturbance - mortality relations in Balkans
#  Preliminary annealing runs for Fagus
#  Fit separate models for each focal species:
#      "ABIALB","FAGSYL","PICABI",and "ACEPSE"
#  12.02.25
###########################################################
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
Dir.Anneal <- file.path(Dir.Base,"Annealing_results")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim, Dir.Dist, Dir.Anneal), function(x) if(!dir.exists(x)) dir.create(x))

library(likelihood)
library(tidyverse)
library(ggplot2)

# get model fitting dataset (all species)
filename <- "9_Full_model_fitting_dataset"
load(file=file.path(Dir.Data,paste0(filename,".Rdata")))


# list of focal species
spp_list <- c("ABIALB","FAGSYL","PICABI","ACEPSE")

# select focal species for model runs
index <- 2
spp <- spp_list[index]

targets <- subset(trees, SPCD == spp )   
targets <- droplevels(targets)

# ggplot(targets) +
#   geom_point(aes(x = precip_mean_ann_int, y = temp_k_mean_ann_int, color = lord), alpha = 0.6)+
#   theme(legend.position = "none")

# targets <- targets %>%
#   filter(!location %in% c("Stara Planina", "Seminic"))

summary(targets)
str(targets)
check <- subset(targets,select=c(treeid))
which(duplicated(check))


# number of stands and strata
targets$stand <- as.factor(targets$stand)
ns <- length(levels(targets$stand))

targets$stratum <- as.factor(targets$stratum)
n.strata <- length(levels(targets$stratum))


# likelihood function
loglikelihood <- function(pred,observed)
{ ifelse(observed == 1, log(pred), log(1-pred)) }


# define separate results directories for each species 
Dir.Out <- file.path(Dir.Anneal, spp)
if(!dir.exists(Dir.Out)) { dir.create(Dir.Out) }


################################################################################
#  Annealing setup
################################################################################
setwd(Dir.Out)
comment.txt <- paste0("Input file: ",spp," ",filename,".Rdata")
initial_temp <- 5
tempred <- 0.9
iter <- 50000  # number of iterations
#
#  Note: All variables other than observed and predicted and stand are accessed 
#    directly from the targets dataset in the global environment
################################################################################

# null means model (separate means for each site)
null_model <- function (PS, stand) { 
					  logit <- PS[stand]
					  logit/(1+logit) }

# number of sites = stands
ns = length(levels(targets$stand))

var <- list(pred = "predicted", observed = "mort", stand="stand")

# set parameter limits
par <- list(PS = rep(10, ns))
par_lo <- list(PS = rep(0.001, ns))
par_hi <- list(PS = rep(1000, ns))

means_site_results <- anneal(model = null_model, par = par, var = var, source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

write_results(means_site_results, paste(spp,"Null means results.txt"), data=F, print_whole_hist=F)
save(means_site_results,file=paste(spp,"Null means results.Rdata"))



##############################
# lognormal size model
model_1_size <- function(PS,sizeX0,sizeXb,sizeXp,stand)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp)/sizeX0)/sizeXb)^2) 
  logit <- (PS[stand] * size.effect)
  logit/(1+logit)
}


var <- list(pred="predicted", observed="mort", stand="stand")

# set parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 217, sizeXb = 0.9, sizeXp = 0.945)

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0)

par_hi <- list(PS = rep(100,ns), sizeX0 = 10000, sizeXb = 15000, sizeXp = 100)

model_1_size_results <- anneal(model = model_1_size, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_1_size_results, file = paste(spp,"model_1_size_results.txt"), data=T, print_whole_hist=T)
save(model_1_size_results, file=paste(spp,"model_1_size_results.Rdata"))

targets <- targets %>% !drop_na()

##############################
# model 2: mean census climate
# lognormal size and climate function
model_2_log_clim <- function(PS, sizeX0, sizeXb, sizeXp, stand, TEMP, PREC,WIND,
                         temp.X0, temp.Xb, temp.Xp, prec.X0, prec.Xb, prec.Xp,wind.Xp, wind.X0,wind.Xb)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp)/sizeX0)/sizeXb)^2) 
  temp.effect <- exp(-0.5*(log((TEMP + temp.Xp)/temp.X0)/temp.Xb)^2)
  prec.effect <- exp(-0.5*(log((PREC + prec.Xp)/prec.X0)/prec.Xb)^2)
  wind.effect <- exp(-0.5*(log((WIND + wind.Xp)/wind.X0)/wind.Xb)^2)
  logit <- PS[stand] * size.effect * temp.effect * prec.effect * wind.effect
  logit/(1+logit)
}


var <- list(pred="predicted", observed="mort", stand="stand",
            TEMP="temp_k_mean_summer_int", PREC="WD_mm_mean_summer_int", WIND = "wind_mean_ann_int")

# set parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 439, sizeXb = 0.598, sizeXp = 98,
			temp.X0 = 336, temp.Xb = 408, temp.Xp = 74.4,
			prec.X0 = 5.42, prec.Xb = 51.4, prec.Xp = 0.0815,
			wind.X0 = 5.42, wind.Xb = 5, wind.Xp = 0.0815)

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0,
				temp.X0 = 0, temp.Xb = 0.1, temp.Xp = 0,
				prec.X0 = 0, prec.Xb = 0.1, prec.Xp = 0,
				wind.X0 = 0, wind.Xb = 0.1, wind.Xp = 0)

par_hi <- list(PS = rep(100,ns), sizeX0 = 10000, sizeXb = 5000, sizeXp = 100,
				temp.X0 = 500, temp.Xb = 10000, temp.Xp = 100,
				prec.X0 = 500, prec.Xb = 10000, prec.Xp = 100,
				wind.X0 = 500, wind.Xb = 10000, wind.Xp = 100)

model_2_clim_results <- anneal(model = model_2_log_clim, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_2_clim_results, file = paste(spp,"model_2_clim_results.txt"), data=T, print_whole_hist=T)
save(model_2_clim_results, file=paste(spp,"model_2_clim_results.Rdata"))
warnings()


##############################
# Gaussian climate terms
model_3_G_clim <- function(PS, sizeX0, sizeXb, sizeXp, stand, TEMP, PREC,
                         temp.a, temp.b, temp.c, prec.a, prec.b, prec.c)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp)/sizeX0)/sizeXb)^2) 
  temp.effect <- temp.a * exp(-0.5*((TEMP - temp.c)/temp.b)^2)
  prec.effect <- prec.a * exp(-0.5*((PREC - prec.c)/prec.b)^2)

  logit <- PS[stand] * size.effect * temp.effect * prec.effect
  logit/(1+logit)
}

targets <-  targets %>%
filter(!is.na(WD_mm_mean_summer_int))

var <- list(pred="predicted", observed="mort", stand="stand",
            TEMP="temp_k_mean_summer_int", PREC="WD_mm_mean_summer_int")

# set parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 903, sizeXb = 0.313, sizeXp = 500,
			temp.a = 1, temp.b = 5793, temp.c = mean(targets$temp_k_mean_summer_int),
			prec.a = 0.473, prec.b = 19941, prec.c = mean(targets$WD_mm_mean_summer_int)  )

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0,
				temp.a = 0, temp.b = 0.5, temp.c = min(targets$temp_k_mean_summer_int),
				prec.a = 0, prec.b = 0.5, prec.c = min(targets$WD_mm_mean_summer_int) )

par_hi <- list(PS = rep(100,ns), sizeX0 = 10000, sizeXb = 15000, sizeXp = 1000,
				temp.a = 1, temp.b = 10000, temp.c = max(targets$temp_k_mean_summer_int),
				prec.a = 1, prec.b = 20000, prec.c = max(targets$WD_mm_mean_summer_int) )

model_3_G_clim_results <- anneal(model = model_3_G_clim, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_3_G_clim_results, file = paste(spp,"model_3_G_clim_results.txt"), data=T, print_whole_hist=T)
save(model_3_G_clim_results, file=paste(spp,"model_3_G_clim_results.Rdata"))



##############################
# Asymmetric exponential climate terms  -- 
# this model is not running yet - there is an error in the model code somewhere
# Canham & Murphy 2018
model_4_asy_exp_climate <- function(PS, stand, sizeX0,sizeXb,sizeXp, 
									tempX0,tempXbl,tempXbh,tempXc,
									precX0,precXbl,precXbh,precXc )
{
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp)/sizeX0)/sizeXb)^2)  

  temp.effect <- ifelse(targets$temp_k_mean_summer_int < tempXc, 
					    tempX0 * (tempXbl/100)^((targets$temp_k_mean_summer_int - tempXc)^2), 
                        tempX0 * (tempXbh/100)^((targets$temp_k_mean_summer_int - tempXc)^2) )
  
  prec.effect <- ifelse(targets$WB_mm_mean_seas_int < precXc, 
					    precX0 * (precXbl/100)^((targets$WB_mm_mean_seas_int - precXc)^2), 
                        precX0 * (precXbh/100)^((targets$WB_mm_mean_seas_int - precXc)^2) )
  
  logit <- PS[stand] * size.effect * temp.effect * prec.effect
  logit/(1+logit)
}

var <- list(pred="predicted", observed="mort", stand="stand")

# get previous model output
res <- load(file=file.path(Dir.Out,paste(spp,"model_1_size_results.Rdata")))
res <- get(res)

par <- list(PS = res$best_pars$PS,
			sizeX0 = res$best_pars$sizeX0,
			sizeXb = res$best_pars$sizeXb,
			sizeXp = res$best_pars$sizeXp,
			tempX0 = 0.5, tempXbl = 0.5, tempXbh = 0.5, tempXc = 5,
			precX0 = 0.5, precXbl = 0.5, precXbh = 0.5, precXc = 1 )

par_lo <- list(PS = rep(0, ns), 
              sizeX0 = 0.001, sizeXb = 0.1, sizeXp = 0, 
			  tempX0 = 0, tempXbl = 0.1, tempXbh = 0.1, tempXc = 0.1, 
			  precX0 = 0, precXbl = 0.1, precXbh = 0.1, precXc = 0.1 )
			  
par_hi <- list(PS = rep(100, ns), 
			  sizeX0 = 1000, sizeXb = 15000, sizeXp = 1000, 
			  tempX0 = 1000, tempXbl = 1, tempXbh = 1, tempXc = 1000, 
			  precX0 = 1000, precXbl = 1, precXbh = 1, precXc = 1000 )

model_4_asy_exp_climate_results <- anneal(model = model_4_asy_exp_climate, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_4_asy_exp_climate_results, file = paste(spp,"model_4_asy_exp_climate_results.txt"), data=T, print_whole_hist=T)
save(model_4_asy_exp_climate_results, file=paste(spp,"model_4_asy_exp_climate_results.Rdata"))




##############################
# N deposition and climate
model_5_Ndep <- function(PS, sizeX0, sizeXb, sizeXp, stand, TEMP, PREC,
                         temp.a, temp.b, temp.c, prec.a, prec.b, prec.c,
                         N0, Nb, Np)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp)/sizeX0)/sizeXb)^2) 
  temp.effect <- temp.a * exp(-0.5*((TEMP - temp.c)/temp.b)^2)
  prec.effect <- prec.a * exp(-0.5*((PREC - prec.c)/prec.b)^2)
  Ndep.effect <- exp(-0.5*(log((targets$Ndep_tot_kg_ha_yr_int + Np)/N0)/Nb)^2) 
  
  logit <- PS[stand] * size.effect * temp.effect * prec.effect * Ndep.effect
  logit/(1+logit)
}


var <- list(pred="predicted", observed="mort", stand="stand",
            TEMP="temp_k_mean_summer_int", PREC="WD_mm_mean_summer_int")

# set parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 0.1, sizeXb = 0.5, sizeXp = 0.5,
			temp.a = 0.5, temp.b = 10, temp.c = mean(targets$temp_k_mean_summer_int),
			prec.a = 0.5, prec.b = 10, prec.c = mean(targets$WD_mm_mean_summer_int),  
			N0 = 1, Nb = 1, Np = 1)

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0,
				temp.a = 0, temp.b = 0.5, temp.c = min(targets$temp_k_mean_summer_int),
				prec.a = 0, prec.b = 0.5, prec.c = min(targets$WD_mm_mean_summer_int),
				N0 = 1, Nb = 0.1, Np = 0 )

par_hi <- list(PS = rep(100,ns), sizeX0 = 5000, sizeXb = 10000, sizeXp = 500,
				temp.a = 1, temp.b = 15000, temp.c = max(targets$temp_k_mean_summer_int),
				prec.a = 1, prec.b = 15000, prec.c = max(targets$WD_mm_mean_summer_int),
				N0 = 1000, Nb = 5000, Np = 500)

model_5_Ndep_results <- anneal(model = model_5_Ndep, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_5_Ndep_results, file = paste(spp,"model_5_Ndep_results.txt"), data=T, print_whole_hist=T)
save(model_5_Ndep_results, file=paste(spp,"model_5_Ndep_results.Rdata"))



#################################
# disturbance and Gaussian climate
model_6_dist <- function(PS, sizeX0, sizeXb, sizeXp, stand, TEMP, PREC,
                         temp.a, temp.b, temp.c, prec.a, prec.b, prec.c,
                         Da, Db, Dc)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp)/sizeX0)/sizeXb)^2) 
  temp.effect <- temp.a * exp(-0.5*((TEMP - temp.c)/temp.b)^2)
  prec.effect <- prec.a * exp(-0.5*((PREC - prec.c)/prec.b)^2)
  dist.effect <- exp(-Da/100 * targets$time_since_dist_max^Db * targets$my_max_sev^Dc )
  
  logit <- PS[stand] * size.effect * temp.effect * prec.effect * dist.effect
  logit/(1+logit)
}


var <- list(pred="predicted", observed="mort", stand="stand",
            TEMP="temp_k_mean_summer_int", PREC="WD_mm_mean_summer_int")

# set parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 0.1, sizeXb = 0.5, sizeXp = 0.5,
			temp.a = 0.5, temp.b = 10, temp.c = mean(targets$temp_k_mean_summer_int),
			prec.a = 0.5, prec.b = 10, prec.c = mean(targets$WD_mm_mean_summer_int),  
			Da = 1, Db = 1, Dc = 1)

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0,
				temp.a = 0, temp.b = 0.5, temp.c = min(targets$temp_k_mean_summer_int),
				prec.a = 0, prec.b = 0.5, prec.c = min(targets$WD_mm_mean_summer_int),
				Da = 1, Db = -4, Dc = -4 )

par_hi <- list(PS = rep(100,ns), sizeX0 = 5000, sizeXb = 10000, sizeXp = 500,
				temp.a = 1, temp.b = 15000, temp.c = max(targets$temp_k_mean_summer_int),
				prec.a = 1, prec.b = 15000, prec.c = max(targets$WD_mm_mean_summer_int),
				Da = 1000, Db = 4, Dc = 4)

model_6_dist_results <- anneal(model = model_6_dist, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_6_dist_results, file = paste(spp,"model_6_dist_results.txt"), data=T, print_whole_hist=T)
save(model_6_dist_results, file=paste(spp,"model_6_dist_results.Rdata"))




##########################################################
# Climate and size effects differentiated by stratum
model_7_G_clim_strata <- function(PS, sizeX0, sizeXb, sizeXp, stand, TEMP, PREC, stratum,
                         temp.a, temp.b, temp.c, prec.a, prec.b, prec.c)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_mm + sizeXp[stratum])/sizeX0[stratum])/sizeXb[stratum])^2) 
  temp.effect <- temp.a[stratum] * exp(-0.5*((TEMP - temp.c[stratum])/temp.b[stratum])^2)
  prec.effect <- prec.a[stratum] * exp(-0.5*((PREC - prec.c[stratum])/prec.b[stratum])^2)

  logit <- PS[stand] * size.effect * temp.effect * prec.effect
  logit/(1+logit)
}

ns <- length(levels(targets$stand))
n.strata <- length(levels(targets$stratum))

var <- list(pred="predicted", observed="mort", stand="stand", stratum="growth",
            TEMP="temp_k_mean_summer_int", PREC="WD_mm_mean_summer_int")

# set parameter limits
par <- list(PS = rep(50,ns), sizeX0 = rep(1,n.strata), sizeXb = rep(1,n.strata), sizeXp = rep(0.5,n.strata),
			temp.a = rep(0.5,n.strata), temp.b = rep(5000,n.strata), temp.c = rep(mean(targets$temp_k_mean_summer_int), n.strata),
			prec.a = rep(0.7,n.strata), prec.b = rep(1000,n.strata), prec.c = rep(mean(targets$WD_mm_mean_summer_int), n.strata) )

par_lo <- list(PS = rep(10,ns), sizeX0 = rep(0.1,n.strata), sizeXb = rep(0.1,n.strata), sizeXp = rep(0,n.strata),
				temp.a = rep(0,n.strata), temp.b = rep(0.5,n.strata), temp.c = rep(min(targets$temp_k_mean_summer_int), n.strata),
				prec.a = rep(0,n.strata), prec.b = rep(0.5,n.strata), prec.c = rep(min(targets$WD_mm_mean_summer_int), n.strata) )

par_hi <- list(PS = rep(100,ns), sizeX0 = rep(1000,n.strata), sizeXb = rep(5000,n.strata), sizeXp = rep(500,n.strata),
				temp.a = rep(1,n.strata), temp.b = rep(10000,n.strata), temp.c = rep(max(targets$temp_k_mean_summer_int),n.strata),
				prec.a = rep(1,n.strata), prec.b = rep(10000,n.strata), prec.c = rep(max(targets$WD_mm_mean_summer_int),n.strata) )

model_7_G_clim_strata_results <- anneal(model = model_7_G_clim_strata, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)


warnings()
# save
write_results(model_7_G_clim_strata_results, file = paste(spp,"model_7_G_clim_with_strata_results.txt"), data=T, print_whole_hist=T)
save(model_7_G_clim_strata_results, file=paste(spp,"model_7_G_clim_with_strata_results.Rdata"))


mort <- function(PS,sizeI,sizeXa,sizeXb,Da,Db,Dc,Ca,Cb,Cc)
  
{ 
  size <- exp(-0.5*(log(data$dbh_mm+sizeI/sizeXa[status])/sizeXb[status])^2) 
  dist <- Da*(data$time_since^Db*data$severity^Dc)
  clim <- Ca*Cb^((data$climate-Cc)^2)
  
  logit <- PS[stand] * size * dist * clim
  logit/(1+logit)
}









