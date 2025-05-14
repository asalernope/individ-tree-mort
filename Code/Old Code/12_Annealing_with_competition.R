###########################################################
#  Disturbance - mortality relations in Balkans
#  Preliminary annealing runs for Fagus
#  Fit separate models for each focal species:
#      "ABIALB","FAGSYL","PICABI" and "ACEPSE"
#  Test competition effect 
#  21.04.25
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
Dir.Ndep <- file.path(Dir.Data,"Ndep")
Dir.Analysis.Files <- file.path(Dir.Data,"Full_data_files")
dirs <- sapply(c(Dir.Base,Dir.Data,Dir.Code,Dir.Clim,Dir.Dist,Dir.Ndep,Dir.Anneal,Dir.Analysis.Files), 
               function(x) if(!dir.exists(x)) dir.create(x))

library(likelihood)

# focal species
spp_list <- c("ABIALB","FAGSYL","PICABI","ACEPSE")

# select focal species
index <- 2
spname <- spp_list[index]


# get model fitting dataset
filename <- "Full_working_dataset_with_all_climate_variables.Rdata"
load(file=file.path(Dir.Analysis.Files, paste(spname,filename,sep="_")))
targets$stand <- as.factor(targets$stand)

# likelihood function
loglikelihood <- function(pred,observed)
{ ifelse(observed == 1, log(pred), log(1-pred)) }

# define separate results directories for each species 
Dir.Out <- file.path(Dir.Anneal, spname)
if(!dir.exists(Dir.Out)) { dir.create(Dir.Out) }


################################################################################
#  Annealing setup
################################################################################
setwd(Dir.Out)
comment.txt <- paste0("Input file: ",spname," ",filename,".Rdata")
initial_temp <- 5
tempred <- 0.9
iter <- 100  # number of iterations
################################################################################

#
# lognormal size and nci (based on taxonomic similarity)
# index lambda based on relatedness classes
# 
m_size_nci_spp <- function(PS,sizeX0,sizeXb,sizeXp,stand,lambda,alpha,beta,C,gamma,D)
{ 
	size.effect <- exp(-0.5*(log((targets$dbh_cm + sizeXp)/sizeX0)/sizeXb)^2) 
	nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
	comp.effect <- exp(-(C/100) * targets$size_ratio^gamma * nci^D)
	logit <- PS[stand] * size.effect * comp.effect
	logit/(1+logit)
}
# stand & lambda levels
ns <- length(levels(targets$stand))
nspp <- length(levels(neighbours[,names(neighbours) %in% "relatedness"]))

var <- list(pred="predicted", observed="mort", stand="stand")

# parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 0.1, sizeXb = 0.5, sizeXp = 0.5,
			alpha = 2, beta = 1, C=1, gamma = 1, D = 2, lambda = rep(0.5,nspp) )

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0,
               alpha = 0, beta = 0, C=0, gamma = -1, D = 1, lambda = rep(0,nspp) )

par_hi <- list(PS = rep(1000,ns), sizeX0 = 10000, sizeXb = 15000, sizeXp = 100,
				alpha = 5, beta = 5, C=1000, gamma = 5, D = 10, lambda = rep(1.0,nspp) )

m_size_nci_results <- anneal(model = m_size_nci_spp, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(m_size_nci_results, file = paste(spname,"model_size_nci_spp_relatedness_results.txt"), data=T, print_whole_hist=T)
save(m_size_nci_results, file=paste(spname,"model_size_nci_spp_relatedness_results.Rdata"))



#############################
# lognormal size and double logistic precip function
model_dblog_clim <- function(PS, sizeX0, sizeXb, sizeXp, stand, PREC,
                               al,bl,cl,ah,bh,ch,lambda,alpha,beta,C,gamma,D)
{ 
  size.effect <- exp(-0.5*(log((targets$dbh_cm + sizeXp)/sizeX0)/sizeXb)^2) 
  prec.effect <- (al + ( (1-al)/(1+(bl/PREC)^cl))) * (ah + ( (1-ah)/(1+(PREC/bh)^ch)))
  nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
  comp.effect <- exp(-(C/100) * targets$size_ratio^gamma * nci^D)
  logit <- PS[stand] * size.effect * prec.effect * comp.effect
  logit/(1+logit)
}
# stand levels
ns <- length(levels(targets$stand))
nspp <- length(levels(neighbours[,names(neighbours) %in% "relatedness"]))

var <- list(pred="predicted", observed="mort", stand="stand",
            PREC="WD_mm_mean_ann_int")

# parameter limits
par <- list(PS = rep(50,ns), sizeX0 = 9000, sizeXb = 960, sizeXp = 99,
			al = 0.8, bl = 1040, cl = 19,
			ah = 0.3, bh = 3600, ch = 360,
			alpha = 0.3, beta = 0.32, C=0.30, gamma = 1.14, D = 1.13, lambda = rep(0.5,nspp))

par_lo <- list(PS = rep(10,ns), sizeX0 = 0, sizeXb = 0, sizeXp = 0,
				al = 0, bl = 0, cl = 0,
				ah = 0.3, bh = 1000, ch = 0,
				alpha = 0, beta = 0, C=0, gamma = -1, D = 1, lambda = rep(0,nspp))

par_hi <- list(PS = rep(100,ns), sizeX0 = 12000, sizeXb = 5000, sizeXp = 200,
				al = 1.0, bl = 2000, cl = 500,
				ah = 0.3, bh = 5000, ch = 500,
				alpha = 5, beta = 5, C=1000, gamma = 5, D = 10, lambda = rep(1.0,nspp))

model_dblog_clim_results <- anneal(model = model_dblog_clim, par = par, var = var,source_data = targets,
                             par_lo = par_lo, par_hi = par_hi, loglikelihood,
                             dep_var = "mort", max_iter = iter, note = comment.txt)

# save
write_results(model_dblog_clim_results, file = paste(spname,"model_dblog_clim_results.txt"), data=T, print_whole_hist=T)
save(model_dblog_clim_results, file=paste(spname,"model_dblog_clim_results.Rdata"))


##### plot results
load(file=file.path(Dir.Out, paste(spname,"model_dblog_clim_results.Rdata")))

pars <- model_dblog_clim_results$best_par
attach(pars)
targets <- model_dblog_clim_results$source_data

# water deficit observed range
prec_min <- floor(min(targets$WD_mm_mean_ann_int))
prec_max <- round(ceiling(max(targets$WD_mm_mean_ann_int)),-1)
prec_vector <- seq(prec_min,prec_max,.1)

# effect on survival
r <- (pars$al + ( (1-pars$al)/(1+(pars$bl/prec_vector)^pars$cl))) * (pars$ah + ( (1-pars$ah)/(1+(prec_vector/pars$bh)^pars$ch)))

# plot
#quartz(height=6,width=6)
plot(prec_vector,r, type="l",las=1, lwd=2, col="blue", ylim=c(0,1))

warnings()




# water deficit observed range
prec_min <- floor(min(targets$WD_mm_mean_ann_int))
prec_max <- round(ceiling(max(targets$WD_mm_mean_ann_int)),-1)
prec_vector <- seq(prec_min,prec_max,.1)

# effect on survival
r <- (pars$al + ( (1-pars$al)/(1+(pars$bl/prec_vector)^pars$cl))) * (pars$ah + ( (1-pars$ah)/(1+(prec_vector/pars$bh)^pars$ch)))

# plot
#quartz(height=6,width=6)
plot(prec_vector,r, type="l",las=1, lwd=2, col="blue", ylim=c(0,1))

warnings()



































