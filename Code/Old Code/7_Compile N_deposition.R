############################################################################
# Disturbance - mortality modeling
# Compile nitrogen and sulphur deposition values for each plot
# Calculate total N (oxidized + reduced forms)
# Calculate census interval mean annual values
# Data from: Engardt, M, and others (2017) Deposition of sulphur and nitrogen in Europe 1900–2050. Model calculations and comparison
#         to historical observations, Tellus B: Chemical and Physical Meteorology, 69:1, 1328945, DOI: 10.1080/16000889.2017.1328945
# 12.02.25
############################################################################

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
Dir.Out <- file.path(Dir.Base,"Annealing_output")
Dir.Ndep <- file.path(Dir.Data,"Ndep")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim, Dir.Dist, Dir.Ndep, Dir.Out), function(x) if(!dir.exists(x)) dir.create(x))



### read tree dataset
targets <- load("~/Desktop/Prague/mort-likelihood/Data/full_trees_sp.Rdata")
targets<-get(targets)
head(targets)
summary(targets)


# rename census date columns
colnames(targets)[names(targets) %in% "date"] <- "date_census_1"
colnames(targets)[names(targets) %in% "date_re"] <- "date_census_2"

# unique plot_ids
plot_ids <- unique(subset(targets,select=c(plotid,plot_id)) )
nrow(plot_ids)
which(duplicated(plot_ids))

# combine census dates
survey_dates <- subset(targets,select=c(plotid,date_census_1,date_census_2))
survey_dates <- unique(survey_dates)
survey_dates[which(duplicated(survey_dates$plotid)),]

# remove NA
survey_dates <- survey_dates[!is.na(survey_dates$date_census_2),]
survey_dates <- survey_dates[order(survey_dates$plotid),]
nrow(survey_dates)


# get N deposition data (from Engardt et al. 2017)
setwd(Dir.Ndep)
Ndep <- load(file="Engardt 2017 Unfiltered raw model estimates for N and S deposition rates 1900 to 2050 all REMFOR plots revised 1.Rdata")
Ndep <- get(Ndep)
colnames(Ndep)[names(Ndep) %in% "Plot_No"] <- "plot_id"
colnames(Ndep)[names(Ndep) %in% "Year"] <- "year"
head(Ndep)


# Convert from mg/m2 to kg/ha
Ndep[,c(2,3,4)] <- Ndep[,c(2,3,4)] * 100*100/(1000*1000)
colnames(Ndep)[names(Ndep) %in% "SOX_S"] <- "SOx_kg_ha_yr"
colnames(Ndep)[names(Ndep) %in% "NOY_N"] <- "NOx_Ndep_dep_kg_ha_yr"  # oxidized N
colnames(Ndep)[names(Ndep) %in% "NHX_N"] <- "NHx_Ndep_dep_kg_ha_yr"  # reduced N
Ndep <- Ndep[,!names(Ndep) %in% "PREC"]


# keep Balkan plots
Ndep_Balkan <- merge(Ndep,plot_ids,by="plot_id")
Ndep_Balkan <- droplevels(Ndep_Balkan)
Ndep_Balkan <- Ndep_Balkan[order(Ndep_Balkan$plotid,Ndep_Balkan$year),]

# calculate total N deposition (oxidized + reduced forms)
Ndep_Balkan$Ndep_tot_kg_ha_yr <- Ndep_Balkan$NOx_Ndep_dep_kg_ha_yr + Ndep_Balkan$NHx_Ndep_dep_kg_ha_yr

head(Ndep_Balkan)
length(levels((Ndep_Balkan$plot_id)))
setdiff(survey_dates$plotid, Ndep_Balkan$plotid)


###########
# loop through plots and calculate mean values per census interval
Ndep_array <- array(data=list(NULL),dim=nrow(survey_dates))

# set length of antecedent climate interval
period_length <- 30

# loop through plot and calculate quantiles
pb <- txtProgressBar(min=0, max=length(survey_dates$plotid), style=3)  # progress bar
tind=0

# loop through plots and compute mean value for each census interval
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
	dat.plot.int <- subset(Ndep_Balkan, year >= start_yr &
								        year <= end_yr &
								        plotid %in% p )
	
	SOx_kg_ha_yr_int <- mean(dat.plot.int$SOx_kg_ha_yr,na.rm=T)
	Ndep_tot_kg_ha_yr_int <- mean(dat.plot.int$Ndep_tot_kg_ha_yr, na.rm=T)
	
	Ndep_array[[i]] <- data.frame("plotid"=p, SOx_kg_ha_yr_int, Ndep_tot_kg_ha_yr_int)
    
    rm(start_yr,end_yr)
}

# convert to dataframe
plot_Ndep <- do.call(rbind,Ndep_array)
plot_Ndep$plotid <- as.factor(plot_Ndep$plotid)
str(plot_Ndep)
summary(plot_Ndep)

# save census interval climate
save(plot_Ndep,file=file.path(Dir.Data,"7_plot_mean_census_nitrogen_deposition.Rdata"))











