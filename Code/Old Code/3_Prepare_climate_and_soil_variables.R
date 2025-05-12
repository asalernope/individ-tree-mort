###############################################################################################
# Prepare climate and soil data for water balance input
# temperature in Kelvin
# calculate soil moisture index following -- Gao et al., 2016. 
#      Assessing various drought indicators in representing summer drought in boreal 
#      forests in Finland. Hydrology and Earth System Sciences, 20(1), pp.175-191.
# Dec 2024
# Feb 2025 - revised error in code - merge of df and awc corrected
###############################################################################################

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


install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
sapply(c("reshape2"), install.load.package)


### function to compute soil moisture index 
smi <- function(qsoil,wilt,sat) { (qsoil-wilt) / (sat-wilt) }  # soil moisture index (SMI)

### load plot soil AWC
load(file=file.path(Dir.Data,"A3_Plots_with_soil_AWC.Rdata"))
summary(plot_soil)


### load climate data
plot_clim <- load(file=file.path(Dir.Clim,"2_Plot_terra_climate_compiled_by_month_1970_to_2023.Rdata"))
plot_clim <- get(plot_clim)
rm(plot_clim_array)

summary(plot_clim)
which(duplicated(plot_clim))

# get monthly mean temp in kelvin
plot_clim$temp_mean_K <- ((plot_clim$tmax + plot_clim$tmin) / 2) + 273.15

# normalize soil moisture data using smi function
plot_clim$SMI <- smi(plot_clim$soil_mm, min(plot_clim$soil_mm), max(plot_clim$soil_mm) )


# Convert numeric month to abbreviated name
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
plot_clim$Month <- months[plot_clim$month]

# climate variables to process
Vars <- c("pet","ppt","temp_mean_K","SMI", "PDSI", "ws", "def")


### convert from long to wide format
index <- 1
for(index in 1:length(Vars))
{
  clim_var <- Vars[index]
  
  # Cast to wide format
  plot_clim_wide <- subset(plot_clim, select=c("plotid","year","Month",clim_var))
  plot_clim_wide <- dcast(plotid+year ~ Month, data=plot_clim_wide, value.var=clim_var)
  
  
  # Order columns
  plot_clim_wide <- subset(plot_clim_wide, select=c("plotid","year",months))
  
  # Rename columns
  new_nam <- paste(months,clim_var,sep="_")
  colnames(plot_clim_wide)[3:length(plot_clim_wide)] <- new_nam
  head(plot_clim_wide)
  
  # combine variables
  if(index == 1) 
  {
    df <- plot_clim_wide
  }else{
    df <- merge(df, plot_clim_wide, by=c("plotid","year"))
  }
}

# add soil AWC
plot_soil$plotid <- as.factor(plot_soil$plotid)

# add awc values
awc <-  plot_soil[,names(plot_soil) %in% c("plotid","AWC_mm_per_m2")]

df <- merge(df, awc, by="plotid", all.x=TRUE)
which(duplicated(df))


### save
save(df, file=file.path(Dir.Clim,"3_Working_climate_and_soil_variables.Rdata"))







