###############################################################
#   Get soil water holding capacity (AWC) for plot locations
#   Download from http://globalchange.bnu.edu.cn/research/soilw
#   Download date = 4 June 2019
#	netCDF format
# 	30 arc-seconds (~1.0 km)
#	Global extent
#	Categorical - 6 ordinal classes with midpoints (FAO 1995)
#	Ref: Shangguan, W., Dai, Y., Duan, Q., Liu, B. and Yuan, H., 
#        2014. A global soil data set for earth system modeling. 
#        Journal of Advances in Modeling Earth Systems, 6(1)
#   Dec 2024
#   Feb 2025 rerun with revised tree dataset
###############################################################

rm(list=ls())
#  set user computer
computer <- "~"

#  directores
#  directores
if(computer == "~"){
  Dir.Base <- "~/Desktop/Prague/mort-likelihood/individ-tree-mort"
}else{
  Dir.Base <- "insert path"
}

Dir.Data <- file.path(Dir.Base, "Data")
Dir.Code <- file.path(Dir.Base, "Code")
Dir.Clim <- file.path(Dir.Data,"Climate")
Dir.Soil <- file.path(Dir.Data, "Soil")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim, Dir.Soil), function(x) if(!dir.exists(x)) dir.create(x) )


# need raster function
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
sapply(c("stringr","raster","ncdf4"), install.load.package)


# load tree data
targets <- load("~/Desktop/Prague/mort-likelihood/Data/full_trees_sp.Rdata")
targets<-get(targets)
# str(plots)
summary(targets)
str(targets)

### subset coordinates and remove NAs
plots <- unique(targets[,names(targets) %in% c("plotid","plot_id","lng","lat")] )
plots <- as.data.frame(plots[!is.na(plots$lng),])
plots_copy <- plots


# Set directory for downloaded gridded soil data
#Dir.Soil <- file.path("/Users/arnebuechling/Documents/Environmental_data/Soil_EU")
setwd(Dir.Soil)

# Create list of model file names
soil_list <- list.files(pattern=c('\\.nc$'), recursive=F)

# Create raster stack
rasStack <- stack(soil_list)  

# Get coordinate system of climate grids 
soil_crs <- crs(rasStack)

# Define as spatial object
coordinates(plots) <- ~ lng + lat
proj4string(plots) <- crs(rasStack)

# Get values using bilinear interpolation (4 nearest cells)
rasValue = extract(rasStack, plots, method='bilinear')

# Combine with plotid
plot_soil <- cbind(plots_copy,rasValue)
sort(unique(plot_soil$available.water.capacity))

# classes defined as follows:
plot_soil$AWC_mm_per_m2 <- ifelse(plot_soil$available.water.capacity <= 1.5, 150,
				 ifelse(plot_soil$available.water.capacity > 1.5 & plot_soil$available.water.capacity <= 2.5 , 125,
				   ifelse(plot_soil$available.water.capacity > 2.5 & plot_soil$available.water.capacity <= 3.5, 100,
				     ifelse(plot_soil$available.water.capacity > 3.5 & plot_soil$available.water.capacity <= 4.5, 75,
				       ifelse(plot_soil$available.water.capacity > 4.5 & plot_soil$available.water.capacity <= 5.5, 50,
				         ifelse(plot_soil$available.water.capacity > 5.5 & plot_soil$available.water.capacity <= 6.5, 15, 
				           10))))))

str(plot_soil)
unique(plot_soil$AWC_mm_per_m)

#### Save
save(plot_soil, file=file.path(Dir.Data,"A3_Plots_with_soil_AWC.Rdata"))










