##############################################################
# Download TerraClimate
# http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
# netcdf file format
# monthly resolution
# 0.04° (~4 km) grid
# period = 1958-2023
# download global files and clip to study area
# select variables = "pet","ppt","tmax","tmin","soil" and ...
# Dec 2024
# Feb 2015: revised tree dataset comprised of same plotids
#           no need to rerun this step
##############################################################

rm(list=ls())

# user computer
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
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Clim), function(x) if(!dir.exists(x)) dir.create(x) )


# install packages
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

#update.packages("raster",repos = "http://cran.us.r-project.org")
#update.packages("terra",repos = "http://cran.us.r-project.org")
sapply(c ("stringr","raster","terra"), install.load.package)
	

# get Remote plot dataset
#plots <- load(file=file.path(Dir.Data,"A1_Balkan_trees_no_cores.Rdata"))

#plots <- load(file=file.path(Dir.Data,"full_trees_sp.Rdata")) # revised expanded dataset
plots <- load("~/Desktop/Prague/mort-likelihood/Data/full_trees_sp.Rdata")
plots<-get(plots)
str(plots)



# get latitude / longitude limits for selected plots
xmin <- floor(min(plots$lng, na.rm = TRUE) ) - 1 
xmax <- ceiling(max(plots$lng, na.rm = TRUE)) + 1
ymin <- floor(min(plots$lat,na.rm = TRUE)) - 1
ymax <- ceiling(max(plots$lat,na.rm = TRUE)) + 1

# define extent
Extent <- extent(xmin,xmax,ymin,ymax)
# raster
rast <- raster(Extent, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# shapefile
region_shp <- rasterToPolygons(rast,dissolve=T)


#################
# download
dataset <- "TerraClimate"

# climate variables (soil=soil moisture)
Vars <- c("PDSI", "ws") #x"ppt","tmax","tmin","soil", "pet",

# define dates
dates <- seq(1970,2023,1)

mo_year_names <- paste( str_pad(string = 1:12, width = 2, "left", 0),rep(1970:2023, each = 12), sep="_")


## Download from web
options(timeout = max(1000, getOption("timeout"))) # increase max allowable download period
getOption("timeout")

Start <- Sys.time()

i <- 2
for(i in 1:length(Vars))
{
  Dir.Iter <- file.path(Dir.Clim, Vars[i])
  if(!dir.exists(Dir.Iter)) {dir.create(Dir.Iter)}
  k <- 1
  for(k in 1:length(dates))
  {
    URL <- paste0("https://climate.northwestknowledge.net/TERRACLIMATE-DATA/",dataset,"_", Vars[i],"_", dates[k],".nc")

    if(!file.exists(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".nc")))){
       download.file(URL, destfile = file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".nc")), mode="wb")
    }
  }
  # stack monthly data
  setwd(Dir.Iter); dir()
  BRICK <- stack(list.files(Dir.Iter, pattern = ".nc"))
  names(BRICK) <- mo_year_names
  
  # crop to extent of Remote plots
  setwd(Dir.Clim)
  region_subset <- crop(BRICK, region_shp)
  region_subset <- mask(region_subset, region_shp)
  
  # write
  writeRaster(x = region_subset, filename = file.path(Dir.Clim, paste("1",dataset, Vars[i], "all_years_region", sep="_")), format = "CDF")
  
  # delete folder and global monthly raster files
  unlink(Dir.Iter, recursive = TRUE) 

}
  
End <- Sys.time()
Process.Time <- End - Start








