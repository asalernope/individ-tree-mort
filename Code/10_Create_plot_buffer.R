########################################################################################
# Create plot buffer with simulated trees for competition analyses
# Buffer area = 3x actual plot area
# Select trees from actual plot and generate random x and y coordinates within the buffer
# Mirror composition and size distribution of the measured within-plot trees 
# 15.04.25
########################################################################################


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

# packages
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec <- c("RPostgreSQL","pool","tidyverse")
sapply(package_vec, install.load.package)


# Distance function
neighdist<-function(targetx, targety, neighborx, neighbory) {
                    sqrt((targetx-neighborx)^2 + (targety-neighbory)^2) }

# load neighbour tree data
neighbours_df <-load("~/Desktop/Prague/mort-likelihood/Data/full_trees_comp.Rdata")
neighbours_df <- get(neighbours_df)
neighbours_df$plotid <- as.factor(neighbours_df$plotid)


# read species table

spp_table <- read_csv("~/Desktop/Prague/mort-likelihood/Data/Species_table.csv")
spp_table <- spp_table[,!names(spp_table) %in% c("Genus_nam","Species_nam","functional.grp","sp_type")]


# add species code
neighbours_df <- merge(neighbours_df,spp_table,by="species",all.x=T)
neighbours_df$SPCD <- as.factor(neighbours_df$SPCD)
levels(neighbours_df$SPCD)

# Calculate neighbour distance to plot centre
neighbours_df$DIST_m <- neighdist(neighbours_df$x_m, neighbours_df$y_m, 0, 0)

# dbh in cm
neighbours_df$dbh_cm <- neighbours_df$dbh_mm / 10
neighbours_df <- neighbours_df[,!names(neighbours_df) %in% "dbh_mm"]


# get plot size from database
KELuser <- dbPool(RPostgreSQL::PostgreSQL(),dbname = 'remoteforestsorg',
                  host = '91.239.201.14',port = 4010,user = 'remoteforests002', 
                  password = 'COVBtxQ5')

# list of plots with missing size data
plotids <- unique(neighbours_df$plotid)

plots <- tbl(KELuser, "plot") %>%
			filter(census %in% 1,           # first survey
			       plotid %in% plotids) %>% 
			collect() %>%
  dplyr::select(date, plotid, plotsize) 

plots$plotid <- as.factor(plots$plotid)

# add plot size to neighbours dataframe
neighbours_df <- merge(neighbours_df, plots[,-1], by="plotid",all.x=T)


###########
# generate buffer trees
pb <- txtProgressBar(min=0, max=length(unique(neighbours_df$plotid)), style=3)  # progress bar
tind=0
seed_no <- 0

plot_list <- levels(neighbours_df$plotid)
dat_array <- array(data=list(NULL),dim=length(plot_list))
  
i <- 1
for(i in 1:length(unique(neighbours_df$plotid))){
    seed_no <- seed_no + 1
	tind=tind+1
	setTxtProgressBar(pb, tind)
	
	plot <- plot_list[i]

	# select plot trees 
	trees_sim <- neighbours_df[neighbours_df$plotid %in% plot,]
	# rename treeid
	trees_sim$treeid <- paste(trees_sim$treeid,"_sim",sep="")
	
	# determine min and max radius for simulated zone
	Rmin <- unique(round(sqrt(trees_sim$plotsize/pi),2))  # radius of actual plot
	Rmax <- unique(round(sqrt(trees_sim$plotsize*3/pi),2)) # radius of an area equal to 3 X actual plot area

	# Simulate random tree positions where radius = distance from plot center and angle_rad = angle in radians
	set.seed(20 + seed_no)
	trees_sim$radius <- runif(nrow(trees_sim), min=Rmin, max=Rmax)
	trees_sim$angle_rad <- runif(nrow(trees_sim), min = 0.001, max = 2*pi)

	# calculate x and y coordinates
	trees_sim$x_m_ran <- trees_sim$radius * sin(trees_sim$angle_rad)
	trees_sim$y_m_ran <- trees_sim$radius * cos(trees_sim$angle_rad)

	# Calculate tree distance to plot centre based on new x and y tree coordinates
	trees_sim$DIST_m <- neighdist(trees_sim$x_m_ran, trees_sim$y_m_ran, 0, 0)

	dat_array[[i]] <- trees_sim		
}

sim_neighbours <- do.call(rbind,dat_array)


##### check
summary(sim_neighbours)
length(unique(sim_neighbours$plotid))  # 212
length(unique(neighbours_df$plotid))   # 212
length(unique(sim_neighbours$treeid))  # 15808
length(unique(neighbours_df$treeid))   # 15808

levels(neighbours_df$SPCD)
#  [1] "ABIALB"   "ACEHEL"   "ACEOBT"   "ACEPLA"   "ACEPSE"   "ACER"     "ANGIO"   
#  [8] "FAGSYL"   "FRAEXC"   "FRAORN"   "FRAXINUS" "LABANA"   "LIANS"    "PICABI"  
# [15] "RHAMNUS"  "SALCAP"   "SAMNIG"   "SORARI"   "SORAUC"   "TAXBAC"   "TILCOR"  
# [22] "TILIA"    "ULMGLA"   "ULMLAE"   "ULMUS"  

### check buffers for a random sample of plots
set.seed(42)
plot_list <- sample(sim_neighbours$plotid, 16)

if(computer=="~")
{ quartz(height=8,width=10)
}else{
  windows(height=8,width=10)
}
par(mfrow=c(4,4),mar=c(2,2,2,1),oma=c(1,2,2,1) ) 

for(index in 1:16)
{  plotid <- plot_list[index]
   plot(y_m_ran ~ x_m_ran, data=sim_neighbours[sim_neighbours$plotid %in% plotid,],col="gray")  # buffer trees
   points(y_m ~ x_m, data=neighbours_df[neighbours_df$plotid %in% plotid,], col="blue")       # measured trees
}


# replace x and y with buffer coordinates
sim_neighbours$x_m <- sim_neighbours$x_m_ran 
sim_neighbours$y_m <- sim_neighbours$y_m_ran
sim_neighbours <- sim_neighbours[,!names(sim_neighbours) %in% "x_m_ran" ]
sim_neighbours <- sim_neighbours[,!names(sim_neighbours) %in% "y_m_ran" ]

# drop radius and angle columns
sim_neighbours <- sim_neighbours[,!names(sim_neighbours) %in% c("radius","angle_rad")]


# Combine within-plot and buffer trees
sim_neighbours$tree_type <- "buffer"
neighbours_df$tree_type <- neighbours_df$treetype
neighbours_df <- rbind(neighbours_df,sim_neighbours)
neighbours_df <- neighbours_df[order(neighbours_df$plotid,neighbours_df$treeid),]


#########
# save
save(neighbours_df,file=file.path(Dir.Data,"10_Neighbour_tree_files.Rdata"))
















