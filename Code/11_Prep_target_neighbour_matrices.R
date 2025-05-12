########################################################################################
# Prepare target and neighbour matrices for regression analyses
# For neighbours: generate matrices of species, distances and dbhs for competition analyses
#                 use neighbour trees >= 6cm dbh
#                 create neighbour groups (via clustering of wood density and LAI values)
# Canham et al 2004. A neighborhood analysis of canopy tree competition: effects of 
#                    shading versus crowding. Canadian Journal of Forest Research, 34(4), pp.778-787.
# Uriarte et al 2004. A neighborhood analysis of tree growth and survival in a hurricane‐driven 
#                     tropical forest. Ecological Monographs
# This step will need to be rerun if filtering criteria change
# Date. 15.04.2025
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
Dir.Analysis.Files <- file.path(Dir.Data,"Full_data_files")
dirs <- sapply(c(Dir.Base,Dir.Data,Dir.Code,Dir.Clim,Dir.Dist,Dir.Ndep,Dir.Out,Dir.Analysis.Files), 
               function(x) if(!dir.exists(x)) dir.create(x))

library("stringr")

##########  data filters  ##########
# Minimum target tree DBH
min.target.dbh.cm <- 6
# Minimum neighbour tree DBH
min.neighbour.dbh.cm <- 6
# Maximum radius (m) for neighbours
neighbour_dist_max <- 15
# min sample size of trees for growth analyses
Nmin <- 100 
# number of classes of wood density
k <- 3
####################################

# function to calculate distance
neighdist <- function(targetx, targety, neighborx, neighbory) { 
	                  sqrt((targetx-neighborx)^2 + (targety-neighbory)^2) }

# Status bar
status_bar <- function(title.txt, n, nmax)
{ m <- paste(title.txt, ":", n,"of",nmax,"Completed")
  barplot(height=c(n), width = c(0.25), main = m, xlim = c(0,nmax), ylim = c(0,1), horiz=T)
}

###################################
# get neighbour tree data
load(file=file.path(Dir.Data,"10_Neighbour_tree_files.Rdata"))

# drop trees with missing trait values
unique(neighbours_df[,names(neighbours_df) %in% c("SPCD","wood_density_gcm3","SLA_cm2g_mean")] )
neighbours_df <- subset(neighbours_df, !SPCD %in% c("ANGIO","LIANS","ACEHEL") )  # missing trait values


# load target tree data
load("~/Desktop/Prague/mort-likelihood/Data/9_Trees_with_derived_climate_variables.Rdata")

# Calculate tree distance to plot centre
trees$DIST_m <- neighdist(trees$x_m, trees$y_m, 0, 0)
trees$dbh_cm <- trees$dbh_mm/10


# select focal species with minimum number of samples
spp_list <- sort(table(trees$SPCD),decreasing = TRUE)
top_species <- spp_list[which(spp_list > Nmin)]


###################################
# loop through selected species and collect neighbour data matrices

for(index in 1:length(top_species)) {
	
	spname <- names(top_species[index])
	targets <- subset(trees, SPCD == spname &
	                  dbh_cm >= min.target.dbh.cm &
					  !is.na(time_since_dist_recent) )
	targets <- droplevels(targets)
	
	# list of plots for this species
	target.plots  <- levels(targets$plotid)
	
	# subset neighbours that are in target plots
	neighbours <- subset(neighbours_df, plotid %in% target.plots)

	# Exclude neighbours under the minimum defined neighbour DBH
	neighbours <- subset(neighbours, dbh_cm >= min.neighbour.dbh.cm )
	neighbours <- droplevels(neighbours)
	

	#############
	# Add neighbor classification schemes 
	# relatedness (conspecific, congeneric, other); add family?
	
	neighbours$genus <- str_split(neighbours$species, boundary("word"), simplify = T)[,1]  # extract genus
    targets$genus <- str_split(targets$species, boundary("word"), simplify = T)[,1]  
	
	neighbours$relatedness <- ifelse(neighbours$SPCD %in% targets$SPCD, "conspec", 
	                            ifelse(neighbours$genus %in% targets$genus, "congen", "other"))
	neighbours$relatedness <- as.factor(neighbours$relatedness)
	
	# classify neighbours species by clustering trait values
	d <- unique(neighbours[,names(neighbours) %in% c("SPCD","wood_density_gcm3","SLA_cm2g_mean")] )
	d$wood_density_kgm3 <- d$wood_density_gcm3 * 1000
	kmeans.out <- kmeans(d[,c(3:4)], centers=k, nstart=20)

	trait_class <- cbind(d,"trait_class"=kmeans.out$cluster)
	trait_class <- trait_class[order(trait_class$trait_class),]

	# Check clusters
	quartz(height=6, width=6)
	plot(SLA_cm2g_mean ~ wood_density_kgm3, data=d, col=(kmeans.out$cluster+1), xlab="Wood density (kg/m3)", ylab="SLA (cm2/g)", 
		 pch=20, cex=2, las=1, main=paste0("focal = ",spname))

	# add clusters to neighbours file
	neighbours <- merge(neighbours, trait_class, by=c("SPCD", "wood_density_gcm3", "SLA_cm2g_mean"),all.x=T)
	neighbours$trait_class <- as.factor(neighbours$trait_class)

	#############
	# generate matrices of species, distances and dbhs for neighbours
	quartz(height=2.5, width=9)
	par(mfrow=c(1,1))

	# determine maximum number of neighbours for a target tree
	max.neighbours <- 0
	for (m in 1:nrow(targets))
	{ 
	  max.neighbours <- max(max.neighbours, 
							nrow(subset(neighbours, plotid %in% targets$plotid[m])) )
	  status_bar(paste("Calculating max neighbours for species", spname), m, nrow(targets))
	}
	# initialize matices
	distances <- matrix(0, nrow=nrow(targets), ncol=max.neighbours)
	dbhs <- matrix(0, nrow=nrow(targets), ncol=max.neighbours)
	species <- matrix(0, nrow=nrow(targets), ncol=max.neighbours)
    trait_clust <- matrix(0, nrow=nrow(targets), ncol=max.neighbours)
	
	# populate matrices
	quartz(height=2.5, width=9)
	for(i in 1:nrow(targets)) 
	{ neighbours.for.tree <- subset(neighbours, neighbours$plotid %in% targets$plotid[i] )

	  if(nrow(neighbours.for.tree) > 0) {
		distances[i,1:nrow(neighbours.for.tree)] <- neighdist(targets$x_m[i], targets$y_m[i], # using distance function                 
															 neighbours.for.tree$x_m, neighbours.for.tree$y_m)                                                
		dbhs[i,1:nrow(neighbours.for.tree)] <- neighbours.for.tree$dbh_cm     
		species[i,1:nrow(neighbours.for.tree)] <- neighbours.for.tree$relatedness
	    trait_clust[i,1:nrow(neighbours.for.tree)] <- neighbours.for.tree$trait_class
	  }
	  status_bar(paste("Populating neighbour matrices for", spname), i, nrow(targets))
	}

	# convert dbhs to meters for derivation of competition index
	dbhs <- dbhs/100.0

	# Replace zeros with NA
	distances <- ifelse(distances == 0, NA, distances)
	dbhs[which(distances %in% NA)] <- NA
	distances[which(dbhs %in% NA)] <- NA
	species[which(dbhs %in% NA)] <- NA
    trait_clust[which(dbhs %in% NA)] <- NA

	# Exclude neighbours beyond a maximum distance of a target
	for (i in 1:nrow(distances))
	{ too_far <- which(distances[i,] > neighbour_dist_max)
	  distances[i,too_far] <- NA
	  dbhs[i,too_far] <- NA
	  species[i,too_far] <- NA
	  trait_clust[i,too_far] <- NA
	}
	colnames(dbhs) <- rep("NDBH", ncol(dbhs))
	colnames(distances) <- rep("NDIST", ncol(dbhs))
	colnames(species) <- rep("NSPP", ncol(dbhs))
    colnames(trait_clust) <- rep("NTR", ncol(dbhs))
	
	# get number of larger neighbours relative to focal tree (as a simple competition index)
	nd <- rowSums(dbhs > (targets$dbh_cm/100),na.rm=T)

	# size ratio = average dbh of all neighbours / dbh target (Canham and Murphy 2017)
	targets$size_ratio <- rowMeans((dbhs/2)^2, na.rm=T) / (((targets$dbh_cm)/200)^2)  # dbhs are in meters
	
	# save
	
	setwd(Dir.Analysis.Files)
	save(targets,neighbours,distances,dbhs,species,trait_clust,nd,trait_class,kmeans.out,
		 file=paste(spname,"Full_working_dataset_with_all_climate_variables.Rdata",sep="_"))

} # species loop



# lognormal size and nci (based on taxonomic similarity)
# index lambda based on relatedness classes
# 
nci_spp <- function(lambda,alpha,beta)
{ 
  nci <- rowSums(lambda[species] * (dbhs ^ alpha)/(distances ^ beta), na.rm=T)
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















