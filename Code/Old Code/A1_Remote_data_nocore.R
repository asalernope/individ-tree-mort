#################################################################
#  Disturbance - mortality modeling Balkans
#  Download Remote tree data with ring widths
#  Canopy layer will need further filtering... 
#  Get x and y tree coordinates from re-survey year if available
#  December 2024
#  Feb 2025: Remote tree filtering was revised by Audrey. This code
#            is no longer used
#################################################################

rm(list=ls(all=TRUE))

# User computer
computer <- "arne"

# Directores
if(computer == "arne"){
	Dir.Base <- "/Users/arnebuechling/Documents/CULS/Projects_miro_lab_group/Audrey/Mortality"
}else{
	Dir.Base <- "insert path"
}

Dir.Data <- file.path(Dir.Base, "Data")
Dir.Code <- file.path(Dir.Base, "Code")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code), function(x) if(!dir.exists(x)) dir.create(x) )


# packages
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}
package_vec <- c("RPostgreSQL","pool","tidyverse")
sapply(package_vec, install.load.package)

# load Balkan tree data
balkan_trees <- read.csv(file=file.path(Dir.Data,"balkan_trees.csv"))
which(duplicated(balkan_trees$treeid))

plotids <- unique(balkan_trees$plotid)
treeids <- unique(balkan_trees$treeid)  # n=9202
head(balkan_trees)

# database
KELuser <- dbPool(RPostgreSQL::PostgreSQL(),dbname = 'remoteforestsorg',
                  host = '91.239.201.14',port = 4010,user = 'remoteforests002', 
                  password = 'COVBtxQ5')


#########################################
# Part 1 -- get live tree data with cores
# all stands 
# 1st census
#########################################

trees <- tbl(KELuser, "plot") %>%
			filter(census %in% 1,           # plot establishment, first survey
			       #foresttype == "beech",
			       plotid %in% plotids) %>% 
			
			select(date, census, plotid, plot_id = id, plotsize, lng, lat, stand, country, foresttype, altitude_m,slope,aspect) %>%
			
			inner_join(.,
			tbl(KELuser, "tree") %>% 
			filter(#!onplot %in% 0,
					#!species %in% "99",
					#!is.na(x_m), 
					#treeid %in% treeids,
					!is.na(dbh_mm),
					status %in% c(1)     # alive undamaged, include 2 ?    
					#growth == 1,        # not shaded
					#!layer %in% 13      # remove trees in lower canopy ?
					) %>%              
			
			select(plot_id, treeid, tree_id=id, species, x_m, y_m, dbh_mm, height_m,crownht_m, crowndiam1_m, crowndiam2_m,
				   treetype, tree_status=status, shading=growth, can_layer=layer, onplot),
			by = "plot_id") %>%
						
			collect()

trees <- data.frame(trees)

trees$foresttype <- as.factor(trees$foresttype)
trees$stand <- as.factor(trees$stand)
trees$plotid <- as.factor(trees$plotid)
trees$treeid <- as.factor(trees$treeid)

summary(trees)
str(trees)  # 157 plots 

length(levels(trees$treeid))  # 12006 ?
length(levels(trees$stand))   # 11 stands
nrow(trees[trees$shading==99,])  # 22 trees with no growth class

# reduced tree list
treeids_2 <- levels(trees$treeid)


###########################################################
# Part 2: get tree data from most recent census 
###########################################################

trees_resurvey <- tbl(KELuser, "plot") %>%
				  filter(!census == 1,  # exclude 8?
						 #plottype %in% c(3,4),
						 plotid %in% plotids) %>%
				  arrange(plotid, desc(date)) %>%
					
				  group_by(plotid) %>% 
				  filter(row_number() == 1) %>%
				  ungroup() %>%
					
				  select(date_re=date, lng, lat, plotid, stand, plot_id = id, plotsize_re = plotsize, census_re = census) %>%
				  
				  inner_join(.,
				  tbl(KELuser, "tree") %>% 
				  filter(#!onplot %in% 0,
				         #!status %in% '99'),   # 3 trees with no status
						 treeid %in% treeids_2) %>% 
						 
				  select(plot_id, treeid,tree_id_re = id, species, x_m_re = x_m, y_m_re = y_m, dbh_mm_re = dbh_mm, status_re = status),	
				  by = "plot_id") %>%
				  
				  collect()

trees_resurvey <- data.frame(trees_resurvey)
trees_resurvey$stand <- as.factor(trees_resurvey$stand)
trees_resurvey$plotid <- as.factor(trees_resurvey$plotid)
trees_resurvey$treeid <- as.factor(trees_resurvey$treeid)
colnames(trees_resurvey)[names(trees_resurvey) %in% "plot_id"] <- "plot_id_re"

trees_resurvey <- trees_resurvey[order(trees_resurvey$plotid,trees_resurvey$treeid),]

levels(trees_resurvey$stand)  # 11
length(levels(trees_resurvey$treeid) )  # 11937 


#########  combine datasets  #########
# join common records from the full 1st-census and re-surveyed datasets
# some trees have an uncertain species designation (species ids do not match in the two censuses)
# some re-surveyed trees are missing in the second census (fallen dead trees not included in tree table in some cases) 
# use inner_join (merge) to exclude missing trees and trees with uncertain species id 

trees_df <- merge(trees,trees_resurvey,by=c("stand","plotid","treeid","lng","lat","species"),all=F)
trees_df <- droplevels(trees_df)

# use x and y tree coordinates from re-survey (remove coordinates from first survey)
trees_df <- trees_df[, !names(trees_df) %in% c("x_m","y_m")]
colnames(trees_df)[names(trees_df) %in% "x_m_re"] <- "x_m"
colnames(trees_df)[names(trees_df) %in% "y_m_re"] <- "y_m"

colnames(trees_df)[names(trees_df) %in% "date"] <- "date_census_1"
colnames(trees_df)[names(trees_df) %in% "date_re"] <- "date_census_2"


# dbh in cm
trees_df$dbh_census1_cm <- trees_df$dbh_mm / 10
trees_df$dbh_census2_cm <- trees_df$dbh_mm_re / 10

summary(trees_df)  
str(trees_df)  # 157 plots and 11811 trees

##########
save(trees_df, file=file.path(Dir.Data,"A1_Balkan_trees_no_cores.Rdata"))

poolClose(KELuser)









