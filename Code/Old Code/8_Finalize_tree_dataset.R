#######################################################################
#  Finalize tree dataset
#  add integer model response variable (mort)
#  assume response = prob of survival
#  live trees = 1 at 2nd census
#  dead trees = 0 at 2nd census
#  edit & reclassify canopy layer data into two classes (upper vs lower)
#  add variable for length of census interval
#  12.02.25
#
#######################################################################

rm(list=ls())

# user computer
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
Dir.Subsample <- file.path(Dir.Data, "Subsampled_data_files")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Subsample), function(x) if(!dir.exists(x)) dir.create(x) )


### get tree plot data (filtered and expanded by Audrey)

targets <- load("~/Desktop/Prague/mort-likelihood/Data/full_trees.Rdata")
targets<-get(targets)

targets$SPCD <- as.factor(targets$SPCD)
colnames(targets)[names(targets) %in% "status"] <- "tree_status_census_1"
colnames(targets)[names(targets) %in% "status_re"] <- "tree_status_census_2"
colnames(targets)[names(targets) %in% "date"] <- "date_census_1"
colnames(targets)[names(targets) %in% "date_re"] <- "date_census_2"


# length of census interval
targets$census_yrs <- targets$date_census_2 - targets$date_census_1

# dbh in cm
# targets$dbh_cm <- targets$dbh_mm / 10
# targets <- targets[,!names(targets) %in% "dbh_mm"]


#### add model response variable (mort)
# 0 = dead
# 1 = alive
targets$mort <- ifelse(targets$tree_status_census_2 %in% c(1,2,3,4),1,0)


#### reclassify canopy layer data (canopy position or strata)
table(targets$layer,targets$growth)
#         0    1   99
#   -1    0    0    1
#   11  269 5716    0
#   12 1995 1324    0
#   13 4065  765    0
#   99    0    0   11

filter(targets,targets$layer == 99)
filter(targets,targets$layer == -1)

# replace -1 or 99 layer values (1st census)
#targets$layer_corr <- ifelse(targets$layer %in% c(-1,99), targets$layer_re, targets$layer)

targets$layer_corr <- ifelse(targets$layer == 99 & targets$dbh_mm < 200, 13,  # assume small trees are shaded
                             ifelse(targets$layer == 99 & targets$dbh_mm >= 200, 11,
                                    targets$layer))
targets$layer_corr <- ifelse(targets$layer == -1 & targets$growth == 99, 11, targets$layer_corr)

table(targets$layer_corr,targets$growth)

# reclassify trees into two strata (upper vs lower canopy)
targets$stratum <- ifelse(targets$layer %in% c(11,12),"upper","lower")
targets$stratum <- as.factor(targets$stratum)


#### save
save(targets, file=file.path(Dir.Data,"8_Trees_with_mortality.Rdata"))







