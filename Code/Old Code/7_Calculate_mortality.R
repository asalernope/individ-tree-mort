###########################################################
#  Calculate mortality at end of census period
#  Dec 2024
###########################################################

rm(list=ls())

# user computer
computer <- "~"

# Directores
if(computer == "~"){
  Dir.Base <- "~/Desktop/Prague/mort-likelihood"
}else{
  Dir.Base <- "insert path"
}


Dir.Data <- file.path(Dir.Base, "Data")
Dir.Code <- file.path(Dir.Base, "Code")
Dir.Subsample <- file.path(Dir.Data, "Subsampled_data_files")
dirs <- sapply(c(Dir.Base, Dir.Data, Dir.Code, Dir.Subsample), function(x) if(!dir.exists(x)) dir.create(x) )


# get tree plot data
load(file=file.path(Dir.Data,"full_trees.Rdata"))

targets <- full_trees

targets$SPCD <- as.factor(targets$SPCD)
colnames(targets)[names(targets) %in% "status"] <- "tree_status_census_1"
colnames(targets)[names(targets) %in% "status_re"] <- "tree_status_census_2"

# live or dead class at end of census period by tree
# 0 = alive
# 1 = dead
targets$lord <- ifelse(!targets$tree_status_census_2 %in% c(1,2),1,0)
unique(targets$lord)
head(targets)


table(targets$SPCD,targets$lord)
#               0    1
#   ABIALB   2318  121
#   ACEOBT     48    5
#   ACEPLA      1    0
#   ACEPSE    183   14
#   ACER        5    0
#   FAGSYL   8037  428
#   FRAEXC      9    2
#   FRAORN      1    0
#   FRAXINUS    1    0
#   LABANA      1    1
#   PICABI    546   54
#   RHAMNUS     3    0
#   SAMNIG      5    0
#   SORARI      3    0
#   SORAUC      4    1
#   TILCOR     11    0
#   TILIA       1    0
#   ULMGLA      6    1
#   ULMUS       1    0


### save
save(targets, file=file.path(Dir.Data,"7_Trees_with_mortality.Rdata"))











