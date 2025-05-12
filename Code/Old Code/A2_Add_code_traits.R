####################################################################
#  Balkan disturbance - survival
#  Add species codes to tree data
#  Dec 2024
#  Revised 10.01.25: use revised expanded tree dataset (157 plotids)
####################################################################

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


# get tree tree data (original filtered dataset)
#load(file=file.path(Dir.Data,"A1_Balkan_trees_no_cores.Rdata"))
#str(trees_df) # n=11811

# load tree dataset revised by Audrey
name <- "Filltered_Balkan_trees_for_mortality_analysis_v2"
targets <- load(file=file.path(Dir.Data,paste0(name,".Rdata")))
targets <- get(targets)
targets$plotid <- as.factor(targets$plotid)
str(targets)   # 157 plotids
length(targets$treeid)

# read species table
setwd(Dir.Data)
spp_table <- read.table("Species_table.txt", header=T)
spp_table <- spp_table[,!names(spp_table) %in% c("Genus_nam","Species_nam","functional.grp","sp_type")]

str(spp_table)

# attach
targets <- merge(targets,spp_table,by="species",all.x=T)
unique(targets$SPCD)
unique(targets$species)
head(targets)

#####
# remove unnecessary columns
#names <- c("id","date_re","plotsize_re","tree_id_re","dbh_mm_re","plot_id_re")
#targets <- targets[,!names(targets) %in% names]


#####
# save
save(targets,file=file.path(Dir.Data,"A2_Balkan_trees_no_cores_code.Rdata"))





