#########################################################
#  Compile species table
#  wood density and specific leaf area values from:
#    1. Forrester, et al. 2017. Generalized biomass and leaf area 
#       allometric equations for European tree species FEM 396,160-175.
#    2. Zanne, Amy E. et al. (2009). Data from: Towards a worldwide wood 
#       economics spectrum [Dataset]. Dryad. https://doi.org/10.5061/dryad.234
#  December 2024
#########################################################

rm(list=ls(all=TRUE))

# User computer
computer <- "arne"

# Directores
if(computer == "arne"){
	Dir.Root <- "/Users/arnebuechling/Documents/CULS/Projects_miro_lab_group/Audrey"
}else{
	Dir.Root <- "insert path"
}
Dir.Base <- file.path(Dir.Root, "Mortality")
Dir.Data <- file.path(Dir.Base, "Data")
Dir.Code <- file.path(Dir.Base, "Code")
Dir.Forr.Traits <- file.path(Dir.Data, "Forrester_trait_tables")
Dir.GWDD <- file.path(Dir.Data, "GlobalWoodDensityDatabase")
dirs <- sapply(c(Dir.Root, Dir.Base, Dir.Data, Dir.Code, Dir.Forr.Traits, Dir.GWDD), function(x) if(!dir.exists(x)) dir.create(x) )


# packages
library("stringr")
library(RPostgreSQL)
library(pool)

# Connect to database
KELuser <- dbPool(RPostgreSQL::PostgreSQL(),
                  dbname = 'remoteforestsorg',
                  host = '91.239.201.14',
                  port = 4010,
                  user = 'remoteforests002', 
                  password = 'COVBtxQ5')

# Get list of tables in database
sort(dbListTables(KELuser))

# Get species list
dbListFields(KELuser, "species_fk")
spec_table <- dbReadTable(KELuser, "species_fk")
spec_table <- spec_table[order(spec_table$id),]
save(spec_table,file=file.path(Dir.Data,"Remote_species_table.Rdata"))


# Get wood density data
#density <- dbReadTable(KELuser, "wood_density")

poolClose(KELuser) 


##########################
# Extract Genus and Species
spec <- spec_table
spec$Genus_nam <- str_split(spec$id, boundary("word"), simplify = T)[,1]  # extract second word
spec$Species_nam <- str_split(spec$id, boundary("word"), simplify = T)[,2]  # extract second word
spec$Genus_nam <- ifelse(spec$Genus_nam %in% '99', NA, spec$Genus)

# Add 6 letter species code based on taxanomic name
spec$SPCD <- ifelse(spec$Species_nam %in% "", spec$Genus, paste(substr(spec$Genus,1,3), substr(as.factor(spec$Species),1,3),sep=""))
spec$SPCD <- toupper(spec$SPCD)
spec$SPCD <- ifelse(spec$SPCD %in% "BROADLEAVES","ANGIO",spec$SPCD) 
spec$SPCD <- ifelse(spec$SPCD %in% "CONIFEROUS","GYMNO",spec$SPCD) 
spec$SPCD <- ifelse(spec$Species_nam %in% "obtusifolium", "ACEOBFO", spec$SPCD)

colnames(spec)[names(spec) %in% "id"] <- "species"


##########################
# read wood density data from Forrester 2017 and edit density values from above
#setwd("/Users/arnebuechling/Documents/Environmental_data/Traits/trait_database_Forrester2017/")

WD <- read.csv(file = file.path(Dir.Forr.Traits,"traits_A3_traits_Forrester_2017.csv"),header=T)
colnames(WD)[names(WD) %in% "Wood.density..g.per.cm3."] <- "wood_density_gcm3"

# calculate mean density per species
WD_Forr_mean <- aggregate(wood_density_gcm3~Species,data=WD, function(x) mean(x,na.rm=T))

# edit names to match spec file from database
WD_Forr_mean$id <- as.character(WD_Forr_mean$Species)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "A.alba", "Abies alba",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "A.pseudoplatanus", "Acer pseudoplatanus",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "A.incana", "Alnus incana",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "A.glutinosa", "Alnus glutinosa",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "B.pendula", "Betula pendula",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "B.pubescens", "Betula pubescens",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "B.pendulaOrpubescens", "Betula",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "C.betulus", "Carpinus betulus",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "C.sativa", "Castania sativa",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "F.sylvatica", "Fagus sylvatica",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "F.excelsior", "Fraxinus excelsior",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "L.decidua", "Larix decidua",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.abies", "Picea abies",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.cembra", "Pinus cembra",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.nigra", "Pinus nigra",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.sylvestris", "Pinus sylvestris",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.tremula", "Populus tremula",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.tremulaXtremuloides", "Populus",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.avium", "Prunus avium",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.serotina", "Prunus serotina",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "P.menziesii", "Pseudotsuga menziesii",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "Q.petraea", "Quercus petraea",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "Q.ilex", "Quercus ilex",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "Q.roburOrpetraea", "Quercus",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "Q.robur", "Quercus robur",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "R.pseudoacacia", "Robinia pseudoacacia",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "S.aucuparia", "Sorbus aucuparia",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "T.cordata", "Tilia cordata",WD_Forr_mean$id)
WD_Forr_mean$id <- ifelse(WD_Forr_mean$Species %in% "T.cordataOrplatyphyllos", "Tilia",WD_Forr_mean$id)

WD_Forr_mean <- WD_Forr_mean[,!names(WD_Forr_mean) %in% "Species"]
colnames(WD_Forr_mean)[names(WD_Forr_mean) %in% "wood_density_gcm3"] <- "wood_density_gcm3_Forrester"
colnames(WD_Forr_mean)[names(WD_Forr_mean) %in% "id"] <- "species"


##########################
# read wood density data from Zanne et al 2009 (GlobalWoodDensityDatabase)

GWDD <- read.csv(file=file.path(Dir.GWDD,"GWDD_zanne_2009.csv"))
unique(GWDD$Region)

# subset Europe
GWDD_europe <- subset(GWDD,Region=="Europe")

colnames(GWDD_europe)[names(GWDD_europe) %in% "Wood.density..g.cm.3...oven.dry.mass.fresh.volume"] <- "wood_density_zanne_2009_gcm3"
colnames(GWDD_europe)[names(GWDD_europe) %in% "Binomial"] <- "species"
GWDD_europe <- GWDD_europe[,!names(GWDD_europe) %in% "Number"]
GWDD_europe <- GWDD_europe[,!names(GWDD_europe) %in% "X"]
GWDD_europe <- GWDD_europe[,!names(GWDD_europe) %in% "X.1"]
GWDD_europe <- GWDD_europe[,!names(GWDD_europe) %in% "Reference.Number"]
GWDD_europe <- GWDD_europe[,!names(GWDD_europe) %in% "Region"]

# compute mean by species
GWDD_europe_mean <- aggregate(wood_density_zanne_2009_gcm3 ~ species, data=GWDD_europe, function(x) mean(x,na.rm=T))


##########################
# merge species info and two wood density datasets
species_table <- merge(spec, WD_Forr_mean, by="species", all.x=T)
species_table <- merge(species_table, GWDD_europe_mean, by="species", all.x=T)


# integrate wood density values
species_table$wood_density_gcm3 <- ifelse(is.na(species_table$wood_density_zanne_2009_gcm3), species_table$wood_density_gcm3_Forrester, 
                                                species_table$wood_density_zanne_2009_gcm3)

species_table$wood_density_gcm3 <- round(species_table$wood_density_gcm3, 3)

# remove columns
#species_table <- species_table[,!(names(species_table) %in% c("wood_density_gcm3_Forrester"))]
species_table <- species_table[,!names(species_table) %in% "sp_group_dist"]
species_table <- species_table[,!names(species_table) %in% "sp_code"]
species_table <- species_table[,!names(species_table) %in% "functional.grp"]


##########################
# get specific leaf area from Table A1 in Forrester et al. 2017
SLA <-read.csv(file=file.path(Dir.Forr.Traits,"specific_leaf_area_Table_A1.csv"))
names(SLA)
# calculate species level means
SLA_spp_means <- aggregate(cbind("SLA_m2_kg_mean"=SLA_m2_kg)~species,data=SLA,function(x) mean(x,na.rm=T))
SLA_spp_means$SLA_m2_kg_mean <- round(SLA_spp_means$SLA_m2_kg_mean,3)

#SLA_spp_means$SLA_cm2g_mean <- round(SLA_spp_means$SLA_m2_kg_mean * 1000 / 100^2, 3)
SLA_spp_means$SLA_cm2g_mean <- round(SLA_spp_means$SLA_m2_kg_mean * 100^2 / 1000, 3)


# combine
species_table <- merge(species_table,SLA_spp_means,by="species",all.x=T)
species_table <- species_table[,!names(species_table) %in% "SLA_m2_kg_mean"]
species_table

# check for duplicates
which(duplicated(species_table$SPCD))


##########################
# Save 
setwd(Dir.Data)
write.table(species_table,"Species_table.txt", sep="\t",row.names=FALSE)












