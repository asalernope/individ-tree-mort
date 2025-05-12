###################################################################################
#  Computing a growing season moisture supply variable ("water_balance_seasonal_mm")
#  from Canham, C.D. and Murphy, L., 2016. 
#       The demography of tree species response to climate: 
#       sapling and canopy tree growth. Ecosphere, 7(10)
#  Dec 2024
###################################################################################

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


# load
load(file=file.path(Dir.Clim,"4_plot_monthly_water_budget_variables.Rdata"))
summary(monthly_sm_df)


#################################################################
#   Effective seasonal water supply ("water_balance_seasonal_mm")
#   defined as soil water storage, plus precip during 
#   growing season
#   - use balance between PET and PPT to define growing season
#   i.e. sum precip during months when PET >= PPT,
#   plus PET when PET < PPT,  and add to
#   that the available water storage
#################################################################

#  columns of monthly_sm_df containing
#    PET: 3:14
#    precip: 15:26

monthly_sm_df$water_balance_seasonal_mm <- monthly_sm_df$AWC_mm_per_m

# sum precip & AWC during months when PET >= PPT
for (i in 1:12)
{  monthly_sm_df$water_balance_seasonal_mm <- ifelse(monthly_sm_df[,2+i] >= monthly_sm_df[,14+i],
                                                     monthly_sm_df$water_balance_seasonal_mm + monthly_sm_df[,14+i],
                                                     monthly_sm_df$water_balance_seasonal_mm)
}

#  add in PET during months when PET < PPT (on assumption that precip will replenish PET)
for (i in 1:12)
{  monthly_sm_df$water_balance_seasonal_mm <- ifelse(monthly_sm_df[,2+i] < monthly_sm_df[,14+i],
                                                     monthly_sm_df$water_balance_seasonal_mm + monthly_sm_df[,2+i],
                                                     monthly_sm_df$water_balance_seasonal_mm)
}


####
save(monthly_sm_df,file=file.path(Dir.Clim,"5_plot_monthly_water_budget_variables.Rdata"))






