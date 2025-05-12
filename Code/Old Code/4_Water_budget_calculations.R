###########################################################
#  Water budget calculations
#  Base on Canham, C.D. and Murphy, L., 2016. 
#       The demography of tree species response to climate: 
#       sapling and canopy tree growth. Ecosphere, 7(10)
#  Dec 2024
###########################################################

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

# month vector
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# load
load(file=file.path(Dir.Clim,"3_Working_climate_and_soil_variables.Rdata"))

#----------------------------------------------------------------------------
# Calculate plot soil moisture for the end of each month
# sms.t1 = sms.t0 + precip.t1 - pet.t1
# ppt and pet are in units of mm
#----------------------------------------------------------------------------

# Define matrix for storage of monthly soil moisture calculations
monthly_sms <- matrix(NA, nrow=nrow(df), ncol=12)
colnames(monthly_sms) <- paste(months, "sms", sep=".")

# index for precip
precip_base <- 14  # January = column 15

# index for columns containing PET
pet_base <- 2  # January = column 3

# Get available water storage for all plots
aws <- df$AWC_mm_per_m   

# January soil moisture: aws + precip[i] - pet[i], since soil is saturated on last day of the year
# pmin and pmax bound the answer between 0 and AWS
monthly_sms[,1] <- pmax(pmin(aws + df[,(precip_base+1)] - df[,(pet_base+1)], aws),0)


### Feb - Nov soil moisture: sms[i-1] + precip[i] - pet[i]
for (i in 2:11) 
{  monthly_sms[,i] <- pmax(pmin(monthly_sms[,i-1] + df[,(precip_base+i)] - df[,(pet_base+i)], aws),0)
}

# December soil moisture: AWS
monthly_sms[,12] <- aws 

# combine variables
df <- cbind(df, monthly_sms)

#
#----------------------------------------------------------------------------
# Calculate AET summed across months
#    SMt1 = SMt0 + precip.t1 - PETt1
#    if SMt1 > 0, AETt1 = PETt1 
#    if SMt1 = 0, AETt1 = SMt0 + precip.t1
#
#  since soil is saturated at end of year, jan and dec AET = PET by definition
#----------------------------------------------------------------------------

# Calculate monthly AET
df$Jan_aet <- df[,pet_base+1]
df$Feb_aet <- ifelse(monthly_sms[,2]> 0,df[,pet_base+2],
                     monthly_sms[,1] + df[,precip_base+2])
df$Mar_aet <- ifelse(monthly_sms[,3]> 0,df[,pet_base+3],
                     monthly_sms[,2] + df[,precip_base+3])
df$Apr_aet <- ifelse(monthly_sms[,4]> 0,df[,pet_base+4],
                     monthly_sms[,3] + df[,precip_base+4])
df$May_aet <- ifelse(monthly_sms[,5]> 0,df[,pet_base+5],
                     monthly_sms[,4] + df[,precip_base+5])
df$Jun_aet <- ifelse(monthly_sms[,6]> 0,df[,pet_base+6],
                     monthly_sms[,5] + df[,precip_base+6])
df$Jul_aet <- ifelse(monthly_sms[,7]> 0,df[,pet_base+7],
                     monthly_sms[,6] + df[,precip_base+7])
df$Aug_aet <- ifelse(monthly_sms[,8]> 0,df[,pet_base+8],
                     monthly_sms[,7] + df[,precip_base+8])
df$Sep_aet <- ifelse(monthly_sms[,9]> 0,df[,pet_base+9],
                     monthly_sms[,8] + df[,precip_base+9])
df$Oct_aet <- ifelse(monthly_sms[,10]> 0,df[,pet_base+10],
                     monthly_sms[,9] + df[,precip_base+10])
df$Nov_aet <- ifelse(monthly_sms[,11]> 0,df[,pet_base+11],
                     monthly_sms[,10] + df[,precip_base+11])
df$Dec_aet <- df[,pet_base+12]


### calculate annual and seasonal variables
### water deficit (WD = PET - AET)
ind1 <- which(colnames(df) == "Jan_aet")
ind2 <- which(colnames(df) == "Dec_aet")
df$aet_ann_mm <- rowSums(df[,ind1:ind2])

ind1 <- which(colnames(df) == "Jan_pet")
ind2 <- which(colnames(df) == "Dec_pet")
df$pet_ann_mm <- rowSums(df[,ind1:ind2])

df$WD_mm_ann <- round((df$pet_ann_mm - df$aet_ann_mm), 3)


### growing season water deficit
ind1 <- which(colnames(df) == "Apr_aet")
ind2 <- which(colnames(df) == "Aug_aet")
df$aet_seasonal_mm <- rowSums(df[,ind1:ind2])

ind1 <- which(colnames(df) == "Apr_pet")
ind2 <- which(colnames(df) == "Aug_pet")
df$pet_seasonal_mm <- rowSums(df[,ind1:ind2])

df$WD_mm_seasonal <- round((df$pet_seasonal_mm - df$aet_seasonal_mm),3)

### calculate annual and seasonal variables - terra
### water deficit (WD = PET - AET)

ind1 <- which(colnames(df) == "Jan_def")
ind2 <- which(colnames(df) == "Dec_def")

df$terra_WD_mm_ann <- round(rowSums(df[,ind1:ind2]), 3)

ind1 <- which(colnames(df) == "Apr_def")
ind2 <- which(colnames(df) == "Aug_def")

df$terra_WD_mm_seasonal <- round(rowSums(df[,ind1:ind2]), 3)


### annual total precip
ind1 <- which(colnames(df) == "Jan_ppt")
ind2 <- which(colnames(df) == "Dec_ppt")
df$ppt_mm_ann_sum <- round(rowSums(df[,ind1:ind2]),3)

### growing season total precip
ind1 <- which(colnames(df) == "Apr_ppt")
ind2 <- which(colnames(df) == "Aug_ppt")
df$ppt_mm_seasonal_sum <- round(rowSums(df[,ind1:ind2]),3)

# min monthly precip during growing season
df$ppt_mm_seasonal_month_min <- round(apply(df[,ind1:ind2],1,function(x) min(x)), 3)



### annual mean temperature
ind1 <- which(colnames(df) == "Jan_temp_mean_K")
ind2 <- which(colnames(df) == "Dec_temp_mean_K")
df$temp_K_ann_ave <- round(rowMeans(df[,ind1:ind2]),3)


### growing season mean temperature
ind1 <- which(colnames(df) == "Apr_temp_mean_K")
ind2 <- which(colnames(df) == "Aug_temp_mean_K")
df$temp_K_seasonal_ave <- round(rowMeans(df[,ind1:ind2]),3)

# maximum monthly temp during growing season
df$temp_K_seasonal_month_max <- round(apply(df[,ind1:ind2],1,function(x)max(x)), 3)


### growing season mean soil moisture index
ind1 <- which(colnames(df) == "Apr_SMI")
ind2 <- which(colnames(df) == "Aug_SMI")
df$SMI_seasonal_ave <- round(rowMeans(df[,ind1:ind2]),3)

# min monthly moisture during growing season
df$SMI_seasonal_month_min <- round(apply(df[,ind1:ind2],1,function(x)min(x)), 3)




summary(df)



##### save
monthly_sm_df <- df
save(monthly_sm_df,file=file.path(Dir.Clim,"4_plot_monthly_water_budget_variables.Rdata"))















