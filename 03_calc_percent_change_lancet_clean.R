#########################################################################################################################
## Child script to calculate percent change between two specified years
## Output: csv files of percent change for all-age counts and age-standardized rates (prevalence and YLDs)
#########################################################################################################################

## Set up ----------------------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")
library(purrr)

#Arguments
args <- commandArgs(trailingOnly = T)

group_level <- args[1]
name <- args[2]

print(paste0(group_level, ", ", name, ", ", year_start, ", ", year_end))

#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id,22)

loc_df <- read.csv(paste0(save_dir, "locations.csv"))
loc_df <- loc_df %>% filter(level %in% loc_level)
loc_ids <- loc_df %>% filter(level %in% loc_level) %>% pull(location_id)

pop_df <- read.csv(paste0(save_dir, "population.csv"))

## Functions
#Summarize draws function
summaries <- function(dt, draw_vars, count_vars){
  sum <- as.data.table(copy(dt))
  sum[, mean_rate := rowMeans(.SD), .SDcols = draw_vars]
  sum[, standard_error_rate := apply(.SD, 1, sd), .SDcols = draw_vars]  
  sum[, lower_rate := apply(.SD, 1, quantile, probs= 0.025), .SDcols = draw_vars]
  sum[, upper_rate := apply(.SD, 1, quantile, probs=0.975), .SDcols = draw_vars]
  sum[, mean_count := rowMeans(.SD), .SDcols = count_vars]
  sum[, standard_error_count := apply(.SD, 1, sd), .SDcols = count_vars]  
  sum[, lower_count := apply(.SD, 1, quantile, probs= 0.025), .SDcols = count_vars]
  sum[, upper_count := apply(.SD, 1, quantile, probs=0.975), .SDcols = count_vars]
  sum[, c(draw_vars, count_vars) := NULL]
  return(sum)
}


## Calculate % change between two specified years for all-age cases and age-standardised rates -------------------------------------
#Use this equation: [(new val - orig val)/orig val]*100

#First get relevant files for all age and age-standardized draws and create one dt, subset to specified years for comparison
if(group_level!="group2"){
  dt_orig <- as.data.table(readRDS(paste0(FILEPATH))) #get count and rate draws
}
if(group_level=="group2"){
  dt_orig <- as.data.table(readRDS(paste0(FILEPATH))) #get count and rate draws
}

dt_agestd <- dt_orig[age_group_id==ID & (year_id==year_start | year_id==year_end), ] #subset to age-standardized and desired years
dt_agestd <- dt_agestd[, c(draw_vars, "population", "age_group_id") := NULL] #get rid of unnecessary variables
dt_count <- dt_orig[age_group_id==ID & (year_id==year_start | year_id==year_end), ] #subset to all age and desired years
dt_count <- dt_count[, c(rate_vars, "population", "age_group_id") := NULL] #get rid of unnecessary variables

dt <- merge(dt_count, dt_agestd, by=c("measure_id", "location_id", "sex_id", "year_id", "short_name")) #combine
for (j in 1:ncol(dt)) set(dt, which(is.infinite(dt[[j]])), j, 0)

#Make wide format by year
#dt_wide <- dcast(dt, measure_id + location_id  + sex_id + short_name ~ year_id, value.var = c(draw_vars, rate_vars))
dt_wide <- reshape(dt, idvar = c("measure_id", "location_id", "sex_id", "short_name"), timevar = "year_id", direction = "wide")


#Divide year_end by year_start cases, then subtract 1 and multiply by 100 to get % change
dt_wide[, paste0("change_count_", 0:nf_draws) := lapply(0:nf_draws, function(x) ((get(paste0("draw_", x, ".", year_end)) / get(paste0("draw_", x, ".", year_start))-1)*100))]
dt_wide[, paste0("change_rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) ((get(paste0("rate_", x, ".", year_end)) / get(paste0("rate_", x, ".", year_start))-1)*100))]
dt_wide <- dt_wide[, c(paste0("rate_", 0:nf_draws,"_1990"), paste0("rate_", 0:nf_draws,".2021"), paste0("draw_", 0:nf_draws,".1990"), paste0("draw_", 0:nf_draws,".2021")) := NULL]

#Collapse
if(name=="klinefelter" | name=="cysticercosis" | name=="malaria" | name=="fetal_alcohol" | name=="intellectual_disability"){
  dt_wide[is.na(dt_wide), ] <- 0  #Klinefelter is male only, cysticercosis and malaria are geographically restricted so set NaN to 0% change
}

if(name=="zika" | name=="covid"){
  for (j in 1:ncol(dt_wide)) set(dt_wide, which(is.infinite(dt_wide[[j]])), j, 100)
  for (j in 1:ncol(dt_wide)) set(dt_wide, which(is.na(dt_wide[[j]])), j, 100)
}

dt_change <- summaries(dt_wide, count_vars = paste0("change_count_", 0:nf_draws), draw_vars = paste0("change_rate_",0:nf_draws))


#Finalize and save
dt_change <- as.data.table(merge(dt_change, loc_df[,c("location_id", "location_name")], by = c("location_id")))
dt_change[, sex := ifelse(sex_id==1, "Male", ifelse(sex_id==2, "Female", "Both"))]
dt_change[, measure := ifelse(measure_id==5, "Prevalence", "YLDs")]
dt_change[, description := paste0("Percent change (out of 100) in counts and age-standardized rates from ", year_start, " to ", year_end)]

write.csv(dt_change, paste0(FILEPATH), row.names = F)


## Calculate female to male ratio for all-age cases and age-standardised rates -----------------------------------------------------------------------
dt <- merge(dt_count[sex_id!=3], dt_agestd[sex_id!=3], 
            by=c("measure_id", "location_id", "sex_id", "year_id", "short_name")) #combine

#Make wide format by sex
dt_wide <- reshape(dt, idvar = c("measure_id", "location_id", "year_id", "short_name"), timevar = "sex_id", direction = "wide")

#Divide female by male cases
dt_wide[, paste0("ratio_count_", 0:nf_draws) := lapply(0:nf_draws, function(x) (get(paste0("draw_", x, ".2")) / get(paste0("draw_", x, ".1"))))]
dt_wide[, paste0("ratio_rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) (get(paste0("rate_", x, ".2")) / get(paste0("rate_", x, ".1"))))]
dt_wide <- dt_wide[, c(paste0("rate_", 0:nf_draws,".1"), paste0("rate_", 0:nf_draws,".2"), paste0("draw_", 0:nf_draws,".1"), paste0("draw_", 0:nf_draws,".2")) := NULL]

#Collapse
if(name=="klinefelter" | name=="cysticercosis" | name=="malaria" | name=="covid" | name=="zika" | name=="fetal_alcohol" | name=="intellectual_disability"){
  dt_wide[is.na(dt_wide), ] <- 0  #Klinefelter is male only, cysticercosis and malaria are geographically restricted, zika and covid are temporally restricted, so set NaN to 0 
}

if(name=="covid" | name=="fetal_alcohol" | name=="intellectual_disability"){
  temp <- dt_wide %>% keep(~ all(. < 10))
  dt_wide <- cbind(dt_wide[,c("measure_id", "location_id", "short_name", "year_id")], temp)
  dt_ratio <- summaries(dt_wide, count_vars = grep("ratio_count_", names(dt_wide)), draw_vars = grep("ratio_rate_", names(dt_wide)))
  }


if(name!="covid" & name!="fetal_alcohol" | name=="intellectual_disability"){
dt_ratio <- summaries(dt_wide, count_vars = paste0("ratio_count_", 0:nf_draws), draw_vars = paste0("ratio_rate_",0:nf_draws))
}

#Finalize and save
dt_ratio[, location_name := "Global"]
dt_ratio[, measure := ifelse(measure_id==5, "Prevalence", "YLDs")]
dt_ratio[,c("measure_id") := NULL]
dt_ratio[, description := "Female to male ratio in counts and age-standardized rates"]

write.csv(dt_ratio, paste0(FILEPATH), row.names = F)

