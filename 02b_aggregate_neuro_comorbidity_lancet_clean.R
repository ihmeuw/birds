###############################################################################################################################################################
## Child script to produce neuro aggregate groups by age, sex, location, accounting for comorbidity between conditions
## Output: csv files of neuro aggregate groupings prevalence and ylds - all-age, age-specific and age-standardized rates, counts (draw files and summary files)
##         apply co-morbidity correction for prevalence
###############################################################################################################################################################

## Set up ------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")

#Arguments
args <- commandArgs(trailingOnly = T)

agg_group <- args[1]
group_level <- args[2]
name <- args[3]

print(paste0(agg_group, ", ", group_level, ", ", name))

#define name variable for all categories because not in map
if(group_level == "all"){
  name="all_neuro"
}

#confirm arguments
print(paste0(agg_group, ", ", group_level, ", ", name))

#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id,ID)

loc_df <- read.csv(paste0(save_dir, "locations.csv"))
loc_ids <- loc_ids$location_id

pop_df <- read.csv(paste0(save_dir, "population.csv"))

## ID map and subsetted IDs -----------------------------------------------------------------------
neuro_map <- as.data.table(read.xlsx(paste0(FILEPATH)))
neuro_map[, grouping0 := gsub(" ", "_", grouping0)] 
neuro_map[, grouping1 := gsub(" ", "_", grouping1)]

if(group_level=="all"){ #remove mental disorders category but include everything else
  neuro_map <- neuro_map[grouping1!="Mental_disorders", ] 
  gbd_id <- unique(neuro_map$short_name)
  gbd_id <- gbd_id[gbd_id!="cp"]
}
if(group_level=="group0"){
  neuro_map <- neuro_map[grouping0==agg_group, ]
  gbd_id <- unique(neuro_map$short_name)
}
if(group_level=="group1"){
  neuro_map <- neuro_map[grouping1==agg_group, ]
  gbd_id <- unique(neuro_map$short_name)
}

group_ids <- as.data.frame(gbd_id)


## Create aggregates ----------------------------------------------------------------------------
#Summarize draws function
summaries <- function(dt, rate_vars, draw_vars){
  sum <- as.data.table(copy(dt))
  sum[, mean_rate := rowMeans(.SD), .SDcols = rate_vars]
  sum[, standard_error_rate := apply(.SD, 1, sd), .SDcols = rate_vars]  
  sum[, lower_rate := apply(.SD, 1, quantile, probs= 0.025), .SDcols = rate_vars]
  sum[, upper_rate := apply(.SD, 1, quantile, probs=0.975), .SDcols = rate_vars]
  sum[, mean_count := rowMeans(.SD), .SDcols = draw_vars]
  sum[, standard_error_count := apply(.SD, 1, sd), .SDcols = draw_vars]  
  sum[, lower_count := apply(.SD, 1, quantile, probs= 0.025), .SDcols = draw_vars]
  sum[, upper_count := apply(.SD, 1, quantile, probs=0.975), .SDcols = draw_vars]
  sum[, c(rate_vars, draw_vars) := NULL]
  return(sum)
}

## Age-standardized rates function
calc_agestd_rates <- function(data_dt){
  rate_vars <- paste0("rate_", 0:nf_draws)
  draw_vars <- paste0("draw_", 0:nf_draws)
  dt_tmp <- copy(data_dt)
  dt_tmp <- dt_tmp[age_group_id!=ID, ]
  age_weights <- get_age_metadata(ID, gbd_round_id=ID)
  dt_tmp[, c(draw_vars) := NULL]
  dt_tmp <- merge(dt_tmp, age_weights[, c("age_group_id", "age_group_weight_value")], by = "age_group_id", all.x=T)
  dt_tmp[, (rate_vars) := lapply(.SD, function(x) sum(x * age_group_weight_value)), 
         by = c("location_id", "year_id", "sex_id", "measure_id"), .SDcols = rate_vars]
  dt_tmp <- unique(dt_tmp, by = c("location_id", "year_id", "sex_id", "measure_id"))
  dt_tmp[,`:=` (age_group_weight_value = NULL, age_group_id = ID)]
  dt_tmp[, mean_rate := rowMeans(.SD), .SDcols = rate_vars]
  dt_tmp[, standard_error_rate := apply(.SD, 1, sd), .SDcols = rate_vars]  
  dt_tmp[, lower_rate := apply(.SD, 1, quantile, probs= 0.025), .SDcols = rate_vars]
  dt_tmp[, upper_rate := apply(.SD, 1, quantile, probs= 0.975), .SDcols = rate_vars]
  return(dt_tmp)
}

#Function to aggregate each specified set of IDs and account for condition co-morbidity from this equation: 
# prev_total = 1 - (1-prev1)*(1-prev2)*(1-prev3)...
aggregate_groups <- function(group_ids, name, como_id){ 
  #Pull desired set of causes (prevalence and YLDs), bind into dataframe and subset to rates
  files <- paste(agg_dir_grp2, unique(group_ids$gbd_id), paste0("_como", como_id, "_nonfatal.rds"), sep = "")
  temp_list = lapply(files, function (x) data.table(readRDS(x)))
  dt = rbindlist(temp_list, fill = TRUE)
  
  #Create copy
  dt <- as.data.table(copy(dt))

  #Remove and re-append population as safety check, remove age-standardized rates
  dt[, c("population") := NULL]
  dt <- as.data.table(merge(dt, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                            by = c("location_id", "sex_id", "year_id", "age_group_id")))
  dt <- dt[age_group_id!=ID, ]
  
  #Recalculate rates accounting for comorbidity of conditions by age, sex, year, location  
    #Subset to prevalence (keep rates) and YLDs (keep counts)
    dt_sub <- dt[measure_id==ID, ]
    dt_sub[, c(paste0("draw_",0:nf_draws)) := NULL]
    dt_oth <- dt[measure_id==ID, ]
    dt_oth[, c(paste0("rate_",0:nf_draws)) := NULL]
    
    #For prevalence co--morbidity correction, first produce 1-prevalence for each draw
    dt_sub[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) 1 - get(paste0("rate_", x)))]
    #Then, multiply all rows by each other by draw
    dt_sub <- dt_sub[, lapply(.SD, prod, na.rm=TRUE), 
             by=c("location_id", "year_id", "sex_id", "age_group_id", "population", "measure_id"), .SDcols = c(rate_vars)]  
    #Then, subtract the product from 1 to get final co-morbidity rates
    dt_sub[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) 1 - get(paste0("rate_", x)))]
    #To get final co-morbidity prevalent cases, multiply by population  
    dt_sub[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("rate_", x)) * population)]
    
    #For YLDs, aggregate across causes and population by age, year, sex, location
    dt_oth <- dt_oth[, lapply(.SD, sum, na.rm=TRUE), 
             by=c("location_id", "year_id", "sex_id", "age_group_id", "population", "measure_id"), .SDcols = c(draw_vars)]
    #To get final YLDs rates, divide by population  
    dt_oth[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    
    #Recombine prevalence with YLDs
    dt <- rbind(dt_sub,dt_oth, fill=T)
    
    #Create summary dts
    dt_summary <- summaries(dt, rate_vars = rate_vars, draw_vars = draw_vars)

  #Produce age-standardized rates - create draws and summary dts
  dt_agestd <- calc_agestd_rates(dt)
  dt_agestd_summary <- copy(dt_agestd)
  dt_agestd_summary[ , c(rate_vars) := NULL]
  dt_agestd[, c("mean_rate", "lower_rate", "upper_rate", "standard_error_rate") := NULL]
  
  #Combine age-standardized with age-specific and all age
  dt <- rbind(dt, dt_agestd, fill=TRUE)
  dt_summary <- rbind(dt_summary, dt_agestd_summary, fill=TRUE)
  
  #Add short name to each data set to be saved
  dt[, short_name := name]
  dt_summary[ ,short_name := name]
  
  #Add age group names and location names
  dt = merge(dt, loc_df[,c("location_id", "lancet_label")], by = c("location_id"), all.x=T)
  dt_summary = merge(dt_summary, loc_df[,c("location_id", "lancet_label", "level")], by = c("location_id"), all.x=T)
  
  dt = merge(dt, age_df[,c("age_group_id", "age_group_name")], by=c("age_group_id"), all.x=T)
  dt_summary = merge(dt_summary, age_df[,c("age_group_id", "age_group_name")], by=c("age_group_id"), all.x=T)
  
  #Write to file as csv
  saveRDS(dt, paste0(FILEPATH))
  write.csv(dt_summary, paste0(FILEPATH), row.names=F)
}

aggregate_groups(group_ids = group_ids, name = name, como_id = como_id)
