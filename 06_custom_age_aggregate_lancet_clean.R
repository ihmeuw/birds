#########################################################################################################################
## Collapse age categories (under 5, 5-19, 20-39, etc.)
## Output: csv files of PAFs, plots of attributable burden
#########################################################################################################################

## Age categories --------------------------------------------------------------------------------------------------------
# Under 5
# 5 to 19 years
# 20-59 years
# 60-79 years
# 80+ years

## Set up ----------------------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")

#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id,22)

loc_df <- read.csv(paste0(save_dir, "locations.csv"))
loc_df <- loc_df %>% filter(level %in% loc_level)
loc_ids <- loc_df %>% filter(level %in% loc_level) %>% pull(location_id)

pop_df <- read.csv(paste0(save_dir, "population.csv"))

## Age aggregation function ----------------------------------------------------------------------------------------------

age_aggregate <- function(data_dt){
  temp_dt <- copy(data_dt)
  
  #Merge population
  temp_dt[, c("population")] <- NULL
  temp_dt <- as.data.table(merge(temp_dt, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                                 by = c("location_id", "sex_id", "year_id", "age_group_id"), all.x=T))
  
  setnames(temp_dt, "cause_name", "short_name")

  #subset to custom age groupings and and population 
  under5_dt <- temp_dt[age_group_id== ID | age_group_id== ID | age_group_id==ID  | age_group_id>ID, ] #under 5
  child_dt <- temp_dt[age_group_id>ID & age_group_id<ID, ] #5-19 years
  adult_dt <- temp_dt[age_group_id>ID & age_group_id<ID, ] #20-59 years
  older_dt <- temp_dt[age_group_id>ID & age_group_id<ID, ] #60-79 years
  over80_dt <- temp_dt[age_group_id>ID & age_group_id<ID | age_group_id==ID, ] #80+ years
  
  #collapse ages and population and combine into one data table
  under5_dt <- under5_dt[, lapply(.SD, sum, na.rm=TRUE), 
                   by=c("location_id", "year_id", "sex_id", "measure_id", "short_name"), .SDcols = c(draw_vars, "population")]
  under5_dt[, age_group_id := ID]
  under5_dt[, age_group_name := "Under 5"]
  
  child_dt <- child_dt[, lapply(.SD, sum, na.rm=TRUE), 
                         by=c("location_id", "year_id", "sex_id", "measure_id", "short_name"), .SDcols = c(draw_vars, "population")]
  child_dt[, age_group_id := ID]
  child_dt[, age_group_name := "5 to 19 years"]
  
  adult_dt <- adult_dt[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "sex_id", "measure_id", "short_name"), .SDcols = c(draw_vars, "population")]
  adult_dt[, age_group_id := ID]
  adult_dt[, age_group_name := "20 to 59 years"]
  
  older_dt <- older_dt[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "sex_id", "measure_id", "short_name"), .SDcols = c(draw_vars, "population")]
  older_dt[, age_group_id := ID]
  older_dt[, age_group_name := "60 to 79 years"]
  
  over80_dt <- over80_dt[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "sex_id", "measure_id", "short_name"), .SDcols = c(draw_vars, "population")]
  over80_dt[, age_group_id := ID]
  over80_dt[, age_group_name := "80+ years"]
  
  combo_dt <- rbind(under5_dt, child_dt, adult_dt, older_dt, over80_dt, fill=T)
  
  #produce rates by dividing by population
  combo_dt[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
  
  #create summaries
  summary_dt <- summaries(combo_dt, rate_vars = rate_vars, draw_vars = draw_vars)
  
  #save files
  saveRDS(combo_dt, paste0(FILEPATH))
  write.csv(summary_dt, paste0(FILEPATH), row.names=F)
  
  return(combo_dt)
}

## Summarize function ----------------------------------------------------------------------------------------------------
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
  sum[, c(draw_vars, rate_vars) := NULL]
  return(sum)
}


## Get all files and paste together --------------------------------------------------------------------------------------
#DALYS for all causes and total neuro, all sexes
dt = readRDS(paste0(FILEPATH))
subset_dt <- dt[sex_id==ID & age_group_id!=ID & age_group_id!=ID, ]


combo_dt = age_aggregate(data_dt = subset_dt)
