################################################################################
## Child script to produce neuro aggregate groups by age, sex, location
## Output: csv files of cause-level neuro aggregates for prevalence and ylds - 
## age-specific, all-age, age-standardized rates, counts (draw files and summary files)
################################################################################

## Set up ------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")

#Arguments
args <- commandArgs(trailingOnly = T)

agg_group <- args[1]
group_level <- args[2]
name <- args[3]
type <- args[4]

print(paste0(agg_group, ", ", group_level, ", ", name, ", ", type))

#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id,22)

loc_df <- read.csv(paste0(save_dir, "locations.csv"))
loc_df <- loc_df %>% filter(level %in% loc_level)
loc_ids <- loc_df %>% filter(level %in% loc_level) %>% pull(location_id)

pop_df <- read.csv(paste0(save_dir, "population.csv"))

## Get ID map and subset IDs for group2 cause specified by agg_group ----------------------------
neuro_ids <- read.xlsx(paste0(FILEPATH))
neuro_ids$grouping2 <- gsub(" ", "_", neuro_ids$grouping2)
neuro_ids$grouping2 <- gsub("'", "", neuro_ids$grouping2)

group_ids <- neuro_ids %>% filter(grouping2 == agg_group) %>% dplyr::select(gbd_id)
group_ids <- unique(group_ids)
gbd_ids = paste0("^", group_ids$gbd_id, "_")

## Create aggregates ----------------------------------------------------------------------------
#Summarize draws function
summaries <- function(dt, draw_vars, rate_vars){
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

## Age-standardized rates function
calc_agestd_rates <- function(data_dt){
  rate_vars <- paste0("rate_", 0:nf_draws)
  count_vars <- paste0("draw_", 0:nf_draws)
  dt_tmp <- copy(data_dt)
  dt_tmp <- dt_tmp[age_group_id!=ID & age_group_id!=ID & age_group_id!=ID,]
  age_weights <- get_age_metadata(ID, gbd_round_id=ID)
  dt_tmp[, c(count_vars)] <- NULL
  dt_tmp <- merge(dt_tmp, age_weights[, c("age_group_id", "age_group_weight_value")], by = "age_group_id", all.x=T)
  dt_tmp[, (rate_vars) := lapply(.SD, function(x) sum(x * age_group_weight_value)), by = c("location_id", "year_id", "sex_id", "measure_id"), .SDcols = rate_vars]
  dt_tmp <- unique(dt_tmp, by = c("location_id", "year_id", "sex_id", "measure_id"))
  dt_tmp[,`:=` (age_group_weight_value = NULL, age_group_id = ID)]
  dt_tmp[, mean_rate := rowMeans(.SD), .SDcols = rate_vars]
  dt_tmp[, standard_error_rate := apply(.SD, 1, sd), .SDcols = rate_vars]  
  dt_tmp[, lower_rate := apply(.SD, 1, quantile, probs= 0.025), .SDcols = rate_vars]
  dt_tmp[, upper_rate := apply(.SD, 1, quantile, probs= 0.975), .SDcols = rate_vars]
  return(dt_tmp)
}

#Function to aggregate each specified set of IDs
aggregate_groups <- function(group_ids, name, save_dir, como_id){ 
  #Pull desired set of IDs, bind into dataframe and subset to counts for prevalence and YLDs
  if(type=="rei" | type=="sequela"){
    files <- list.files(prev_dir) 
    files_prev <- grep(paste(gbd_ids, collapse="|"), files, value=TRUE)
    files_prev <- grep("prev_como", files_prev, value=TRUE)
    files_prev <- paste(prev_dir, files_prev, sep = "")
    temp_list_prev = lapply(files_prev, function (x) data.table(readRDS(x)))
    dt_prev = rbindlist(temp_list_prev, fill = TRUE)
    dt_prev <- dt_prev[measure_id==ID, ] #safety check
  }
  
  if(type=="cause"){
    files <- list.files(prev_dir) 
    files_prev <- grep(paste(gbd_ids, collapse="|"), files, value=TRUE)
    files_prev <- grep("prev_cause_como", files_prev, value=TRUE)
    files_prev <- paste(prev_dir, files_prev, sep = "")
    temp_list_prev = lapply(files_prev, function (x) data.table(readRDS(x)))
    dt_prev = rbindlist(temp_list_prev, fill = TRUE)  
    dt_prev <- dt_prev[measure_id==ID, ] #safety check
  }
  
  if(name=="cp"){
    files <- list.files(yld_dir) 
    files_yld <- unique(grep(paste(gbd_ids, collapse="|"), files, value=TRUE))
    files_yld <- grep("ylds_cp_como", files_yld, value=TRUE)
    files_yld <- paste(yld_dir, files_yld, sep = "")
    temp_list_yld = lapply(files_yld, function (x) data.table(readRDS(x)))
    dt_yld = rbindlist(temp_list_yld, fill = TRUE)
    dt_yld <- dt_yld[measure_id==ID, ] #safety check
  } 
  
  if(name!="cp" & type!="cause"){
    files <- list.files(yld_dir) 
    files_yld <- unique(grep(paste(gbd_ids, collapse="|"), files, value=TRUE))
    files_yld <- grep("ylds_como", files_yld, value=TRUE)
    files_yld <- paste(yld_dir, files_yld, sep = "")
    temp_list_yld = lapply(files_yld, function (x) data.table(readRDS(x)))
    dt_yld = rbindlist(temp_list_yld, fill = TRUE)
    dt_yld <- dt_yld[measure_id==ID, ] #safety check
  }
  
  if(type=="cause"){
    files <- list.files(yld_dir) 
    files_yld <- grep(paste(gbd_ids, collapse="|"), files, value=TRUE)
    files_yld <- grep("ylds_cause_como", files_yld, value=TRUE)
    files_yld <- paste(yld_dir, files_yld, sep = "")    
    temp_list_yld = lapply(files_yld, function (x) data.table(readRDS(x)))
    dt_yld = rbindlist(temp_list_yld, fill = TRUE)
    dt_yld <- dt_yld[measure_id==ID, ] #safety check
  }
  
  dt <- rbind(dt_prev,dt_yld, fill=T)
  
  #Change Inf to 0
  for (j in 1:ncol(dt)) set(dt, which(is.infinite(dt[[j]])), j, 0)
  
  #Only counts
  dt <- dt[metric_id==ID, ]
  
  #Reappend population as safety check
  dt[,c("population", "version_id") := NULL]
  dt <- as.data.table(merge(dt, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                            by = c("location_id", "sex_id", "year_id", "age_group_id")))

  #Aggregate across sequela/rei and population by age, year, sex, location (can still run this for cause and injury IDs, but nothing will aggregate)
  draw_vars <- paste0("draw_", 0:nf_draws)
  dt <- dt[, lapply(.SD, sum, na.rm=TRUE), 
           by=c("location_id", "year_id", "sex_id", "age_group_id", "population", "measure_id", "metric_id"), .SDcols = c(draw_vars)]
  
  #Aggregate early + late neonatal to "neonatal", aggregate 1-5 months + 6-11 months to "post neonatal"
  dt_neo <- dt[age_group_id==ID | age_group_id==ID, ]
  dt_post <- dt[age_group_id==ID | age_group_id==ID, ]

  dt_neo <- dt_neo[, lapply(.SD, sum, na.rm=TRUE), 
           by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
  dt_post <- dt_post[, lapply(.SD, sum, na.rm=TRUE), 
                   by=c("location_id", "year_id", "sex_id", "measure_id"), .SDcols = c(draw_vars, "population")]
  dt_neo[, age_group_id := ID]
  dt_post[, age_group_id := ID]
  
  dt <- rbind(dt, dt_neo, dt_post, fill=TRUE)
  
  #Revert to rates
  dt[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
  setnafill(dt, fill = 0)
  
  #Create summary dts
  dt_summary <- summaries(dt, draw_vars=draw_vars, rate_vars=rate_vars)  
  
  #Age-standardize and create draws and summary dts
  dt_agestd <- calc_agestd_rates(dt)
  dt_agestd_summary <- copy(dt_agestd)
  dt_agestd_summary[ , c(rate_vars)] <- NULL
  dt_agestd[, c("mean_rate", "lower_rate", "upper_rate", "standard_error_rate")] <- NULL
  
  #Combine age-standardized dts with age-specific and all-age dts
  dt <- rbind(dt, dt_agestd, fill=TRUE)
  dt_summary <- rbind(dt_summary, dt_agestd_summary, fill=TRUE)
  
  #Add short name to each data set to be saved
  dt[,`:=` (short_name = name, version_id = NULL, metric_id = NULL)]
  dt_summary[,`:=` (short_name = name, version_id = NULL, metric_id = NULL)]
  

  #Write to file as csv
  saveRDS(dt, paste0(FILEPATH))
  write.csv(dt_summary, paste0(FILEPATH), row.names=F)
}

aggregate_groups(group_ids = group_ids, name = name, save_dir = save_dir, como_id = como_id)
