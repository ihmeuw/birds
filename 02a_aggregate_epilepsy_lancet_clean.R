################################################################################
## Child script to produce all epilepsy minus other conditions in analysis, grouped by age, sex, location
## Output: csv files of cause-level epilepsy aggregate for prevalence and ylds - age-specific, 
## all-age, age-standardized rates, counts (draw files and summary files)
################################################################################

## Set up ------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")


#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id)

loc_df <- read.csv(paste0(save_dir, "locations.csv"))
loc_df <- loc_df %>% filter(level %in% loc_level)
loc_ids <- loc_df %>% filter(level %in% loc_level) %>% pull(location_id)

pop_df <- read.csv(paste0(save_dir, "population.csv"))

name = "epilepsy"

#Functions
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
  dt_tmp <- dt_tmp[age_group_id!=ID &age_group_id!=ID & age_group_id!=ID,]
  age_weights <- get_age_metadata(ID, gbd_round_id=ID)
  dt_tmp[, c(count_vars, "population")] <- NULL
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


#Get epilepsy impairment draws
ep_imp_us <- get_draws(gbd_id_type = "rei_id",
                    gbd_id = ID,
                    source = "como",
                    release_id = RELEASE_ID,
                    version_id = como_id,
                    location_id = loc_ids,
                    year_id = years,
                    age_group_id = ages,
                    sex_id = c(ID),
                    measure_id = c(ID),
                    metric_id = ID,
                    num_workers = 10)
 
setDT(ep_imp) 

##Convert to cases
dt <- copy(ep_imp)
dt[,c("version_id", "rei_id")] <- NULL
dt <- merge(dt, pop_df[c("population", "year_id", "sex_id", "location_id", "age_group_id")])
dt[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) * population)]

##Produce both sex and all age cases
both_sex <- dt[, lapply(.SD, sum, na.rm=TRUE), 
                   by=c("location_id", "year_id", "age_group_id", "measure_id", "metric_id", "cause_id"), .SDcols = c(draw_vars, "population")]    
both_sex[, sex_id := ID]
dt <- rbind(dt, both_sex, fill=T)

all_age <- dt[, lapply(.SD, sum, na.rm=TRUE), 
                   by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id", "cause_id"), .SDcols = c(draw_vars, "population")]
all_age[, age_group_id := ID]
dt <- rbind(dt, all_age, fill=T)

##Subtract prevalence and YLDs from ten causes with epilepsy: meningitis, encephalitis, tetanus, 
# sepsis, jaundice, preterm birth, NCC, neonatal encephalopathy, malaria, echino
#Subset to causes and total epilepsy, aggregate for the causes
id_list = c(IDS)
total_dt <- dt[cause_id==ID, ]
cause_dt <- dt[cause_id %in% id_list, ]

cause_dt <- cause_dt[, lapply(.SD, sum, na.rm=TRUE), 
                              by=c("location_id", "year_id", "sex_id", "measure_id", "age_group_id", "population"), .SDcols = c(draw_vars)]

#Subtract cause aggregate from total to get remainder for YLDs and prevalence
total_vars <- paste0("total_", 0:nf_draws)
setnames(total_dt, draw_vars, total_vars)
temp <- merge(total_dt, cause_dt, by=c("location_id", "year_id", "sex_id", "measure_id", "age_group_id", "population"))
temp[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("total_", x)) - get(paste0("draw_", x)))]

#Aggregate early + late neonatal to "neonatal", aggregate 1-5 months + 6-11 months to "post neonatal"
dt_neo <- temp[age_group_id==ID | age_group_id==ID, ]
dt_post <- temp[age_group_id==ID | age_group_id==ID, ]

dt_neo <- dt_neo[, lapply(.SD, sum, na.rm=TRUE), 
                 by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
dt_post <- dt_post[, lapply(.SD, sum, na.rm=TRUE), 
                   by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
dt_neo[, age_group_id := ID]
dt_post[, age_group_id := ID]

temp <- rbind(temp, dt_neo, dt_post, fill=TRUE)

#Divide by population to get rates
temp[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
temp[,c(total_vars, "rei_id", "metric_id", "version_id", "cause_id")] <- NULL
temp[, short_name := "epilepsy"]


##Create age-standardized rates and append
final_agestd <- calc_agestd_rates(data_dt=temp)
final_agestd_summary <- copy(final_agestd)
final_agestd_summary[ , c(rate_vars)] <- NULL
final_agestd[, c("mean_rate", "lower_rate", "upper_rate", "standard_error_rate")] <- NULL


##Summarize and save draws and summary estimates for epilepsy
temp_summary <- summaries(temp, draw_vars=draw_vars, rate_vars=rate_vars)

final <- rbind(temp, final_agestd, fill=TRUE)
final_summary <- rbind(temp_summary, final_agestd_summary, fill=TRUE)
final[,c("population")] <- NULL
final_summary[,c("population")] <- NULL


saveRDS(final, paste0(FILEPATH))
write.csv(final_summary, paste0(FILEPATH), row.names=F)
