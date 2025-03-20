#########################################################################################################################
## Child script to pull draws for neuro causes with fatal outcomes by age, sex, location, year for a given CoDCorrect run
## Output: rds files of DALYs and deaths (rates and counts) draws and csv summary files
#########################################################################################################################


## Set up ------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")

#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id,22)

loc_df <- read.csv(paste0(save_dir, "locations.csv"))
loc_df <- loc_df %>% filter(level %in% loc_level)

pop_df <- read.csv(paste0(save_dir, "population.csv"))

#Cause IDs with fatal component
id_list <- c(ID)

#Summaries function for draws
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
  rate_vars <- paste0("rate_", 0:f_draws)
  case_vars <- paste0("draw_", 0:f_draws)
  dt_tmp <- as.data.table(copy(data_dt))
  dt_tmp <- dt_tmp[dt_tmp$age_group_id!=ID & dt_tmp$age_group_id!=ID & dt_tmp$age_group_id!=ID & dt_tmp$age_group_id!=ID,]
  age_weights <- get_age_metadata(ID, gbd_round_id=ID)
  dt_tmp[, c(case_vars) := NULL]
  dt_tmp <- merge(dt_tmp, age_weights[, c("age_group_id", "age_group_weight_value")], by = "age_group_id", all.x=T)
  dt_tmp[, (rate_vars) := lapply(.SD, function(x) sum(x * age_group_weight_value)), by = c("location_id", "year_id", "sex_id", "measure_id", "cause_id", "cause_name"), .SDcols = rate_vars]
  dt_tmp <- unique(dt_tmp, by = c("location_id", "year_id", "sex_id", "measure_id", "cause_id", "cause_name"))
  dt_tmp[,`:=` (age_group_weight_value = NULL, age_group_id = ID)]
  dt_tmp[, mean_rate := rowMeans(.SD), .SDcols = rate_vars]
  dt_tmp[, lower_rate := apply(.SD, 1, quantile, probs= 0.025), .SDcols = rate_vars]
  dt_tmp[, upper_rate := apply(.SD, 1, quantile, probs= 0.975), .SDcols = rate_vars]
  return(dt_tmp)
}

#Locations for countries and other
locs <- loc_df %>% filter(level %in% loc_level) %>% pull(location_id)

#Function to get death and YLL draws for specified locations and id list, and aggregate to total deaths and YLLs 
fatal_draws <- function(locs, id_list, gbd_round, cod_id, daly_id, ages, pop_df, death_dir){
  #First, get draws of death, YLL, and DALY counts for all fatal neuro causes, for each age where deaths occur, all sexes, all specified locations
  print(paste0("Starting deaths get_draws call"))
  temp_death <- get_draws(gbd_id_type = "cause_id",
                    gbd_id = id_list,
                    source = "codcorrect",
                    gbd_round_id = gbd_round,
                    version_id = cod_id,
                    location_id = locs,
                    year_id = years,
                    age_group_id = ages,
                    sex_id = c(ID),
                    measure_id = c(ID),
                    metric_id = ID,
                    decomp_step = "STEP",
                    downsample = T,
                    n_draws = 500,
                    num_workers = 12)

  print(paste0("Starting dalys get_draws call"))
  temp_daly <- get_draws(gbd_id_type = "cause_id",
                         gbd_id = id_list,
                         source = "dalynator",
                         gbd_round_id = gbd_round,
                         version_id = daly_id,
                         location_id = locs,
                         year_id = years,
                         age_group_id = ages,
                         sex_id = c(ID),
                         measure_id = ID,
                         metric_id = ID,
                         decomp_step = "STEP",
                         downsample = T,
                         n_draws = 500,
                         num_workers = 15)
  
  temp <- plyr::rbind.fill(temp_death, temp_daly)

  #Save intermediary files just in case
  saveRDS(temp, paste0(death_dir, "intermediary_get_draws_dalys_deaths_500_draws.rds"))
  print(paste0("Intermediary drawes saved"))
  temp = readRDS(paste0(death_dir, "intermediary_get_draws_dalys_deaths_500_draws.rds"))
  
  #Second, merge on cause name using cause ID
  causes <- get_ids("cause")
  temp <-merge(temp, causes[,c("cause_id", "cause_name")], by = c("cause_id"))
  
  #Second Part 2, aggregate CNS cancers and neuroblastoma/PNS tumors into an aggregate "nervous system cancer" category
  setDT(temp)
  temp_cancer <- temp[cause_id == ID | cause_id == ID, ]
  temp_no_cancer <- temp[cause_id != ID & cause_id != ID, ]
  temp_cancer[,`:=` (cause_name = "Nervous system cancer", cause_id = 477)] #create aggregate name
  temp_cancer <- temp_cancer[, lapply(.SD, sum, na.rm=TRUE), 
                 by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id", "age_group_id", "cause_name", "cause_id"), .SDcols = c(draw_vars)]
  temp <- rbind(temp_no_cancer, temp_cancer, fill=TRUE)
  temp[,c("metric_id", "version_id")] <- NULL
  
  rm(temp_cancer, temp_no_cancer)
  
  #Third, merge on population, create custom age groups (neonatal and postnatal) and produce rates
  temp <- as.data.table(merge(temp, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                              by = c("location_id", "sex_id", "year_id", "age_group_id"), all.x = T))
  
  temp_neo <- temp[age_group_id==ID | age_group_id==ID, ]
  temp_post <- temp[age_group_id==ID | age_group_id==ID, ]
  temp_neo <- temp_neo[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "sex_id", "measure_id", "cause_id", "cause_name"), .SDcols = c(draw_vars, "population")]
  temp_post <- temp_post[, lapply(.SD, sum, na.rm=TRUE), 
                         by=c("location_id", "year_id", "sex_id", "measure_id", "cause_id", "cause_name"), .SDcols = c(draw_vars, "population")]
  temp_neo[, paste0("rate_", 0:f_draws) := lapply(0:f_draws, function(x) get(paste0("draw_", x)) / population)]
  temp_post[, paste0("rate_", 0:f_draws) := lapply(0:f_draws, function(x) get(paste0("draw_", x)) / population)]
  temp_neo[, age_group_id := ID]
  temp_post[, age_group_id := ID]
  temp <- rbind(temp, temp_neo, temp_post, fill=TRUE)
  
  temp[, paste0("rate_", 0:f_draws) := lapply(0:f_draws, function(x) get(paste0("draw_", x)) / population)]
  
  rm(temp_neo, temp_post)
  
  #Fourth, aggregate to total YLLs and deaths for all neuro causes using counts, and produce rates
  total_fatal_dt <- as.data.table(copy(temp))
  #retain "stroke" cause ID but not individual types for creating aggregate; do not aggregate DALYs (done later in script)
  total_fatal_dt <- total_fatal_dt[cause_id!=ID & cause_id!=ID & cause_id!=ID & measure_id!=ID,] 
  total_fatal_dt[,c(rate_vars) := NULL] 
  total_fatal_dt <- total_fatal_dt[, lapply(.SD, sum, na.rm=TRUE), 
           by=c("location_id", "year_id", "sex_id", "age_group_id", "population", "measure_id"), .SDcols = c(draw_vars)]
  total_fatal_dt[, paste0("rate_", 0:f_draws) := lapply(0:f_draws, function(x) get(paste0("draw_", x)) / population)]
  total_fatal_dt[, cause_name := "Total neuro"]
  total_fatal_dt[, cause_id := ID]
  
  #Fourth Part 2, create DALYs aggregate using DALYs for fatal causes and YLDs for non-fatal only causes (where YLDs=DALYs)
  daly_dt <- as.data.table(copy(temp))
  daly_dt <- daly_dt[measure_id==ID & cause_id!=ID & cause_id!=ID & cause_id!=ID & cause_id!=ID, ]

  files_yld <- list.files(file.path(agg_dir_grp2))
  files_yld <- grep(".rds", files_yld, value=TRUE)
  temp_list = lapply(paste0(agg_dir_grp2, files_yld), function (x) data.table(readRDS(x)))
  yld_dt = rbindlist(temp_list, fill = TRUE)
  yld_list = c("adhd", "autism", "epilepsy", "fetal_alcohol", "gbs", "tt_headache", "migraine", 
               "intellectual_disability", "sci", "tbi", "congenital_grp2", "covid", "down",
               "echinococcosis", "klinefelter", "malaria", "neo_jaundice", "neonatal_sepsis", 
               "other_chrom", "neo_preterm", "syphilis", "zika", "diabetes")
  yld_dt <- yld_dt[short_name %in% yld_list,  ]
  yld_dt <- yld_dt[measure_id==ID, ]
  
  rep_str = c("epilepsy" = "Epilepsy", "migraine" = "Migraine", "tt_headache" = "Tension-type headache", "fetal_alcohol" = "Fetal alcohol syndrome",
              "autism"="Autism spectrum disorder", "neonatal_sepsis"="Neonatal sepsis", "tbi" = "Traumatic brain injury", "sci" = "Spinal cord injury", 
              "intellectual_disability" = "Idiopathic intellectual disability", "adhd" = "AD/HD", "covid"="Covid-19", 
              "malaria"="Malaria", "gbs"="Guillain-Barre Syndrome", "syphilis"="Syphilis", "echinococcosis"="Echinococcosis", "zika"="Zika virus", 
              "down"="Down syndrome", "klinefelter"="Klinefelter", "congenital_grp2"="Congenital birth defects", "neo_jaundice"="Neonatal jaundice",
              "neo_preterm" = "Preterm birth", "other_chrom"="Chromosomal abnormalities", "diabetes" = "Diabetes")
  yld_dt[, short_name := str_replace_all(short_name, rep_str)]
  setnames(yld_dt, "short_name", "cause_name") 
  
  #Fourth Part 2.5, aggregate YLDS for epilepsy impairment and YLLs for idiopathic epilepsy to get total epilepsy DALYs 
  epilepsy_yld = yld_dt[cause_name == "Epilepsy" & measure_id==ID & age_group_id!=ID, ]
  epilepsy_yld[,c(rate_vars, "population")] <- NULL
  epilepsy_yll = temp[temp$cause_id == ID & measure_id==ID, ]
  epilepsy_yll[,c(rate_vars, "cause_id", "cause_name", "measure_id")] <- NULL
  setnames(epilepsy_yld, draw_vars, paste0("yld_", 0:nf_draws))
  
  epilepsy_daly <- merge(epilepsy_yld, epilepsy_yll, by=c("location_id", "year_id", "sex_id", "age_group_id"))
  epilepsy_daly[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("yld_", x)) + get(paste0("draw_", x)))]
  epilepsy_daly[,c(paste0("yld_", 0:nf_draws))] <- NULL
  epilepsy_daly[, paste0("rate_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
  epilepsy_daly[, measure_id :=ID]
  epilepsy_daly[, cause_id := ID]
  epilepsy_daly[, cause_name := "Epilepsy"]
  
  yld_dt = yld_dt[cause_name != "Epilepsy", ]
  
  #Add updated epilepsy DALYs to temp
  temp <- temp[ !(measure_id == ID & cause_id == ID),] 
  temp = rbind(temp, epilepsy_daly, fill=T)
  temp[cause_name=="Idiopathic epilepsy", cause_name := "Epilepsy"]
    
  #Fourth Part 2 final, combine YLDs and DALYs datasets to get DALYs for each of the 36 conditions 
  daly_dt <- rbind(daly_dt, yld_dt, epilepsy_daly, fill = TRUE) 
  cause_daly <- copy(daly_dt) #to save all dalys by individual cause across all conditions
  cause_daly[, measure_id := ID]
  
  daly_dt[,c(rate_vars, "cause_id", "cause_name", "measure_id", "short_name") := NULL]
  daly_dt <- daly_dt[age_group_id!=ID, ]
  daly_dt <- daly_dt[, lapply(.SD, sum, na.rm=TRUE), 
                    by=c("location_id", "year_id", "sex_id", "age_group_id", "population"), .SDcols = c(draw_vars)]
  daly_dt[, paste0("rate_", 0:f_draws) := lapply(0:f_draws, function(x) get(paste0("draw_", x)) / population)]
  daly_dt[, cause_name := "Total neuro"]
  daly_dt[,cause_id := ID]
  daly_dt[,measure_id := ID]
  
  #Fourth, Part 3, save draws for total neuro DALYs and each cause (combined nonfatal + fatal and nonfatal only) for all locs and age groups
  combo_daly <- rbind(cause_daly, total_fatal_dt, daly_dt, fill=T)
  combo_daly <- merge(combo_daly, loc_df[,c("location_id", "level", "lancet_label")], by = c("location_id"))
  combo_daly[,c("cause_id", "level")] <- NULL
  saveRDS(combo_daly, paste0(agg_dir, "granular_daly_draws_cod", cod_id, "_como", como_id, ".rds"))
  
  combo_daly_total_neuro <- combo_daly[cause_name=="Total neuro"]
  combo_daly_total_neuro[, cause_id := ID]
  combo_daly_agestd <- calc_agestd_rates(data_dt=combo_daly_total_neuro)
  combo_daly_agestd[, paste0("rate_", 0:nf_draws) := NULL]
  
  
  #Fifth, create age-standardised death rates for individual causes and total neuro
  agestd_dt <- rbind(temp, total_fatal_dt, daly_dt, fill=T)
  agestd_draw_dt <- calc_agestd_rates(data_dt=agestd_dt)
  agestd_sum_dt <- copy(agestd_draw_dt)
  agestd_sum_dt[,c(rate_vars, paste0("yld_", 0:nf_draws))] <- NULL
  agestd_sum_dt <- merge(agestd_sum_dt, loc_df[,c("location_id", "lancet_label", "level")], by = c("location_id"), all.x=T)
  
  saveRDS(agestd_draw_dt, paste0(FILEPATH))
  print(paste0("Draws of age-standardized deaths saved for id list"))
  
  write.csv(agestd_sum_dt, paste0(FILEPATH), row.names = F)
  print(paste0("Summary of age-standardized deaths saved for id list"))
  
  #Sixth, combine aggregate and individual causes into one data set, collapse draws, and save draws and summaries
  final_draw_dt <- plyr::rbind.fill(temp, total_fatal_dt, daly_dt)
  final_draw_dt <- as.data.table(merge(final_draw_dt, loc_df[,c("location_id", "lancet_label", "level")], by = c("location_id"), all.x=T))
  final_draw_country <- final_draw_dt[level==3 & age_group_id==ID, ]
  final_draw_other <- final_draw_dt[level!=3 & age_group_id!=ID & age_group_id!=ID, ]
  final_sum_country <- summaries(dt = final_draw_country, draw_vars = draw_vars, rate_vars = rate_vars)
  final_sum_other <- summaries(dt = final_draw_other, draw_vars = draw_vars, rate_vars = rate_vars)
  
  rm(final_draw_country, final_draw_other)
  final_sum_country[, paste0("yld_", 0:nf_draws) := NULL]
  final_sum_other[, paste0("yld_", 0:nf_draws) := NULL]
  
  saveRDS(final_draw_dt, paste0(FILEPATH))
  print(paste0("Draws of deaths and YLLs saved for id list"))
  
  write.csv(final_sum_country, paste0(FILEPATH), row.names = F)
  write.csv(final_sum_other, paste0(FILEPATH), row.names = F)
  print(paste0("Summary of deaths and YLLs saved for id list"))
  
  #Seventh, calculate percent change from draw files and save summaries
  final_draw_chg <- copy(final_draw_dt)
  final_draw_chg <- final_draw_chg[final_draw_chg$age_group_id==ID, ]
  final_draw_chg[,c(rate_vars, "age_group_id", "population")] <- NULL
  agestd_draw_dt[,c("age_group_id", "population", "mean_rate", "lower_rate", "upper_rate")] <- NULL
  change_dt <- as.data.table(merge(final_draw_chg, agestd_draw_dt), by=c("location_id", "sex_id", "year_id", "cause_id", "cause_name", "measure_id"))
  
  #make wide format by year
  change_wide <- change_dt[(year_id==year_start | year_id==year_end), ]
  change_wide[, c("cause_id", "lancet_label")] <- NULL
  change_wide <- reshape(change_wide, idvar = c("measure_id", "location_id", "sex_id", "cause_name"), timevar = "year_id", direction = "wide")
  
  #divide year_end by year_start cases, then subtract 1 and multiply by 100 to get % change
  change_wide[, paste0("change_count_", 0:f_draws) := lapply(0:f_draws, function(x) ((get(paste0("draw_", x, ".", year_end)) / get(paste0("draw_", x, ".", year_start))-1)*100))]
  change_wide[, paste0("change_rate_", 0:f_draws) := lapply(0:f_draws, function(x) ((get(paste0("rate_", x, ".", year_end)) / get(paste0("rate_", x, ".", year_start))-1)*100))]
  change_wide_keep <- change_wide[, c(paste0("change_count_", 0:f_draws), paste0("change_rate_", 0:f_draws), "measure_id", "location_id", "sex_id", "cause_name"), with=FALSE]
  
  #Change NA to 0 for neurocysticercosis because of geographic restrictions
  change_wide_keep[is.na(change_wide_keep), ] <- 0 
  
  #collapse
  change_final <- summaries(change_wide_keep, draw_vars = paste0("change_count_", 0:f_draws), rate_vars = paste0("change_rate_",0:f_draws))
  
  #finalize and save
  change_final[, measure := ifelse(measure_id==1, "Deaths", ifelse(measure_id==4, "YLLs", "DALYs"))]
  change_final[, description := paste0("Percent change (out of 100) in counts and age-standardized rates from ", year_start, " to ", year_end)]
  
  write.csv(change_final, paste0(FILEPATH))
  print(paste0("Percent change saved"))
    
  #Eighth, calculate female to male ratio from draw files and save summaries
  rate_vars = paste0("rate_", 0:f_draws)
  sex_dt <- as.data.table(readRDS(paste0(FILEPATH)))
  sex_dt[,c("mean_rate", "lower_rate", "upper_rate") := NULL]
  sex_dt <- sex_dt[sex_id!=ID, ]
  
  #Make wide format by sex
  sex_dt <- reshape(sex_dt, idvar = c("measure_id", "location_id", "year_id", "cause_name"), timevar = "sex_id", direction = "wide")
  
  #Divide female by male rates
  sex_dt[, paste0("ratio_rate_", 0:f_draws) := lapply(0:f_draws, function(x) (get(paste0("rate_", x, ".2")) / get(paste0("rate_", x, ".1"))))]
  sex_dt <- sex_dt[,c(paste0("ratio_rate_", 0:f_draws), "measure_id", "location_id", "year_id", "cause_name"), with=FALSE]
  
  #Change NA to 0 for cysticercosis because of geographic restrictions
  sex_dt[is.na(sex_dt), ] <- 0 
  
  #Collapse
  rate_vars <- paste0("ratio_rate_", 0:f_draws)
  sex_dt[, mean_rate := rowMeans(.SD), .SDcols = rate_vars]
  sex_dt[, lower_rate := apply(.SD, 1, quantile, probs= 0.025), .SDcols = rate_vars]
  sex_dt[, upper_rate := apply(.SD, 1, quantile, probs=0.975), .SDcols = rate_vars]
  sex_dt[, c(rate_vars) := NULL]
  
  #Finalize and save
  sex_dt <- sex_dt[location_id==1,]
  sex_dt[, location_name := "Global"]
  sex_dt[, measure := ifelse(measure_id==ID, "Deaths", ifelse(measure_id==ID, "DALYs", "YLLs"))]
  sex_dt[, description := "Female to male ratio for age-standardized rates"]
  
  write.csv(sex_dt, paste0(death_dir, "sex_ratio_fatal.csv"), row.names = F)
  print("Sex ratios saved for deaths and YLLs")
  
  #Ninth, save DALY summaries together in aggregate summary folder
  #causes with fatal component plus total neuro dalys
  daly_sum <- rbind(final_sum_other[measure_id==ID], final_sum_country[measure_id==ID], fill=TRUE)
  daly_agestd <- agestd_sum_dt[measure_id==ID]
  daly_chg <- change_final[measure_id==ID]
  daly_chg <- merge(daly_chg, loc_df[,c("location_id", "lancet_label")], by=c("location_id"))
  daly_sex <- sex_dt[measure=="DALYs"]
  daly_sex <- merge(daly_sex, loc_df[,c("location_id", "lancet_label")], by=c("location_id"))
  
  
  #causes with only nonfatal
  yld_list = c("adhd", "autism", "fetal_alcohol", "gbs", "tt_headache", "migraine", "intellectual_disability", 
               "sci", "tbi", "congenital_grp2", "covid", "down",
               "echinococcosis", "klinefelter", "malaria", "neo_jaundice", "neonatal_sepsis", 
               "other_chrom", "neo_preterm", "syphilis", "zika", "diabetes")
  rep_str = c("migraine" = "Migraine", "tt_headache" = "Tension-type headache", 
              "fetal_alcohol" = "Fetal alcohol syndrome","autism"="Autism spectrum disorder", 
              "neonatal_sepsis"="Neonatal sepsis", "tbi" = "Traumatic brain injury", "sci" = "Spinal cord injury", 
              "intellectual_disability" = "Idiopathic intellectual disability", 
              "adhd" = "AD/HD", "covid"="Covid-19", "malaria"="Malaria", 
              "gbs"="Guillain-Barre Syndrome", "syphilis"="Syphilis", 
              "echinococcosis"="Echinococcosis", "zika"="Zika virus", "down"="Down syndrome", 
              "klinefelter"="Klinefelter", "congenital_grp2"="Congenital birth defects", 
              "neo_jaundice"="Neonatal jaundice",
              "neo_preterm" = "Preterm birth", "other_chrom"="Chromosomal abnormalities", "diabetes"="Diabetes")
  
  files_yld_sum <- list.files(file.path(aggsum_dir_grp2))
  temp_list_sum = lapply(paste0(aggsum_dir_grp2, files_yld_sum), function (x) data.table(read.csv(x)))
  yld_sum = rbindlist(temp_list_sum, fill = TRUE)
  yld_sum <- yld_sum[short_name %in% yld_list, ]
  yld_sum <- merge(yld_sum, loc_df[c("location_id", "level")])
  yld_sum_country <- yld_sum[level==3 & measure_id==ID & (age_group_id==ID | age_group_id==ID), ] #ONLY AGE GROUP 22 and 27 FOR COUNTRIES
  yld_sum_oth <- yld_sum[level!=3 & measure_id==ID, ]
  yld_sum <- rbind(yld_sum_oth, yld_sum_country, fill=TRUE)
  yld_sum[, short_name := str_replace_all(short_name, rep_str)]
  setnames(yld_sum, "short_name", "cause_name") 
  
  files_yld_chg <- list.files(file.path(change_dir))
  temp_list_chg = lapply(paste0(change_dir, files_yld_chg), function (x) data.table(read.csv(x)))
  yld_oth = rbindlist(temp_list_chg, fill = TRUE)
  yld_oth <- yld_oth[short_name %in% yld_list & measure=="YLDs", ]
  setnames(yld_oth, "short_name", "cause_name")
  yld_oth[, cause_name := str_replace_all(cause_name, rep_str)]
  yld_oth <- merge(yld_oth, loc_df[,c("location_id", "lancet_label")], by=c("location_id"))
  
  #combine into all causes
  cause_daly <- rbind(daly_sum, daly_agestd, yld_sum, fill=TRUE)
  cause_daly[, sex := ifelse(sex_id==ID, "Male", ifelse(sex_id==ID, "Female", "Both"))]
  cause_daly[, measure := "DALYs"]
  cause_daly[,c("cause_id", "measure_id", "population", "lancet_label") := NULL]
  cause_daly <- merge(cause_daly, loc_df[,c("location_id", "lancet_label")], by=c("location_id"))
  cause_daly_country <- cause_daly[level==3, ]
  cause_daly_other <- cause_daly[level!=3, ]
  write.csv(cause_daly_other, paste0(FILEPATH))
  write.csv(cause_daly_country, paste0(FILEPATH))
  print(paste0("Dalys saved for id list")) 
  
  oth_daly <- rbind(daly_chg, daly_sex, yld_oth, fill=TRUE)
  oth_daly[, sex := ifelse(sex_id==ID, "Male", ifelse(sex_id==ID, "Female", "Both"))]
  oth_daly[, measure := "DALYs"]
  oth_daly[, age_group_name := "All-ages/Age-standardised"]
  oth_daly[,c("measure_id", "location_name") := NULL]
  oth_daly <- merge(oth_daly, loc_df[c("location_id", "lancet_label")], by = c("location_id"))
  write.csv(oth_daly, paste0(FILEPATH))
  print(paste0("Dalys % change and sex ratio saved for id list")) 
    
  return(final_draw_dt)
}

death_dt <- fatal_draws(locs=locs, id_list=id_list, gbd_round=gbd_round, cod_id=cod_id, 
                        daly_id=daly_id, ages=ages, pop_df=pop_df, death_dir=death_dir)  
  