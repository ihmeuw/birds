##############################################################################################################################
## Child script to pull draws for neuro IDs by age, sex, location, year - computes custom YLDs as needed
## Output: rds files of prevalence, incidence, and YLD rate and count draws 
##############################################################################################################################

## Set up ------------------------------------------------------------------------------------------------
rm(list=ls())

source("FILEPATH")


#Arguments
args <- commandArgs(trailingOnly = T)

id <- args[1]
id_type <- args[2]
region_id <- args[3]

print(paste0(id, ", ", id_type, ", ", region_id))


#Directories and shared functions
functions_dir <- paste0(central_lib, "FILEPATH")
neuro_ids <- as.data.table(read.xlsx(paste0(neuro_dir, "FILEPATH")))
functs <- c("get_draws.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

#Set up demographics
age_df <- read.csv(paste0(save_dir, "ages.csv"))
ages <- c(age_df$age_group_id)
pop_df <- read.csv(paste0(save_dir, "population.csv"))
loc_df <- read.csv(paste0(save_dir, "locations.csv"))

#Set up country ID list for region specified by region_id (argument passed from parent script)
if(region_id!=ID) {
  loc_ids <- loc_df[loc_df$parent_id==region_id, ]
  loc_ids <- c(region_id, loc_ids$location_id)
}
if(region_id==ID) {
  loc_ids <- loc_df[loc_df$parent_id==region_id, ]
  loc_ids <- loc_ids$location_id
}

#Function to get draws for specified locations saved by id
if(id_type=="sequela" & id!=ID & id!=ID & id!=ID & id!=ID){
  neuro_draws <- function(region_id, loc_ids, id, gbd_round, como_id, ages, decomp_step, prev_dir){
    temp <- get_draws(gbd_id_type = "sequela_id",
                      gbd_id = id,
                      source = "como",
                      gbd_round_id = gbd_round,
                      version_id = como_id,
                      location_id = loc_ids,
                      year_id = years,
                      age_group_id = ages,
                      sex_id = c(ID),
                      measure_id = c(ID),
                      metric_id = ID,
                      decomp_step = decomp_step,
                      num_workers = 10)

    #merge on population, produce case draws, append to original rate draws and save; subset to ylds and return ylds
    temp <- as.data.table(merge(temp, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                                by = c("location_id", "sex_id", "year_id", "age_group_id")))
    temp[, c("version_id", "sequela_id") := NULL]
    counts <- copy(temp)
    counts[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) * population)]
    counts[, metric_id := ID]
    
    #aggregate to both sexes and append to rate and count dts
    both_sex <- counts[, lapply(.SD, sum, na.rm=TRUE), 
                      by=c("location_id", "year_id", "age_group_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]    
    both_sex[, sex_id := ID]
    both_sex_rates <- copy(both_sex)  
    both_sex_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    both_sex_rates[, metric_id := ID]
    
    temp2 <- rbind(temp, both_sex_rates, fill=T)
    counts2 <- rbind(counts, both_sex, fill=T)
    
    #aggregate to all ages
    all_age <- counts2[, lapply(.SD, sum, na.rm=TRUE), 
                      by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
    all_age[, age_group_id := ID]
    all_age_rates <- copy(all_age)  
    all_age_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    all_age_rates[, metric_id := ID]
    
    temp_final <- rbind(temp2,counts2, all_age, all_age_rates, fill=T)
    ylds <- temp_final[measure_id==ID, ]
    temp_final <- temp_final[measure_id==ID, ]
    saveRDS(temp_final, paste0(FILEPATH))
    print(paste0("Prevalence saved for ", id, " for region ", region_id))
    return(ylds)
  }
}

if(id_type=="sequela" & (id==ID | id==ID | id==ID | id==ID)){ 
  neuro_draws <- function(region_id, loc_ids, id, gbd_round, como_id, ages, decomp_step, prev_dir){
    temp <- get_draws(gbd_id_type = "sequela_id",
                      gbd_id = id,
                      source = "como",
                      gbd_round_id = gbd_round,
                      version_id = como_id,
                      location_id = loc_ids,
                      year_id = years,
                      age_group_id = ages,
                      sex_id = ID,
                      measure_id = c(ID),
                      metric_id = ID,
                      decomp_step = decomp_step,
                      num_workers = 10)
    
    #merge on population, produce male case draws, 
    temp <- as.data.table(merge(temp, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                                by = c("location_id", "sex_id", "year_id", "age_group_id")))
    temp[, version_id := NULL]
    m_counts <- copy(temp)
    m_counts[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) * population)]
    m_counts[,`:=` ( metric_id = ID, population = NULL)]
    
    #create female and both sex counts, where female counts = 0 and both sex counts = male counts, merge together and re-add population 
    f_counts = copy(m_counts)
    f_counts[, paste0("draw_", 0:nf_draws) := 0]
    f_counts[, sex_id := ID]
    b_counts = copy(m_counts)
    b_counts[, sex_id := ID]
    counts = rbind(m_counts, f_counts, b_counts, fill=T)
    counts[, metric_id := ID]
    counts <- as.data.table(merge(counts, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                                by = c("location_id", "sex_id", "year_id", "age_group_id")))
    
    all_age <- counts[, lapply(.SD, sum, na.rm=TRUE), 
                      by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
    all_age[, age_group_id := ID]
    all_age_rates <- copy(all_age)  
    all_age_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    all_age_rates[, metric_id := ID]
    
    #re-add population and produce rates
    rates <- copy(counts)
    rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    rates[, metric_id := ID]
    
    #append case and rate draws and save; subset to ylds and return ylds
    combo <- rbind(counts, rates, all_age, all_age_rates, fill=T)
    combo[, population := NULL]
    ylds <- combo[measure_id==ID, ]
    temp <- combo[measure_id==ID, ]
    saveRDS(temp, paste0(FILEPATH))
    print(paste0("Prevalence and incidence saved for ", id, " for region ", region_id))
    return(ylds)
  }
}

if(id_type=="cause"){
  neuro_draws <- function(region_id, loc_ids, id, gbd_round, como_id, ages, decomp_step, prev_dir){
    temp <- get_draws(gbd_id_type = "cause_id",
                      gbd_id = id,
                      source = "como",
                      gbd_round_id = gbd_round,
                      version_id = como_id,
                      location_id = loc_ids,
                      year_id = years,
                      age_group_id = ages,
                      sex_id = c(ID),
                      measure_id = c(ID),
                      metric_id = ID,
                      decomp_step = decomp_step,
                      num_workers = 10)
    
    #merge on population, produce case draws, aggregate to all ages, append to original rate draws and save; subset to ylds and return ylds
    temp <- as.data.table(merge(temp, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                                by = c("location_id", "sex_id", "year_id", "age_group_id")))
    temp[, c("version_id", "sequela_id") := NULL]
    counts <- copy(temp)
    counts[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) * population)]
    counts[, metric_id := ID]
    
    #aggregate to both sexes and append to rate and count dts
    both_sex <- counts[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "age_group_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]    
    both_sex[, sex_id := 3]
    both_sex_rates <- copy(both_sex)  
    both_sex_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    both_sex_rates[, metric_id := ID]
    
    temp2 <- rbind(temp, both_sex_rates, fill=T)
    counts2 <- rbind(counts, both_sex, fill=T)
    
    #aggregate to all ages
    all_age <- counts2[, lapply(.SD, sum, na.rm=TRUE), 
                     by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
    all_age[, age_group_id := ID]
    all_age_rates <- copy(all_age)  
    all_age_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    all_age_rates[, metric_id := ID]

    temp_final <- rbind(temp2, counts2, all_age, all_age_rates, fill=T)
    ylds <- temp_final[measure_id==ID, ]
    temp_final <- temp_final[measure_id==ID, ]
    saveRDS(temp_final, paste0(FILEPATH))
    print(paste0("Prevalence saved for ", id, " for region ", region_id))
    return(ylds)
  }
}

if(id_type=="rei"){
  #Function to get draws for specified locations saved by REI id
  neuro_draws <- function(region_id, loc_ids, id_list, gbd_round, como_id, ages, decomp_step, prev_dir){
    temp <- get_draws(gbd_id_type = "rei_id",
                      gbd_id = id,
                      source = "como",
                      gbd_round_id = gbd_round,
                      version_id = como_id_injuries,
                      location_id = loc_ids,
                      year_id = years,
                      age_group_id = ages,
                      sex_id = c(ID),
                      measure_id = c(ID),
                      metric_id = ID,
                      decomp_step = decomp_step,
                      num_workers = 20)
    temp <- temp[temp$cause_id==ID, ]

    #merge on population, produce case draws, append to original rate draws and save; subset to ylds and return ylds
    temp <- as.data.table(merge(temp, pop_df[,c("location_id", "sex_id", "year_id", "age_group_id", "population")], 
                                by = c("location_id", "sex_id", "year_id", "age_group_id")))
    temp[, c("version_id", "sequela_id") := NULL]
    counts <- copy(temp)
    counts[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) * population)]
    counts[, metric_id := ID]
    
    #aggregate to both sexes and append to rate and count dts
    both_sex <- counts[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "age_group_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]    
    both_sex[, sex_id := ID]
    both_sex_rates <- copy(both_sex)  
    both_sex_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    both_sex_rates[, metric_id := ID]
    
    temp2 <- rbind(temp, both_sex_rates, fill=T)
    counts2 <- rbind(counts, both_sex, fill=T)
    
    #aggregate to all ages
    all_age <- counts2[, lapply(.SD, sum, na.rm=TRUE), 
                       by=c("location_id", "year_id", "sex_id", "measure_id", "metric_id"), .SDcols = c(draw_vars, "population")]
    all_age[, age_group_id := ID]
    all_age_rates <- copy(all_age)  
    all_age_rates[, paste0("draw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / population)]
    all_age_rates[, metric_id := ID]
    
    temp_final <- rbind(temp2, counts2, all_age, all_age_rates, fill=T)
    ylds <- temp_final[measure_id==ID, ]
    temp <- temp_final[measure_id==ID, ]
    saveRDS(temp_final, paste0(FILEPATH))
    print(paste0("Prevalence saved for ", id, " for region ", region_id))
    return(ylds)
  }
}

ylds <- neuro_draws(region_id=region_id, loc_ids=loc_ids, id=id, gbd_round=gbd_round, decomp_step=decomp_step, como_id=como_id, 
                    ages=ages, prev_dir=prev_dir)

#Pull or calculate YLDs
# If no custom YLD needed because health state has only neurological components and not for cerebral palsy
neuro_id_subset <- neuro_ids[short_name != "cp" & gbd_id==id, ]
neuro_id_subset <- unique(neuro_id_subset[,c("detailed_cause", "sequela_name", "gbd_id", "diff_dw", "gbd_match_id", "type")])
neuro_id_subset <- neuro_id_subset[type==id_type, ]

if(length(neuro_id_subset!=0)){
  if(neuro_id_subset$diff_dw==0 & (id_type=="sequela" | id_type=="rei")){
    saveRDS(ylds, paste0(yld_dir, id, "_loc", region_id, "_ylds_como", como_id, ".rds"))
    print(paste0("YLDs saved for ", id, " for region ", region_id))
  }

  if(neuro_id_subset$diff_dw==0 & id_type=="cause"){
    saveRDS(ylds, paste0(FILEPATH))
    print(paste0("YLDs saved for ", id, " for region ", region_id))
  }
}

#Function for custom ylds if combined health state includes non-neurological components
#pull prevalence and ylds for id_adj - use separate get_draws calls so easier to get wide format
dw_adj_draws <- function(region_id, loc_ids, id, id_adj, gbd_round, como_id, ages, decomp_step, prev_dir, pop_df){
  if(id!=ID & id!=ID){
    temp_yld <- get_draws(gbd_id_type = "sequela_id",
                        gbd_id = id_adj,
                        source = "como",
                        gbd_round_id = gbd_round,
                        version_id = como_id,
                        location_id = loc_ids,
                        year_id = years,
                        age_group_id = c(ages,ID),
                        sex_id = c(ID),
                        measure_id = ID,
                        metric_id = ID,
                        decomp_step = decomp_step,
                        num_workers = 10)
  
    temp_prev <- get_draws(gbd_id_type = "sequela_id",
                         gbd_id = id_adj,
                         source = "como",
                         gbd_round_id = gbd_round,
                         version_id = como_id,
                         location_id = loc_ids,
                         year_id = years,
                         age_group_id = c(ages,ID),
                         sex_id = c(ID),
                         measure_id = ID,
                         metric_id = ID,
                         decomp_step=decomp_step,
                         num_workers = 10)
  }
  
  if(id==ID | id==ID){
    temp_yld = as.data.table(readRDS(paste0(FILEPATH)))
    temp_prev = as.data.table(readRDS(paste0(FILEPATH))) 
    temp_yld = temp_yld[metric_id==ID, ]
    temp_prev = temp_prev[metric_id==ID, ]
  }
  
  setnames(temp_prev, draw_vars, prev_vars)
  temp_prev[,c("metric_id", "version_id", "measure_id", "sequela_id", "population") := NULL]
  temp_yld[,c("metric_id", "version_id", "measure_id", "sequela_id", "population") := NULL]
  
  #merge and divide ylds by prevalence to get disability weight draws
  temp <- as.data.table(merge(temp_yld, temp_prev, by = c("location_id", "year_id", "sex_id", "age_group_id")))
  temp[, paste0("dw_", 0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) / get(paste0("prev_", x)))]
  temp[, c(prev_vars, draw_vars) := NULL]
  setnafill(temp, fill = 0)
  
  #get prevalence rates for original sequela ID (saved above) 
  prev <- as.data.table(readRDS(paste0(prev_dir, id, "_loc", region_id, "_prev_como", como_id, ".rds")))
  prev <- prev[metric_id==ID & measure_id==ID, ]
  prev[,c("population", "metric_id", "version_id", "measure_id") := NULL]
  setnames(prev, draw_vars, prev_vars)
  
  #combine and multiply prevalence draws by disability weight draws to get custom yld rate draws
  dt_combo <- as.data.table(merge(temp, prev, by = c("location_id", "year_id", "sex_id", "age_group_id")))
  dt_combo[, paste0("draw_",0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("dw_", x)) * get(paste0("prev_", x)))]
  dt_combo[, c(prev_vars, paste0("dw_", 0:nf_draws), "population") := NULL]
  dt_combo[, metric_id := ID]
  
  #merge on population and multiply to get counts
  dt_counts <- as.data.table(merge(dt_combo, pop_df[,c("location_id", "year_id", "sex_id", "age_group_id", "population")], 
                                  by = c("location_id", "year_id", "sex_id", "age_group_id")))
  dt_counts[, paste0("draw_",0:nf_draws) := lapply(0:nf_draws, function(x) get(paste0("draw_", x)) * population)]
  dt_counts[, c("population") := NULL]
  dt_counts[, metric_id := ID]
  
  #combine rates and counts and return
  dt <- rbind(dt_combo, dt_counts, fill=T)
  dt[, measure_id := ID]
  return(dt)
}    

# If custom YLD is needed because combined health state includes non-neurological states - determine correct ID and use function defined above
if(length(neuro_id_subset!=0)){
  if(neuro_id_subset$diff_dw==1){
  
    #get sequela_id with only neuro state and use to create custom ylds 
    id_adj <- neuro_id_subset$gbd_match_id
  
    dt <- dw_adj_draws(region_id=region_id, loc_ids=loc_ids, id=id, id_adj=id_adj, gbd_round=gbd_round, decomp_step=decomp_step, como_id=como_id, 
                   ages=ages, prev_dir=prev_dir, pop_df=pop_df)

  saveRDS(dt, paste0(FILEPATH))
  print(paste0("YLDs saved for ", id, " for region ", region_id))  
  }
}

# If custom YLD is needed for cerebral palsy where only motor component is included in YLD calculation - determine correct ID, pull disability weight draws, merge to prevalence rates and multiply
neuro_id_subset <- neuro_ids[short_name == "cp" & gbd_id==id, ]
if(length(neuro_id_subset!=0)){
  if(neuro_id_subset$diff_dw==1){
    #get sequela_id with only motor state and use to create custom ylds
    id_adj <- neuro_id_subset$gbd_match_id   
  
    dt <- dw_adj_draws(region_id=region_id, loc_ids=loc_ids, id=id, id_adj=id_adj, gbd_round=gbd_round, decomp_step=decomp_step, como_id=como_id, 
               ages=ages, prev_dir=prev_dir, pop_df=pop_df)

    saveRDS(dt, paste0(FILEPATH))
    print(paste0("YLDs saved for ", id, " for region ", region_id))  
  }
}
