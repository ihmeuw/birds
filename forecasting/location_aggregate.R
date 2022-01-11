##########################################################################
### Author: USERNAME
### Date: 05/11/2019
### Project: GBD Nonfatal Estimation
### Purpose: Dementia Forecasting - Create Location Aggregates
##########################################################################

# LOAD DATA ---------------------------------------------------------------

message("loading data")

files <- list.files(paste0("FILEPATH"), full.names = T)[!grepl("location", list.files(paste0("FILEPATH")))]
draw_dt <- rbindlist(lapply(files, function(x) readr::read_rds(x)))

# FIND LOCS WE WANT -------------------------------------------------------

message("getting location data")

reportloc_dt <- get_location_metadata(location_set_id = ID, gbd_round_id = ID)

reportloc_dt <- reportloc_dt[location_type_id %in% IDS]
reportloc_dt <- reportloc_dt[!(location_type_id == ID & level == ID)]
normalaggregate_dt <- reportloc_dt[location_id %in% setdiff(reportloc_dt[, unique(location_id)], draw_dt[, unique(location_id)]) & !location_type_id %in% IDS]

fullloc_dt <- get_location_metadata(location_set_id = ID, gbd_round_id = ID)

# GET SDI LOCS - BUT FOR OUR LOCATIONS... ---------------------------------

sdiloc_dt <- get_location_metadata(location_set_id = ID, gbd_round_id = ID)
sdilocs <- sdiloc_dt[level == 0, location_id]

sdiloc_dt[!location_id %in% sdilocs, pttp := as.numeric(stringr::str_extract(path_to_top_parent, "[0-9]*$"))]
fixsdi_locdt <- copy(fullloc_dt)
sdiloc_dt <- merge(sdiloc_dt, fixsdi_locdt[, .(hier = path_to_top_parent, pttp = location_id)], by = "pttp")

## TRY ONE LEVEL UP - MULTIPLE TIMES
sdiloc_dt[location_id %in% draw_dt[, unique(location_id)], match := as.numeric(location_id)]
sdiloc_dt[is.na(match), upone := gsub("\\,[0-9]*$", "", hier)][, upone := stringr::str_extract(upone, "[0-9]*$")]
sdiloc_dt[upone %in% draw_dt[, unique(location_id)], match := as.numeric(upone)]
sdiloc_dt[is.na(match), uptwo := gsub("\\,[0-9]*\\,[0-9]*$", "", hier)][, uptwo := stringr::str_extract(uptwo, "[0-9]*$")]
sdiloc_dt[uptwo %in% draw_dt[, unique(location_id)], match := as.numeric(uptwo)]
sdiloc_dt[is.na(match), upthree := gsub("\\,[0-9]*\\,[0-9]*$", "", hier)][, upthree := gsub("\\,[0-9]*$", "", upthree)]
sdiloc_dt[is.na(match), upthree := ifelse(grepl("\\,", upthree), as.numeric(stringr::str_extract(upthree, "[0-9]*$")), as.numeric(upthree))]
sdiloc_dt[upthree %in% draw_dt[, unique(location_id)], match := as.numeric(upthree)]
sdiloc_dt <- sdiloc_dt[!is.na(match)] 

sdiloc_dt <- unique(sdiloc_dt, by = "match")
sdiloc_dt <- sdiloc_dt[, .(location_id = match, path_to_top_parent = paste0(parent_id, ",", match))]

# GET POPULATION ----------------------------------------------------------

message("getting population data")

pop_dt <- readr::read_rds(paste0(forecast_dir, "FILEPATH"))
pop_dt <- pop_dt[age_group_id %in% IDS, c("year_id", "location_id", "age_group_id", "sex_id", draws), with = F]
pop_dt <- melt(pop_dt[year_id %in% forecast_years], id.vars = c("year_id", "age_group_id", "sex_id", "location_id"), measure.vars = draws,
                 variable.name = "draw", value.name = "pop")
pop_dt[, draw := as.numeric(gsub("^draw_", "", draw)) + 1]
setkeyv(pop_dt, c("year_id", "age_group_id", "sex_id", "location_id", "draw"))

# AGGREGATE TO SDI POPS ---------------------------------------------------

aggregate_pop <- function(agg_loc, loc_dt, data_dt){
  print(agg_loc)
  agg_locs <- loc_dt[grepl(paste0(",", agg_loc, ","), path_to_top_parent) | grepl(paste0(",", agg_loc, "$"), path_to_top_parent) 
                     | grepl(paste0("^", agg_loc, ","), path_to_top_parent), location_id]
  dt <- copy(data_dt[location_id %in% agg_locs, .(location_id, year_id, age_group_id, sex_id, draw, pop)])
  dt[, pop := sum(pop), by = c("year_id", "sex_id", "age_group_id", "draw")]
  dt <- unique(dt, by = c("year_id", "sex_id", "age_group_id", "draw"))
  dt[, location_id := agg_loc]
  return(dt)
}

sdi_pops <- rbindlist(lapply(sdilocs, 
                               function(x) aggregate_pop(agg_loc = x, loc_dt = sdiloc_dt, data_dt = pop_dt)))

pop_dt <- rbind(sdi_pops, pop_dt)

# AGGREGATE TO ALL AGE POPS -----------------------------------------------

allage_popdt <- copy(pop_dt)
allage_popdt[, pop := sum(pop), by = c("location_id", "year_id", "sex_id", "draw")]
allage_popdt <- unique(allage_popdt, by = c("location_id", "year_id", "sex_id", "draw"))
allage_popdt[, age_group_id := 22]

# AGGREGATE COUNTRIES/REGIONS ---------------------------------------------

aggregate_rate <- function(agg_loc, loc_dt, data_dt, population_dt = pop_dt){
  print(agg_loc)
  agg_locs <- loc_dt[grepl(paste0(",", agg_loc, ","), path_to_top_parent) | grepl(paste0(",", agg_loc, "$"), path_to_top_parent) 
                     | grepl(paste0("^", agg_loc, ","), path_to_top_parent), location_id]
  dt <- copy(data_dt[location_id %in% agg_locs & !age_group_id %in% IDS, .(location_id, year_id, age_group_id, sex_id, draw, prev)])
  setkeyv(dt, c("year_id", "age_group_id", "sex_id", "location_id", "draw"))
  dt <- merge(dt, population_dt)
  dt[, total_pop := sum(pop), by = c("year_id", "sex_id", "age_group_id", "draw")]
  dt[, prev := sum(prev * pop / total_pop), by = c("year_id", "sex_id", "age_group_id", "draw")]
  dt <- unique(dt, by = c("year_id", "sex_id", "age_group_id", "draw"))
  dt[,`:=` (pop = NULL, total_pop = NULL, location_id = agg_loc)]
  return(dt)
}

aggregated_rates <- rbindlist(lapply(normalaggregate_dt[, location_id], 
                                                 function(x) aggregate_rate(agg_loc = x, loc_dt = fullloc_dt, data_dt = draw_dt)))
sdi_rates <- rbindlist(lapply(sdilocs, 
                              function(x) aggregate_rate(agg_loc = x, loc_dt = sdiloc_dt, data_dt = draw_dt)))
aggregated_rates <- rbind(aggregated_rates, sdi_rates)

# CALC ALL AGE AND AGE STANDARDIZED RATES ---------------------------------

message("calculating all age and age standardized")

allage <- copy(aggregated_rates)
setkeyv(allage, c("year_id", "age_group_id", "sex_id", "location_id", "draw"))
aggpop_dt <- copy(pop_dt[location_id %in% allage[, unique(location_id)]])
allage <- merge(allage, aggpop_dt, all = T)
allage[is.na(prev), prev := 0]
allage[, `:=` (total_pop = sum(pop)), by = c("year_id", "sex_id", "draw", "location_id")]
allage[, prev := sum(prev * pop / total_pop), by = c("year_id", "sex_id", "draw", "location_id")]
allage <- unique(allage, by = c("year_id", "sex_id", "draw", "location_id"))
allage[, `:=` (age_group_id = ID, pop = NULL, total_pop = NULL)]

source(paste0(functions_dir, "get_age_metadata.R"))
age_weights <- get_age_metadata(ID, gbd_round_id = ID)

agestd <- copy(aggregated_rates)
agestd <- merge(agestd, age_weights[, .(age_group_id, weight = age_group_weight_value)], by = "age_group_id")
agestd[, prev := sum(prev * weight), by = c("year_id", "sex_id", "draw", "location_id")]
agestd <- unique(agestd, by = c("year_id", "sex_id", "draw", "location_id"))
agestd[, `:=` (age_group_id = ID, weight = NULL)]

aggregated_rates <- rbind(aggregated_rates, allage, agestd)

# CALCULATE NUMBERS FROM ALL-AGE RATES AND POPS ---------------------------

num_dt <- merge(aggregated_rates[!age_group_id == ID], rbind(pop_dt, allage_popdt), by = c("location_id", "year_id", "sex_id", "draw", "age_group_id"))
num_dt[, num := prev * pop]
num_dt[, c("prev", "pop") := NULL]
aggregated_counts <- copy(num_dt)

# COMBINE -----------------------------------------------------------------

message("combining and calculating means")

setkeyv(aggregated_rates, c("year_id", "age_group_id", "sex_id", "location_id", "draw"))
setkeyv(aggregated_counts, c("year_id", "age_group_id", "sex_id", "location_id", "draw"))
region_dt <- merge(aggregated_rates, aggregated_counts, all = T)

# CALCULATE MEANS ---------------------------------------------------------

mean_dt <- copy(region_dt)
mean_dt[, `:=` (mean_prev = mean(prev), lower_prev = quantile(prev, probs= 0.025), upper_prev = quantile(prev, probs = 0.975)),
        by = c("location_id", "year_id", "age_group_id", "sex_id")]
mean_dt[!age_group_id == ID, `:=` (mean_num = mean(num), lower_num = quantile(num, probs= 0.025), upper_num = quantile(num, probs = 0.975)),
        by = c("location_id", "year_id", "age_group_id", "sex_id")]
mean_dt <- unique(mean_dt, by = c("location_id", "year_id", "age_group_id", "sex_id"))
mean_dt[, c("prev", "num", "draw") := NULL]

# SAVE RESULTS ------------------------------------------------------------

message("saving results")

readr::write_rds(mean_dt, paste0("FILEPATH"))
readr::write_rds(region_dt, paste0("FILEPATH"))
