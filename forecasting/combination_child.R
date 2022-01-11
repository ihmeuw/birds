##########################################################################
### Author: USERNAME
### Date: 05/11/2019
### Project: GBD Nonfatal Estimation
### Purpose: Dementia Forecasting - Child Combination
##########################################################################

rm(list=ls())

Sys.setenv(MKL_VERBOSE = 0)

if (Sys.info()["sysname"] == "Linux") {
  j_root <- "FILEPATH" 
  h_root <- "FILEPATH"
  l_root <- "FILEPATH"
  functions_dir <- "FILEPATH"
} else { 
  j_root <- "FILEPATH"
  h_root <- "FILEPATH"
  l_root <- "FILEPATH"
  functions_dir <- "FILEPATH"
}

pacman::p_load(data.table, openxlsx, ggplot2, parallel, magrittr, lme4, boot, arm, forecast)
library(mortcore, lib = "FILEPATH")
forecast_dir <- "FILEPATH"
functions_dir <- paste0(functions_dir, "FILEPATH")
draws <- paste0("draw_", 0:999)
forecast_years <- c(2019, 2020, 2030, 2040, 2050)

date <- gsub("-", "_", Sys.Date())

## GET ARGS AND ITEMS
args<-commandArgs(trailingOnly = TRUE)
map_path <-args[1]
save_dir <-args[2]
date_dir <-args[3]
sim_dir <- paste0(date_dir, "regression_sims/")
resid_dir <- paste0(date_dir, "residuals/")
params <- fread(map_path)
task_id <- Sys.getenv("SGE_TASK_ID")
loc_id <- params[, location][as.numeric(task_id)]
print(loc_id)
print(sim_dir)

# USER FUNCTION -----------------------------------------------------------

dummy_cols <- function(.data,
                       select_columns = NULL,
                       remove_first_dummy = TRUE) {
  
  # Grabs column names that are character or factor class -------------------
  if (!is.null(select_columns)) {
    char_cols <- select_columns
    cols_not_in_data <- char_cols[!char_cols %in% names(.data)]
    char_cols <- char_cols[!char_cols %in% cols_not_in_data]
    if (length(char_cols) == 0) {
      stop("select_columns is/are not in data. Please check data and spelling.")
    }
  } else if (ncol(.data) == 1) {
    char_cols <- names(.data)
  } else {
    char_cols <- sapply(.data, class)
    char_cols <- char_cols[char_cols %in% c("factor", "character")]
    char_cols <- names(char_cols)
  }
  
  for (col_name in char_cols) {
    # If factor type, order by assigned levels
    if (is.factor(.data[[col_name]])) {
      unique_vals <- levels(.data[[col_name]])
      if (any(is.na(.data[[col_name]]))) {
        unique_vals <- c(unique_vals, NA)
      }
      # Else by alphabetical order.
    } else {
      unique_vals <- unique(.data[[col_name]])
      unique_vals <- stringi::stri_sort(unique_vals, na_last = TRUE, locale = "en_US")
    }
    unique_vals <- as.character(unique_vals)
    
    if (remove_first_dummy) {
      unique_vals <- unique_vals[-1]
    }
    
    data.table::alloc.col(.data, ncol(.data) + length(unique_vals))
    data.table::set(.data, j = paste0(col_name, "_", unique_vals), value = 0L)
    for (unique_value in unique_vals) {
      data.table::set(.data, i =
                        which(data.table::chmatch(
                          as.character(.data[[col_name]]),
                          unique_value, nomatch = 0) == 1L),
                      j = paste0(col_name, "_", unique_value), value = 1L)
      
      
      # Sets NA values to NA, only for columns that are not the NA columns
      if (!is.na(unique_value)) {
        data.table::set(.data, i =
                          which(is.na(.data[[col_name]])),
                        j = paste0(col_name, "_", unique_value), value = NA)
      }
      if (is.na(unique_value)) {
        .data[[paste0(col_name, "_", unique_value)]][which(!is.na(.data[[col_name]]))] <- 0
      }
    }
  }
  return(.data)
}

# SLEEP FOR RANDOM TIME  -------------------------------------------------

set.seed(loc_id)
rtime <- runif(1, min = 0, max = 60)
Sys.sleep(rtime)

# LOAD ALL DATA -----------------------------------------------------------

model_dt <- readr::read_rds(paste0(date_dir, "FILEPATH"))
cov_dt <- readr::read_rds(paste0(date_dir, "FILEPATH"))
scalar_dt <- readr::read_rds(paste0(date_dir, "FILEPATH"))
oldprev_dt <- readr::read_rds(paste0(date_dir, "FILEPATH"))

# READ IN THE MODEL BETA SIMS ---------------------------------------------

file_dt <- data.table(file = list.files(sim_dir, full.names = T)[grepl("^draw", list.files(sim_dir))])
file_dt[, draw_num := gsub("_[0-9]\\.rds$", "", file)]
file_dt[, draw_num := as.numeric(stringr::str_extract(draw_num, pattern = "[0-9]*$"))]
file_dt[, sex_id := ifelse(grepl("2.rds$", file), 2, 1)]
file_dt <- file_dt[order(draw_num)]

get_matrices <- function(row, data){
  file_name <- data[, file][row]
  draw_num <- data[, draw_num][row]
  sim <- readr::read_rds(file_name)
  mat <- cbind(sim@coef, draw = rep(draw_num, nrow(sim@coef)), sigma = sim@sigma)
  return(mat)
}

female_sims <- do.call(rbind, lapply(1:1000, function(x) get_matrices(row = x, data = file_dt[sex_id == 2])))
male_sims <- do.call(rbind, lapply(1:1000, function(x) get_matrices(row = x, data = file_dt[sex_id == 1])))

# SAMPLE MATRICES DOWN TO 1000 ROWS ---------------------------------------

small_fsims <- female_sims[sample(nrow(female_sims), 1000, replace = T),]
small_msims <- male_sims[sample(nrow(male_sims), 1000, replace = T),]

# GET PREDICTIONS ---------------------------------------------------------

get_preds <- function(row, sims, sid){
  print(row)
  sim_row <- sims[row, ]
  
  ## GET PREDICTION MATRIX
  pred_mat <- unique(model_dt[sex_id == sid, .(age_group_id, region_name, location_id)])
  pred_mat[, region_name := factor(region_name, levels = pred_mat[, sort(unique(region_name))])]
  cov_mat <- copy(cov_dt[year_id %in% forecast_years & sex_id == sid, c("age_group_id", "location_id", 
                  paste0("education_", sim_row[["draw"]]),  
                  "year_id"), with = F])
  cov_mat <- cov_mat[!location_id %in% setdiff(cov_mat[, location_id], pred_mat[, location_id])]
  mat <- merge(pred_mat, cov_mat, by = c("age_group_id", "location_id"))
  mat <- mat[order(region_name, age_group_id)]
  mat[, age_group_id := factor(age_group_id, levels = unique(age_group_id))]
  mat <- dummy_cols(mat, select_columns = c("age_group_id", "region_name"), remove_first_dummy = T)
  
  region_cols <- names(mat)[grepl("^region_name", names(mat)) & !names(mat) == "region_name"]
  interaction_mat <- copy(mat[, c("location_id", region_cols, "year_id"), with = F])
  interaction_mat[, (region_cols) := lapply(.SD, function(x) x * year_id), .SDcols = region_cols]
  interaction_mat[, year_id := NULL]
  
  final_mat <- cbind(rep(1, nrow(mat[location_id == loc_id])), mat[location_id == loc_id, c(names(mat)[grepl("^age_group", names(mat)) & !names(mat) == "age_group_id"], 
                     paste0("education_", sim_row[["draw"]]), region_cols), with = F])
  
  ## CALCULATE PREDICTION
  beta_row <- sim_row[-match("draw", names(sim_row)):-match("sigma", names(sim_row))]
  draw <- rnorm(nrow(final_mat), as.matrix(final_mat) %*% beta_row, 0) 
  draw_row <- data.table(location_id = loc_id, sex_id = sid, year_id = mat[location_id == loc_id, year_id], age_group_id = mat[location_id == loc_id, age_group_id],
                         value = draw, draw = sim_row['draw'], draw_id = row)
  return(draw_row)
}

fpreds <- rbindlist(lapply(1:nrow(small_fsims), function(x) get_preds(row = x, sims = small_fsims, sid = 2)))
mpreds <- rbindlist(lapply(1:nrow(small_msims), function(x) get_preds(row = x, sims = small_msims, sid = 1)))

allpreds <- rbind(fpreds, mpreds)
allpreds[, age_group_id := as.numeric(as.character(age_group_id))]

# GET RESIDUALS -----------------------------------------------------------

rfile_dt <- data.table(file = list.files(resid_dir, full.names = T)[grepl("^draw", list.files(resid_dir))])
rfile_dt[, draw_num := gsub("\\.rds$", "", file)]
rfile_dt[, draw_num := as.numeric(stringr::str_extract(draw_num, pattern = "[0-9]*$"))]
rfile_dt <- rfile_dt[order(draw_num)]

get_resids <- function(row, data){
  file_name <- data[, file][row]
  draw_num <- data[, draw_num][row]
  resid_dt <- readr::read_rds(file_name)
  resid_dt[, draw := draw_num]
  resid_dt <- resid_dt[location_id == loc_id]
  return(resid_dt)
}

all_resids <- rbindlist(lapply(1:nrow(rfile_dt), function(x) get_resids(row = x, data = rfile_dt)))

# GET RANDOM WALK  --------------------------------------------------------

walk <- function(id_var, data){
  walk_dt <- copy(data[id == id_var & !year_id == 2019])
  rw <- Arima(y = as.matrix(walk_dt[, resid]), order = c(0,1,0)) 
  rw_dt <- data.table(year_id = 2019:2050, walk_value = as.numeric(simulate(rw, nsim = 32)), id = rep(id_var, 32))
  rw_dt <- rw_dt[year_id %in% forecast_years]
  return(rw_dt)
}

get_walks <- function(n){
  print(n)
  key <- allpreds_map[n]
  baby_resid <- copy(all_resids[draw == key[, draw] & sex_id == key[, sex_id]])
  baby_resid[, id := .GRP, by = c("age_group_id")]
  walks <- rbindlist(lapply(1:baby_resid[, max(id)], function(x) walk(x, data = baby_resid)))
  walks <- merge(walks, unique(baby_resid[, .(id, age_group_id, sex_id)]), by = "id")
  walks[, `:=` (draw = key[, draw], draw_id = key[, draw_id], id = NULL)]
  return(walks)
}

allpreds_map <- unique(allpreds[, .(draw, draw_id, sex_id)])
allwalk_dt <- rbindlist(mclapply(1:nrow(allpreds_map), function(x) get_walks(x), mc.cores = 9))

# ADD RANDOM WALK ON TO RESULTS -------------------------------------------

pred_dt <- merge(allpreds, allwalk_dt, by = c("age_group_id", "year_id", "sex_id", "draw", "draw_id"))
pred_dt[, value := value + walk_value]
pred_dt[, c("walk_value", "draw_id", "draw") := NULL] 
setorderv(pred_dt, c("age_group_id", "year_id", "sex_id"))
pred_dt[, draw := rep(draws, nrow(pred_dt)/1000)]
pred_dt[, value := inv.logit(value)]

# APPLY FUTURE SCALARS ----------------------------------------------------

scalar_dt <- melt(scalar_dt, id.vars = c("location_id", "year_id", "age_group_id", "sex_id"), measure.vars = draws, variable.name = "draw", value.name = "scalar")
pred_dt <- merge(pred_dt, scalar_dt, by = c("age_group_id", "year_id", "sex_id", "location_id", "draw"))
pred_dt[, value := value * scalar]
pred_dt[, scalar := NULL]

# CALCULATE SHIFT ---------------------------------------------------------

## SET UP DATA TO COMBINE
newshift_dt <- copy(pred_dt[year_id == 2019])
newshift_dt[, age_group_id := as.numeric(as.character(age_group_id))]
oldshift_dt <- copy(oldprev_dt[year_id == 2019, c("location_id", "year_id", "age_group_id", "sex_id", draws), with = F])
oldshift_dt <- melt(oldshift_dt, id.vars = c("location_id", "year_id", "age_group_id", "sex_id"), measure.vars = draws, variable.name = "draw")

setnames(oldshift_dt, "value", "old")

## COMBINE AND CALCULATE
shift_dt <- merge(oldshift_dt, newshift_dt, by = c("location_id", "year_id", "age_group_id", "sex_id", "draw"), sort = F)
shift_dt[, shift := value/old]
shift_dt[, c("old", "value", "year_id") := NULL]

# COMBINE AND APPLY SHIFT -------------------------------------------------

postshift_dt <- merge(pred_dt, shift_dt, by = c("sex_id", "age_group_id", "location_id", "draw"))
postshift_dt[, value := value/shift]
postshift_dt[, `:=` (shift = NULL)]
results_dt <- copy(postshift_dt)
setnames(results_dt, "value", "prev")
results_dt[, draw := as.numeric(stringi::stri_extract(draw, regex = "[0-9]*$")) + 1]

pop_dt <- as.data.table(readr::read_rds(paste0(forecast_dir, "FILEPATH"))) 
pop_dt <- pop_dt[, c("year_id", "location_id", "age_group_id", "sex_id", draws), with = F]
pop_dt_long <- melt(pop_dt, id.vars = c("year_id", "location_id", "age_group_id", "sex_id"), measure.vars = draws)
pop_dt_long[, draw := as.numeric(gsub("draw_", "", variable))+1]
pop_dt_long[, variable := NULL]

setkeyv(pop_dt_long, c("year_id", "location_id", "age_group_id", "sex_id", "draw"))
setkeyv(results_dt, c("year_id", "location_id", "age_group_id", "sex_id", "draw"))
num_dt <- merge(results_dt, pop_dt_long)
num_dt[, num := prev * value]
num_dt[, c("prev", "value") := NULL]

draw_dt <- merge(results_dt, num_dt, by = c("location_id", "year_id", "age_group_id", "sex_id", "draw"))

# CALCULATE ALL-AGE AND AGE-STANDARDIZED ----------------------------------

pop_long <- readr::read_rds(paste0(forecast_dir, "FILEPATH")) 
pop_long <- pop_long[age_group_id %in% IDS, c("year_id", "location_id", "age_group_id", "sex_id", draws), with = F]
pop_long <- melt(pop_long[location_id == loc_id & year_id %in% forecast_years], id.vars = c("year_id", "age_group_id", "sex_id"), measure.vars = draws,
                 variable.name = "draw", value.name = "pop")
pop_long[, draw := as.numeric(gsub("^draw_", "", draw)) + 1]

allage <- copy(draw_dt)
allage <- merge(allage, pop_long[, .(year_id, draw, age_group_id, sex_id, pop)], by = c("year_id", "draw", "age_group_id", "sex_id"), all = T)
allage[is.na(location_id), `:=` (prev = 0, num = 0)]
allage[, `:=` (total_pop = sum(pop)), by = c("year_id", "sex_id", "draw")]
allage[, prev := sum(prev * pop / total_pop), by = c("year_id", "sex_id", "draw")]
allage[, num := prev * total_pop]
allage <- unique(allage, by = c("year_id", "sex_id", "draw"))
allage[, `:=` (age_group_id = ID, pop = NULL, total_pop = NULL, location_id = loc_id)]

source(paste0(functions_dir, "get_age_metadata.R"))
age_weights <- get_age_metadata(ID, gbd_round_id = ID)

agestd <- copy(draw_dt)
agestd <- merge(agestd, age_weights[, .(age_group_id, weight = age_group_weight_value)], by = "age_group_id")
agestd[, prev := sum(prev * weight), by = c("year_id", "sex_id", "draw")]
agestd <- unique(agestd, by = c("year_id", "sex_id", "draw"))
agestd[, `:=` (age_group_id = ID, weight = NULL, num = NA)]

draw_dt <- rbind(draw_dt, allage, agestd)

# CALC MEANS --------------------------------------------------------------

mean_dt <- copy(draw_dt)
mean_dt[, `:=` (mean_prev = mean(prev), lower_prev = quantile(prev, probs= 0.025), upper_prev = quantile(prev, probs = 0.975)),
        by = c("location_id", "year_id", "age_group_id", "sex_id")]
mean_dt[!age_group_id == 27, `:=` (mean_num = mean(num), lower_num = quantile(num, probs= 0.025), upper_num = quantile(num, probs = 0.975)),
        by = c("location_id", "year_id", "age_group_id", "sex_id")]
mean_dt <- unique(mean_dt, by = c("location_id", "year_id", "age_group_id", "sex_id"))

# SAVE RESULTS ------------------------------------------------------------

dir.create(paste0(save_dir, "means"))
dir.create(paste0(save_dir, "draws"))

readr::write_rds(mean_dt, paste0(save_dir, "FILEPATH"))
readr::write_rds(draw_dt, paste0(save_dir, "FILEPATH"))
