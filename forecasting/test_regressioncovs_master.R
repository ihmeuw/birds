##########################################################################
### Author: USERNAME
### Date: 05/11/2019
### Project: GBD Nonfatal Estimation
### Purpose: Dementia Forecasting Test Covariates
##########################################################################

rm(list=ls())

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

pacman::p_load(data.table, openxlsx, ggplot2, parallel, magrittr, lme4, ncdf4)
library(mortcore, lib = "FILEPATH")
date <- gsub("-", "_", Sys.Date())

# SET OBJECTS -------------------------------------------------------------

model_version_fp <- "FILEPATH"
forecast_dir <- "FILEPATH"
plot_dir <- paste0(j_root, "PFILEPATH")
code_dir <- paste0(h_root, "FILEPATH")
meid <- ID
functions_dir <- paste0(functions_dir, "FILEPATH")
draws <- paste0("draw_", 0:999)
forecast_years <- c(2019, 2020, 2030, 2040, 2050)
covs <- c("activity", "education", "drugs_alcohol", "metab_bmi", "metab_ldl", "metab_sbp", "smoking_direct", "air_pm", "metab_fpg")

# SOURCE FUNCTIONS --------------------------------------------------------

functs <- c("get_model_results.R", "get_location_metadata.R", "get_age_metadata.R", "get_draws.R",
            "get_population.R", "get_ids.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

# GET INFO ----------------------------------------------------------------

age_dt <- get_age_metadata(ID, gbd_round_id = ID)
age_dt <- age_dt[age_group_id >=ID]

loc_dt <- get_location_metadata(location_set_id = ID, gbd_round_id = ID)

for (cov in covs){
  ## GET COVARIATE
  cov_dt <- as.data.table(readr::read_rds(paste0(forecast_dir, cov, "FILEPATH")))
  setnames(cov_dt, draws, paste0(cov, "_", 0:999))
  
  ## GET EDUCATION
  education_dt <- as.data.table(readr::read_rds(paste0(forecast_dir, "FILEPATH")))
  setnames(education_dt, draws, paste0("education_", 0:999))
  
  ## GET PREVALENCE
  prev <- as.data.table(readr::read_rds(paste0(forecast_dir, "FILEPATH")))
  
  ## GET SCALARS
  scalar <- as.data.table(readr::read_rds(paste0(forecast_dir, "FILEPATH")))
  setnames(scalar, draws, paste0("scalar_", 0:999))
  
  ## APPLY SCALARS
  model_dt <- merge(prev, cov_dt, by = c("location_id", "age_group_id", "sex_id", "year_id"))
  model_dt <- merge(model_dt, education_dt, by = c("location_id", "age_group_id", "sex_id", "year_id"))
  model_dt <- merge(model_dt, scalar, by = c("location_id", "age_group_id", "sex_id", "year_id"))
  model_dt[, (draws) := lapply(0:999, function(x) get(paste0("draw_", x)) / get(paste0("scalar_", x)))]
  model_dt[, paste0("scalar_", 0:999) := NULL]
  
  ## TRANSFORM AND FORMAT
  model_dt[, (draws) := lapply(.SD, gtools::logit), .SDcols = draws]
  model_dt <- merge(model_dt, loc_dt[, .(location_id, region_name, super_region_name, location_name)], by = "location_id")
  
  ## SAVE DATA FOR REGRESSION JOBS
  dir.create(paste0("FILEPATH"), recursive = T)
  dir.create(paste0("FILEPATH"))
  save_dir <- paste0(forecast_dir, date, "_", cov, "/")
  readr::write_rds(model_dt, paste0(save_dir, "FILEPATH"))
  params <- data.table(draw = 0:999)
  params[, task_num := 1:.N]
  map_path <- paste0(save_dir, "FILEPATH")
  write.csv(params, map_path, row.names = F)
  
  ## RUN JOBS
  array_qsub(jobname = paste0("forecasts_", cov),
             shell = paste0(forecast_dir, "FILEPATH"),
             code = paste0(code_dir, "FILEPATH"),
             pass = list(map_path, save_dir, cov),
             proj = "proj_yld", 
             num_tasks = nrow(params), archive_node = F,
             cores = 2, mem = 10, log = T, submit = F)
  
  ## CHECK FILES
  all_files <- paste0(apply(expand.grid(paste0("draw_", 0:999), c(1,2)), 1, paste, collapse="_"), ".rds")
  assertable::check_files(all_files, paste0(save_dir, "FILEPATH"))
}


