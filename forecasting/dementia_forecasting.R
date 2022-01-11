##########################################################################
### Author: USERNAME
### Date: 05/11/2019
### Project: GBD Nonfatal Estimation
### Purpose: Dementia Forecasting
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
plot_dir <- paste0(j_root, "FILEPATH")
code_dir <- paste0(h_root, "FILEPATH")
meid <- ID
cid <- ID
codcorrect_version <- ID
functions_dir <- paste0(functions_dir, "FILEPATH")
draws <- paste0("draw_", 0:999)
forecast_years <- c(2019, 2020, 2030, 2040, 2050)

# SOURCE FUNCTIONS --------------------------------------------------------

functs <- c("get_model_results.R", "get_location_metadata.R", "get_age_metadata.R", "get_draws.R",
            "get_population.R", "get_ids.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

# GET INFO ----------------------------------------------------------------

age_dt <- get_age_metadata(12, gbd_round_id = 6)
age_dt <- age_dt[age_group_id >=13]

loc_dt <- get_location_metadata(location_set_id = 1, gbd_round_id = 6)

# GET COVARIATES ----------------------------------------------------------

## GBD Risks = FPG, BMI, smoking

## Include covs where male/female estimates are both in the same (correct) direction
covs <- c("education")

get_cov <- function(cov){
  cov_dt <- readr::read_rds(paste0(forecast_dir, cov, "FILEPATH")) 
  setnames(cov_dt, draws, paste0(cov, "_", 0:999))
  setkeyv(cov_dt, cols = c("location_id", "age_group_id", "sex_id", "year_id"))
  return(cov_dt)
}

allcov_dt <- Reduce(function(d1, d2) merge(d1, d2),
                    lapply(covs, get_cov))

# GET SCALARS -------------------------------------------------------------

scalar <- as.data.table(readr::read_rds(paste0(forecast_dir, "FILEPATH"))) 
setnames(scalar, draws, paste0("scalar_", 0:999))

# GET PREVALENCE ----------------------------------------------------------

prev <- as.data.table(readr::read_rds(paste0(forecast_dir, "FILEPATH")))

# GET REGRESSION DATA -----------------------------------------------------

model_dt <- merge(prev, allcov_dt, by = c("location_id", "age_group_id", "sex_id", "year_id"))
model_dt <- merge(model_dt, scalar, by = c("location_id", "age_group_id", "sex_id", "year_id"))
model_dt[, (draws) := lapply(0:999, function(x) get(paste0("draw_", x)) / get(paste0("scalar_", x)))]
model_dt[, paste0("scalar_", 0:999) := NULL]
model_dt[, (draws) := lapply(.SD, gtools::logit), .SDcols = draws]
model_dt <- merge(model_dt, loc_dt[, .(location_id, region_name, super_region_name, location_name)], by = "location_id")

# SAVE DATA FOR JOBS AND FORMAT -------------------------------------------

dir.create(paste0(forecast_dir, date, "/regression_sims"), recursive = T)
dir.create(paste0(forecast_dir, date, "/residuals"))
save_dir <- paste0(forecast_dir, date, "/")

readr::write_rds(model_dt, paste0(save_dir, "model_data.rds"))

params <- data.table(draw = 0:999)
params[, task_num := 1:.N]
map_path <- paste0(save_dir, "task_map.csv")
write.csv(params, map_path, row.names = F)

# RUN JOBS ----------------------------------------------------------------

array_qsub(jobname = "forecasts",
           shell = paste0(forecast_dir, "FILEPATH"),
           code = paste0(code_dir, "FILEPATH"),
           pass = list(map_path, save_dir),
           proj = "proj_yld",
           num_tasks = nrow(params), archive_node = F,
           cores = 2, mem = 30, log = T, submit = F)

# CHECK FILES -------------------------------------------------------------

all_files <- paste0(apply(expand.grid(paste0("draw_", 0:999), c(1,2)), 1, paste, collapse="_"), ".rds")
assertable::check_files(all_files, paste0(save_dir, "FILEPATH"))

# READ FILES IN BY YEAR AND COLLAPSE --------------------------------------

newscalar <- copy(scalar[year_id >= 2019])
setnames(newscalar, paste0("scalar_", 0:999), draws)

readr::write_rds(allcov_dt, paste0(save_dir, "FILEPATH"))
readr::write_rds(newscalar, paste0(save_dir, "FILEPATH"))
readr::write_rds(prev, paste0(save_dir, "FILEPATH"))

dir.create(paste0(forecast_dir, date, "FILEPATH"))
combined_dir <- paste0(forecast_dir, date, "FILEPATH")

comb_params <- data.table(location = model_dt[, unique(location_id)])
comb_params[, task_num := 1:.N]
comb_mappath <- paste0(combined_dir, "task_map.csv")
write.csv(comb_params, comb_mappath, row.names = F)

array_qsub(jobname = "forecast_combination",
           shell = paste0(forecast_dir, "FILEPATH"),
           code = paste0(code_dir, "FILEPATH"),
           pass = list(comb_mappath, combined_dir, save_dir),
           proj = "proj_yld",
           num_tasks = nrow(comb_params), archive_node = F,
           cores = 10, mem = 40, log = T, submit = F)

# CHECK FILES -------------------------------------------------------------

pop_dt <- readr::read_rds(paste0(forecast_dir, "FILEPATH"))
locs <- pop_dt[, unique(location_id)]
locs <- locs[locs %in% model_dt[, location_id]]
assertable::check_files(paste0("FILEPATH"), paste0(combined_dir, "FILEPATH"))

# CREATE REGION AGGREGATES ------------------------------------------------

source(paste0(code_dir, "FILEPATH"))
