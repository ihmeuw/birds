##########################################################################
### Author: USERNAME
### Date: 05/11/2019
### Project: GBD Nonfatal Estimation
### Purpose: Dementia Forecasting - Child
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
date <- gsub("-", "_", Sys.Date())

# SET OBJECTS -------------------------------------------------------------

model_version_fp <- "FILEPATH"
forecast_dir <- "FILEPATH"
meid <- ID
functions_dir <- paste0(functions_dir, "FILEPATH")
draws <- paste0("draw_", 0:999)
forecast_years <- c(2019, 2020, 2030, 2040, 2050)

## GET ARGS AND ITEMS
args<-commandArgs(trailingOnly = TRUE)
map_path <-args[1]
save_dir <-args[2]
params <- fread(map_path)
task_id <- Sys.getenv("SGE_TASK_ID")
draw <- params[, draw][as.numeric(task_id)] 

sim_dir <- paste0(save_dir, "regression_sims/")
resid_dir <- paste0(save_dir, "residuals/")

# SLEEP FOR RANDOM TIME TO HELP FILE SYSTEM -------------------------------

print(task_id)
print(draw)
set.seed(draw)
rtime <- runif(1, min = 0, max = 60)
Sys.sleep(rtime)

# LOAD ALL DATA -----------------------------------------------------------

model_dt <- readr::read_rds(paste0(save_dir, "FILEPATH"))

# RUN MODELS --------------------------------------------------------------

male <- lm(get(paste0("draw_", draw)) ~ as.factor(age_group_id) + get(paste0("education_", draw)) + 
           region_name, 
           data = model_dt[sex_id == 1])

female <- lm(get(paste0("draw_", draw)) ~ as.factor(age_group_id) + get(paste0("education_", draw)) + 
             region_name, 
             data = model_dt[sex_id == 2])

# GET BETA SIMS -----------------------------------------------------------

male_draws <- sim(male, 1000)
female_draws <- sim(female, 1000)

readr::write_rds(male_draws, paste0("FILEPATH"))
readr::write_rds(female_draws, paste0("FILEPATH"))

# GET RESIDUALS -----------------------------------------------------------

male_resid <- cbind(model_dt[sex_id == 1, .(age_group_id, year_id, super_region_name, location_id)], resid = residuals(male))
male_resid[, sex_id := 1]
female_resid <- cbind(model_dt[sex_id == 2, .(age_group_id, year_id, super_region_name, location_id)], resid = residuals(female))
female_resid[, sex_id := 2]
resid_dt <- rbind(male_resid, female_resid)

readr::write_rds(resid_dt, paste0("FILEPATH"))



