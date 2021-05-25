##########################################################################
### Author: USERNAME
### Date: 07/22/2019
### Project: GBD Nonfatal Estimation
### Purpose: Save Custom COD Results to COD
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


pacman::p_load(data.table, openxlsx, ggplot2)
date <- gsub("-", "_", Sys.Date())
date <- "2019_07_24"

# SET OBJECTS -------------------------------------------------------------

functions_dir <- paste0("FILEPATH")
save_dir <- paste0("FILEPATH", date, "/interpolated/")
cid <- ID

# SOURCE FUNCTIONS --------------------------------------------------------

functs <- c("save_results_cod.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

# RUN SAVE RESULTS --------------------------------------------------------

save_results_cod(input_dir = save_dir, input_file_pattern = "{location_id}.csv", 
                 cause_id = cid, description = "Custom Results Dementia Mortality",
                 metric_id = 3, decomp_step = "step4", mark_best = T, gbd_round_id = 6)
