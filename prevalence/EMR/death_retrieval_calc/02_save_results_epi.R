##########################################################################
### Author: USERNAME
### Date: 07/22/2019
### Project: GBD Nonfatal Estimation
### Purpose: Save Custom COD Results to MEID
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
date <- "2019_04_26"

# SET OBJECTS -------------------------------------------------------------

functions_dir <- paste0("FILEPATH")
save_dir <- paste0("FILEPATH", date, "/rate/")
meid <- ID

# SOURCE FUNCTIONS --------------------------------------------------------

functs <- c("save_results_epi.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

# RUN SAVE RESULTS --------------------------------------------------------

save_results_epi(input_dir = save_dir, input_file_pattern = "{location_id}.csv", 
                 modelable_entity_id = meid, description = "COD Results Estimation Years",
                 measure_id = 6, decomp_step = "step3", mark_best = T)

