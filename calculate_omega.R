##########################################################################
### Author: USERNAME
### Date: 04/01/2020
### Project: GBD Nonfatal Estimation
### Purpose: IRT Paper - get Omega
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

pacman::p_load(data.table, openxlsx, ggplot2, readr, gridExtra, grid)
library("mirt", lib.loc = paste0("FILEPATH"))
date <- gsub("-", "_", Sys.Date())

# SET OBJECTS -------------------------------------------------------------

harmonization_dir <- paste0("FILEPATH")
model_dir <- paste0("FILEPATH")
figure_dir <- paste0("FILEPATH")
adams_dir <- paste0("FILEPATH")

model_name <- "model_2factorsSE_hrs_adams_2020_08_26.rds"

# LOAD MODEL --------------------------------------------------------------

model <- read_rds(paste0(model_dir, model_name))

# LOAD DATA ---------------------------------------------------------------

adams_dt <- readr::read_rds("FILEPATH")
hrs_dt <- readr::read_rds("FILEPATH")
categories <- as.data.table(read.xlsx("FILEPATH"))

communalities <- extract.mirt(model, "h2")

F_matrix <- extract.mirt(model, "F")
standardized_loadings <- ifelse(F_matrix[,1] == 0, F_matrix[,2],F_matrix[,1])

get_dt <- function(item){
  data <- c()
  if (item %in% names(adams_dt)) data <- adams_dt[, get(item)]
  if (length(data) == 0) data <- hrs_dt[, get(item)]
  if (item %in% names(hrs_dt) & !length(data) == 0) data <- c(data, hrs_dt[, get(item)])
  row <- data.table(item_name = item, item_variance = var(data, na.rm = T), 
                    communality = communalities[(item)], loading = standardized_loadings[(item)])
  return(row)
}

omega_dt <- rbindlist(lapply(names(communalities), get_dt))
omega_dt[, residual_variance := 1-communality]
omega_dt <- merge(omega_dt, categories, by.x = "item_name", by.y = "item")
omega_dt[, .(omega = sum(loading^2)/(sum(loading^2+residual_variance^2))), by = "factor"]
