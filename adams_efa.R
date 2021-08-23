##########################################################################
### Author: USERNAME
### Date: 02/03/2020
### Project: GBD Nonfatal Estimation
### Purpose: EFA models for ADAMS
##########################################################################

rm(list=ls())

if (Sys.info()["sysname"] == "Linux") {
  j_root <- "FILEPATH" 
  h_root <- "FILEPATH"
  l_root <- "FILEPATH"
} else { 
  j_root <- "FILEPATH"
  h_root <- "FILEPATH"
  l_root <- "FILEPATH"
}

pacman::p_load(data.table, openxlsx, ggplot2, ggrepel, RColorBrewer, MASS, lmtest, nnet, gridExtra, grid, lattice)
package_dir <- paste0(j_root, "FILEPATH")
library("mirt", lib.loc = paste0(j_root, "FILEPATH"))
date <- gsub("-", "_", Sys.Date())

# SET UP OBJECTS ----------------------------------------------------------

repo_dir <- paste0(h_root, "FILEPATH")
graph_dir <- paste0(j_root, "FILEPATH")
model_dir <- paste0(j_root, "FILEPATH")
appendix_dir <- paste0(j_root, "FILEPATH")
draws <- paste0("draw_", 0:999)

# GET DATA ----------------------------------------------------------------

source(paste0(repo_dir, "FILEPATH"))

# SET UP DATA FOR MODEL ---------------------------------------------------

set_up_model <- function(test_dt){
  ## SET UP MODEL
  id_vars <- c("hh_id", "pn", "unique_id", "pweight")
  include_vars <- names(test_dt)[!names(test_dt) %in% id_vars]
  test_dt <- test_dt[!rowSums(is.na(dplyr::select(test_dt, include_vars))) == ncol(dplyr::select(test_dt, include_vars))] ## MAKE SURE NO ROWS ARE ALL NA
  test_dt <- test_dt[!pweight == 0]
  pweights <- test_dt[, pweight]
  model_dt <- dplyr::select(test_dt, include_vars)
  
  ## GET MODEL TYPES AND ANCHOR ITEMS
  types_dt <- data.table(question = names(model_dt), 
                         class = unlist(lapply(names(model_dt), function(x) return(model_dt[, class(get(x))]))))
  model_types <- rep(NA, nrow(types_dt))
  model_types[types_dt[class %in% c("numeric", "integer"), which = T]] <- "graded"
  model_types[types_dt[class == "logical", which = T]] <- "2PL"
  
  return(list(dt = model_dt, types = model_types, weights = pweights))
}

test_objects <- set_up_model(test_dt_adams)


# RUN MODELS --------------------------------------------------------------

run_efa <- function(objects, factors){
  print(factors)
  model <- mirt(test_objects$dt, factors, itemtype = test_objects$types)
  readr::write_rds(model, paste0(model_dir, "adams_efa_", factors, "factor_", date, ".rds"))
}

fs <- 1:9

lapply(fs, function(x) run_efa(objects = test_objects, factors = x))

# GET SUMMARIES -----------------------------------------------------------

## GET MODELS
read_file <- function(d, f){
  model <- readr::read_rds(paste0(model_dir, "adams_efa_", f, "factor_", d, ".rds"))
  return(model)
}

all_models <- function(date){
  all_files <- list.files(model_dir)
  files <- all_files[grepl("adams_efa_", all_files) & grepl(date, all_files)]
  factors <- gsub("factor.*$", "", files)
  factors <- gsub("^.*_", "", factors)
  models <- lapply(factors, function(x) read_file(d = date, f = x))
  return(models)
}

models <- all_models("2021_05_03")

## GET SUMMARIES
get_summaries <- function(model){
  sum <- summary(model, "oblimin", verbose = F)$rotF
  sum_dt <- as.data.table(sum)
  sum_dt[, item := rownames(sum)]
  factors <- colnames(sum)
  factors <- as.numeric(gsub("^[A-Z]*", "", factors))
  factors <- max(factors)
  sum_dt[, num_factors := factors]
  return(sum_dt)
}

summaries <- lapply(models, get_summaries)

## GET EIGENVALUES
get_eigen <- function(sum_dt){
  f_cols <- names(sum_dt)[grepl("^F", names(sum_dt))]
  eigen <- sum_dt[, apply(.SD, 2, function(x) sum(x^2)), .SDcols = f_cols]
  return(eigen)
}

eigenvalues <- lapply(summaries, get_eigen)

# VIZUALIZE ---------------------------------------------------------------

scree_plot <- function(eigen){
  plot_dt <- data.table(eigenvalue = eigen)
  plot_dt <- plot_dt[order(-eigenvalue)]
  plot_dt[, factor := 1:.N]
  gg <- ggplot(plot_dt, aes(x = factor, y = eigenvalue)) +
    geom_point() +
    geom_line() +
    scale_x_continuous(breaks = seq(1,plot_dt[, max(factor)])) +
    labs(x = "Factor", y = "Eigenvalues") +
    ggtitle(paste0(plot_dt[, max(factor)], " Factors")) +
    theme_classic()
  return(gg)
}

pdf(paste0("FILEPATH"))
lapply(eigenvalues, scree_plot)
dev.off()

plots <- lapply(eigenvalues[-1], scree_plot)
graph <- grid.arrange(grobs = plots, top = textGrob("Scree Plots using ADAMS data", gp=gpar(fontface="bold", fontsize=20)))
ggsave(filename = paste0("FILEPATH"), plot = graph, width = 8, height = 8)