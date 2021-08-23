##########################################################################
### Author: USERNAME
### Date: 11/8/2019
### Project: GBD Nonfatal Estimation
### Purpose: General IRT Equating Script
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

pacman::p_load(data.table, openxlsx, ggplot2, ggrepel, RColorBrewer, MASS, lmtest, nnet)
library("mirt", lib.loc = paste0("FILEPATH"))
library("arules", lib.loc = paste0("FILEPATH"))
library("gplots", lib.loc = paste0("FILEPATH"))
library("carData", lib.loc = paste0("FILEPATH"))
library("car", lib.loc = paste0("FILEPATH"))
date <- gsub("-", "_", Sys.Date())

# SET UP OBJECTS ----------------------------------------------------------

repo_dir <- paste0("FILEPATH")
harmonization_dir <- paste0("FILEPATH")
graph_dir <- paste0("FILEPATH")
model_dir <- paste0("FILEPATH")
draws <- paste0("draw_", 0:999)

sample <- "hrs" 

# SOURCE DATA PREP CODE ---------------------------------------------------

doc_dt <- as.data.table(read.xlsx("FILEPATH"))

source(paste0("FILEPATH")) ## SOURCE ADAMS DATA PREP CODE
source(paste0("FILEPATH")) ## SOURCE HRS DATA PREP CODE

# DISCRETIZE FUNCTION -----------------------------------------------------

create_interval_var <- function(dt, variable, new_variable, groups = 10){
  dt[, c(new_variable) := as.numeric(discretize(get(variable), method = "interval", breaks = groups,
                                                labels = 1:groups))]
  n <- nrow(dt[!is.na(get(new_variable))])
  cats <- sort(dt[!is.na(get(new_variable)), unique(get(new_variable))])
  for (cat in cats){
    cats <- sort(dt[!is.na(get(new_variable)), unique(get(new_variable))]) ## recalculate categories based on what we have
    if (nrow(dt[get(new_variable) == cat])/n < 0.05 & cat == min(cats)){ ## if less than 5% of data and not the max, subtract 1
      dt[get(new_variable) == cat, c(new_variable) := min(cats[!cats == cat])]
    }
    if (nrow(dt[get(new_variable) == cat])/n < 0.05 & !cat %in% c(max(cats), min(cats))){
      dt[get(new_variable) == cat, c(new_variable) := cat + 1]
    }
    if (nrow(dt[get(new_variable) == cat])/n < 0.05 & cat == max(cats)){
      dt[get(new_variable) == cat, c(new_variable) := max(cats[!cats == cat])]
    }
  }
  dt[, c(new_variable) := as.numeric(as.factor(get(new_variable)))]
  return(dt)
}

# COMBINE AND SET UP ANALYSIS ---------------------------------------------

dts <- ls(pattern = "test_dt")
dt_adams <- get(dts[grepl("adams", dts)])
dt_other <- get(dts[!grepl("adams", dts)])
other_name <- gsub("test_dt_", "", dts[!grepl("adams", dts)])

## DROP INDIVIDUALS IN BOTH ADAMS AND HRS FROM HRS
if (other_name == "hrs"){
  dt_other[, unique_id := paste0(hh_id, "_", pn)]
  dt_adams[, unique_id := paste0(gsub("^0", "", hh_id), "_", gsub("^0", "", pn))]
  dt_other <- dt_other[!unique_id %in% intersect(dt_other[, unique_id], dt_adams[, unique_id])]
}

## SET GROUP AND COMBINE
dt_adams[, group := "adams"]
dt_other[, group := other_name]
test_dt <- rbind(dt_adams, dt_other, fill = T, use.names = T)

## SET UP MODEL (no indiviudals who are all NA on either factor, no 0 pweights)
id_vars <- c("hh_id", "pn", "pweight", "group")
if (other_name == "hrs") id_vars <- c(id_vars, "unique_id")
include_vars <- names(test_dt)[!names(test_dt) %in% id_vars]
categories <- as.data.table(read.xlsx(paste0(j_root, "Project/Dementia/harmonization/item_categories.xlsx")))
func_qs <- include_vars[include_vars %in% categories[category == "function", item]]
cog_qs <- include_vars[include_vars %in% categories[category == "cognition", item]]
test_dt <- test_dt[!rowSums(is.na(dplyr::select(test_dt, cog_qs))) == ncol(dplyr::select(test_dt, cog_qs))] ## MAKE SURE NO ROWS ARE ALL NA
test_dt <- test_dt[!rowSums(is.na(dplyr::select(test_dt, func_qs))) == ncol(dplyr::select(test_dt, func_qs))]
test_dt <- test_dt[!pweight == 0]

## DISCRETIZE CONTINUOUS VARIABLES 
for (var in include_vars){
  if (test_dt[!is.na(get(var)), length(unique(get(var)))] > 12){ 
    test_dt <- create_interval_var(test_dt, var, var)
  } 
}
model_dt <- dplyr::select(test_dt, include_vars)

## GET MODEL TYPES AND INVARIANT ITEMS
types_dt <- data.table(question = names(model_dt), 
                       class  = unlist(lapply(names(model_dt), function(x) return(model_dt[, class(get(x))]))))
model_types <- rep(NA, nrow(types_dt))
model_types[types_dt[class %in% c("numeric", "integer"), which = T]] <- "graded"
model_types[types_dt[class == "logical", which = T]] <- "2PL"

## SEPEARTE FACTORS AND DESIGN
cat_mergecheck <- merge(data.table(item = names(model_dt)),
                        categories, all.x = T, by = "item")
if (nrow(cat_mergecheck[is.na(category)]) > 0) stop("Need to add items to item categories")

factor2_qs <- match(categories[category == "function" & item %in% names(model_dt), item], names(model_dt))
factor2_qs <- paste0(sort(factor2_qs), ", ", collapse = '')
factor2_qs <- gsub(", $", "", factor2_qs)
factor1_qs <- match(categories[category == "cognition" & item %in% names(model_dt), item], names(model_dt))
factor1_qs <- paste0(sort(factor1_qs), ", ", collapse = '')
factor1_qs <- gsub(", $", "", factor1_qs)

modeldesign2factor <- paste0("F1 = ", factor1_qs, "
                             F2 = ", factor2_qs, "
                             COV = F1*F2") 

# # RUN MODEL ---------------------------------------------------------------

model_equate2f <- multipleGroup(model_dt, modeldesign2factor, group = test_dt[, group],
                    itemtype = model_types, invariance = c(invariant, "free_means", "free_varcov"),
                    TOL = 0.001, method = 'QMCEM', SE = T) 

readr::write_rds(model_equate2f, paste0("FILEPATH"))

# DIAGNOSTIC PLOTS --------------------------------------------------------

get_thresholds <- function(n, thresh){
  print(n)
  par_matrix <- thresh[[n]]
  par_matrix <- par_matrix['par',]
  i <- names(thresh[n])
  threshold_vec <- par_matrix[names(par_matrix)[grepl("^d", names(par_matrix))]]
  disc <- par_matrix[names(par_matrix)[grepl("^a", names(par_matrix))]]
  disc <- disc[!disc == 0]
  threshold_vec <- -threshold_vec / disc
  tdt <- data.table(item = i, threshold = threshold_vec, label_name = c(i, rep('', length(threshold_vec) - 1)))
  return(tdt)
}

## FACTOR SCORE GRAPH
score_dt <- as.data.table(fscores(model_equate2f))
hist_dt <- data.table(factor = c(rep("Cognition", nrow(score_dt)), rep("Function", nrow(score_dt))),
                      score = c(score_dt[, F1], score_dt[, F2]),
                      group = rep(test_dt[, group], 2))

gg_density <- ggplot(hist_dt, aes(score, fill = as.factor(group))) +
  geom_density(alpha = 0.4) +
  labs(x = "Factor Score", y = "Density") +
  scale_fill_manual(values = c("red", "blue")) +
  facet_wrap(~factor) +
  theme_classic()

if (length(dif_items) == 0){
  ## THRESHOLD/LOADING GRAPH
  loadings_m <- summary(model_equate2f, verbose = F)[[1]]$rotF ## just picking one group
  loadings <- data.table(item = row.names(loadings_m), 
                         cog_f = loadings_m[, 'F1'],
                         func_f = loadings_m[, 'F2'])
  thresholds <- coef(model_equate2f)$adams ## just picking one group, thresholds are currently constrained to be the same for both groups
  thresholds <- thresholds[!grepl("GroupPars", names(thresholds))]
  
  thresh_dt <- rbindlist(lapply(1:length(thresholds), function(x) get_thresholds(n = x, thresh = thresholds)))
  graph_dt <- merge(thresh_dt, loadings, by = "item")
  graph_dt[cog_f == 0, cog_f := NA][func_f == 0, func_f := NA]
  graph_dt <- melt(graph_dt, id.vars = c("item", "threshold", "label_name"), measure.vars = c("cog_f", "func_f"), 
                   variable.name = "factor", value.name = "loading")
  graph_dt[factor == "cog_f", factor := "Cognition"][factor == "func_f", factor := "Function"]
  
  getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
  gg_thresh <- ggplot(graph_dt, aes(x = threshold, y = loading, color = item, label = label_name)) +
    geom_point() +
    geom_path(aes(group = item)) + 
    geom_text_repel(size = 2, segment.alpha = 0.3) +
    facet_wrap(~factor) +
    theme_classic() +
    scale_color_manual(values = getPalette(graph_dt[, length(unique(item))])) +
    labs(y = "Loading", x = "Ability") +
    theme(legend.position = "none")
  
 } else if (length(dif_items) > 0){
  loadings_m_a <- summary(model_equate2f, verbose = F)[[1]]$rotF ## first group is adams
  loadings_a <- data.table(item = row.names(loadings_m_a), 
                         cog_f = loadings_m_a[, 'F1'],
                         func_f = loadings_m_a[, 'F2'])
  thresholds_a <- coef(model_equate2f)$adams 
  thresholds_a <- thresholds_a[!grepl("GroupPars", names(thresholds_a))]
  thresh_dt_a <- rbindlist(lapply(1:length(thresholds_a), function(x) get_thresholds(n = x, thresh = thresholds_a)))
  thresh_dt_a[, sample := "adams"]
  both_dt_a <- merge(thresh_dt_a, loadings_a, by = "item")
  
  loadings_m_o <- summary(model_equate2f, verbose = F)[[2]]$rotF ## first group is adams
  loadings_o <- data.table(item = row.names(loadings_m_o), 
                           cog_f = loadings_m_o[, 'F1'],
                           func_f = loadings_m_o[, 'F2'])
  thresholds_o <- coef(model_equate2f)[[other_name]]
  thresholds_o <- thresholds_o[!grepl("GroupPars", names(thresholds_o))]
  thresh_dt_o <- rbindlist(lapply(1:length(thresholds_o), function(x) get_thresholds(n = x, thresh = thresholds_o)))
  thresh_dt_o[, sample := other_name]
  both_dt_o <- merge(thresh_dt_o, loadings_o, by = "item")
  
  graph_dt <- rbind(both_dt_a, both_dt_o)
  graph_dt[cog_f == 0, cog_f := NA][func_f == 0, func_f := NA]
  graph_dt <- melt(graph_dt, id.vars = c("item", "threshold", "label_name", "sample"), measure.vars = c("cog_f", "func_f"), 
                   variable.name = "factor", value.name = "loading")
  graph_dt[factor == "cog_f", factor := "Cognition"][factor == "func_f", factor := "Function"]
  
  getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
  gg_thresh <- ggplot(graph_dt, aes(x = threshold, y = loading, color = item, label = label_name)) +
    geom_point() +
    geom_path(aes(group = item)) + 
    geom_text_repel(size = 2, segment.alpha = 0.3) +
    facet_wrap(~factor+sample) +
    theme_classic() +
    scale_color_manual(values = getPalette(graph_dt[, length(unique(item))])) +
    labs(y = "Loading", x = "Ability") +
    theme(legend.position = "none")
}



