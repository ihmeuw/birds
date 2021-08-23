##########################################################################
### Author: USERNAME
### Date: 11/18/19
### Project: GBD Nonfatal Estimation
### Purpose: Item Response Theory Prevalence Calculation ADAMS General
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

pacman::p_load(data.table, openxlsx, ggplot2, ggrepel, RColorBrewer, MASS, lmtest, nnet, 
               gridExtra, grid, lattice, survey, stringr)
package_lib <- paste0("FILEPATH")
library("mirt", lib.loc = package_lib)
library("gplots", lib.loc = package_lib)
library("pROC", lib.loc = package_lib)
library("carData", lib.loc = package_lib)
library("car", lib.loc = package_lib)
library("ResourceSelection", lib.loc = package_lib)
library("caret", lib.loc = package_lib)
date <- gsub("-", "_", Sys.Date())

# SET UP OBJECTS ----------------------------------------------------------

repo_dir <- paste0("FILEPATH")
harmonization_dir <- paste0("FILEPATH")
prev_dir <- paste0("FILEPATH")
graph_dir <- paste0("FILEPATH")
results_dir <- paste0("FILEPATH")
model_dir <- paste0("FILEPATH")
functions_dir <- paste0("FILEPATH")
draws <- paste0("draw_", 0:999)

sample <- "hrs"

# FUNCTIONS ---------------------------------------------------------------

functs <- c("get_age_metadata.R", "get_ids.R", "get_outputs.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

setup_design <- function(df, var){
  options(survey.lonely.psu = 'adjust')
  options(survey.adjust.domain.lonely = TRUE)
  for (i in c("strata", "psu", "pweight")) {
    ## Assign to *_formula the variable if it exists and nonmissing, else NULL
    assign(paste0(i, "_formula"),
           as.formula(ifelse(i %in% names(df) & nrow(df[!is.na(i)]) > 0, paste("~", i), "~1")))
  }
  ## Set svydesign
  if (strata_formula == ~1) strata_formula <- NULL
  return(svydesign(id = psu_formula, weight = pweight_formula, strat = strata_formula, data = df, nest = TRUE))
}

sims <- function(object, n.sims){
  object.class <- class(object)[[1]]
  summ <- summary (object, correlation=TRUE, dispersion = object$dispersion)
  coef <- summ$coef[,1:2,drop=FALSE]
  dimnames(coef)[[2]] <- c("coef.est","coef.sd")
  beta.hat <- coef[,1,drop=FALSE]
  sd.beta <- coef[,2,drop=FALSE]
  corr.beta <- summ$corr
  n <- summ$df[1] + summ$df[2]
  k <- summ$df[1]
  V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
  beta <- array (NA, c(n.sims,k))
  dimnames(beta) <- list (NULL, dimnames(beta.hat)[[1]])
  for (s in 1:n.sims){
    beta[s,] <- MASS::mvrnorm (1, beta.hat, V.beta)
  }
  
  beta2 <- array (0, c(n.sims,length(coefficients(object))))
  dimnames(beta2) <- list (NULL, names(coefficients(object)))
  beta2[,dimnames(beta2)[[2]]%in%dimnames(beta)[[2]]] <- beta
  
  sigma <- rep (sqrt(summ$dispersion), n.sims)
  
  ans <- list(coef = beta2, sigma = sigma)
  return(ans)
}

mround <- function(x,base){ 
  base*floor(x/base) 
} 

summaries <- function(dt, draw_vars){
  sum <- copy(dt)
  sum[, mean := rowMeans(.SD), .SDcols = draw_vars]
  sum[, lower := apply(.SD, 1, quantile, probs= 0.025, type = 8), .SDcols = draw_vars]
  sum[, upper := apply(.SD, 1, quantile, probs=0.975, type = 8), .SDcols = draw_vars]
  sum[, c(draw_vars) := NULL]
  return(sum)
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# SOURCE DATA PREP CODE ---------------------------------------------------

source(paste0("FILEPATH")) ## SOURCE ADAMS DATA PREP CODE
source(paste0("FILEPATH")) ## SOURCE HRS DATA PREP CODE

categories <- as.data.table(read.xlsx(paste0("FILEPATH")))
setnames(categories, "factor", "category")

# COMBINE AND SET UP ANALYSIS ---------------------------------------------

dts <- ls(pattern = "test_dt")
dt_adams <- get(dts[grepl("adams", dts)])
dt_other <- get(dts[!grepl("adams", dts)])
other_name <- gsub("test_dt_", "", dts[!grepl("adams", dts)])

## DROP INDIVIDUALS IN BOTH ADAMS AND HRS FROM HRS
dt_other[, unique_id := paste0(hh_id, "_", pn)]
dt_adams[, unique_id := paste0(gsub("^0", "", hh_id), "_", gsub("^0", "", pn))]
if (other_name == "hrs"){
    dt_other <- dt_other[!unique_id %in% intersect(dt_other[, unique_id], dt_adams[, unique_id])]
}

## GET OTHER INFO
adams_demog <- copy(adamsa[, .(hhid, pn, age, sex_id)])
setnames(adams_demog, c("hhid", "sex_id"), c("hh_id", "sex"))
adams_demog[, unique_id := paste0(gsub("^0", "", hh_id), "_", gsub("^0", "", pn))]
adams_demog[, c("hh_id", "pn") := NULL]
other_demog <- get(other_name)[, .(hh_id, pn, age_year, sex_id)]
other_demog[, sex := ifelse(sex_id == 1, "male", "female")]
other_demog[, sex_id := NULL]
setnames(other_demog, "age_year", "age")
other_demog[, unique_id := paste0(hh_id, "_", pn)]
other_demog[, c("hh_id", "pn") := NULL]
dt_other <- merge(dt_other, other_demog, by = c("unique_id"), sort = F)
dt_adams <- merge(dt_adams, adams_demog, by = c("unique_id"), sort = F)

## SET GROUP AND COMBINE
dt_adams[, group := "adams"]
dt_adams[, unique_id := paste0(gsub("^0", "", hh_id), "_", gsub("^0", "", pn))]
dt_other[, group := other_name]
dt_other[, unique_id := paste0(hh_id, "_", pn)]
test_dt <- rbind(dt_adams, dt_other, fill = T, use.names = T)

## SET UP MODEL (no indiviudals who are all NA on either factor, no 0 sample weights)
id_vars <- c("hh_id", "pn", "pweight", "group", "age", "sex")
if (other_name == "hrs") id_vars <- c(id_vars, "unique_id")
include_vars <- names(test_dt)[!names(test_dt) %in% id_vars]
func_qs <- include_vars[include_vars %in% categories[category == "function", item]]
cog_qs <- include_vars[include_vars %in% categories[category == "cognition", item]]
test_dt <- test_dt[!rowSums(is.na(dplyr::select(test_dt, cog_qs))) == ncol(dplyr::select(test_dt, cog_qs))] ## MAKE SURE NO ROWS ARE ALL NA
test_dt <- test_dt[!rowSums(is.na(dplyr::select(test_dt, func_qs))) == ncol(dplyr::select(test_dt, func_qs))]
test_dt <- test_dt[!pweight == 0]
model_dt <- dplyr::select(test_dt, include_vars)

# GET MOST RECENT MODELS --------------------------------------------------

all_models <- list.files(model_dir)
models <- all_models[grepl(paste0("model_2factorsSE_", other_name), all_models)]
dates <- gsub(paste0("model_2factorsSE_", other_name, "_adams_"), "", models)
dates <- gsub(".rds$", "", dates)
dates <- gsub("_", "-", dates)

model_equate2f <- readr::read_rds(paste0("FILEPATH"))

# GET FACTOR SCORES AND OTHER INFORMATION ---------------------------------

scores2f <- as.data.table(fscores(model_equate2f))
score_dt <- data.table(group = test_dt[, group], score1 = scores2f$F1, score2 = scores2f$F2, 
                       unique_id = test_dt[, unique_id])

# SET UP AND RUN REGRESSION -----------------------------------------------

## SET UP DATA
prev_values <- c(1,3,10:11,13)
new_values <- c(1:4, 8, 10:11, 13, 15:19)

dementia <- c("probable AD", "possible AD", "probable vascular dementia", "possible vascular dementia", "normal pressure hydroencephalus", "dementia of undermined etiology", "frontal lobe dementia",
              "alcoholic dementia", "hypoperfusion dementia", "probable lewy body dementia")
regress_dt <- dplyr::select(adamsa, c("hhid", "pn", "strata", "sweight_cross", "cluster", "sex_id", "age", "consensus_diagnosis_final", "consensus_2diagnosis_final", "consensus_3diagnosis_final"))
regress_dt[, unique_id := paste0(gsub("^0", "", hhid), "_", gsub("^0", "", pn))]
setnames(regress_dt, c("sex_id", "cluster", "sweight_cross"), c("sex", "psu", "pweight"))
regress_dt[consensus_diagnosis_final %in% dementia | consensus_2diagnosis_final %in% dementia | consensus_3diagnosis_final %in% dementia, dementia := 1]
regress_dt[is.na(dementia), dementia := 0]
regress_dt <- merge(regress_dt, score_dt[group == "adams"], by = "unique_id", all.x = T)
regress_dt[, sex := factor(sex, levels = c("male", "female"))]
regress_dt[, age_cat := cut(age, c(seq(mround(regress_dt[, min(age)], 5), 95, by = 5), 120), right = F)]
regress_dt <- regress_dt[!is.na(score1)] ## only people that have scores

## DO REGRESSION
adams_design <- setup_design(regress_dt, var = "dementia")
regression <- survey::svyglm(dementia ~ score1 + score2 + age + sex, design = adams_design, data = regress_dt, family = binomial(link = "logit")) ## omits missing values on purpose

# HOSMER-LEMESHOW GOODNESS OF FIT -----------------------------------------

fit <- hoslem.test(regression$y, fitted(regression), g = 10)

# CROSS-VALIDATION --------------------------------------------------------

set.seed(12494)
train <- createFolds(regress_dt$dementia, k = 10, list = T, returnTrain = T)

prediction_draws <- function(X_pred, model){
  n_pred <- dim(X_pred)[1]
  sim <- sims(model, 1)
  p_pred <- boot::inv.logit(X_pred %*% t(sim$coef))
  y_pred <- rbinom(n_pred, 1, p_pred)
}

get_preds <- function(fold){
  print(fold)
  dt <- copy(regress_dt[train[[fold]]])
  regression <- survey::svyglm(dementia ~ score1 + score2 + age + sex, design = adams_design, data = dt, family = binomial(link = "logit")) ## omits missing values on purpose
  test_dt <- copy(regress_dt[-train[[fold]]])
  mean_preds <- predict(regression, type = "response")[-train[[fold]]]
  
  X_apred <- as.matrix(cbind(rep(1, nrow(test_dt)), test_dt[, .(score1, score2, age)], as.numeric(test_dt$sex == "female")))
  apred_draws <- as.data.table(replicate(1000, prediction_draws(X_apred, regression)))
  setnames(apred_draws, names(apred_draws), draws)
  apred_draws <- cbind(test_dt, apred_draws)
  apred_draws[, mean := mean_preds]
  return(apred_draws)
}

cv_predicts <- rbindlist(lapply(1:10, get_preds))

# ROC Curve ---------------------------------------------------------------

roc_curve <- plot.roc(cv_predicts$dementia, cv_predicts$mean, percent = T, ci = T, print.auc = T)
ci_object <- ci.se(roc_curve, specificities = seq(0,100,2))
plot(ci_object, type = "shape", col = "#1c61b6AA")

# SENSITIVITY AND SPECIFICITY ---------------------------------------------

sens_funct <- function(draw, dt){
  tp <- nrow(dt[dementia == 1 & get(paste0("draw_", draw)) == 1])
  fn <- nrow(dt[dementia == 1 & get(paste0("draw_", draw)) == 0])
  sens <- tp/(tp+fn)
  return(sens)
}
spec_funct <- function(draw, dt){
  tn <- nrow(dt[dementia == 0 & get(paste0("draw_", draw)) == 0])
  fp <- nrow(dt[dementia == 0 & get(paste0("draw_", draw)) == 1])
  spec <- tn/(tn+fp)
  return(spec)
}
accur_funct <- function(draw, dt){
  correct <- nrow(dt[(dementia == 0 & get(paste0("draw_", draw)) == 0) | (dementia == 1 & get(paste0("draw_", draw)) == 1)])
  accuracy <- correct/nrow(dt)
  return(accuracy)
}
draw_num <- 0:999
sens <- unlist(parallel::mclapply(1:length(draw_num), function(x) sens_funct(draw = draw_num[x], dt = cv_predicts), mc.cores = 9))
spec <- unlist(parallel::mclapply(1:length(draw_num), function(x) spec_funct(draw = draw_num[x], dt = cv_predicts), mc.cores = 9))
accuracy <- unlist(parallel::mclapply(1:length(draw_num), function(x) accur_funct(draw = draw_num[x], dt = cv_predicts), mc.cores = 9))

vector_sum <- function(x){
  print(paste0(round(mean(x), 2), " (", round(quantile(x, 0.025), 2), " - ", round(quantile(x, 0.975), 2), ")"))
}
vector_sum(sens)
vector_sum(spec)
vector_sum(accuracy)

# PREDICT IN OTHER SAMPLE -------------------------------------------------

predict_dt <- get(paste0(other_name, "_all"))
vars <- c("pn", "hh_id", "pweight", "age_year", "sex_id")
if ("strata" %in% names(predict_dt)) vars <- c(vars, "strata")
if ("psu" %in% names(predict_dt)) vars <- c(vars, "psu")
predict_dt[is.na(pweight), pweight := 0]
predict_dt <- predict_dt[, c(vars), with = F]
predict_dt[, unique_id := paste0(hh_id, "_", pn)]
predict_dt[, sex := ifelse(sex_id == 1, "male", "female")]
predict_dt[, sex := factor(sex, levels = c("male", "female"))]
predict_dt[, sex_id := NULL]
setnames(predict_dt, "age_year", "age")
predict_dt <- merge(predict_dt, score_dt, by = "unique_id", all.x = T)
predict_dt[, age_cat := cut(age, c(seq(mround(predict_dt[, min(age, na.rm = T)], 5), 95, by = 5), 120), right = F)]

X_pred <- as.matrix(cbind(rep(1, nrow(predict_dt)), predict_dt[, .(score1, score2, age)], as.numeric(predict_dt$sex == "female")))
prediction_draws <- function(X_pred, model){
  n_pred <- dim(X_pred)[1]
  sim <- sims(model, 1)
  p_pred <- boot::inv.logit(X_pred %*% t(sim$coef))
  y_pred <- rbinom(n_pred, 1, p_pred)
}
pred_draws <- as.data.table(replicate(1000, prediction_draws(X_pred, regression)))
setnames(pred_draws, names(pred_draws), draws)
pred_draws <- cbind(predict_dt, pred_draws)

# GET MEANS USING SURVEY WEIGHTS ------------------------------------------

other_design <- setup_design(pred_draws, var = "draw_0")
other_design_subset <- subset(other_design, subset = (!is.na(score1) & age >= 70))

## INCORPORATE UNCERTAINTY FROM SURVEY DESIGN
predicted_means <- as.data.table(survey::svyby(~draw_0, by = ~sex+age_cat, design = other_design_subset, svymean))[, 1:3]
new_mean_col <- function(x){
  new_col <- c(survey::svyby(~get(paste0("draw_", x)),
                             by = ~sex+age_cat, 
                             design = other_design_subset, 
                             svymean)[3])
  return(new_col)
}
new_cols <- parallel::mclapply(1:999, new_mean_col, mc.cores = 9) 
for (n in 1:999){
  new_col <- new_cols[[n]]
  predicted_means[, paste0("draw_", n) := new_col]
}
predicted_means[, se := apply(.SD, 1, sd), .SDcols = draws]

# GET SUMMARIES -----------------------------------------------------------

predicted_summaries <- summaries(predicted_means, draws)
predicted_summaries[, cv := se/mean]
predicted_summaries[, outlier := ifelse(cv < 1 & !is.na(cv), 0, 1)]

# COLLAPSE IF NOT ENOUGH PEOPLE/TOO UNCERTAIN -----------------------------

collapse_ages <- function(dt, s){
  dt <- copy(dt[sex == s])
  while (nrow(dt[outlier == 1]) > 0){
    dt[, sample_size := (mean*(1-mean)/se^2)] 
    dt[, cases := mean*sample_size]
    dt[, `:=` (age_start = gsub("^.", "", age_cat), age_end = gsub("^.*\\,", "", age_cat))]
    dt[, `:=` (age_start = as.numeric(gsub("\\,.*$", "", age_start)), age_end = as.numeric(gsub("\\)$", "", age_end)))]
    ages <- sort(dt[, age_start])
    for (age in ages){
      ages <- sort(dt[, age_start])
      if (!age %in% ages | nrow(dt[age_start == age & outlier == 1]) == 0){ ## skip age if already aggregated up into age grouop or if no outliers
        next
      }
      if (dt[age_start == age, outlier] == 1 & !age == max(ages)){ ## if outlier and is not the max
        dt[age_start %in% c(age, ages[which(ages == age) + 1]), collapse := 1]
      }
      if (dt[age_start == age, outlier] == 1 & age == max(ages)){
        dt[age_start %in% c(age, ages[which(ages == age) - 1]), collapse := 1]
      }
      dt[collapse == 1, `:=` (cases = sum(cases), sample_size = sum(sample_size),
                              age_start = min(age_start), age_end = max(age_end))]
      dt[collapse == 1, `:=` (age_cat = paste0("[", age_start, ",", age_end, ")"),
                              mean = cases/sample_size)]
      z <- qnorm(0.975)
      dt[collapse == 1, se := sqrt(mean*(1-mean)/sample_size + z^2/(4*sample_size^2))]
      dt[, cv := se/mean]
      dt[, outlier := ifelse(cv < 1 & !is.na(cv), 0, 1)]
      dt[collapse == 1, lower := ifelse(cases <= 5, mean - z*se, 
                                        1/(1+z^2/sample_size)*(mean + z^2/(2*sample_size) - z*se))]
      dt[lower < 0, lower := 0]
      dt[collapse == 1, upper := ifelse(cases > 5, mean + z*se,
                                        1/(1+z^2/sample_size)*(mean + z^2/(2*sample_size) + z*se))]
      dt[upper > 1, upper := 0]
      dt <- unique(dt, by = c("sex", "age_cat"))
      dt[, collapse := NULL]
    }
  }
  dt <- dt[, .(sex, age_cat, se, mean, lower, upper, cv, outlier)]
  return(dt)
}

collapsed_summaries <- rbindlist(lapply(c("male", "female"), function(x) collapse_ages(dt = predicted_summaries, s = x)))

# SAVE FILES AND GRAPHS ---------------------------------------------------

## SAVE FILES
readr::write_rds(list(main_preds = pred_draws, adams_preds = cv_predicts, model = regression),
                paste0("FILEPATH"))

