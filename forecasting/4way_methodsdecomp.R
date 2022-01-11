##########################################################################
### Author: USERNAME
### Date: 05/11/2019
### Project: GBD Nonfatal Estimation
### Purpose: Dementia Decomp graph
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

pacman::p_load(data.table, openxlsx, ggplot2, parallel, magrittr, grid, gridExtra, cowplot)
package_lib <- paste0(j_root, "FILEPATH")
library("gridGraphics", lib.loc = package_lib)
today_date <- gsub("-", "_", Sys.Date())

# SET OBJECTS -------------------------------------------------------------

forecast_dir <- "FILEPATH"
plot_dir <- paste0(j_root, "FILEPATH")
functions_dir <- paste0(functions_dir, "FILEPATH")
draws <- paste0("draw_", 0:999)

nocov_forecast_date <- "2021_02_19"
forecast_date <- "2021_03_04"

# GET FUNCITONS -----------------------------------------------------------

functs <- c("get_model_results.R", "get_location_metadata.R", "get_age_metadata.R", "get_draws.R",
            "get_population.R", "get_ids.R")
invisible(lapply(functs, function(x) source(paste0(functions_dir, x))))

# GET DEMOGRAPHICS --------------------------------------------------------

age_dt_full <- get_age_metadata(ID, gbd_round_id = ID)
age_dt <- copy(age_dt_full[age_group_id >= ID])

loc_dt <- get_location_metadata(ID, gbd_round_id = ID)

pop_dt <- readr::read_rds(paste0(forecast_dir, "population_draws_allage.rds"))
pop_dt[, population := rowMeans(.SD), .SDcols = draws]
pop_dt[, (draws) := NULL]
oldpop_dt <- readr::read_rds(paste0(forecast_dir, "oldpops.rds"))
pop_dt <- rbind(pop_dt[year_id == 2050], oldpop_dt[year_id == 2019], fill = T)

# GET DATA ----------------------------------------------------------------

result_dt <- readr::read_rds(paste0(forecast_dir, forecast_date, "FILEPATH"))
nocov_result_dt <- readr::read_rds(paste0(forecast_dir, nocov_forecast_date, "FILEPATH"))
setnames(nocov_result_dt, "mean_prev", "nocov_prev")
result_dt <- merge(result_dt[year_id == 2050], nocov_result_dt[year_id == 2050, .(age_group_id, sex_id, location_id, nocov_prev)],
                   by = c("age_group_id", "sex_id", "location_id"))
result_dt[, population := mean_num/mean_prev]
result_dt[, mean_prev := mean_prev/nocov_prev]
pre_dt <- readr::read_rds(paste0(forecast_dir, "FILEPATH"))
pre_dt[, population := mean_num/mean_prev]
pre_dt[, nocov_prev := mean_prev]
pre_dt[, mean_prev := 1]

## MERGE TO CREATE DECOMP DT WITH ALL AGE GROUPS - AGE RESTIRCTIONS = 0
decomp_dt <- rbind(result_dt[year_id == 2050], pre_dt[year_id == 2019])
decomp_dt <- decomp_dt[location_id %in% loc_dt[location_type_id %in% IDS, location_id] & !age_group_id %in% IDS]

# PREP AND CALCULATE ------------------------------------------------------

decomp_dt[, age_structure := population/sum(population), by = c("location_id", "sex_id", "year_id")]
decomp_dt[, population := sum(population), by = c("location_id", "sex_id", "year_id")]
decomp_dt[, year_id := ifelse(year_id == 2019, 1, 2)]
decomp_dt <- dcast(decomp_dt, location_id + age_group_id + sex_id ~ year_id, value.var = c("population", "age_structure", "mean_prev", "mean_num", "nocov_prev"))

## CALCULATE 
decomp_dt[, population_effect := ((age_structure_1 * mean_prev_1 * nocov_prev_1 + age_structure_2 * mean_prev_2 * nocov_prev_2) / 4 +
                                    (age_structure_1 * mean_prev_1 * nocov_prev_2 + age_structure_1 * mean_prev_2 * nocov_prev_1 +
                                     age_structure_2 * mean_prev_1 * nocov_prev_1 + age_structure_2 * mean_prev_2 * nocov_prev_1 +
                                     age_structure_2 * mean_prev_1 * nocov_prev_2 + age_structure_1 * mean_prev_2 * nocov_prev_2) / 12) * 
                                     (population_2 - population_1)]
decomp_dt[, age_structure_effect := ((population_1 * mean_prev_1 * nocov_prev_1 + population_2 * mean_prev_2 * nocov_prev_2) / 4 +
                                       (population_1 * mean_prev_1 * nocov_prev_2 + population_1 * mean_prev_2 * nocov_prev_1 +
                                          population_2 * mean_prev_1 * nocov_prev_1 + population_2 * mean_prev_2 * nocov_prev_1 +
                                          population_2 * mean_prev_1 * nocov_prev_2 + population_1 * mean_prev_2 * nocov_prev_2) / 12) * 
                                         (age_structure_2 - age_structure_1)]
decomp_dt[, xtracov_effect := ((population_1 * age_structure_1 * nocov_prev_1 + population_2 * age_structure_2 * nocov_prev_2) / 4 +
                             (population_1 * age_structure_1 * nocov_prev_2 + population_1 * age_structure_2 * nocov_prev_1 +
                              population_2 * age_structure_1 * nocov_prev_1 + population_2 * age_structure_2 * nocov_prev_1 +
                              population_2 * age_structure_1 * nocov_prev_2 + population_1 * age_structure_2 * nocov_prev_2) / 12) * 
                                (mean_prev_2 - mean_prev_1)]
decomp_dt[, paf_effect := ((population_1 * age_structure_1 * mean_prev_1 + population_2 * age_structure_2 * mean_prev_2) / 4 +
                             (population_1 * age_structure_1 * mean_prev_2 + population_1 * age_structure_2 * mean_prev_1 +
                              population_2 * age_structure_1 * mean_prev_1 + population_2 * age_structure_2 * mean_prev_1 +
                              population_2 * age_structure_1 * mean_prev_2 + population_1 * age_structure_2 * mean_prev_2) / 12) * 
                                (nocov_prev_2 - nocov_prev_1)]


decomp_dt <- decomp_dt[, lapply(.SD, sum), by = "location_id",
                       .SDcols = c("mean_num_1", "mean_num_2", "population_effect", "age_structure_effect", "paf_effect", "xtracov_effect")]
aaic_dt <- copy(decomp_dt)

# PREP GRAPH --------------------------------------------------------------

decomp_dt[, total_pct := (mean_num_2 - mean_num_1) / mean_num_1] 

decomp_dt[, c("population_effect", "paf_effect", "xtracov_effect", "age_structure_effect") :=
            .(population_effect / mean_num_1,
              paf_effect / mean_num_1,
              xtracov_effect / mean_num_1,
              age_structure_effect / mean_num_1)]
decomp_dt[, c("mean_num_1", "mean_num_2") := .(NULL, NULL)]

decomp_dt <- melt(decomp_dt,
                  id.vars = c("location_id"),
                  value.vars = c(grep("effect", names(decomp_dt), value = TRUE), "total_pct"),
                  variable.name = "factor_format",
                  value.name = "change")
total_dt <- decomp_dt[factor_format=="total_pct",list(location_id,change)]
setnames(total_dt,"change","total")
decomp_dt <- merge(decomp_dt[factor_format!="total_pct",],total_dt,by=c("location_id"))

decomp_dt <- merge(decomp_dt, loc_dt[, .(location_id, lancet_label, sort_order)])
decomp_dt <- decomp_dt[order(-sort_order)]
decomp_dt[, lancet_label := factor(lancet_label, levels = unique(lancet_label))]
decomp_dt[, factor_format := factor(factor_format,
                                    level = c("age_structure_effect", "population_effect", "paf_effect", "xtracov_effect"),
                                    label = c("Change due to population ageing",
                                              "Change due to population growth",
                                              "Change due to GBD risk factors",
                                              "Change due to education"))]

# MAKE GRAPH --------------------------------------------------------------

colors <- c("Change due to population ageing" = "#A2C851",
            "Change due to population growth" = "#218380",
            "Change due to GBD risk factors" = "#ffbc42",
            "Change due to education" = "navyblue")
plot <- ggplot() +
  geom_point(data=decomp_dt,aes(x=lancet_label,y=total),
             size=-1,na.rm=T,color="white")+ 
  geom_bar(data = decomp_dt[change < 0,],
           aes(x = lancet_label, y = change*100, fill = factor_format), stat = "identity", width = .75, na.rm = TRUE) +
  geom_bar(data = decomp_dt[change >= 0,],
           aes(x = lancet_label, y = change*100, fill = factor_format), stat = "identity", width = .75, na.rm = TRUE) +
  geom_point(data = decomp_dt,
             aes(x = lancet_label, y = total*100, color = "Total percent change"),
             size = 3, na.rm = TRUE) +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = "dashed") +
  scale_fill_manual(values = colors, drop = FALSE) +
  scale_color_manual(values = c("Total percent change" = "black")) +
  coord_flip() +
  labs(y = "Percent Change (%)", x = "", fill = "") +
  theme_bw() +
  theme(plot.title = element_text(vjust = 2, hjust = 0, size = 15),
        axis.text.x = element_text(size=7),
        legend.key.size = unit(0.5, "cm"),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow = 2, reverse = FALSE))

ggsave(plot = plot, filename = paste0(plot_dir, "decomp_", today_date, ".jpg"), height = 7, width = 10)

pdf(paste0(plot_dir, "decomp_", today_date, ".pdf"), width = 10)
plot
dev.off()
