# Create catch-all weighted average to use for "Other" diagnosis

rm(list=ls())

library(data.table)
library(dplyr)

# Get IHME shared functions
invisible(sapply(list.files("FILEPATH", full.names = T), source))

ncodes <- get_rei_metadata(rei_set_id = 7, release_id = 9)
detailed <- ncodes[level == 2]

# Get pre-saved flat file of ncode incidence rate draws 
dt <- fread("FILEPATH")
dt <- merge(dt, ncodes[,.(rei_id,rei,rei_name)], by = c("rei_id"))

# Get corresponding population to use to convert to counts
population <- get_population(age_group_id = 22,
                             location_id = 102,
                             year_id = c(2015:2020),
                             sex_id = 3,
                             release_id = 9)

counts <- merge(dt, population, by = c("year_id","location_id","age_group_id","sex_id"))

# Do draw-level computation to get counts of n-code incidence and sum across years
pattern <- "^draw_\\d+$"
draw_cols <- grep(pattern, names(dt), value = TRUE)
counts <- counts[, (draw_cols) := lapply(.SD, function(x) x * population), .SDcols = draw_cols]
counts_sum <- counts[, lapply(.SD, sum), by = .(rei,rei_id,rei_name), .SDcols = draw_cols]

# Get DWs
dws <- fread("FILEPATH", na.strings = "")

# Drop any "UNTRT" dws
dws <- dws[treatment == "TRT" | is.na(treatment)]

# If both ST & LT present, only keep ST
dws <- dws %>%
  group_by(ncode) %>%
  filter(duration == max(duration) | is.na(duration)) %>%
  ungroup()
setDT(dws)

draws <- paste0("draw", 0:499)    

# Make sure DWs and counts are sorted the same way
dws <- dws[order(dws$ncode),]
counts_sum <- counts_sum[order(counts_sum$rei),]

# Subset to only keep draws
dw_draws <- dws[,6:505]
counts_draws <- counts_sum[, 4:503]

# Multiply the counts by the DWs
numerator <- dw_draws*counts_draws
numerator <- as.data.frame(t(colSums(numerator)))
denominator <- as.data.frame(t(colSums(counts_draws)))
weighted_result <- numerator/denominator

# Summarize!
mean_ui <- function (draw_dt) {
  dt <- copy(draw_dt)
  dt[, val := rowMeans(.SD), .SDcols = draws]
  dt[, lower := apply(.SD, 1, quantile, probs = 0.025), .SDcols = draws]
  dt[, upper := apply(.SD, 1, quantile, probs = 0.975), .SDcols = draws]
  dt[, (draws) := NULL]
  return(dt)
} 

setDT(weighted_result)
weighted_dw <- mean_ui(weighted_result)

fwrite(weighted_dw, "FILEPATH")