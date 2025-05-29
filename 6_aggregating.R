# pulling in diabetes population to calculate:
# 1. absolute numbers
# 2. proportion diagnosed among all diabetes, proportion treated among diagnosed & proportion controlled among treated
# 3. age-standardized
# 4. all ages (15+)
# 5. aggregate 15-39, 40-64, 65+
# 6. aggregate 10 yr age groups
# 7. sex ratio
# 8. absolute change over time 

library(dplyr)
library(tidyr)
library(ggplot2)

invisible(sapply(list.files("FILEPATH", full.names = T), source))

draws <- paste0("draw_", 0:249)

start_year <- 2000
end_year <- 2023


## summarize function
summaries <- function(dt, draw_vars){
  sum <- as.data.table(copy(dt))
  sum[, mean := rowMeans(.SD), .SDcols = draw_vars]
  sum[, standard_error := apply(.SD, 1, sd), .SDcols = draw_vars]  
  sum[, lower := apply(.SD, 1, quantile, probs= 0.025), .SDcols = draw_vars]
  sum[, upper := apply(.SD, 1, quantile, probs=0.975), .SDcols = draw_vars]
  sum[, c(draw_vars) := NULL]
  return(sum)
}

## metadata 
loc_meta <- as.data.table(get_location_metadata(location_set_id = 35, release_id = 16))
loc_meta2 <- loc_meta[level %in% c(0, 1, 2, 3), ]
loc_meta_merge <- subset(loc_meta, select = c(location_id, location_name, region_name, super_region_name))

## reading in post-scale draws
dir <- "FILEPATH"

undiagnosed <- as.data.table(read.csv(paste0(dir, "undiagnosed_post_scaled_draws.csv")))
diagnosed_untx <- as.data.table(read.csv(paste0(dir, "diagnosed_untx_post_scaled_draws.csv")))
diagnosed_tx_uncont <- as.data.table(read.csv(paste0(dir, "diagnosed_tx_uncont_post_scaled_draws.csv")))
diagnosed_tx_cont <- as.data.table(read.csv(paste0(dir, "diagnosed_tx_cont_post_scaled_draws.csv")))

undiagnosed <- undiagnosed[, measure := "undiagnosed"]
diagnosed_untx <- diagnosed_untx[, measure := "diagnosed_untx"]
diagnosed_tx_uncont <- diagnosed_tx_uncont[, measure := "diagnosed_tx_uncont"]
diagnosed_tx_cont <- diagnosed_tx_cont[, measure := "diagnosed_tx_cont"]

list <- list(undiagnosed, diagnosed_untx, diagnosed_tx_uncont, diagnosed_tx_cont)
df <- rbindlist(list)

df <- subset(df, select = -c(X))

df[, denom := "diabetes"]
df[, metric := "proportion"]

## pulling dm results and converting to counts
dm <- as.data.table(read.csv(paste0(dir, "dm_count_draws.csv")))
dm <- subset(dm, select = -c(X))


# 1. absolute numbers -----------------------------------------------------------------------------
## multiplying by dm prevalence counts
df <- pivot_longer(df, cols = draws, names_to = "draw", values_to = "val")

df <- merge(df, dm, all.x=TRUE, by = c("age_group_id", "year_id", "sex_id", "location_id", "draw"))
df <- as.data.table(df)

df_count <- copy(df)
df_count[, val := val*dm_count]
df_count[, metric := "count"]


# 2. proportion diagnosed among all diabetes, proportion treated among diagnosed & proportion controlled among treated -----------------------------------------------------------------------------
## calculating counts of diagnosed and treated for denominators
df_count_diag <- df_count[measure == "undiagnosed", ]
df_count_diag[, diag_count := dm_count-val]
df_count_diag <- subset(df_count_diag, select = -c(measure, val, dm_count))

df_count_tx <- df_count[measure %in% c("diagnosed_tx_uncont", "diagnosed_tx_cont"), ]
df_count_tx <- group_by(df_count_tx, year_id, sex_id, location_id, age_group_id, metric, draw)
df_count_tx <- summarize(df_count_tx, tx_count = sum(val))

## calculating proportions
df_prop_tx_diag <- merge(df_count_diag, df_count_tx)
df_prop_tx_diag[, val := tx_count/diag_count]
df_prop_tx_diag[, metric := "proportion"]
df_prop_tx_diag[, measure := "tx"]
df_prop_tx_diag[, denom := "diagnosed"]
df_prop_tx_diag <- subset(df_prop_tx_diag, select = -c(tx_count))

df_prop_cont_tx <- merge(df_count[measure == "diagnosed_tx_cont",], df_count_tx)
df_prop_cont_tx[, val := val/tx_count]
df_prop_cont_tx[, metric := "proportion"]
df_prop_cont_tx[, denom := "treated"]
df_prop_cont_tx <- subset(df_prop_cont_tx, select = -c(dm_count))

df_prop_diag <- df_count[measure == "undiagnosed", ]
df_prop_diag[, val := (dm_count-val)/dm_count]
df_prop_diag[, metric := "proportion"]
df_prop_diag[, measure := "diagnosed"]
df_prop_diag[, denom := "diabetes"]

df_props <- rbind(df, df_prop_tx_diag, fill = TRUE)
df_props <- rbind(df_props, df_prop_cont_tx, fill = TRUE)
df_props <- rbind(df_props, df_prop_diag, fill = TRUE)

## calculating counts
df_count_tx_diag <- merge(df_count_diag, df_count_tx)
df_count_tx_diag[, val := tx_count]
df_count_tx_diag[, metric := "count"]
df_count_tx_diag[, measure := "tx"]
df_count_tx_diag[, denom := "diagnosed"]
df_count_tx_diag <- subset(df_count_tx_diag, select = -c(tx_count))

df_count_cont_tx <- merge(df_count[measure == "diagnosed_tx_cont",], df_count_tx)
df_count_cont_tx[, metric := "count"]
df_count_cont_tx[, denom := "treated"]
df_count_cont_tx <- subset(df_count_cont_tx, select = -c(dm_count))

df_count_diag <- df_count[measure == "undiagnosed", ]
df_count_diag[, val := dm_count-val]
df_count_diag[, metric := "count"]
df_count_diag[, measure := "diagnosed"]
df_count_diag[, denom := "diabetes"]

df_count <- rbind(df_count, df_count_tx_diag, fill = TRUE)
df_count <- rbind(df_count, df_count_cont_tx, fill = TRUE)
df_count <- rbind(df_count, df_count_diag, fill = TRUE)


# 3. age-standardized -----------------------------------------------------------------------------
## calculating age weights based on global both sex 2010 age pattern
age_weights <- dm[location_id == 1 & sex_id == 3 & year_id == 2010, ]
age_weights <- group_by(age_weights, age_group_id)
age_weights <- summarise(age_weights, mean_draws = mean(dm_count))
age_weights_mean <- sum(age_weights$mean_draws)
age_weights <- as.data.table(age_weights)

age_weights[, weight := mean_draws/age_weights_mean]

## calculating age-stand
df_std <- merge(df_props, age_weights, by = "age_group_id")
df_std[, val_weighted := val*weight]
df_std <- group_by(df_std, year_id, sex_id, location_id, measure, denom, metric, draw)
df_std <- summarize(df_std, val = sum(val_weighted))
df_std <- as.data.table(df_std)
df_std <- df_std[, age_group_id := 27]


# 4. all ages (15+) -----------------------------------------------------------------------------
## aggregating 
df_all_age <- group_by(df_count, year_id, sex_id, location_id, measure, metric, denom, draw)
df_all_age <- summarize(df_all_age, sum_num = sum(val), sum_dm = sum(dm_count), sum_diag = sum(diag_count), sum_tx = sum(tx_count))

df_all_age_count <- as.data.table(df_all_age)
df_all_age_count <- df_all_age_count[, val := sum_num]
df_all_age_count <- df_all_age_count[, dm_count := sum_dm]
df_all_age_count <- df_all_age_count[, diag_count := sum_diag]
df_all_age_count <- df_all_age_count[, tx_count := sum_tx]
df_all_age_count <- df_all_age_count[, age_group_id := 22]
df_all_age_count <- subset(df_all_age_count, select = -c(sum_num, sum_dm, sum_diag, sum_tx))

df_count <- rbind(df_count, df_all_age_count, fill = TRUE)

## redividing by denominator
df_all_age <- as.data.table(df_all_age)
df_all_age <- df_all_age[denom == "diabetes", `:=` (val = (sum_num/sum_dm), metric = "proportion")]
df_all_age <- df_all_age[denom == "diagnosed", `:=` (val = sum_num/sum_diag, metric = "proportion")]
df_all_age <- df_all_age[denom == "treated", `:=` (val = sum_num/sum_tx, metric = "proportion")]
df_all_age <- df_all_age[, age_group_id := 22]
df_all_age <- df_all_age[, dm_count := sum_dm]
df_all_age <- df_all_age[, diag_count := sum_diag]
df_all_age <- df_all_age[, tx_count := sum_tx]
df_all_age <- subset(df_all_age, select = -c(sum_num, sum_dm, sum_diag, sum_tx))

df_props <- rbind(df_props, df_all_age, fill = TRUE)


# 5. aggregate 15-39, 40-64, 65+ -----------------------------------------------------------------------------
## aggregating
df_wide_age <- df_count[!age_group_id %in% c(22, 27), ]
df_wide_age[age_group_id %in% c(8:12), age_group_id := 197]
df_wide_age[age_group_id %in% c(13:17), age_group_id := 220]
df_wide_age[age_group_id %in% c(18:20, 30:32, 235), age_group_id := 154]

df_wide_age <- group_by(df_wide_age, year_id, sex_id, location_id, measure, metric, denom, draw, age_group_id)
df_wide_age <- summarize(df_wide_age, sum_num = sum(val), sum_dm = sum(dm_count), sum_diag = sum(diag_count), sum_tx = sum(tx_count))

df_wide_age_count <- as.data.table(df_wide_age)
df_wide_age_count <- df_wide_age_count[, val := sum_num]
df_wide_age_count <- df_wide_age_count[, dm_count := sum_dm]
df_wide_age_count <- df_wide_age_count[, diag_count := sum_diag]
df_wide_age_count <- df_wide_age_count[, tx_count := sum_tx]
df_wide_age_count <- subset(df_wide_age_count, select = -c(sum_num, sum_dm, sum_diag, sum_tx))

df_count <- rbind(df_count, df_wide_age_count, fill = TRUE)

## redividing by denominator
df_wide_age <- as.data.table(df_wide_age)
df_wide_age <- df_wide_age[denom == "diabetes", `:=` (val = (sum_num/sum_dm), metric = "proportion")]
df_wide_age <- df_wide_age[denom == "diagnosed", `:=` (val = sum_num/sum_diag, metric = "proportion")]
df_wide_age <- df_wide_age[denom == "treated", `:=` (val = sum_num/sum_tx, metric = "proportion")]
df_wide_age <- df_wide_age[, dm_count := sum_dm]
df_wide_age <- df_wide_age[, diag_count := sum_diag]
df_wide_age <- df_wide_age[, tx_count := sum_tx]
df_wide_age <- subset(df_wide_age, select = -c(sum_num, sum_dm, sum_diag, sum_tx))

df_props <- rbind(df_props, df_wide_age, fill = TRUE)


# 6. aggregate 10 yr age groups -----------------------------------------------------------------------------
## aggregating
df_wide_age2 <- df_count[!age_group_id %in% c(22, 27, 197, 220, 154, 235), ]
df_wide_age2[age_group_id %in% c(8:9), age_group_id := 149]
df_wide_age2[age_group_id %in% c(10:11), age_group_id := 150]
df_wide_age2[age_group_id %in% c(12:13), age_group_id := 151]
df_wide_age2[age_group_id %in% c(14:15), age_group_id := 152]
df_wide_age2[age_group_id %in% c(16:17), age_group_id := 153]
df_wide_age2[age_group_id %in% c(18:19), age_group_id := 232]
df_wide_age2[age_group_id %in% c(20:30), age_group_id := 243]
df_wide_age2[age_group_id %in% c(31:32), age_group_id := 269]

df_wide_age2 <- group_by(df_wide_age2, year_id, sex_id, location_id, measure, metric, denom, draw, age_group_id)
df_wide_age2 <- summarize(df_wide_age2, sum_num = sum(val), sum_dm = sum(dm_count), sum_diag = sum(diag_count), sum_tx = sum(tx_count))

df_wide_age2_count <- as.data.table(df_wide_age2)
df_wide_age2_count <- df_wide_age2_count[, val := sum_num]
df_wide_age2_count <- df_wide_age2_count[, dm_count := sum_dm]
df_wide_age2_count <- df_wide_age2_count[, diag_count := sum_diag]
df_wide_age2_count <- df_wide_age2_count[, tx_count := sum_tx]
df_wide_age2_count <- subset(df_wide_age2_count, select = -c(sum_num, sum_dm, sum_diag, sum_tx))

df_count <- rbind(df_count, df_wide_age2_count, fill = TRUE)

## redividing by denominator
df_wide_age2 <- as.data.table(df_wide_age2)
df_wide_age2 <- df_wide_age2[denom == "diabetes", `:=` (val = (sum_num/sum_dm), metric = "proportion")]
df_wide_age2 <- df_wide_age2[denom == "diagnosed", `:=` (val = sum_num/sum_diag, metric = "proportion")]
df_wide_age2 <- df_wide_age2[denom == "treated", `:=` (val = sum_num/sum_tx, metric = "proportion")]
df_wide_age2 <- df_wide_age2[, dm_count := sum_dm]
df_wide_age2 <- df_wide_age2[, diag_count := sum_diag]
df_wide_age2 <- df_wide_age2[, tx_count := sum_tx]
df_wide_age2 <- subset(df_wide_age2, select = -c(sum_num, sum_dm, sum_diag, sum_tx))

df_props <- rbind(df_props, df_wide_age2, fill = TRUE)


# 7. sex ratio -----------------------------------------------------------------------------
## calculating sex ratio
male <- df_props[sex_id == 1, ]
female <- df_props[sex_id == 2, ]

male <- subset(male, select = -c(sex_id, dm_count, diag_count, tx_count))
female <- subset(female, select = -c(sex_id, dm_count, diag_count, tx_count))

sex_ratio <- merge(female, male, by = c("year_id", "location_id", "measure", "metric", "denom", "draw", "age_group_id"))
sex_ratio <- sex_ratio[, val := val.x/val.y]
sex_ratio <- subset(sex_ratio, select = -c(val.x, val.y))


# 8. absolute change over time -----------------------------------------------------------------------------
## calculating absolute change - proportion
df_start_year <- df_props[year_id == start_year, ]
df_end_year <- df_props[year_id == end_year, ]

df_start_year <- subset(df_start_year, select = -c(year_id, dm_count, diag_count, tx_count))
df_end_year <- subset(df_end_year, select = -c(year_id, dm_count, diag_count, tx_count))

abs_change_props <- merge(df_start_year, df_end_year, by = c("sex_id", "location_id", "measure", "metric", "denom", "draw", "age_group_id"))
abs_change_props <- abs_change_props[, val := val.y-val.x]
abs_change_props <- subset(abs_change_props, select = -c(val.x, val.y))

## calculating percent change - count
df_start_year <- df_count[year_id == start_year, ]
df_end_year <- df_count[year_id == end_year, ]

df_start_year <- subset(df_start_year, select = -c(year_id, dm_count, diag_count, tx_count))
df_end_year <- subset(df_end_year, select = -c(year_id, dm_count, diag_count, tx_count))

abs_change_count <- merge(df_start_year, df_end_year, by = c("sex_id", "location_id", "measure", "metric", "denom", "draw", "age_group_id"))
abs_change_count <- abs_change_count[, val := (val.y-val.x)/val.x]
abs_change_count <- subset(abs_change_count, select = -c(val.x, val.y))

abs_change <- rbind(abs_change_props, abs_change_count)


# summarizing draws ------------------------------------------------------------------------------------------
df_count2 <- subset(df_count, select = -c(dm_count, diag_count, tx_count))
df_props2 <- subset(df_props, select = -c(dm_count, diag_count, tx_count))

df_count2 <- pivot_wider(df_count2, names_from = "draw", values_from = "val")
df_props2 <- pivot_wider(df_props2, names_from = "draw", values_from = "val")
df_std2 <- pivot_wider(df_std, names_from = "draw", values_from = "val")
sex_ratio2 <- pivot_wider(sex_ratio, names_from = "draw", values_from = "val")
abs_change2 <- pivot_wider(abs_change, names_from = "draw", values_from = "val")

df_count_summary <- summaries(df_count2, draws)
df_props_summary <- summaries(df_props2, draws)
df_std_summary <- summaries(df_std2, draws)
sex_ratio_summary <- summaries(sex_ratio2, draws)
abs_change_summary <- summaries(abs_change2, draws)

## writing out 
write.csv(df_count_summary, paste0(dir, "count_summary.csv"))
write.csv(df_props_summary, paste0(dir, "props_summary.csv"))
write.csv(df_std_summary, paste0(dir, "agestd_summary.csv"))
write.csv(sex_ratio_summary, paste0(dir, "sex_ratio_summary.csv"))
write.csv(abs_change_summary, paste0(dir, "abs_change_summary.csv"))

write.csv(df_count, paste0(dir, "count_draws.csv"))
write.csv(df_props, paste0(dir, "props_draws.csv"))
write.csv(df_std, paste0(dir, "agestd_draws.csv"))
write.csv(sex_ratio, paste0(dir, "sex_ratio_draws.csv"))
write.csv(abs_change, paste0(dir, "abs_change_draws.csv"))


# plotting ------------------------------------------------------------------------------------------
outdir <- "FILEPATH"

## all age - proportion diabetes
temp <- df_props_summary

temp <- as.data.table(temp)
temp <- temp[measure %in% c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed') & denom == "diabetes" & age_group_id == 22 & !sex_id == 3, ]
temp <- merge(temp, loc_meta_merge, by = "location_id")

temp <- temp[sex_id == 1, sex := "Male"]
temp <- temp[sex_id == 2, sex := "Female"]

temp$measure <- factor(temp$measure, levels = c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed'))

temp1 <- temp[location_id %in% loc_meta[level == 3, ]$location_id, ]
pdf(paste0(outdir, "15plus_over_year.pdf"), width = 16, height = 9)
for (x in unique(temp1$region_name)) {
  temp2<-temp1[region_name==x,]
  region_name <- unique(temp2$region_name)
  plot <- ggplot(data = temp2, aes(x = year_id, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
    xlab("Year") + ylab("Proportion") + ggtitle(paste0("Region Name: ", region_name)) +
    facet_grid(vars(sex), vars(location_name)) + theme_bw() + theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Greens")
  print(plot)
}
dev.off()

pdf(paste0(outdir, "15plus_over_year_regions.pdf"), width = 16, height = 9)
temp2<-temp[location_id %in% loc_meta2[level == 2,]$location_id,]
plot <- ggplot(data = temp2, aes(x = year_id, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
  xlab("Year") + ylab("Proportion") +
  facet_grid(vars(sex), vars(location_name)) + theme_bw() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens")
print(plot)
dev.off()

pdf(paste0(outdir, "15plus_over_year_global_sr.pdf"), width = 16, height = 9)
temp2<-temp[location_id %in% loc_meta2[level %in% c(0,1),]$location_id,]
plot <- ggplot(data = temp2, aes(x = year_id, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
  xlab("Year") + ylab("Proportion") +
  facet_grid(vars(sex), vars(location_name)) + theme_bw() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens")
print(plot)
dev.off()


## age standardized - proportion diabetes
temp <- df_std_summary

temp <- as.data.table(temp)
temp <- temp[measure %in% c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed') & denom == "diabetes" & !sex_id == 3, ]
temp <- merge(temp, loc_meta_merge, by = "location_id")

temp <- temp[sex_id == 1, sex := "Male"]
temp <- temp[sex_id == 2, sex := "Female"]

temp$measure <- factor(temp$measure, levels = c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed'))

temp1 <- temp[location_id %in% loc_meta[level == 3, ]$location_id, ]
pdf(paste0(outdir, "15plus_agestd_over_year.pdf"), width = 16, height = 9)
for (x in unique(temp1$region_name)) {
  temp2<-temp1[region_name==x,]
  region_name <- unique(temp2$region_name)
  plot <- ggplot(data = temp2, aes(x = year_id, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
    xlab("Year") + ylab("Proportion") + ggtitle(paste0("Region Name: ", region_name)) +
    facet_grid(vars(sex), vars(location_name)) + theme_bw() + theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Greens")
  print(plot)
}
dev.off()

pdf(paste0(outdir, "15plus_agestd_over_year_regions.pdf"), width = 16, height = 9)
temp2<-temp[location_id %in% loc_meta2[level == 2,]$location_id,]
plot <- ggplot(data = temp2, aes(x = year_id, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
  xlab("Year") + ylab("Proportion") +
  facet_grid(vars(sex), vars(location_name)) + theme_bw() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens")
print(plot)
dev.off()

pdf(paste0(outdir, "15plus_agestd_over_year_global_sr.pdf"), width = 16, height = 9)
temp2<-temp[location_id %in% loc_meta2[level %in% c(0,1),]$location_id,]
plot <- ggplot(data = temp2, aes(x = year_id, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
  xlab("Year") + ylab("Proportion") +
  facet_grid(vars(sex), vars(location_name)) + theme_bw() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens")
print(plot)
dev.off()

