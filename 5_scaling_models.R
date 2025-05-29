# aggregating to region/super region/global & income groups, scaling models to sum to 100%
# scale cont&uncont to alt_tx model, and undiag&diag_untx to 1-alt_tx model

library(dplyr)
library(tidyr)
library(ggplot2)

invisible(sapply(list.files("FILEPATH", full.names = T), source))

num_draws <- 250
draws <- paste0("draw_", 0:249)

years <- c(1990, 1995, 2000, 2005, 2010, 2015, 2020, 2023)

plot_year <- 2023

dir <- "FILEPATH"

## post processing codebook - update model version IDs here first!
cb <- as.data.table(read.csv(paste0(dir, "inputs/codebooks/post_processing_codebook.csv")))

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
loc_meta_merge <- subset(loc_meta, select = c(location_id, region_id))
loc_meta_merge2 <- subset(loc_meta, select = c(location_id, super_region_id))

loc_income <- as.data.table(get_location_metadata(location_set_id = 26, release_id = 16))
loc_merge_level <- subset(loc_meta, select = c(location_id, level, parent_id))
loc_income_agg <- loc_income[most_detailed == 1, ]
loc_income_agg <- subset(loc_income_agg, select = c(location_id, location_name, region_id, region_name))
loc_income_agg <- left_join(loc_income_agg, loc_merge_level, by = "location_id")
loc_income_agg <- as.data.table(loc_income_agg)
loc_income_agg <- loc_income_agg[level == 4, location_id := parent_id]
loc_income_agg <- loc_income_agg[level == 5 & region_id == 44575, location_id := 95]
loc_income_agg <- loc_income_agg[level == 5 & region_id == 44577, location_id := 163]
loc_income_agg <- subset(loc_income_agg, select = c(location_id, region_id))
loc_income_agg <- unique(loc_income_agg)

## regional scalars
scalar_df <- get_regional_scalars(release_id = 16)
scalar_df <- subset(scalar_df, select = c(location_id, year_id, mean))

## reading in draws
undiagnosed <- as.data.table(get_draws(gbd_id_type = "modelable_entity_id", gbd_id = cb[measure_upload=="undiagnosed",]$me_id, source = "epi", 
                                       version_id = cb[measure_upload=="undiagnosed",]$model_version_id, release_id = 16, 
                                       age_group_id = c(8:20, 30:32, 235), sex_id = c(1,2,3), location_id = loc_meta[level ==3,]$location_id,
                                       downsample = TRUE, n_draws = num_draws,
                                       year_id = years))
diagnosed_untx <- as.data.table(get_draws(gbd_id_type = "modelable_entity_id", gbd_id = cb[measure_upload=="diagnosed_untx",]$me_id, source = "epi", 
                                          version_id = cb[measure_upload=="diagnosed_untx",]$model_version_id, release_id = 16, 
                                          age_group_id = c(8:20, 30:32, 235), sex_id = c(1,2,3), location_id = loc_meta[level ==3,]$location_id,
                                          downsample = TRUE, n_draws = num_draws,
                                          year_id = years))
diagnosed_tx_uncont <- as.data.table(get_draws(gbd_id_type = "modelable_entity_id", gbd_id = cb[measure_upload=="diagnosed_tx_uncont",]$me_id, source = "epi", 
                                               version_id = cb[measure_upload=="diagnosed_tx_uncont",]$model_version_id, release_id = 16, 
                                               age_group_id = c(8:20, 30:32, 235), sex_id = c(1,2,3), location_id = loc_meta[level ==3,]$location_id,
                                               downsample = TRUE, n_draws = num_draws,
                                               year_id = years))
diagnosed_tx_cont <- as.data.table(get_draws(gbd_id_type = "modelable_entity_id", gbd_id = cb[measure_upload=="diagnosed_tx_cont",]$me_id, source = "epi", 
                                             version_id = cb[measure_upload=="diagnosed_tx_cont",]$model_version_id, release_id = 16, 
                                             age_group_id = c(8:20, 30:32, 235), sex_id = c(1,2,3), location_id = loc_meta[level ==3,]$location_id,
                                             downsample = TRUE, n_draws = num_draws,
                                             year_id = years))
alt_tx <- as.data.table(get_draws(gbd_id_type = "modelable_entity_id", gbd_id = cb[measure_upload=="alt_tx",]$me_id, source = "epi", 
                                  version_id = cb[measure_upload=="alt_tx",]$model_version_id, release_id = 16, 
                                  age_group_id = c(8:20, 30:32, 235), sex_id = c(1,2,3), location_id = loc_meta[level ==3,]$location_id,
                                  downsample = TRUE, n_draws = num_draws,
                                  year_id = years))

undiagnosed <- undiagnosed[, measure := "undiagnosed"]
diagnosed_untx <- diagnosed_untx[, measure := "diagnosed_untx"]
diagnosed_tx_uncont <- diagnosed_tx_uncont[, measure := "diagnosed_tx_uncont"]
diagnosed_tx_cont <- diagnosed_tx_cont[, measure := "diagnosed_tx_cont"]
alt_tx <- alt_tx[, measure := "tx"]

list <- list(undiagnosed, diagnosed_untx, diagnosed_tx_uncont, diagnosed_tx_cont, alt_tx)
df <- rbindlist(list)

## pulling dm results and converting to counts
dm <- as.data.table(get_draws(gbd_id_type = "cause_id", gbd_id = NUM, source = "como", 
                              release_id = 16, 
                              age_group_id = c(8:20, 30:32, 235), 
                              sex_id = c(1,2,3), location_id = loc_meta[level %in% c(0:3),]$location_id,
                              metric_id = 3, measure_id = 5, downsample = TRUE, n_draws = num_draws, 
                              year_id = years, version_id = 1591)) 

dm <- pivot_longer(dm, cols = draws, names_to = "draw", values_to = "val")

pop <- get_population(release_id = 16, year_id= years, location_id = loc_meta[level %in% c(0:3),]$location_id,
                      sex_id = c(1,2,3), age_group_id = c(8:20, 30:32, 235))

dm <- merge(dm, pop, all.x = TRUE, by = c("age_group_id", "year_id", "sex_id", "location_id"))
dm <- as.data.table(dm)

dm <- dm[, dm_count := val*population]

dm <- subset(dm, select = c(age_group_id, year_id, sex_id, location_id, draw, dm_count))


# 1. GBD regions, super regions, global -----------------------------------------------------------------------------
## multiplying by dm prevalence counts
df2 <- pivot_longer(df, cols = draws, names_to = "draw", values_to = "val")

df2 <- merge(df2, dm, all.x=TRUE, by = c("age_group_id", "year_id", "sex_id", "location_id", "draw"))
df2 <- as.data.table(df2)

df_region <- copy(df2)
df_region[, val := val*dm_count]
df_region <- subset(df_region, select = -c(dm_count))

## aggregating to regions + multiplying by scalars
df_region <- merge(df_region, loc_meta_merge, by = "location_id")
df_region <- group_by(df_region, year_id, sex_id, region_id, age_group_id, measure, draw)
df_region <- summarize(df_region, val = sum(val))

df_region <- as.data.table(df_region)
df_region[, location_id := region_id]
df_region <- merge(df_region, scalar_df, by = c("location_id", "year_id"), all.x=TRUE)
df_region[, val := val*mean]
df_region <- subset(df_region, select = -c(mean, region_id))

## aggregating to super regions
df_superregion <- merge(df_region, loc_meta_merge2, by = "location_id")
df_superregion <- group_by(df_superregion, year_id, sex_id, super_region_id, age_group_id, measure, draw)
df_superregion <- summarize(df_superregion, val = sum(val))

df_superregion <- as.data.table(df_superregion)
df_superregion[, location_id := super_region_id]
df_superregion <- subset(df_superregion, select = -c(super_region_id))

## aggregating to global
df_global <- group_by(df_superregion, year_id, sex_id, age_group_id, measure, draw)
df_global <- summarize(df_global, val = sum(val))

df_global <- as.data.table(df_global)
df_global[, location_id := 1]

## combining and dividing by dm pop for proportion 
list <- list(df_region, df_superregion, df_global)
df_combine <- rbindlist(list, use.names = TRUE)

df_combine <- merge(df_combine, dm, all.x=TRUE, by = c("age_group_id", "year_id", "sex_id", "location_id", "draw"))
df_combine <- as.data.table(df_combine)

df_combine[, val := val/dm_count]
df_combine <- subset(df_combine, select = -c(dm_count))


# 2. World Bank income groups -----------------------------------------------------------------------------
## multiplying by dm prevalence counts
df_income <- copy(df2)
df_income[, val := val*dm_count]

## aggregating to income groups
df_income <- merge(df_income, loc_income_agg, by = "location_id")
df_income <- group_by(df_income, year_id, sex_id, region_id, age_group_id, measure, draw)
df_income <- summarize(df_income, val = sum(val), dm_count = sum(dm_count))

df_income <- as.data.table(df_income)
df_income[, location_id := region_id]
df_income <- subset(df_income, select = -c(region_id))

## separating out dm counts to save
df_income_dm <- copy(df_income)
df_income_dm <- subset(df_income_dm, select = c(year_id, sex_id, location_id, age_group_id, draw, dm_count))
df_income_dm <- unique(df_income_dm)
dm <- rbind(df_income_dm, dm)
write.csv(dm, paste0(dir, "outputs/dm_count_draws.csv"))

## dividing by dm pop for proportion 
df_income[, val := val/dm_count]
df_income <- subset(df_income, select = -c(dm_count))


# combining and writing out pre scaled results -----------------------------------------------------------------------------
df_combine <- pivot_wider(df_combine, names_from = "draw", values_from = "val")
df_income <- pivot_wider(df_income, names_from = "draw", values_from = "val")

df <- subset(df, select = -c(measure_id, metric_id, model_version_id, modelable_entity_id))

list <- list(df_combine, df_income, df)
df_write <- rbindlist(list, use.names = TRUE)

write.csv(df_write, paste0(dir, "outputs/loc_agg_pre_scaled_draws.csv"))

summary <- summaries(df_write, draws)

write.csv(summary, paste0(dir, "outputs/loc_agg_pre_scaled.csv"))


# scaling to sum to 100% -----------------------------------------------------------------------------
## split out
undiagnosed <- df_write[measure == "undiagnosed", ]
diagnosed_untx <- df_write[measure == "diagnosed_untx",]
diagnosed_tx_uncont <- df_write[measure == "diagnosed_tx_uncont",]
diagnosed_tx_cont <- df_write[measure == "diagnosed_tx_cont",]
alt_tx <- df_write[measure == "tx", ]

undiagnosed <- subset(undiagnosed, select = -c(measure))
diagnosed_untx <- subset(diagnosed_untx, select = -c(measure))
diagnosed_tx_uncont <- subset(diagnosed_tx_uncont, select = -c(measure))
diagnosed_tx_cont <- subset(diagnosed_tx_cont, select = -c(measure))
alt_tx <- subset(alt_tx, select = -c(measure))

## merge
undiagnosed <- pivot_longer(undiagnosed, cols = draws, names_to = "draw", values_to = "val")
diagnosed_untx <- pivot_longer(diagnosed_untx, cols = draws, names_to = "draw", values_to = "val")
diagnosed_tx_uncont <- pivot_longer(diagnosed_tx_uncont, cols = draws, names_to = "draw", values_to = "val")
diagnosed_tx_cont <- pivot_longer(diagnosed_tx_cont, cols = draws, names_to = "draw", values_to = "val")
alt_tx <- pivot_longer(alt_tx, cols = draws, names_to = "draw", values_to = "val")

undiagnosed <- rename(undiagnosed, val_undiagnosed = val)
diagnosed_untx <- rename(diagnosed_untx, val_diagnosed_untx = val)
diagnosed_tx_uncont <- rename(diagnosed_tx_uncont, val_diagnosed_tx_uncont = val)
diagnosed_tx_cont <- rename(diagnosed_tx_cont, val_diagnosed_tx_cont = val)
alt_tx <- rename(alt_tx, val_alt_tx = val)

alt_tx <- as.data.table(alt_tx)
alt_untx <- alt_tx[, val_alt_untx := 1-val_alt_tx]
alt_untx <- subset(alt_untx, select = -c(val_alt_tx))

df_combined_tx <- merge(diagnosed_tx_uncont, diagnosed_tx_cont, by = c("age_group_id", "location_id", "sex_id", "year_id", "draw"))
df_combined_tx <- merge(df_combined_tx, alt_tx, by = c("age_group_id", "location_id", "sex_id", "year_id", "draw"))
df_combined_untx <- merge(undiagnosed, diagnosed_untx, by = c("age_group_id", "location_id", "sex_id", "year_id", "draw"))
df_combined_untx <- merge(df_combined_untx, alt_untx, by = c("age_group_id", "location_id", "sex_id", "year_id", "draw"))

## sum
df_combined_tx <- as.data.table(df_combined_tx)
df_combined_tx <- df_combined_tx[, val_sum := val_diagnosed_tx_uncont+val_diagnosed_tx_cont]
df_combined_tx <- df_combined_tx[, val_prop_diagnosed_tx_uncont := (val_diagnosed_tx_uncont/val_sum)*(val_alt_tx)]
df_combined_tx <- df_combined_tx[, val_prop_diagnosed_tx_cont := (val_diagnosed_tx_cont/val_sum)*(val_alt_tx)]

df_combined_untx <- as.data.table(df_combined_untx)
df_combined_untx <- df_combined_untx[, val_sum := val_undiagnosed+val_diagnosed_untx]
df_combined_untx <- df_combined_untx[, val_prop_undiagnosed := (val_undiagnosed/val_sum)*(val_alt_untx)]
df_combined_untx <- df_combined_untx[, val_prop_diagnosed_untx := (val_diagnosed_untx/val_sum)*(val_alt_untx)]

## split out
undiagnosed <- subset(df_combined_untx, select = c(age_group_id, location_id, sex_id, year_id, draw, val_prop_undiagnosed))
diagnosed_untx <- subset(df_combined_untx, select = c(age_group_id, location_id, sex_id, year_id, draw, val_prop_diagnosed_untx))
diagnosed_tx_uncont <- subset(df_combined_tx, select = c(age_group_id, location_id, sex_id, year_id, draw, val_prop_diagnosed_tx_uncont))
diagnosed_tx_cont <- subset(df_combined_tx, select = c(age_group_id, location_id, sex_id, year_id, draw, val_prop_diagnosed_tx_cont))

undiagnosed <- rename(undiagnosed, val = val_prop_undiagnosed)
diagnosed_untx <- rename(diagnosed_untx, val = val_prop_diagnosed_untx)
diagnosed_tx_uncont <- rename(diagnosed_tx_uncont, val = val_prop_diagnosed_tx_uncont)
diagnosed_tx_cont <- rename(diagnosed_tx_cont, val = val_prop_diagnosed_tx_cont)

undiagnosed <- pivot_wider(undiagnosed, names_from = "draw", values_from = "val")
diagnosed_untx <- pivot_wider(diagnosed_untx, names_from = "draw", values_from = "val")
diagnosed_tx_uncont <- pivot_wider(diagnosed_tx_uncont, names_from = "draw", values_from = "val")
diagnosed_tx_cont <- pivot_wider(diagnosed_tx_cont, names_from = "draw", values_from = "val")

## write out
write.csv(undiagnosed, paste0(dir, "outputs/undiagnosed_post_scaled_draws.csv"))
write.csv(diagnosed_untx, paste0(dir, "outputs/diagnosed_untx_post_scaled_draws.csv"))
write.csv(diagnosed_tx_uncont, paste0(dir, "outputs/diagnosed_tx_uncont_post_scaled_draws.csv"))
write.csv(diagnosed_tx_cont, paste0(dir, "outputs/diagnosed_tx_cont_post_scaled_draws.csv"))

undiagnosed_summary <- summaries(undiagnosed, draws)
diagnosed_untx_summary <- summaries(diagnosed_untx, draws)
diagnosed_tx_uncont_summary <- summaries(diagnosed_tx_uncont, draws)
diagnosed_tx_cont_summary <- summaries(diagnosed_tx_cont, draws)

write.csv(undiagnosed_summary, paste0(dir, "outputs/undiagnosed_post_scaled.csv"))
write.csv(diagnosed_untx_summary, paste0(dir, "outputs/diagnosed_untx_post_scaled.csv"))
write.csv(diagnosed_tx_uncont_summary, paste0(dir, "outputs/diagnosed_tx_uncont_post_scaled.csv"))
write.csv(diagnosed_tx_cont_summary, paste0(dir, "outputs/diagnosed_tx_cont_post_scaled.csv"))


## vetting scaling ---------------------------------------------------------------------------------
## metadata 
age_meta <- as.data.table(get_age_metadata(release_id = 16))

sex_meta <- get_ids("sex")

## post scale model results 
undiagnosed <- as.data.table(undiagnosed_summary)
diagnosed_untx <- as.data.table(diagnosed_untx_summary)
diagnosed_tx_uncont <- as.data.table(diagnosed_tx_uncont_summary)
diagnosed_tx_cont <- as.data.table(diagnosed_tx_cont_summary)

undiagnosed <- undiagnosed[, measure := "undiagnosed"]
diagnosed_untx <- diagnosed_untx[, measure := "diagnosed_untx"]
diagnosed_tx_uncont <- diagnosed_tx_uncont[, measure := "diagnosed_tx_uncont"]
diagnosed_tx_cont <- diagnosed_tx_cont[, measure := "diagnosed_tx_cont"]

undiagnosed <- subset(undiagnosed, select = c(location_id, year_id, age_group_id, sex_id, mean, measure))
diagnosed_untx <- subset(diagnosed_untx, select = c(location_id, year_id, age_group_id, sex_id, mean, measure))
diagnosed_tx_uncont <- subset(diagnosed_tx_uncont, select = c(location_id, year_id, age_group_id, sex_id, mean, measure))
diagnosed_tx_cont <- subset(diagnosed_tx_cont, select = c(location_id, year_id, age_group_id, sex_id, mean, measure))

df <- rbind(undiagnosed, diagnosed_untx)
df <- rbind(df, diagnosed_tx_uncont)
df <- rbind(df, diagnosed_tx_cont)

df <- merge(df, age_meta, by = "age_group_id", all.x=TRUE)
df <- merge(df, loc_meta, by = "location_id", all.x=TRUE)
df <- merge(df, sex_meta, by = "sex_id", all.x=TRUE)

df[, type := "Post-Scaled"]

## reading in pre scale model results 
df_prescale <- as.data.table(read.csv(paste0(dir, "outputs/loc_agg_pre_scaled.csv")))
df_prescale <- df_prescale[!measure == "tx", ]
df_prescale <- subset(df_prescale, select = c(location_id, year_id, age_group_id, sex_id, mean, measure))

df_prescale <- merge(df_prescale, age_meta, by = "age_group_id", all.x=TRUE)
df_prescale <- merge(df_prescale, loc_meta, by = "location_id", all.x=TRUE)
df_prescale <- merge(df_prescale, sex_meta, by = "sex_id", all.x=TRUE)

df_prescale[, type := "Pre-Scaled"]

## stacked bar plots for each location, by sex and year
temp <- df

temp$measure <- factor(temp$measure, levels = c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed'))

outdir <- "FILEPATH"
pdf(paste0(outdir, "model_results_bar_plots_post_scale.pdf"), width = 16, height = 9)
for (x in unique(temp$location_id)) {
  temp2<-temp[location_id==x,]
  location_name <- unique(temp2$location_name)
  plot <- ggplot(data = temp2, aes(x = age_group_name_short, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
    xlab("Age") + ylab("Proportion") + ggtitle(paste0("Location Name: ", location_name)) +
    facet_grid(vars(sex), vars(year_id)) + theme_bw() + theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Greens")
  print(plot)
}
dev.off()


# stacked bar plots, pre and post scaling, both sex over age 
temp <- rbind(df, df_prescale)
temp <- temp[sex_id == 3, ]

temp$measure <- factor(temp$measure, levels = c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed'))
temp$type <- factor(temp$type, levels = c("Pre-Scaled", "Post-Scaled"))

pdf(paste0(outdir, "model_results_bar_plots_pre_post_scale.pdf"), width = 16, height = 9)
for (x in unique(temp$location_id)) {
  temp2<-temp[location_id==x,]
  location_name <- unique(temp2$location_name)
  plot <- ggplot(data = temp2, aes(x = age_group_name_short, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
    xlab("Age") + ylab("Proportion") + ggtitle(paste0("Location Name: ", location_name)) +
    facet_grid(vars(type), vars(year_id)) + theme_bw() + theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Greens")
  print(plot)
}
dev.off()


# stacked bar plots for final reporting year, both sex
## super region
temp <- df[location_id %in% loc_meta[level == 1, ]$location_id & year_id == plot_year & sex_id == 3, ]

temp$measure <- factor(temp$measure, levels = c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed'))

pdf(paste0(outdir, "model_results_bar_plots_sr_post_scale.pdf"), width = 16, height = 9)
plot <- ggplot(data = temp, aes(x = age_group_name_short, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
  xlab("Age") + ylab("Proportion") + ggtitle(paste0("Both Sex, ", plot_year)) +
  facet_wrap(~location_name) + theme_bw() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens")
print(plot)
dev.off()

## region 
temp <- df[location_id %in% loc_meta[level == 2, ]$location_id & year_id == plot_year & sex_id == 3, ]

temp$measure <- factor(temp$measure, levels = c('diagnosed_tx_cont', 'diagnosed_tx_uncont', 'diagnosed_untx', 'undiagnosed'))

pdf(paste0(outdir, "model_results_bar_plots_region_post_scale.pdf"), width = 16, height = 9)
plot <- ggplot(data = temp, aes(x = age_group_name_short, y = mean, fill = measure)) + geom_bar(position="stack", stat="identity") +
  xlab("Age") + ylab("Proportion") + ggtitle(paste0("Both Sex, ", plot_year)) +
  facet_wrap(~location_name) + theme_bw() + theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Greens")
print(plot)
dev.off()


