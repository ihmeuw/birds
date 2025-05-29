rm(list = ls())

library(xlsx)
library(dplyr)
invisible(sapply(list.files("FILEPATH", full.names = T), source))
library(general.utilities, lib.loc = "FILEPATH")

gbd_id_type <- "sequela_id"
source <- "como"

como_draw_version_id <- NUMBER

location_set_id <- NUMBER
age_group_set_id <- NUMBER

year_id <- c(1990, 2023)
sex_id <- c(3)
metric_id <- c(3)
measure_id <- c(5)
release_id <- 16

seqs <- get_sequela_metadata(sequela_set_id = 2, release_id = release_id)
esrd <- seqs %>% filter(sequela_name %like% "End-stage renal")
sq_ids <- esrd$sequela_id

locs <- get_location_metadata(location_set_id = location_set_id, release_id = release_id)
super_regions <- locs[level == 1, location_id]
regions <- locs[level == 2, location_id]
country <- locs[level == 3, location_id]

ages <- get_age_metadata(age_group_set_id = age_group_set_id, release_id = release_id)

age_group_id <- c(22, 27, unique(ages$age_group_id))

location_id <- c(1, super_regions, regions, country)

etio <- get_draws(
  gbd_id_type = gbd_id_type,
  gbd_id = sq_ids,
  source = source,
  measure_id = measure_id,
  location_id = location_id,
  year_id = year_id,
  age_group_id = age_group_id,
  sex_id = sex_id,
  metric_id = metric_id,
  release_id = release_id,
  version_id = como_draw_version_id
)

etio2 <- left_join(etio, seqs)

etio2$etiology <- ifelse(etio2$sequela_name %like% "type 1", "Diabetes Mellitus Type 1",
                         ifelse(etio2$sequela_name %like% "type 2", "Diabetes Mellitus Type 2",
                                ifelse(etio2$sequela_name %like% "hypertension", "Hypertension",
                                       ifelse(etio2$sequela_name %like% "glomerulonephritis", "Glomerulonephritis", "Other and unspecified causes")
                                )
                         )
)

etio3 <- etio2 %>% select(-c("level", "most_detailed", "parent_id", "path_to_top_parent", "sort_order", "modelable_entity_id", "cause_id", "cause_name", "rei_id", "rei_name", "healthstate_id", "healthstate_name", "version_id", "sequela_set_version_id", "sequela_set_id", "sequela_set_name"))

etio3 <- reshape_draws(etio3)

pop <- get_population(
  age_group_id = unique(etio3$age_group_id), location_id = unique(etio3$location_id),
  year_id = unique(etio3$year_id), sex_id = unique(etio3$sex_id), with_ui = TRUE,
  release_id = release_id
) %>% rename("lower_pop" = "lower", "upper_pop" = "upper")

etio4 <- left_join(etio3, pop)

etio4$sequela_id <- NULL
etio4$sequela_name <- NULL

etio5 <- etio4 %>%
  group_by(year_id, age_group_id, location_id, sex_id, etiology, draw) %>%
  reframe(
    num = sum(values * population),
    num_lower = sum(values * lower_pop),
    num_upper = sum(values * upper_pop)
  )

etio6 <- etio5 %>% group_by (year_id, age_group_id, location_id, sex_id,draw) %>% 
  mutate(
    num_total = sum(num),
    num_lower_total = sum(num_lower),
    num_upper_total = sum(num_upper)
  ) 

etio7 <- etio6 %>% group_by (year_id, age_group_id, location_id, sex_id,draw) %>% 
  mutate(prop = num/num_total,
         prop_lower = num_lower/num_lower_total,
         prop_upper = num_upper/num_upper_total)

etio8 <- etio7 %>% group_by (year_id, age_group_id, location_id,etiology, sex_id) %>% 
  na.omit() %>% 
  reframe(proportion = (mean(prop)*100),
          proportion_lower = (quantile(prop, 0.025)*100),
          proportion_upper = (quantile(prop, 0.975)*100))

final <- merge(x = etio8, y = locs)
final < - merge(final, ages)

final <- final %>% mutate(proportion = round(proportion, digits = 1),
                          proportion_lower = round(proportion_lower,digits = 1),
                          proportion_upper = round(proportion_upper, digits = 1))

write.csv(final, paste0( "FILEPATH"), row.names = FALSE)
