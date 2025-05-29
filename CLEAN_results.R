rm(list = ls())

invisible(sapply(list.files("FILEPATH", full.names = TRUE), source))
library(dplyr)
library(writexl)
library(purrr)
library(general.utilities, lib.loc = "FILEPATH")

age_type <- c("age_groups", "age_standardized", "all_age")

for (a in age_type) {
  output_file <- ifelse(a == "age_groups",
    paste0("FILEPATH", Sys.Date(), "/age_group/"),
    ifelse(a == "all_age", paste0("FILEPATH", Sys.Date(), "/all_age/"),
      paste0("FILEPATH", Sys.Date(), "/age_std/")
    )
  )

  dir.create(output_file, recursive = TRUE)

  como_draw_version_id <- NUMBER

  gbd_id_type <- "sequela_id"
  source <- "como"

  location_set_id <- 35
  age_group_set_id <- 24

  year_id <- c(1990, 2023)
  sex_id <- c(1, 2, 3)
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

  if (a == "age_groups") {
    age_group_id <- unique(ages$age_group_id)
  } else if (a == "all_age") {
    age_group_id <- 22
  } else {
    age_group_id <- 27
  }

  location_id <- c(1, super_regions, regions, country)

  df_draws <- get_draws(
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

  df_seqs <- left_join(df_draws, seqs)

  df_seqs$rrt <- ifelse(df_seqs$sequela_name %like% "transplant", "Transplant", "Dialysis")

  df_seqs$etio <- ifelse(df_seqs$sequela_name %like% "type 1", "Diabetes Mellitus Type 1",
    ifelse(df_seqs$sequela_name %like% "type 2", "Diabetes Mellitus Type 2",
      ifelse(df_seqs$sequela_name %like% "hypertension", "Hypertension",
        ifelse(df_seqs$sequela_name %like% "glomerulonephritis", "Glomerulonephritis", "Other and unspecified causes")
      )
    )
  )

  df_seqs2 <- reshape_draws(df_seqs)

  pop <- get_population(
    age_group_id = unique(df_seqs2$age_group_id), location_id = unique(df_seqs2$location_id),
    year_id = unique(df_seqs2$year_id), sex_id = unique(df_seqs2$sex_id), with_ui = TRUE,
    release_id = release_id
  ) %>% rename("lower_pop" = "lower", "upper_pop" = "upper")

  if (a == "age_standardized") {
    pop$age_group_id <- 27
  }

  df_seqs3 <- left_join(df_seqs2, pop)

  df_seqs3$sequela_id <- NULL
  df_seqs3$sequela_name <- NULL

  df_seqs3 <- merge(x = df_seqs3, y = locs[, c("location_id", "location_name", "level", "lancet_label")], by = "location_id")

  df_seqs3 <- merge(x = df_seqs3, y = ages[, c("age_group_id", "age_group_name")], by = "age_group_id", all.x = TRUE)

  df_seqs3$age_group_name <- ifelse(df_seqs3$age_group_id == 22, "All Ages",
    ifelse(df_seqs3$age_group_id == 27, "Age-standardized", df_seqs3$age_group_name)
  )

  eskd <- df_seqs3 %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, level.y, sex_id, measure_id, draw) %>%
    reframe(
      cases = sum(values * population),
      cases_lower = sum(values * lower_pop),
      cases_upper = sum(values * upper_pop)
    )
  eskd <- eskd %>%
    na.omit(cases) %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, level.y, sex_id) %>%
    reframe(
      cases = signif(mean(cases), digits = 3),
      cases_lower = signif(quantile(cases_lower, 0.025), digits = 3),
      cases_upper = signif(quantile(cases_upper, 0.975), digits = 3)
    )

  rrt <- df_seqs3 %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, sex_id, rrt, measure_id, level.y, draw) %>%
    reframe(
      cases = sum(values * population),
      cases_lower = sum(values * lower_pop),
      cases_upper = sum(values * upper_pop)
    )

  rrt <- rrt %>%
    na.omit(cases) %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, sex_id, rrt, measure_id, level.y) %>%
    reframe(
      cases = signif(mean(cases), digits = 3),
      cases_lower = signif(quantile(cases_lower, 0.025), digits = 3),
      cases_upper = signif(quantile(cases_upper, 0.975), digits = 3)
    )

  etio <- df_seqs3 %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, sex_id, etio, measure_id, level.y, draw) %>%
    reframe(
      cases = sum(values * population),
      cases_lower = sum(values * lower_pop),
      cases_upper = sum(values * upper_pop)
    )

  etio <- etio %>%
    na.omit(cases) %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, sex_id, measure_id, etio, level.y) %>%
    reframe(
      cases = signif(mean(cases), digits = 3),
      cases_lower = signif(quantile(cases_lower, 0.025), digits = 3),
      cases_upper = signif(quantile(cases_upper, 0.975), digits = 3)
    )

  eskd_rate <- df_seqs3 %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id, level.y, draw) %>%
    reframe(
      rate = sum(values),
      rate_lower = sum(values),
      rate_upper = sum(values)
    )

  eskd_rate <- eskd_rate %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id, level.y) %>%
    reframe(
      rate = signif(mean(rate) * 100000, digits = 3),
      rate_lower = signif(quantile(rate_lower, 0.025) * 100000, digits = 3),
      rate_upper = signif(quantile(rate_upper, 0.975) * 100000, digits = 3)
    )

  etio_rate <- df_seqs3 %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id, etio, level.y, draw) %>%
    reframe(
      rate = sum(values),
      rate_lower = sum(values),
      rate_upper = sum(values)
    )

  etio_rate <- etio_rate %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id, level.y, etio) %>%
    reframe(
      rate = signif(mean(rate) * 100000, digits = 3),
      rate_lower = signif(quantile(rate_lower, 0.025) * 100000, digits = 3),
      rate_upper = signif(quantile(rate_upper, 0.975) * 100000, digits = 3)
    )

  rrt_rate <- df_seqs3 %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id, rrt, level.y, draw) %>%
    reframe(
      rate = sum(values),
      rate_lower = sum(values),
      rate_upper = sum(values)
    )

  rrt_rate <- rrt_rate %>%
    group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id, level.y, rrt) %>%
    reframe(
      rate = signif(mean(rate) * 100000, digits = 3),
      rate_lower = signif(quantile(rate_lower, 0.025) * 100000, digits = 3),
      rate_upper = signif(quantile(rate_upper, 0.975) * 100000, digits = 3)
    )

  write.csv(df_seqs3, paste0(output_file, "draws.csv"), row.names = FALSE)

  write.csv(eskd_rate, paste0(output_file, "eskd_rate.csv"), row.names = FALSE)
  write.csv(rrt_rate, paste0(output_file, "rrt_rate.csv"), row.names = FALSE)
  write.csv(etio_rate, paste0(output_file, "etio_rate.csv"), row.names = FALSE)

  write.csv(eskd, paste0(output_file, "eskd_count.csv"), row.names = FALSE)
  write.csv(rrt, paste0(output_file, "rrt_count.csv"), row.names = FALSE)
  write.csv(etio, paste0(output_file, "etio_count.csv"), row.names = FALSE)
}
