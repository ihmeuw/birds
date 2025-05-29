rm(list = ls())
invisible(sapply(list.files("FILEPATH", full.names = TRUE), source))

release <- NUMBER

output_file <- "FILEPATH"

df <- read.csv("FILEPATH")
df2 <- read.csv("FILEPATH")

combo <- rbind(df, df2)

etio_df <- combo %>%
  group_by(year_id, age_group_id, age_group_name, location_id, location_name, lancet_label, measure_id, sex_id,etio, draw) %>%
  reframe(
    rate = sum(values),
    rate_lower = sum(values),
    rate_upper = sum(values)
  ) %>%
  distinct() %>%
  filter(year_id %in% c(1990, 2023)) %>%
  pivot_wider(
    names_from = year_id, values_from = c(rate, rate_lower, rate_upper),
    names_glue = "{.value}_{year_id}"
  ) %>%
  mutate(
    percentage_change_rate = ifelse(!is.na(rate_1990) & !is.na(rate_2023),
                                    ((rate_2023 - rate_1990) / rate_1990) * 100,
                                    NA
    ),
    percentage_change_lower_rate = ifelse(!is.na(rate_lower_1990) & !is.na(rate_lower_2023),
                                          ((rate_lower_2023 - rate_lower_1990) / rate_lower_1990) * 100,
                                          NA
    ),
    percentage_change_upper_rate = ifelse(!is.na(rate_upper_1990) & !is.na(rate_upper_2023),
                                          ((rate_upper_2023 - rate_upper_1990) / rate_upper_1990) * 100,
                                          NA
    )
  ) %>%
  drop_na(rate_1990, rate_2023)

etio_df_final <- etio_df %>%
  group_by(location_id, location_name, lancet_label, measure_id, age_group_name, age_group_id, sex_id, etio) %>%
  reframe(
    rate_1990 = signif(mean(rate_1990) * 100000, digits = 3),
    rate_lower_1990 = signif(quantile(rate_lower_1990, 0.025) * 100000, digits = 3),
    rate_upper_1990 = signif(quantile(rate_upper_1990, 0.975) * 100000, digits = 3),
    rate_2023 = signif(mean(rate_2023) * 100000, digits = 3),
    rate_lower_2023 = signif(quantile(rate_lower_2023, 0.025) * 100000, digits = 3),
    rate_upper_2023 = signif(quantile(rate_upper_2023, 0.975) * 100000, digits = 3),
    percentage_change_rate = signif(mean(percentage_change_rate), digits = 3),
    percentage_change_lower_rate = signif(quantile(percentage_change_lower_rate, 0.025), digits = 3),
    percentage_change_upper_rate = signif(quantile(percentage_change_upper_rate, 0.975), digits = 3)
  )


locs <- get_location_metadata(location_set_id = 35, release_id = release)

all <- left_join(etio_df_final, locs)

fwrite(all, "FILEPATH")