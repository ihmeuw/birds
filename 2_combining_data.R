# Combining microdata and tabs for upload

rm(list = ls())
library(data.table)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

invisible(sapply(list.files("FILEPATH", full.names = T), source))

# setting objects 
release_id <- 16

# metadata
loc_meta <- get_location_metadata(location_set_id = 35, release_id = release_id)
loc_merge <- subset(loc_meta, select = c(location_id, location_name))

#------------------------
# pulling in data 
df_collapsed <- as.data.table(read.csv("FILEPATH/aggregated_microdata.csv"))

lit <- as.data.table(read.csv("FILEPATH/extractions.csv"))
lit <- lit[!is.na(nid), ]
lit$type <- "tabs"

# tabs manipulations
lit <- lit[measure=="alt_diagnosed"&is.na(mean), mean := cases/sample_size]
lit <- lit[measure=="alt_diagnosed"&is.na(sample_size), sample_size := cases/mean]

lit <- lit %>%
  mutate(mean = ifelse(measure=="alt_diagnosed", 1-mean, mean),
         cases = ifelse(measure=="alt_diagnosed", NA, cases),
         upper = ifelse(measure=="alt_diagnosed", NA, upper),
         lower = ifelse(measure=="alt_diagnosed", NA, lower),
         measure = ifelse(measure=="alt_diagnosed", "undiagnosed", measure))

lit[is.na(upper)&is.na(lower), uncertainty_type_value := NA]

lit <- lit[sample_size <= 10, is_outlier := 1]
lit <- lit[is.na(is_outlier), is_outlier := 0]

# combining 
df_combined <- rbind(df_collapsed, lit, fill = TRUE)


# writing combined file
write.csv(df_combined, "FILEPATH/capstone_combined_data.csv", row.names = FALSE)

