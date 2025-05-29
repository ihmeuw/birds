# Subsetting and aggregating microdata

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

# pulling in cleaned microdata file
df <- as.data.table(read.csv("FILEPATH/cleaned_individual_microdata_gbd2023_2024_09_17.csv"))

# correct standardized_tx for a few nids
df <- df %>%
  mutate(standardized_tx = ifelse(nid == 108841, "survey did not ask about treatment status", standardized_tx)) %>%
  mutate(standardized_tx = ifelse((nid == 112185 | nid == 140031) & self_report_diab_ever=="Told on med", "on treatment", standardized_tx)) %>%
  mutate(standardized_tx = ifelse((nid == 112185 | nid == 140031 | nid==65079) & self_report_diab_ever=="Told not on med", "not on treatment", standardized_tx))

# map all Health Survey for England data sources to England (sample size issues)
df[study_name %like% "Health Survey for England", location_id := 4749]

# assumptions: if negative self report, switch 'unknown' treatment status to 'not on treatment' : likely a skip pattern
# if they are on treatment, should have a positive self report
df <- df %>%
  mutate(standardized_tx= ifelse((standardized_self_report=="negative self report" | standardized_self_report =="positive self report") & standardized_tx=="unknown", "not on treatment", standardized_tx), 
         standardized_self_report = ifelse(standardized_tx=="on treatment" & (standardized_self_report=="negative self report" | standardized_self_report=="unknown"),
                                           "positive self report", standardized_self_report))

df <- df[standardized_case_definition %like% "_tx"] # treatment status needed in order to accurately define people with diabetes
df <- df[!nid == 110018, ] # British virgin islands, not GBD loc 
df <- df[!nid == 22306, ] # Health Survey for England 1994 (extremely high undiagnosed, very low treated)
df <- df[!nid == 327780, ] # too high undiagnosed proportion (75% undiagnosed in Norway)
df <- df[!nid == 112545, ] # too low undiagnosed proportion (7% undiagnosed in Jamaica)
df <- df[!nid == 107337, ] # SS is small and changing age binning has big impacts on estimates - control too high (Togo)

df <- merge(df, loc_merge, by = "location_id", all.x=TRUE)

# dropping ages below 15
df<-df[!age_start<15,]

# assigning sample weight a 1 if missing for nids with no sample weights, dropping missing for surveys with sample weights
missing_weight <- unique(df[is.na(sample_weight) | sample_weight == "",]$nid)
no_missing_weight <- unique(df[!(is.na(sample_weight) | sample_weight == ""),]$nid)
missing_weight2 <- missing_weight[(missing_weight %in% no_missing_weight)] #these nids have NA or "" rows but other rows have sample weight values

df <- df[!(nid %in% missing_weight2 & (is.na(sample_weight) | sample_weight == "")), ]
df <- df[is.na(sample_weight) | sample_weight == "", sample_weight := 1]

df <- df[!sample_weight == 0,]

# updating to 10 year age bins
df$age_start[df$age>=15]<-15
df$age_start[df$age>=25]<-25
df$age_start[df$age>=35]<-35
df$age_start[df$age>=45]<-45
df$age_start[df$age>=55]<-55
df$age_start[df$age>=65]<-65
df$age_start[df$age>=75]<-75
df$age_start[df$age>=85]<-85
df <- as.data.table(df)
df <- df[, age_end := age_start + 9]
df <- df[age_start == 85, age_end := 99]

# updating to 20 year age bins for sources where 20% or more data points are outliered due to small sample size
# skip this chunk to run collapse code first, then check percentage outliers and re-run 
df$age_start[(df$nid %in% check_age_bins_nids) & df$age>=15]<-15
df$age_start[(df$nid %in% check_age_bins_nids) & df$age>=35]<-35
df$age_start[(df$nid %in% check_age_bins_nids) & df$age>=55]<-55
df$age_start[(df$nid %in% check_age_bins_nids) & df$age>=75]<-75
df <- as.data.table(df)
df <- df[nid %in% check_age_bins_nids, age_end := age_start + 19]
df <- df[nid %in% check_age_bins_nids & age_start == 75, age_end := 99]


# collapsing - fpg
temp_fpg <- copy(df)

temp_fpg <- temp_fpg[!is.na(fpg_mmol),]
temp_fpg <- temp_fpg[fpg_mmol<2.5|fpg_mmol>25, is_outlier := 1]
temp_fpg <- temp_fpg[is.na(is_outlier), is_outlier := 0]

temp_fpg <- temp_fpg %>%
  mutate(dm = case_when((fpg_mmol>=7 | standardized_tx == "on treatment") ~ 1, 
                        TRUE ~ 0), 
         undiagnosed = case_when(standardized_self_report=="negative self report" | standardized_self_report=="unknown" ~ 1, #assumes that 'unknown' self report is undiagnosed
                                 standardized_self_report=="positive self report" ~ 0, 
                                 standardized_self_report=="survey did not ask about diagnosis status" ~ NA), 
         diagnosed_untx = case_when(standardized_self_report=="positive self report" & (standardized_tx == "not on treatment" | standardized_tx =="unknown") ~ 1, #assumes that 'unknown' tx is not treated
                                    standardized_self_report=="survey did not ask about diagnosis status" ~ NA, 
                                    TRUE ~ 0), 
         alt_tx = case_when(standardized_tx == "on treatment" ~ 1, 
                            TRUE ~ 0), 
         diagnosed_tx_fpg_cont_5.6 = case_when(standardized_tx == "on treatment" & fpg_mmol<5.6 ~ 1, 
                                               TRUE ~ 0), 
         diagnosed_tx_fpg_cont_7 = case_when(standardized_tx == "on treatment" & fpg_mmol<7 ~ 1, 
                                             TRUE ~ 0),
         diagnosed_tx_fpg_cont_7.2 = case_when(standardized_tx == "on treatment" & fpg_mmol<7.2 ~ 1, 
                                             TRUE ~ 0),
         diagnosed_tx_fpg_cont_7.22 = case_when(standardized_tx == "on treatment" & fpg_mmol<7.22 ~ 1, 
                                               TRUE ~ 0),
         diagnosed_tx_fpg_cont_8.6 = case_when(standardized_tx == "on treatment" & fpg_mmol<8.6 ~ 1, 
                                                TRUE ~ 0),
         diagnosed_tx_fpg_uncont_5.6 = case_when(standardized_tx == "on treatment" & fpg_mmol>=5.6 ~ 1, 
                                                 TRUE ~ 0),
         diagnosed_tx_fpg_uncont_7 = case_when(standardized_tx == "on treatment" & fpg_mmol>=7 ~ 1, 
                                                 TRUE ~ 0),
         diagnosed_tx_fpg_uncont_7.2 = case_when(standardized_tx == "on treatment" & fpg_mmol>=7.2 ~ 1, 
                                               TRUE ~ 0),
         diagnosed_tx_fpg_uncont_8.6 = case_when(standardized_tx == "on treatment" & fpg_mmol>=8.6 ~ 1, 
                                                 TRUE ~ 0)
         )

temp_fpg[, case_definition := "fpg_7_tx"]


# collapsing - a1c
temp_a1c <- copy(df)

temp_a1c <- temp_a1c[!is.na(a1c),]
temp_a1c <- temp_a1c[a1c<2|a1c>20, is_outlier := 1]
temp_a1c <- temp_a1c[is.na(is_outlier), is_outlier := 0]

temp_a1c <- temp_a1c %>%
  mutate(dm = case_when((a1c>=6.5 | standardized_tx == "on treatment") ~ 1, 
                        TRUE ~ 0), 
         undiagnosed = case_when(standardized_self_report=="negative self report" | standardized_self_report=="unknown" ~ 1, #assumes that 'unknown' self report is undiagnosed
                                 standardized_self_report=="positive self report" ~ 0, 
                                 standardized_self_report=="survey did not ask about diagnosis status" ~ NA), 
         diagnosed_untx = case_when(standardized_self_report=="positive self report" & (standardized_tx == "not on treatment" | standardized_tx =="unknown") ~ 1, #assumes that 'unknown' tx is not treated
                                    standardized_self_report=="survey did not ask about diagnosis status" ~ NA, 
                                    TRUE ~ 0), 
         alt_tx = case_when(standardized_tx == "on treatment" ~ 1, 
                            TRUE ~ 0), 
         diagnosed_tx_a1c_cont_6.5 = case_when(standardized_tx == "on treatment" & a1c<6.5 ~ 1, 
                                               TRUE ~ 0),
         diagnosed_tx_a1c_cont_7 = case_when(standardized_tx == "on treatment" & a1c<7 ~ 1, 
                                             TRUE ~ 0),
         diagnosed_tx_a1c_cont_8 = case_when(standardized_tx == "on treatment" & a1c<8 ~ 1, 
                                             TRUE ~ 0),
         diagnosed_tx_a1c_uncont_6.5 = case_when(standardized_tx=="on treatment" & a1c>=6.5 ~ 1, 
                                                 TRUE ~ 0),              
         diagnosed_tx_a1c_uncont_7 = case_when(standardized_tx=="on treatment" & a1c>=7 ~ 1, 
                                               TRUE ~ 0),
         diagnosed_tx_a1c_uncont_8 = case_when(standardized_tx=="on treatment" & a1c>=8 ~ 1, 
                                               TRUE ~ 0)
         )

temp_a1c[, case_definition := "hba1c_6.5_tx"]


#aggregate data for modeling
df_agg <- rbind(temp_fpg, temp_a1c, fill = TRUE)

df_agg<-as.data.table(df_agg[is_outlier == 0, ])

#stratify
stratify<-c('nid','location_name','year_start','year_end','study_name', 'age_start', 'age_end', 'sex', 'location_id', 'case_definition')

#aggregate and calculate prevalence among study pop
df_prev<- df_agg[,j=list(cases=sum(dm),sample_size=.N,effective_sample_size=((sum(sample_weight))^2)/sum(sample_weight^2),
                         mean=sum(dm*sample_weight)/sum(sample_weight)),
                 by=stratify]
df_prev <- as.data.table(df_prev)
df_prev$measure <- "prevalence"

#aggregate and calculate definitions and alternates
df_agg2 <- df_agg[dm == 1, ]

measures <- c("undiagnosed", "diagnosed_untx", "alt_tx", "diagnosed_tx_a1c_cont_6.5", "diagnosed_tx_a1c_cont_7", "diagnosed_tx_a1c_cont_8",
              "diagnosed_tx_a1c_uncont_6.5", "diagnosed_tx_a1c_uncont_7", "diagnosed_tx_a1c_uncont_8", "diagnosed_tx_fpg_cont_5.6", "diagnosed_tx_fpg_cont_7", "diagnosed_tx_fpg_cont_7.2",
              "diagnosed_tx_fpg_cont_7.22", "diagnosed_tx_fpg_cont_8.6", "diagnosed_tx_fpg_uncont_5.6", "diagnosed_tx_fpg_uncont_7", "diagnosed_tx_fpg_uncont_7.2", "diagnosed_tx_fpg_uncont_8.6")
list <- list()

#Loop through all the measures
for (measure in measures) {
  df_agg3 <- df_agg2[!is.na(get(measure))]
  df_result <- df_agg3[, .(cases = sum(get(measure)),
                           sample_size = .N,
                           effective_sample_size = (sum(sample_weight)^2) / sum(sample_weight^2),
                           mean = sum(get(measure) * sample_weight) / sum(sample_weight)),
                       by = stratify]
  df_result[, measure := measure]  # Add the measure as a new column
  list[[measure]] <- df_result  # Add the result to the list
}

df_collapsed <- rbindlist(list, use.names = TRUE, fill = TRUE)


#outlier if sample size is <=10 sample size
df_collapsed <- as.data.table(df_collapsed)
df_collapsed <- df_collapsed[sample_size <= 10, is_outlier := 1]
df_collapsed <- df_collapsed[is.na(is_outlier), is_outlier := 0]

# generating list of nids where 20% or more data points are outliered due to sample size 
check_age_bins <- as.data.table(table(df_collapsed$nid, df_collapsed$is_outlier))
check_age_bins <- pivot_wider(check_age_bins, names_from = V2, values_from = N)
check_age_bins <- as.data.table(check_age_bins)
check_age_bins <- check_age_bins[, percent_outlier := `1`/(`0`+`1`)]
check_age_bins_nids <- check_age_bins[percent_outlier >= 0.2, ]$V1 # Rerun collapse code using 20 year age bins for these nids

# removing effective sample size for NIDs without sample weights
df_collapsed <- df_collapsed[!(nid %in% no_missing_weight), effective_sample_size := NA]

# other outliers 
df_collapsed[nid == 458144 & age_end <35, is_outlier := 1] # age pattern in Afghanistan is opposite of expected - young ages lowest in undiagnosed and highest in control

# joining
df_prev <- df_prev %>%
  dplyr::rename(dm_prev = mean) %>%
  select(-c(cases, sample_size, effective_sample_size, measure))
df_collapsed <- left_join(df_collapsed, df_prev)
df_collapsed$type <- "microdata"

# split off data to upload to bundle
case_defs_list <- unique(subset(df, select = c(nid, standardized_case_definition)))
case_defs_list <- dplyr::rename(case_defs_list, case_definition=standardized_case_definition)
df_collapsed_subset <- left_join(case_defs_list, df_collapsed, by = c("nid", "case_definition"))
df_collapsed_subset <- df_collapsed_subset[measure %in% c("undiagnosed", "diagnosed_untx", "alt_tx", "diagnosed_tx_fpg_cont_7.2", "diagnosed_tx_fpg_uncont_7.2", "diagnosed_tx_a1c_cont_7", "diagnosed_tx_a1c_uncont_7"),]

# writing out all collapsed data for crosswalking
write.csv(df_collapsed, "FILEPATH/aggregated_microdata.csv", row.names = FALSE)

# writing out subset collapsed data for bundle upload
write.csv(df_collapsed_subset, "FILEPATH/aggregated_microdata.csv", row.names = FALSE)

