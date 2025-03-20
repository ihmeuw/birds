#########################################################################################################################
## Standalone child script to visualize risk factors for neuro conditions
## Output: csv files of PAFs
#########################################################################################################################


## Set up ------------------------------------------------------------------------------
rm(list=ls())

#Source config file and functions, set directory
source(FILEPATH)

loc_ids = ID

ages <- c(IDS)

## Attributable burden function -----------------------------------------------------------------
df <- get_outputs(topic = "rei", cause_id = c(IDS), 
                  measure_id = ID, metric_id = ID, 
                  compare_version_id = compare_risk_id, age_group_id = ID, sex_id = c(ID),
                  year_id = year_end, location_id = ID, gbd_round_id = ID, decomp_step = "STEP")

df$cause_name = factor(
  df$cause_name, 
  levels = c("Stroke", "Idiopathic developmental intellectual disability", 
             "Alzheimer's disease and other dementias", "Multiple sclerosis", 
             "Idiopathic epilepsy", "Parkinson's disease", "Meningitis", "Encephalitis"))
write.csv(df, paste0(FILEPATH))