###############################################################################
# Title: Uploading bundle data for cascade of care models
###############################################################################

rm(list = ls())

#Source central functions
invisible(sapply(list.files("FILEPATH", full.names = T), source))

#libraries
library(openxlsx)

#reading in data
df <- as.data.table(read.csv("FILEPATH/capstone_combined_data.csv"))

# reading in codebook 
cb <- as.data.table(read.csv("FILEPATH/bundle_upload_codebook.csv"))

# checking to make sure cb vars is up to date
unique(df$measure)
unique(cb$vars)

# loop for uploading bundle
for (x in unique(cb$bundle_id)) {
  bundle_id <- x
  measure_upload <- cb[bundle_id == x, ]$measure_upload
  vars <- cb[bundle_id == x, ]$vars
  me_id <- cb[bundle_id == x, ]$me_id
  
  new_bun <- df[measure %like% vars, ]
  
  # clearing bundle -----------------------------------------------------------------
  bun<-get_bundle_data(bundle_id =bundle_id)
  
  if (nrow(bun) > 0) {
    empty_dt <- data.table("seq" = bun[, seq])
    
    path<-'FILEPATH'
    filename <- paste0(measure_upload, "_empty_bundle.xlsx")
    
    write.xlsx(empty_dt, paste0(path,filename), sheetName = "extraction")
    
    upload_bundle_data(bundle_id = bundle_id, filepath = paste0(path,filename))
  }
  
  # updating bundle, saving bundle version -----------------------------------------------------------------------
  
  new_bun$seq <- NA
  new_bun$design_effect <- NA
  new_bun$representative_name <- "Unknown"
  new_bun$recall_type <- "Point"
  new_bun$input_type <- NA
  new_bun$sampling_type <- NA
  new_bun$source_type <- "Unidentifiable"
  new_bun$urbanicity_type <- "Unknown"
  new_bun$unit_value_as_published <- 1
  new_bun$recall_type_value <- NA
  new_bun$underlying_nid <- NA
  new_bun$unit_type <- 'Person'
  new_bun$definition <- new_bun$measure
  new_bun$measure <- 'proportion'
  
  new_bun[effective_sample_size<cases, effective_sample_size := NA]
  
  path<-'FILEPATH'
  filename <- paste0(measure_upload, "_bundle.xlsx")
  
  write.xlsx(new_bun, paste0(path,filename), sheetName = "extraction")
  
  upload_bundle_data(bundle_id = bundle_id, filepath = paste0(path,filename))
  
  bun_new<-get_bundle_data(bundle_id =bundle_id)
  
  if (!nrow(bun_new) == nrow(new_bun)) {
    stop(paste0("Bundle ID ", bundle_id, " does not have same number of rows as uploaded data"))
  }
  
  result <- save_bundle_version(bundle_id = bundle_id)
  
  cb_update <- read.csv("FILEPATH/bundle_upload_codebook.csv")
  cb_update$bundle_version_id[cb_update$bundle_id==x] <- result$bundle_version_id
  write.csv(cb_update, file = "FILEPATH/bundle_upload_codebook.csv", row.names = FALSE)

}
