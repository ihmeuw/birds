# data processing - age splitting

rm(list=ls())

library(data.table)
library(dplyr)
library(knitr)
library(kableExtra)
library(purrr)
library(ggplot2)

# Pull in bundle version and crosswalk version id map (must manually update age pattern MVIDs first!)
process_cv_ref_map_path <- "FILEPATH/bundle_upload_codebook.csv"
bvids <- fread(process_cv_ref_map_path)
cvids <- bvids

release_id <- 16

# output path
output_path <- "FILEPATH"

# to save (1) or not (0) a fully processed and age split crosswalk version
save_cx_post_age_split <- 1

for (bid in cvids$bundle_id){
  
  # SET OBJECTS -------------------------------------------------------------
  # set objects
  bvid<-bvids[bundle_id==bid, bundle_version_id]
  cvid<-cvids[bundle_id==bid, cvid_not_age_split]
  measure_upload<-cvids[bundle_id==bid, measure_upload]
  me_id<-bvids[bundle_id==bid, me_id]
  age_pattern_dismod_id<-bvids[bundle_id==bid, age_pattern_dismod_id] # this needs to be manually entered in csv
  diabetes_population_dismod_id<-bvids[bundle_id==bid, diabetes_population_dismod_id] # this needs to be manually entered in csv
  output_file_name<-paste0(measure_upload, "_post_age_split_fully_processed_bid_", bid, "_bvid_",bvid)
  
  location_set_id <- 9
  age_group_set_id <- 32
  
  year_id <- 2010 # year id you want the age pattern from
  
  draws <- paste0("draw_", 0:999)
  
  message(paste("Bundle version id is", bvid))

  # dt for tracking # datapoints involved in each processing step
  tracking_dt <- data.table()
  tracking_dt[, bundle_id := bid]
  tracking_dt[, bundle_version_id := bvid]
  
  # SOURCE FUNCTIONS ------------------------------------------------
  invisible(sapply(list.files("FILEPATH", full.names = T), source))
  
  library(reticulate)
  reticulate::use_python("FILEPATH")
  splitter <- import("pydisagg.ihme.splitter")
  
  # GET DATA, PATTERN, POPULATION ----------------------------------------------------------------
  # pulling crosswalk to age split and subsetting to data that needs age splitting
  message(paste0("getting crosswalk version ", cvid, " for bid ", bid))
  dt <- as.data.table(get_crosswalk_version(crosswalk_version_id = cvid, export = FALSE)) 
  
  dt_tosplit <- dt[((age_end-age_start)>=25) & !(mean == 0 | mean == 1), ] # pydisagg cannot split mean = 1 in logodds models, pydisagg produces se = NA for mean = 0 datapoints to split
  
  dt_nosplit <- setdiff(dt, dt_tosplit)
  
  ## calculate needed columns & subset
  dt_tosplit[, year_id := floor((year_start+year_end)/2)]
  dt_tosplit[, year_id := round(year_id/5)*5]
  dt_tosplit[year_id < 1990, year_id := 1990]
  dt_tosplit[year_id > 2020, year_id := 2020]

  dt_tosplit[sex == "Male", sex_id := 1]
  dt_tosplit[sex == "Female", sex_id := 2]
  
  dt_tosplit_subset <- subset(dt_tosplit, select = c(nid, seq, sex_id, age_start, age_end, location_id,
                                                     year_id, measure, mean, standard_error))
  
  ## add number datapoints to be age split to tracking dt 
  tracking_dt[, age_split_input_datapoints := nrow(dt_tosplit)]
  
  
  # pulling age split dismod model for pattern, merge with metadata
  loc_dt <- as.data.table(get_location_metadata(location_set_id = location_set_id, release_id = release_id))
  
  measure_dt <- as.data.table(get_ids("measure"))
  measure_dt <- measure_dt[measure %in% dt_tosplit$measure, ]
  
  pattern <- get_draws(gbd_id_type = "modelable_entity_id", gbd_id = me_id, source = "epi",
                       measure_id = measure_dt$measure_id, location_id = loc_dt[level == 1, ]$location_id, year_id = year_id,
                       sex_id = c(1,2), metric_id = 3, version_id = age_pattern_dismod_id,
                       release_id = release_id)
  
  loc_dt_merge <- subset(loc_dt, select = c(location_id, super_region_id))
  pattern <- rename(pattern, super_region_id = location_id)
  
  pattern <- merge(pattern, loc_dt_merge, by = "super_region_id", allow.cartesian = TRUE)
  pattern <- subset(pattern, select = -c(super_region_id))

  age_dt <- get_age_metadata(age_group_set_id = age_group_set_id, release_id = release_id)
  age_dt <- subset(age_dt, select = c(age_group_years_start, age_group_years_end, age_group_id))
  
  pattern <- merge(pattern, age_dt, all.x = TRUE)
  
  ## subset to needed columns
  pattern <- subset(pattern, select = -c(metric_id, model_version_id, modelable_entity_id))
  
  ## removing rows where all draws are 0 (pydisagg will error out) & fixing age_start in dt to split
  pattern[, draw_mean := rowMeans(.SD), .SDcols = paste0("draw_", 0:999)]
  pattern <- pattern[!draw_mean == 0, ]
  pattern <- subset(pattern, select = -c(draw_mean))
  
  pattern_min_age <- min(pattern$age_group_years_start)
  
  dt_tosplit_subset[age_start < pattern_min_age, age_start := pattern_min_age]
  
  
  # pulling population
  # use diabetes model for population since diabetes is denominator
  pop_diabetes <- get_model_results(gbd_team = 'epi', gbd_id = NUM, measure_id = 5, sex_id = c(1,2),
                                    release_id = 9, age_group_id = unique(pattern$age_group_id),
                                    model_version_id = diabetes_population_dismod_id)
  
  pop_count <- get_population(location_id = 'all', year_id = 'all', sex_id = c(1,2), 
                              release_id = release_id, age_group_id = unique(pattern$age_group_id))
  
  pop <- as.data.table(merge(pop_diabetes, pop_count, by = c("year_id", "sex_id", "location_id", "age_group_id")))
  pop[, population := mean*population]
  pop <- subset(pop, select = c(population, year_id, sex_id, location_id, age_group_id))
  
  
  # PYDISAGG ----------------------------------------------------------------
  data_config <- splitter$AgeDataConfig(
    index=c("nid","seq", "location_id", "year_id", "sex_id"),
    age_lwr="age_start",
    age_upr="age_end",
    val="mean",
    val_sd="standard_error")
  
  draw_cols <-grep("draw_",names(pattern), value = TRUE)
  
  pattern_config <- splitter$AgePatternConfig(
    by=list("sex_id", "location_id"),
    age_key="age_group_id",
    age_lwr="age_group_years_start",
    age_upr="age_group_years_end",
    draws=draw_cols)
  
  pop_config <- splitter$AgePopulationConfig(
    index=c("age_group_id", "location_id", "year_id", "sex_id"),
    val="population")
  
  age_splitter <- splitter$AgeSplitter(
    data=data_config, 
    pattern=pattern_config, 
    population=pop_config)
  
  age_split_output <- data.table()
  
  dt_tosplit_subset2 <- dt_tosplit_subset[measure == "proportion", ]
  pattern2 <- pattern[measure_id == 18, ]
  
  result <- age_splitter$split(
    data=dt_tosplit_subset2,
    pattern=pattern2,
    population=pop,
    model="logodds",  
    output_type="rate")
  
  result <- as.data.table(result)
  result <- result[, measure := "proportion"]
  
  age_split_output <- rbind(age_split_output, result)

  ## add number datapoints post age split to tracking dt 
  tracking_dt[, age_split_output_datapoints := nrow(age_split_output)]
  
  
  # CLEAN UP / COMBINING ----------------------------------------------------------------
  write.csv(age_split_output, paste0(output_path, "4_age_split/outputs/csvs/bundle_", bid, "_age_split.csv"))
  
  age_split_output[, split_cases := age_split_result*pop_population_aligned]
  age_split_output[, age_start := pat_age_group_years_start]
  age_split_output[, age_end := pat_age_group_years_end]
  age_split_output[, mean := age_split_result]
  age_split_output[, standard_error := age_split_result_se]
  
  age_split_output2 <- subset(age_split_output, select = c(seq, age_start, age_end,
                                                           mean, standard_error))
  
  dt_tosplit2 <- subset(dt_tosplit, select = -c(year_id, sex_id, age_start, age_end, mean, standard_error,
                                                lower, upper, cases, sample_size, uncertainty_type_value, variance, effective_sample_size,
                                                design_effect))
  
  dt_postsplit <- merge(dt_tosplit2, age_split_output2, by= "seq")
  
  dt_new <- rbind(dt_nosplit, dt_postsplit, fill = TRUE)
  
  
  # VETTING ----------------------------------------------------------------
  ## over age, facet pre and post age split for each super region
  age_split_output[, type := "Post split"]
  dt_tosplit_subset[, type := "Pre split"]
  
  age_split_plot <- rbind(age_split_output, dt_tosplit_subset, fill = TRUE)
  
  loc_dt_merge <- subset(loc_dt, select = c(location_id, super_region_name))
  age_split_plot <- merge(age_split_plot, loc_dt_merge, by = "location_id")
  
  for (i in unique(age_split_plot$measure)) {
    age_split_plot2 <- age_split_plot[measure == i, ]
    pdf(paste0(output_path,"4_age_split/outputs/plots/bundle_",bid, "_over_age_pre_post_split_", i, ".pdf"), width = 18, height = 8.5)
    for (x in unique(age_split_plot2$super_region_name)) {
      age_split_plot3 <- age_split_plot2[super_region_name == x, ]
      age_split_plot3[sex_id == 1, sex := "Male"]
      age_split_plot3[sex_id == 2, sex := "Female"]
      age_split_plot3[, age_midpoint := ((age_start+age_end)/2)]
      p <- ggplot(age_split_plot3) +
        geom_segment(aes(x = age_start, xend = age_end, y = mean, yend = mean), alpha = 0.5) +
        facet_wrap(~ sex + type)+
        labs(x="Age",
             y="Mean",
             title=paste0("Pre- and post-age split means over age: ", x, ", ", i),
             caption=paste0("Bundle ID: ", bid, ". Bundle Version: ", bvid)) +
        theme_light() + 
        theme(strip.text = element_text(size = 10))
      print(p)
    }
    dev.off()
  }
  
  
  # MISC OUTLIERS ----------------------------------------------------------------

  # SAVING CROSSWALK ----------------------------------------------------------------
  dt_new[, seq := NA]
  
  # write all data
  message("writing full processed dataset after age splitting")
  writexl::write_xlsx(list(extraction=dt_new),paste0(output_path,output_file_name, ".xlsx"))
  
  # save full crosswalk version
  if(save_cx_post_age_split==1) {
    description <- "Fully processed"
    result <- save_crosswalk_version(bvid,data_filepath = paste0(output_path,output_file_name, ".xlsx"),
                                     description = description)
    
    bvids<-fread(process_cv_ref_map_path)
    bvids$cvid_final[bvids$bundle_id==bid] <- result$crosswalk_version_id
    fwrite(bvids, file = process_cv_ref_map_path, append = FALSE)
  }
  
  # save tracking dt
  tracking_dt_all<-as.data.table(read.csv(paste0(output_path, "codebooks/data_processing_tracking.csv")))
  tracking_dt_all2 <- tracking_dt_all[bundle_id == bid, ]
  if ("age_split_input_datapoints" %in% colnames(tracking_dt_all)) {
    tracking_dt_all2 <- subset(tracking_dt_all2, select = -c(age_split_input_datapoints, age_split_output_datapoints))
  }
  tracking_dt_all2 <- merge(tracking_dt_all2, tracking_dt, by = c("bundle_id", "bundle_version_id"), all.x = TRUE)
  tracking_dt_all <- tracking_dt_all[!bundle_id == bid, ]
  tracking_dt_all <- rbind(tracking_dt_all, tracking_dt_all2, fill = TRUE)
  fwrite(tracking_dt_all, file = paste0(output_path, "codebooks/data_processing_tracking.csv"), append = FALSE)
  
}
