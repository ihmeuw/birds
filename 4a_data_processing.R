# data processing - within-study age sex splitting, sex splitting, applying crosswalks

rm(list = ls())

#Source functions
source("FILEPATH/sex_split.R")
invisible(sapply(list.files("FILEPATH", full.names = T), source))

library(ggplot2)
library(reticulate)
reticulate::use_python("FILEPATH")
cw <- reticulate::import("crosswalk")

# reading in codebook
cb_filepath <- "FILEPATH/bundle_upload_codebook.csv"
cb <- as.data.table(read.csv(cb_filepath))

# output path
output_path <- "FILEPATH"

# set reference control and uncontrol defs for crosswalking
cont_ref_def <- "diagnosed_tx_fpg_cont_7.2"
uncont_ref_def <- "diagnosed_tx_fpg_uncont_7.2"

# save age specific crosswalk version?
save_cx_age_spec <- 1

# save pre age split crosswalk version?
save_cx_pre_age_split <- 1

# loop for data processing (within-study age sex splitting, sex splitting, applying crosswalks)
for (x in unique(cb$bundle_id)) {
  bundle_id <- x
  bundle_version_id <- cb[bundle_id == x, ]$bundle_version_id
  measure_upload <- cb[bundle_id == x, ]$measure_upload
  
  # pulling bundle version
  bv <- as.data.table(get_bundle_version(bundle_version_id = bundle_version_id))
  
  # dt for tracking # datapoints involved in each processing step
  tracking_dt <- data.table()
  tracking_dt[, bundle_id := bundle_id]
  tracking_dt[, bundle_version_id := bundle_version_id]
  tracking_dt[, measure_upload := measure_upload]
  
  
  
  #  within-study age-sex split -------------------------------------
  message("within-study age-sex split ---------------------")
  
  if (1 %in% unique(bv$within_age_sex_split)) {
    # Subset to data needing to be age-sex split
    dt_split<-copy(bv[within_age_sex_split==1])
    dt_nosplit <- setdiff(bv, dt_split)
    
    # Create datatable with just sex specific data
    dt_sex<-dt_split[sex %in% c("Male", "Female"), ]
    
    # Calculate proportion of cases male and female
    dt_sex[, cases_total:= sum(cases), by = c("nid", "location_id", "year_start", "year_end")]
    dt_sex[, prop_cases := cases / cases_total]
    
    # Calculate proportion of sample male and female
    dt_sex[, ss_total:= sum(sample_size), by = c("nid", "location_id", "year_start", "year_end")]
    dt_sex[, prop_ss := sample_size / ss_total]
    
    # Calculate standard error of % cases & sample_size M and F
    dt_sex[, se_cases:= sqrt(prop_cases*(1-prop_cases) / cases_total)]
    dt_sex[, se_ss:= sqrt(prop_ss*(1-prop_ss) / ss_total)]
    
    # Estimate the ratio and standard error of the ratio of %cases/%sample 
    dt_sex[, ratio := prop_cases / prop_ss]
    dt_sex[, se_ratio:= sqrt( (prop_cases^2 / prop_ss^2)  * (se_cases^2/prop_cases^2 + se_ss^2/prop_ss^2) )]
    
    # Subset to necessary cols
    dt_sex <- subset(dt_sex, select = c("nid", "location_id", "year_start", "year_end", "seq", "sex", "ratio", "se_ratio", "prop_cases", "prop_ss"))

    # Rename seq of sex-specific rows
    setnames(dt_sex,"seq","parent_seq_sex")
    
    # Prep age-specific rows for split
    dt_age<-copy(dt_split[sex == "Both"])
    
    # Rename seq of age-specific rows 
    setnames(dt_age,"seq","parent_seq_age")
    
    # Duplicate age dt for both sexes
    dt_age_m<-copy(dt_age[,sex:="Male"])
    dt_age_f<-copy(dt_age[,sex:="Female"])
    dt_age<-rbindlist(list(dt_age_m,dt_age_f),use.names=T)
    
    # Merge on sex ratios
    dt_age<-merge(dt_age,dt_sex,by=c("nid", "location_id", "year_start", "year_end", "sex"))
    
    # Calculate age-sex specific mean, standard_error, cases, sample_size
    dt_age[, mean := mean * ratio]
    dt_age[, standard_error := sqrt(standard_error^2 * se_ratio^2 + standard_error^2 * ratio^2 + se_ratio^2 * mean^2)]
    dt_age[, cases := cases * prop_cases]
    dt_age[, sample_size := sample_size * prop_ss]
    if (!("note_modeler")%in%names(dt_age)) dt_age[,note_modeler:=NA]
    dt_age[, note_modeler := ifelse(is.na(note_modeler),paste("age,sex split using sex ratio", round(ratio, digits = 2)),
                                    paste(note_modeler, "| age,sex split using sex ratio", round(ratio, digits = 2)))]
    dt_age[, c("ratio", "se_ratio", "prop_cases", "prop_ss") := NULL]
    
    # Assign the original age seq as the crosswalk_parent_seq
    dt_age[!is.na(parent_seq_age),crosswalk_parent_seq:=parent_seq_age]
    
    # Clear other uncertainty information
    dt_age[,`:=`(lower=NA, upper=NA, uncertainty_type_value=NA, uncertainty_type=NA, effective_sample_size=NA)]
    
    # Check rows where mean > 1
    print(paste0("any rows with mean>1 after processing? ", any(dt_age$mean>1)))
    dt_age[mean > 1, `:=` (group_review = 0, note_modeler = paste0(note_modeler, " | group reviewed out because age-sex split over 1"))]
    
    # write out csvs with pre and post split data for possible vetting later
    write.csv(dt_split, paste0(output_path,"1_age_sex_lit/csvs/bundle_",bundle_id, "_pre_age_sex_lit_split.csv"))
    write.csv(dt_age, paste0(output_path,"1_age_sex_lit/csvs/bundle_",bundle_id, "_post_age_sex_lit_split.csv"))
    
    # combine
    bv<-rbindlist(list(dt_nosplit,dt_age),use.names=T,fill=T)
    
    # add to tracking dt
    tracking_dt[, age_sex_split_input_datapoints := nrow(dt_split)]
    tracking_dt[, age_sex_split_output_datapoints := nrow(dt_age)]
    
    # run scatter plot function
    vet_agesexsplit <- dt_split[, parent_seq_age := seq]
    vet_agesexsplit <- vet_agesexsplit[sex == "Both", ]
    vet_agesexsplit <- merge(vet_agesexsplit, dt_age, by = c("parent_seq_age", "location_id", "ihme_loc_id"))
    
    merge_sr <- get_location_metadata(release_id = 16, location_set_id = 35)
    merge_sr <- subset(merge_sr, select = c(location_id, super_region_name))
    
    vet_agesexsplit <- merge(vet_agesexsplit, merge_sr, by = "location_id")
    
    pdf(paste0(output_path,"1_age_sex_lit/plots/bundle_",bundle_id, "_scatter.pdf"), width = 18, height = 8.5)
    for (i in unique(vet_agesexsplit$measure.x)) {
      vet_agesexsplit_plot <- vet_agesexsplit[measure.x == i, ]
      p <- ggplot(vet_agesexsplit_plot, aes(x=mean.x, y=mean.y, color=sex.y)) +
        geom_point(size=1.2) +
        geom_abline(alpha = 0.5, color = "lightgrey") +
        geom_errorbar(aes(ymin=mean.y-1.96*standard_error.y, ymax=mean.y+1.96*standard_error.y), alpha=0.5) +
        facet_wrap(~ super_region_name)+
        labs(x="Mean - Before age-sex split",
             y="Mean - After age-sex split",
             title=paste0("Mean before and after within-literature age-sex split: ", i),
             color="Sex",
             caption=paste0("Bundle ID: ", bundle_id, ". Bundle Version: ", bundle_version_id)) +
        theme_light() + 
        theme(strip.text = element_text(size = 10))
      print(p)
    }
    dev.off()

  } else {
    message("no data to within-study age-sex split")
  }
  
  
  
  # sex split --------------------------------------------------------------------------------------------------------------------------
  
  # splitting into both sex and sex specific, removing to be sex split data from bv
  to_split <- bv[sex == "Both", ]
  sex_specific <- bv[!sex == "Both", ]
  
  to_split[, note_modeler := NA]
  sex_specific[, note_modeler := NA]
  bv[, note_modeler := NA]
  
  to_split[, specificity := NA]
  sex_specific[, specificity := NA]
  bv[, specificity := NA]
  
  # writing csvs
  write.csv(to_split, paste0(output_path, "2_sex_split/inputs/to_split_", measure_upload, ".csv"))
  write.csv(sex_specific, paste0(output_path, "2_sex_split/inputs/sex_specific_", measure_upload, ".csv"))
  
  # removing to be sex split data from bdt
  bv <- setdiff(bv, to_split)
  
  # read in sex_split_ref_map for sex split model settings
  sex_split_ref_map <- as.data.table(read.csv(paste0(output_path, "codebooks/sex_split_codebook.csv")))
  sex_split_ref_map <- sex_split_ref_map[bundle_id==x, ]
  
  # BIRDS team adapted sex split function
  if (nrow(to_split) > 0) {
  message("running sex split function")
  
  date <- strtrim(gsub(":", "_", gsub(" ", "__", gsub("-", "_", Sys.time()))), 17)
  sex_split_path <- paste0(output_path,"2_sex_split/outputs/",measure_upload, "/sex_split_", date, "/")
  topic_name <- measure_upload
  
  sex_split(topic_name = topic_name, 
            output_dir = paste0(output_path, "2_sex_split/outputs/",measure_upload, "/"), 
            bundle_version_id = NULL, 
            data_all_csv = NULL, 
            data_to_split_csv = paste0(output_path, "2_sex_split/inputs/to_split_", measure_upload, ".csv"), 
            data_raw_sex_specific_csv = paste0(output_path, "2_sex_split/inputs/sex_specific_", measure_upload, ".csv"), 
            nids_to_drop = c(), 
            cv_drop = c(), 
            mrbrt_model = NULL, 
            mrbrt_model_age_cv = sex_split_ref_map$mrbrt_model_age_cv, 
            mrbrt_model_linear_age_cv = sex_split_ref_map$mrbrt_model_linear_age_cv,
            mrbrt_model_age_spline_knots = as.numeric(unlist(strsplit(sex_split_ref_map$mrbrt_model_age_spline_knots,","))),
            mrbrt_model_age_spline_degree = sex_split_ref_map$mrbrt_model_age_spline_degree,
            mrbrt_model_age_spline_linear_tail_left = sex_split_ref_map$mrbrt_model_age_spline_linear_tail_left,
            mrbrt_model_age_spline_linear_tail_right = sex_split_ref_map$mrbrt_model_age_spline_linear_tail_right,
            release_id = 9, 
            measure = "proportion",
            vetting_plots = TRUE)
  } else {
    message("no data to sex split")
  }
  
  # reading in post sex split data
  postsplit_dt <- read.csv(paste0(sex_split_path, "sex_split_", date, "_", topic_name, "_post_sex_split_data.csv"))
  postsplit_dt <- postsplit_dt[, intersect(colnames(bv), colnames(postsplit_dt))]
  
  postsplit_dt <- as.data.table(postsplit_dt)
  
  # binding onto bv, fixing seqs
  postsplit_dt[, crosswalk_parent_seq := seq]
  
  bv <- rbind(bv, postsplit_dt, fill = TRUE)
  
  # add to tracking dt
  sex_split_vetting_dt <- as.data.table(read.csv(paste0(sex_split_path, "sex_split_", date, "_", topic_name, "_vetting_data.csv")))
  
  tracking_dt[, sex_split_sex_matches_datapoints := sex_split_vetting_dt[Statement == "Sex matches found", ]$Number]
  tracking_dt[, sex_split_input_datapoints := sex_split_vetting_dt[Statement == "Observations to be sex split", ]$Number]
  tracking_dt[, sex_split_output_datapoints := sex_split_vetting_dt[Statement == "Number of sex split observations", ]$Number]
  
  
  
  # crosswalking --------------------------------------------------------------------------------------------------------------------------
  if(measure_upload %in% c("diagnosed_tx_uncont", "diagnosed_tx_cont")) {
    message("apply crosswalks")
    
    if(measure_upload == "diagnosed_tx_uncont") {
      ref_def <- uncont_ref_def
    } else {
      ref_def <- cont_ref_def
    }
    
    defs <- unique(bv$definition)
    alt_defs <- defs[!defs==ref_def]
    
    for (x in alt_defs){
      comparison <- paste0(x, "_to_", ref_def)
      message(comparison)
      
      # split into alt and reference data sets
      toxwalk_dt <- bv[definition == x, ]
      noxwalk_dt <- bv[!definition == x, ]
      
      # calculate missing SE
      toxwalk_dt[is.na(standard_error) & !is.na(lower) & !is.na(upper), standard_error := (upper-lower)/3.92]
      z <- qnorm(0.975)
      toxwalk_dt[is.na(standard_error) & measure == "proportion", standard_error := sqrt(mean*(1-mean)/sample_size + z^2/(4*sample_size^2))]
      
      # pull in mrbrt model
      model_parent_dir<-paste0(output_path,"crosswalks/mrbrt_models/")
      
      fit1 <- py_load_object(filename = paste0(model_parent_dir, comparison, "_mrbrt_object.pkl"), pickle = "dill")
      
      # first deal with 0s and 1s in data to crosswalk - cannot be adjusted for log or logit scale crosswalks
      toxwalk_dt <- toxwalk_dt[!mean == 0 & !mean == 1, ]
      
      # apply crosswalk adjustment using crosswalk package 
      pred_tmp <- fit1$adjust_orig_vals(
        df = toxwalk_dt,
        orig_dorms = "definition",
        orig_vals_mean = "mean", # mean in un-transformed space
        orig_vals_se = "standard_error" # SE in un-transformed space
      )
      
      #convert mean and se back to logit space
      toxwalk_dt[, c(
        "meanvar_adjusted", "sdvar_adjusted", "pred_logit",
        "pred_se_logit", "data_id")] <- pred_tmp
      
      toxwalk_dt[, c("meanvar_adj_logit", "sdvar_adj_logit")]<-as.data.frame(cw$utils$linear_to_logit(
        mean = array(toxwalk_dt$meanvar_adjusted),
        sd = array(toxwalk_dt$sdvar_adjusted)))
      toxwalk_dt$lower_adj_logit<- toxwalk_dt$meanvar_adj_logit - 1.96*toxwalk_dt$sdvar_adj_logit
      toxwalk_dt$upper_adj_logit<- toxwalk_dt$meanvar_adj_logit + 1.96*toxwalk_dt$sdvar_adj_logit
      
      # now inverse logit out upper and lower out
      toxwalk_dt[,lower_adj:=exp(lower_adj_logit)/(1+exp(lower_adj_logit)) ]
      toxwalk_dt[,upper_adj:=exp(upper_adj_logit)/(1+exp(upper_adj_logit)) ]
      
      # write out csvs with pre and post crosswalk data for possible vetting later
      write.csv(toxwalk_dt, paste0(output_path,"3_crosswalk/csvs/bundle_",bundle_id, "_", comparison, "_crosswalk.csv"))
      
      # recode adjusted values 
      toxwalk_dt[,mean:=meanvar_adjusted]
      toxwalk_dt[,standard_error:=sdvar_adjusted]
      toxwalk_dt[,upper:= upper_adj]
      toxwalk_dt[,lower:=lower_adj]
      
      # remove extra columns
      toxwalk_dt <- as.data.frame(toxwalk_dt)
      toxwalk_dt <- toxwalk_dt[, intersect(colnames(noxwalk_dt), colnames(toxwalk_dt))]
      toxwalk_dt <- as.data.table(toxwalk_dt)
      
      # fix uncertainty information
      toxwalk_dt[,`:=`(uncertainty_type_value=95)]
      
      # change seq to crosswalk_parent_seq where crosswalk parent seq is na 
      toxwalk_dt[is.na(crosswalk_parent_seq),crosswalk_parent_seq:=as.integer(seq)]
      toxwalk_dt[,seq:=NULL]
      
      bv <- rbindlist(list(noxwalk_dt, toxwalk_dt), use.names=T, fill = T)
    }
    
    # crosswalk plots
    message("Crosswalk vetting plots")
    
    # read in csvs from apply_mrbrt_xwalk function
    file<-paste0(output_path,"3_crosswalk/csvs/")
    filenames <- list.files(file, pattern=paste0("bundle_",bundle_id, "_*"), full.names=TRUE)
    dt_xwalk <- as.data.table(rbindlist(lapply(filenames, fread), fill=T)) 
    
    merge_sr <- get_location_metadata(release_id = 16, location_set_id = 35)
    merge_sr <- subset(merge_sr, select = c(location_id, super_region_name))
    
    vet_xwalk <- merge(dt_xwalk, merge_sr, by = "location_id")
    
    vet_xwalk[,age_midpoint:=(age_start+age_end)/2]
    
    vet_xwalk[age_midpoint < 20, age_midpoint_label := "Age midpoint less than 20 years"]
    vet_xwalk[age_midpoint >= 20 & age_midpoint < 40, age_midpoint_label := "Age midpoint between 20 and 40 years"]
    vet_xwalk[age_midpoint >= 40 & age_midpoint < 60, age_midpoint_label := "Age midpoint between 40 and 60 years"]
    vet_xwalk[age_midpoint >= 60 & age_midpoint < 80, age_midpoint_label := "Age midpoint between 60 and 80 years"]
    vet_xwalk[age_midpoint >= 80, age_midpoint_label := "Age midpoint greater than 80 years"]
    
    # scatter pre and post means
    pdf(paste0(output_path,"3_crosswalk/plots/",measure_upload, "_scatter_facet_superregion.pdf"), width = 18, height = 8.5)
    for (i in unique(vet_xwalk$measure)) {
      vet_xwalk_plot <- vet_xwalk[measure == i, ]
      p <- ggplot(vet_xwalk_plot, aes(x=mean, y=meanvar_adjusted, color=definition)) +
        geom_point(size=1.2) +
        geom_abline(alpha = 0.5, color = "lightgrey") +
        facet_wrap(~ super_region_name)+
        labs(x="Mean - Before crosswalking",
             y="Mean - After crosswalking",
             title=paste0("Mean before and after crosswalking: ", i),
             color="Original alternative definition")
        theme_light() + 
        theme(strip.text = element_text(size = 10))
      print(p)
    }
    dev.off()
    
    # post means over age with UI
    pdf(paste0(output_path,"3_crosswalk/plots/",measure_upload, "_over_age_post_crosswalk_means.pdf"), width = 18, height = 8.5)
    for (i in unique(vet_xwalk$measure)) {
      vet_xwalk_plot <- vet_xwalk[measure == i, ]
      p <- ggplot(vet_xwalk_plot, aes(group = definition, color=definition)) +
        geom_pointrange(aes(x = age_midpoint, y = meanvar_adjusted, ymin = lower_adj, ymax = upper_adj), alpha = 0.5) +
        facet_wrap(~ super_region_name)+
        labs(x="Midage",
             y="Mean - After crosswalking",
             title=paste0("Mean after crosswalking over age: ", i),
             color="Original alternative definition")
        theme_light() + 
        theme(strip.text = element_text(size = 10))
      print(p)
    }
    dev.off()
    
    # post means and orig reference means over age with UI
    dt_xwalk[,mean:=meanvar_adjusted]
    dt_xwalk[,standard_error:=sdvar_adjusted]
    dt_xwalk[,upper:= upper_adj]
    dt_xwalk[,lower:=lower_adj]
    
    input_dt_plot <- rbind(bv, dt_xwalk, fill = TRUE)
    
    vet_xwalk <- merge(input_dt_plot, merge_sr, by = "location_id")
    vet_xwalk[,age_midpoint:=(age_start+age_end)/2]
    
    pdf(paste0(output_path,"3_crosswalk/plots/",measure_upload, "_over_age_compare_means.pdf"), width = 18, height = 8.5)
    for (i in unique(vet_xwalk$measure)) {
      vet_xwalk_plot <- vet_xwalk[measure == i, ]
      for (x in unique(vet_xwalk_plot$super_region_name)) {
        vet_xwalk_plot2 <- vet_xwalk_plot[super_region_name == x, ]
        p <- ggplot(vet_xwalk_plot2, aes(group = definition, color=definition)) +
          geom_segment(aes(x = age_start, xend = age_end, y = mean, yend = mean), alpha = 0.5) +
          geom_pointrange(aes(x = age_midpoint, y = mean, ymin = lower, ymax = upper), alpha = 0.5) +
          facet_wrap(~ sex)+
          labs(x="Midage",
               y="Mean",
               title=paste0("Post-crosswalk means compared to reference means over age: ", i, ", ", x),
               color="Original case definition")
          theme_light() + 
          theme(strip.text = element_text(size = 10))
        print(p)
      }
    }
    dev.off()
    
    pdf(paste0(output_path,"3_crosswalk/plots/",measure_upload, "_over_age_compare_means_no_ui.pdf"), width = 18, height = 8.5)
    for (i in unique(vet_xwalk$measure)) {
      vet_xwalk_plot <- vet_xwalk[measure == i, ]
      for (x in unique(vet_xwalk_plot$super_region_name)) {
        vet_xwalk_plot2 <- vet_xwalk_plot[super_region_name == x, ]
        p <- ggplot(vet_xwalk_plot2, aes(group = definition, color=definition)) +
          geom_segment(aes(x = age_start, xend = age_end, y = mean, yend = mean), alpha = 0.5) +
          geom_point(aes(x = age_midpoint, y = mean), alpha = 0.5) +
          facet_wrap(~ sex)+
          labs(x="Midage",
               y="Mean",
               title=paste0("Post-crosswalk means compared to reference means over age: ", i, ", ", x),
               color="Original case definition")
          theme_light() + 
          theme(strip.text = element_text(size = 10))
        print(p)
      }
    }
    dev.off()
    
    # add to tracking dt
    file<-paste0(output_path,"3_crosswalk/csvs/")
    filenames <- list.files(file, pattern=paste0("bundle_",bundle_id, "_*"), full.names=TRUE)
    dt_xwalk <- as.data.table(rbindlist(lapply(filenames, fread), fill=T))
    
    for (a in alt_defs){
      tracking_dt[, paste0("crosswalk_", a, "_datapoints") := nrow(dt_xwalk[definition==a])]
    }
  }
  
  
  
  # clean up & saving crosswalk version --------------------------------------------------------------------------------------------------------------------------
  message("clean up")
  
  # outlier based on sample size
  bv <- bv[sample_size <= 10, is_outlier := 1]
  bv <- bv[is.na(is_outlier), is_outlier := 0]
  
  # assign crosswalk_parent_seq, clear seq for crosswalk upload
  bv[is.na(crosswalk_parent_seq), crosswalk_parent_seq := seq]
  bv[, seq := NA]
  
  # fix uncertainty for validations
  bv[is.na(lower)|is.na(upper), `:=` (lower = NA, upper = NA, uncertainty_type_value = NA)]
  
  # save only age-specific data to estimate age pattern for age-splitting
  bdt_agespec <- bv
  bdt_agespec$age_range <- bdt_agespec$age_end-bdt_agespec$age_start
  
  # keep only data with age range <25 years
  bdt_agespec <- bdt_agespec[bdt_agespec$age_range< 25,]
  
  # keep only data with sample size > 50
  bdt_agespec <- bdt_agespec[bdt_agespec$sample_size>30,]
  
  # add to tracking dt
  tracking_dt[, age_split_dismod_input_datapoints := nrow(bdt_agespec)]
  
  # plot age spec data 
  merge_sr <- get_location_metadata(release_id = 16, location_set_id = 35)
  merge_sr <- subset(merge_sr, select = c(location_id, super_region_name))
  
  input_dt <- merge(bdt_agespec, merge_sr, by = "location_id")
  
  input_dt[,age_midpoint:=(age_start+age_end)/2]
  
  pdf(paste0(output_path,"4_age_split/inputs/plots/",measure_upload, "_age_spec_inputs_over_age_facet_superregion.pdf"), width = 18, height = 8.5)
  for (i in unique(input_dt$measure)) {
    input_dt_plot <- input_dt[measure == i, ]
    p <- ggplot(input_dt_plot, aes(group = sex, color=sex)) +
      geom_pointrange(aes(x = age_midpoint, y = mean, ymin = lower, ymax = upper), alpha = 0.5) +
      facet_wrap(~ super_region_name)+
      labs(x="Midage",
           y="Mean",
           title=paste0("Post sex split/crosswalked - age specific data for age splitting DisMod model: ", i),
           color="Sex")
      theme_light() + 
      theme(strip.text = element_text(size = 10))
    print(p)
    p <- ggplot(input_dt_plot, aes(group = sex, color=sex)) +
      geom_point(aes(x = age_midpoint, y = mean), alpha = 0.5) +
      facet_wrap(~ super_region_name)+
      labs(x="Midage",
           y="Mean",
           title=paste0("Post sex split/crosswalked - age specific data for age splitting DisMod model: ", i),
           color="Sex")
      theme_light() + 
      theme(strip.text = element_text(size = 10))
    print(p)
  }
  dev.off()
  
  pdf(paste0(output_path,"4_age_split/inputs/plots/",measure_upload, "_age_spec_inputs_over_age_facet_sex.pdf"), width = 18, height = 8.5)
  for (i in unique(input_dt$measure)) {
    input_dt_plot <- input_dt[measure == i, ]
    p <- ggplot(input_dt_plot, aes(group = super_region_name, color=super_region_name)) +
      geom_pointrange(aes(x = age_midpoint, y = mean, ymin = lower, ymax = upper), alpha = 0.5) +
      facet_wrap(~ sex)+
      labs(x="Midage",
           y="Mean",
           title=paste0("Post sex split/crosswalked - age specific data for age splitting DisMod model: ", i),
           color="Super Region")
      theme_light() + 
      theme(strip.text = element_text(size = 10))
    print(p)
    p <- ggplot(input_dt_plot, aes(group = super_region_name, color=super_region_name)) +
      geom_point(aes(x = age_midpoint, y = mean), alpha = 0.5) +
      facet_wrap(~ sex)+
      labs(x="Midage",
           y="Mean",
           title=paste0("Post sex split/crosswalked - age specific data for age splitting DisMod model: ", i),
           color="Super Region")
      theme_light() + 
      theme(strip.text = element_text(size = 10))
    print(p)
  }
  dev.off()
  
  # write age spec data, save crosswalk if specified
  bundle_id2 <- bundle_id
  
  message("writing age-specific dataset after processing")
  output_file_name2<-paste0(measure_upload, "_age_specific_sex_split_and_crosswalked_bid_", bundle_id, "_bvid_",bundle_version_id)
  writexl::write_xlsx(list(extraction=bdt_agespec),paste0(output_path,"4_age_split/inputs/", output_file_name2, ".xlsx"))

  if(save_cx_age_spec==1) {
    description <- "Age-specific data only to model age pattern, age range <25 years and sample size > 30"
    result <- save_crosswalk_version(bundle_version_id,data_filepath = paste0(output_path,"4_age_split/inputs/", output_file_name2, ".xlsx"),
                                     description = description)
    
    bvids<-fread(cb_filepath)
    bvids$cvid_age_specific[bvids$bundle_id==bundle_id2] <- result$crosswalk_version_id
    fwrite(bvids, file = cb_filepath, append = FALSE)
  }
  
  # write all data
  message("writing full dataset after processing, not age split")
  output_file_name <- paste0(measure_upload, "_sex_split_and_crosswalked_bid_", bundle_id, "_bvid_",bundle_version_id)
  writexl::write_xlsx(list(extraction=bv),paste0(output_path,output_file_name, ".xlsx"))
  
  # save full crosswalk version
  if(save_cx_pre_age_split==1) {
    description <- "Processed but not age split"
    result <- save_crosswalk_version(bundle_version_id,data_filepath = paste0(output_path,output_file_name, ".xlsx"),
                                     description = description)
    
    bvids<-fread(cb_filepath)
    bvids$cvid_not_age_split[bvids$bundle_id==bundle_id2] <- result$crosswalk_version_id
    fwrite(bvids, file = cb_filepath, append = FALSE)
  }
  
  # save tracking dt
  tracking_dt_all<-as.data.table(read.csv(paste0(output_path, "codebooks/data_processing_tracking.csv")))
  tracking_dt_all <- tracking_dt_all[!bundle_id == bundle_id2, ]
  tracking_dt_all <- rbind(tracking_dt_all, tracking_dt, fill = TRUE)
  fwrite(tracking_dt_all, file = paste0(output_path, "codebooks/data_processing_tracking.csv"), append = FALSE)
  
}

