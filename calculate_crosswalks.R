# calculating crosswalk adjustments


# SET UP -------------------------------------------------------
rm(list = ls())

output_dir <- "FILEPATH"
invisible(sapply(list.files("FILEPATH", full.names = T), source))

# Import packages
require(data.table)
require(ggplot2)
require(msm)
library(dplyr)
library(reticulate)
library(RColorBrewer)
library(stringr)
reticulate::use_python("FILEPATH")
cw <- reticulate::import("crosswalk")
plots <- reticulate::import("crosswalk.plots")

# Read in aggregated microdata
dt <- as.data.table(read.csv(paste0(output_dir, "aggregated_microdata.csv")))

# Set output filepaths for prepped data, diagnostic plots, and MRBRT model outputs
out_dir<-paste0(output_dir,"prepped_data")
plot_dir<-paste0(output_dir,"plots")
mrbrt_dir<-paste0(output_dir,"mrbrt_models")

# Set release id 
release_id <- 16

# Set match vars 
match_vars<-c("nid","location_id","year_start","year_end","sex","age_start","age_end")

# Set reference definitions & get alt definitions
defs <- unique(dt$measure)

cont_ref_def <- "diagnosed_tx_fpg_cont_7.2"
uncont_ref_def <- "diagnosed_tx_fpg_uncont_7.2"

cont_alt_defs <- defs[grepl("_cont", defs)]
cont_alt_defs <- cont_alt_defs[!cont_alt_defs==cont_ref_def]

uncont_alt_defs <- defs[grepl("_uncont", defs)]
uncont_alt_defs <- uncont_alt_defs[!uncont_alt_defs==uncont_ref_def]

alt_defs <- c(cont_alt_defs, uncont_alt_defs)


# ----------------------------------------------

# Prep microdata sources ------------------------------------------------------------------
# Loop for each alt definition
set.seed(555)

for (x in alt_defs) {
  alt <- x
  if (alt %in% cont_alt_defs) {
    ref <- cont_ref_def
  } else {
    ref <- uncont_ref_def
  }

  comparison <- paste0(alt, "_to_", ref)
  message(comparison)
  
  # Read in aggregated microdata
  dt <- as.data.table(read.csv(paste0(output_dir, "aggregated_microdata.csv")))
  
  # Dropping outliers
  dt <- dt[!is_outlier == 1, ]
  
  # Calculate se 
  dt[!is.na(effective_sample_size), standard_error:= ifelse(mean!=0 & mean!=1, sqrt(mean*(1-mean) /effective_sample_size), 
                                         sqrt((mean*(1-mean)/effective_sample_size)+((1.96^2)/(4*effective_sample_size^2)))/(1+ (1.96^2)/effective_sample_size))]
  dt[is.na(effective_sample_size), standard_error := ifelse(mean!=0 & mean!=1, sqrt(mean*(1-mean) /sample_size), 
                                                            sqrt((mean*(1-mean)/sample_size)+((1.96^2)/(4*sample_size^2)))/(1+ (1.96^2)/sample_size))]
  
  # Only keep rows for specified comparison
  dt <- dt[measure == alt | measure == ref, ]
  
  # dropping rows where mean = 0 or 1
  message("dropping rows where mean = 0")
  drop_count <- nrow(dt[mean == 0 | mean == 1])
  tot_count <- nrow(dt)
  dt <- dt[!(mean == 0 | mean == 1),]
  message(paste0("dropped ", drop_count, " of ", tot_count, " rows"))
  
  # Subset to a dt where data was extracted using only reference
  message("separating into reference and alterate datasets")
  dt_ref <- dt[measure == ref,]
  dt_alt <- dt[measure == alt,]
  
  # Set names of comparison columns
  comp_cols <- c("mean", "standard_error")
  setnames(dt_ref, comp_cols, paste0(comp_cols, "_ref"))
  setnames(dt_alt, comp_cols, paste0(comp_cols, "_alt"))
  
  # Merge ref and alt data by match variables
  message(paste("merging reference and alternate dataset. matching on", paste(match_vars, collapse = ", ")))
  dt_merge <- merge(dt_ref, dt_alt, by = match_vars)
  
  # Transform alternate and ref into appropriate space
  message("Calculating logit difference and logit difference se")
  # Logit alternative
  dt_merge[, c("mean_alt_logit", "standard_error_alt_logit")]<-as.data.frame(
    cw$utils$linear_to_logit(
      mean = array(dt_merge$mean_alt),
      sd = array(dt_merge$standard_error_alt))
  )
  
  # Logit reference
  dt_merge[, c("mean_ref_logit", "standard_error_ref_logit")]<-as.data.frame(
    cw$utils$linear_to_logit(
      mean = array(dt_merge$mean_ref),
      sd = array(dt_merge$standard_error_ref))
  )
  
  # Logit difference 
  dt_merge<-  dt_merge %>%
    mutate(
      diff_logit_mean = mean_alt_logit - mean_ref_logit,
      diff_logit_se = sqrt(standard_error_alt_logit^2 + standard_error_ref_logit^2)
    )
  
  # Remove .x and .y vars
  merge_vars <- grep("\\.x|\\.y", names(dt_merge), value = T)
  dt_merge[, (merge_vars) := NULL]
  
  # Write outputs
  write.csv(x = dt_merge, file = paste0(out_dir, "/", comparison, ".csv"), row.names = F)
  message(paste0('writing CSV: ',out_dir, "/", comparison, ".csv" ))
  
  # Make vars for plotting
  dt_merge[, age_bin := paste(age_start, "to", age_end)]
  
  # merging loc 
  loc_meta <- get_location_metadata(location_set_id = 35, release_id = 16)
  loc_merge <- subset(loc_meta, select = c(location_id, region_name))
  
  dt_merge <- merge(dt_merge, loc_merge, all.x = TRUE, by = "location_id")
  
  # Write plots
  # Scatters
  message('creating diagnostics ')
  pdf(paste0(plot_dir, "/", comparison, "_scatters_age_facet.pdf"), width = 12)
  gg <- if (nrow(dt_merge) > 0) {
    ggplot(data = dt_merge) +
      geom_point(aes(x = mean_ref, y = mean_alt, col = region_name)) +
      facet_wrap(~age_bin, scales = "free") +
      geom_abline() +
      theme_bw() +
      xlab(paste0("Reference: Mean - ", ref)) +
      ylab(paste0("Alternative: Mean - ", alt))
  }
  print(gg)
  dev.off()
  
  pdf(paste0(plot_dir, "/", comparison, "_scatters_sex_facet.pdf"), width = 12)
  gg <- if (nrow(dt_merge) > 0) {
    ggplot(data = dt_merge) +
      geom_point(aes(x = mean_ref, y = mean_alt, col = region_name)) +
      facet_wrap(~sex, scales = "free") +
      geom_abline() +
      theme_bw() +
      xlab(paste0("Reference: Mean - ", ref)) +
      ylab(paste0("Alternative: Mean - ", alt))
  }
  print(gg)
  dev.off()
  
  pdf(paste0(plot_dir, "/", comparison, "_scatters.pdf"), width = 12)
  gg <- if (nrow(dt_merge) > 0) {
    ggplot(data = dt_merge) +
      geom_point(aes(x = mean_ref, y = mean_alt, col = region_name)) +
      geom_abline() +
      theme_bw() +
      xlab(paste0("Reference: Mean - ", ref)) +
      ylab(paste0("Alternative: Mean - ", alt))
  }
  print(gg)
  dev.off()
  
  #Plot ratios
  pdf(paste0(plot_dir, "/", comparison, "_ratios_by_age_sex.pdf"), width = 12)
  gg <- if (nrow(dt_merge) > 0) {
    ggplot(data = dt_merge, aes(x = age_start, y = diff_logit_mean, col = region_name)) +
      geom_jitter(size = 2, width = 1, alpha = .5) +
      facet_wrap(~sex) +
      theme_bw()
  }
  print(gg)
  dev.off()
  
  #Plot mean vs standard error
  pdf(paste0(plot_dir, "/", comparison, "_mean_se.pdf"), width = 12)
  gg <- ggplot(data = dt_merge, use.names = T, fill = T) +
    aes(x = diff_logit_mean, y = diff_logit_se, col = region_name) +
    geom_point(alpha = 0.5) +
    theme_bw()
  print(gg)
  dev.off()
  
}



#   -----------------------------------------------------------------------
# Run MRBRT model ------------------------------------------------------------------
# Loop for each alt definition
set.seed(555)

for (x in alt_defs) {
  alt <- x
  if (alt %in% cont_alt_defs) {
    ref <- cont_ref_def
  } else {
    ref <- uncont_ref_def
  }
  
  comparison <- paste0(alt, "_to_", ref)
  message(comparison)
  
  # read in prepped data for specified stage
  dt<-fread(paste0(out_dir, "/", comparison, ".csv"))
  
  # drop rows where val or SE are NA
  dt<-dt[!is.na(diff_logit_se)]
  dt<-dt[!is.na(diff_logit_mean)]
  
  # create dorms variables
  dt[, equation_alt := alt]
  dt[, equation_ref := ref]
  
  # CWData(), CovModel() and CWModel() functions
  dt$id <- as.integer(as.factor(dt$nid))
  
  # format data for meta-regression; pass in data.frame and variable names
  dat1 <- cw$CWData(
    df = dt,
    obs = "diff_logit_mean",   # matched differences in log space
    obs_se = "diff_logit_se",  # SE of matched differences in log space
    alt_dorms = "equation_alt",   # var for the alternative def/method
    ref_dorms = "equation_ref",   # var for the reference def/method
    covs = list(),            # list of (potential) covariate column names
    study_id = "id",          # var for random intercepts; i.e. (1|study_id)
    add_intercept = TRUE
  )
  
  fit1 <- cw$CWModel(
    cwdata = dat1,           # result of CWData() function call
    obs_type = "diff_logit",   # must be "diff_logit" or "diff_log"
    cov_models = list(       # specify covariate details
      cw$CovModel("intercept") ),
    gold_dorm = ref  # level of 'ref_dorms' that's the gold standard
  )
  
  fit1$fit(inlier_pct=0.9)
  
  # model predictions
  df_tmp<-as.data.table(fit1$create_result_df())
  
  # writing model results
  cols<-c("dorms", "cov_names", "beta", "beta_sd", "gamma")
  df_tmp2<-df_tmp %>% dplyr::select(cols) %>% filter(dorms==alt)
  df_tmp2$gamma<- unique(df_tmp[!is.na(gamma),]$gamma)
  
  write.csv(x = df_tmp, file = paste0(mrbrt_dir, "/", comparison, "_full_predictions.csv"), row.names = F)
  write.csv(x = df_tmp2, file = paste0(mrbrt_dir, "/", comparison, "_summary_predictions.csv"), row.names = F)
  
  # save model object
  py_save_object(object = fit1, filename = paste0(mrbrt_dir, "/", comparison, "_mrbrt_object.pkl"), pickle = "dill")
  
  # funnel plot
  plots$funnel_plot(
    cwmodel = fit1,
    cwdata = dat1,
    continuous_variables = list(),
    obs_method = alt,
    plot_note = paste0(alt," vs ",ref),
    plots_dir = paste0(plot_dir, "/"),
    file_name = paste0(comparison,"_mrbrt_funnel_plot"),
    write_file = TRUE
  )

}


#   -----------------------------------------------------------------------
# Compare betas plot ------------------------------------------------------------------
## loop to read in all new crosswalk betas
betas <- data.table()

for (x in alt_defs) {
  alt <- x
  if (alt %in% cont_alt_defs) {
    ref <- cont_ref_def
  } else {
    ref <- uncont_ref_def
  }
  
  comparison <- paste0(alt, "_to_", ref)
  message(comparison)
  
  dt<-as.data.table(fread(paste0(mrbrt_dir, "/", comparison, "_summary_predictions.csv")))
  
  if (alt %in% cont_alt_defs) {
    dt[, type := "controlled"]
  } else {
    dt[, type := "uncontrolled"]
  }
  
  if (dt$dorms %like% "fpg") {
    test <- "fpg"
  } else {
    test <- "hba1c"
  }
  
  numbers <- strsplit(dt$dorms, split = "_")
  numbers <- numbers[[1]]
  numbers <- last(numbers)
  
  dt[, threshold := paste(test, numbers, sep = "_")]
  
  betas <- rbind(betas, dt)
}

# calculating lower and upper for new betas
betas[,beta_sd_total:=sqrt(beta_sd^2+gamma)] 
betas[,lower:=beta-1.96*beta_sd_total]
betas[,upper:=beta+1.96*beta_sd_total]

# plotting 
g<-ggplot(betas, aes(x= beta, y=threshold, color = type)) +
  geom_point()+ 
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    width = 0.2,
    alpha= 0.3) +
  theme_bw()+
  geom_vline(xintercept=0, alpha=0.5) +
  scale_color_brewer(palette = "Dark2")

ggsave(paste0(plot_dir, "/betas_compare_by_threshold.pdf"),g, width=12, height = 10)



