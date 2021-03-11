### Forecast to 2050 ---------------------------------------------------------------------------------------------------------------------------
### Purpose: Forecast non-fatal estimates for hearing to the year 2050 by predicting rates and then multiplying by forecasted population
## Set up
rm(list=ls())
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "FILEPATH"
  h_root <- "FILEPATH"
  l_root <- "FILEPATH"
} else {
  j_root <- "FILEPATH"
  h_root <- "FILEPATH"
  l_root <- "FILEPATH"
}
pacman::p_load(data.table, openxlsx, ggplot2, magrittr, dplyr, stringr, splines)
library(janitor, lib.loc = "FILEPATH")

date <- gsub("-", "_", Sys.Date())
date <- Sys.Date()
my_dir <- paste0(j_root, "FILEPATH")
source("FILEPATH")
source("FILEPATH")
## Set variables here (only one severity for hearing loss)
year_forecast_list <- c(2020, 2030, 2040, 2050)
severity <- "Hearing loss 35+ dB"
loss <- "35+"
year_forecast<- 2050
for (year_forecast in year_forecast_list){
  
  ## Pull in age-specific rates for Dismod years 1990-2019
  if(loss == "35+"){  
    df <- read.csv(paste0(my_dir, "FILEPATH"))
  }
  if(loss == "20+"){
    df <- read.csv(paste0(my_dir, "FILEPATH"))
  }
  dt <- as.data.table(df)
  dt <- dt[dt$sex_id!=3, ]
  dt <- dt[age_group_id!=22, ]
  
  
  ## WHO and World Bank regions
  regions <- read.csv(paste0(my_dir, "ihme_world_bank_locs.csv"))
  #Pull in forecasted population numbers
  df_pop <- read.csv("FILEPATH")
  df_pop <- df_pop[df_pop$year_id==year_forecast, ]
  setnames(df_pop, "population", "pop_size")
  df_pop <- merge(df_pop, regions[,c("location_id","region_name")], all.x=T)
  df_pop[,c("quantile","scenario", "income_group_id", "income_group")] <- NULL
  
  #Subset to those with no forecasted population (df_NA) and those with forecasted population (df_pop2)
  df_NA <- df_pop[is.na(df_pop$pop_size), ]
  df_pop2 <- df_pop[df_pop$location_id!=1 & !is.na(df_pop$pop_size), ]
  
  #Aggregate to WHO and World Bank regions
  df_pop_agg <- aggregate(df_pop2$pop_size, by=list(Category=df_pop2$region_name, df_pop2$age_group_id, df_pop2$sex_id), FUN=sum)
  setnames(df_pop_agg, c("Category", "Group.2", "Group.3", "x"), c("region_name", "age_group_id", "sex_id", "pop_size"))
  df_pop_agg$year_id <- as.numeric(year_forecast)
  df_pop_agg$X <- NULL
  
  
  ## Add age variables to prediction matrix and dt
  ages_2019 = get_age_metadata(gbd_round_id=6)
  #Add to pred matrix
  df_pop_agg <- merge(df_pop_agg, ages_2019[,c("age_group_id", "age_group_name", "age_group_years_start", "age_group_years_end")])
  df_pop_agg$midage <- (df_pop_agg$age_group_years_end + df_pop_agg$age_group_years_start)/2
  df_pop_agg[,c("age_group_years_start", "age_group_years_end")] <- NULL
  
  df_pop_agg_m <- df_pop_agg[df_pop_agg$sex_id==1, ]
  df_pop_agg_f <- df_pop_agg[df_pop_agg$sex_id==2, ]
  df_pop_agg_b <- df_pop_agg[df_pop_agg$sex_id==3, ]
  
  pop_m <- df_pop_agg_m$pop_size
  pop_f <- df_pop_agg_f$pop_size
  pop_b <- df_pop_agg_b$pop_size
  
  
  #Add midage variable to dt
  dt <- merge(dt, ages_2019[,c("age_group_id", "age_group_years_start", "age_group_years_end", "age_group_name")], by = "age_group_id")
  dt$midage <- (dt$age_group_years_end + dt$age_group_years_start)/2
  dt[,c("age_group_years_start", "age_group_years_end")] <- NULL
  
  
  ## Run regression models after logit-transforming data
  dt <- as.data.frame(dt)
  cols <- c(paste0("draw_",0:999))
  dt[cols] <- gtools::logit(dt[cols])
  
  dt_m <- dt[dt$sex_id==1, ]
  dt_f <- dt[dt$sex_id==2, ]
  
  fits_m <- lapply(0:999, function(i) {
    #i <- 1 # dev
    df_tmp <- as.data.frame(dt_m)[, c("region_name", "year_id", "midage", paste0("draw_", i))]
    names(df_tmp)[4] <- "draw_i"
    lm(draw_i ~ year_id * region_name + bs(midage,knots = c(27,50,75)), data = df_tmp) #interaction term region and year, cubic spline on age (midage)
  })
  
  fits_f <- lapply(0:999, function(i) {
    #i <- 1 # dev
    df_tmp <- as.data.frame(dt_f)[, c("region_name", "year_id", "midage", paste0("draw_", i))]
    names(df_tmp)[4] <- "draw_i"
    lm(draw_i ~ year_id * region_name + bs(midage,knots = c(27,50,75)), data = df_tmp) #interaction term region and year, cubic spline on age (midage)
  })
  
  ## Predict out rates using fits
  preds_m <- lapply(fits_m, function(x) {
    # x <- fits[[1]] # dev
    predict(x, newdata = df_pop_agg_m)
  })
  
  preds_f <- lapply(fits_f, function(x) {
    # x <- fits[[1]] # dev
    predict(x, newdata = df_pop_agg_f)
  })
  
  
  ## Combine predictions, back to normal space
  dat_pred_m <- do.call("cbind", preds_m)
  dat_pred_m <- gtools::inv.logit(dat_pred_m)
  dat_pred_m <- as.data.table(dat_pred_m)
  dat_pred_f <- do.call("cbind", preds_f)
  dat_pred_f <- gtools::inv.logit(dat_pred_f)
  dat_pred_f <- as.data.table(dat_pred_f)
  
  #Add region and age variables back to df
  dat_pred_m$region_name <- df_pop_agg_m$region_name 
  dat_pred_m$age_group_name <- df_pop_agg_m$age_group_name
  dat_pred_m$age_group_id <- df_pop_agg_m$age_group_id
  dat_pred_f$region_name <- df_pop_agg_f$region_name 
  dat_pred_f$age_group_name <- df_pop_agg_f$age_group_name
  dat_pred_f$age_group_id <- df_pop_agg_f$age_group_id
  
  region_agg_m <- as.data.frame(dat_pred_m)
  region_agg_f <- as.data.frame(dat_pred_f)
  moveMe <- function(data, tomove, where = "last", ba = NULL) {
    temp <- setdiff(names(data), tomove)
    x <- switch(
      where,
      first = data[c(tomove, temp)],
      last = data[c(temp, tomove)],
      before = {
        if (is.null(ba)) stop("must specify ba column")
        if (length(ba) > 1) stop("ba must be a single character string")
        data[append(temp, values = tomove, after = (match(ba, temp)-1))]
      },
      after = {
        if (is.null(ba)) stop("must specify ba column")
        if (length(ba) > 1) stop("ba must be a single character string")
        data[append(temp, values = tomove, after = (match(ba, temp)))]
      })
    x
  }
  
  region_agg_m <- moveMe(region_agg_m, "region_name", "before", "V1")
  region_agg_m <- moveMe(region_agg_m, "age_group_name", "before", "V1")
  region_agg_m <- moveMe(region_agg_m, "age_group_id", "before", "V1")
  region_agg_f <- moveMe(region_agg_f, "region_name", "before", "V1")
  region_agg_f <- moveMe(region_agg_f, "age_group_name", "before", "V1")
  region_agg_f <- moveMe(region_agg_f, "age_group_id", "before", "V1")
  
  #Subset age-specific rate draws as separate df
  rates_m <- as.data.table(region_agg_m)
  rates_f <- as.data.table(region_agg_f)
  ## Calculate age-specific forecasted case number
  region_agg_m <- as.data.table(region_agg_m)
  region_agg_m[, paste0("draw_", 0:999) := lapply(1:1000, function(x) get(paste0("V", x)) *pop_m)]
  region_agg_m[, c(paste0("V", 1:1000))] <- NULL
  region_agg_f <- as.data.table(region_agg_f)
  region_agg_f[, paste0("draw_", 0:999) := lapply(1:1000, function(x) get(paste0("V", x)) *pop_f)]
  region_agg_f[, c(paste0("V", 1:1000))] <- NULL
  
  #Create rates for both sex by sex-aggregating age-specific cases and then dividing by population
  rates_b <- plyr::rbind.fill(region_agg_m, region_agg_f)
  setDT(rates_b)
  rates_b <- rates_b[, lapply(.SD, sum, na.rm=TRUE), by=c("region_name", "age_group_id", "age_group_name") ]
  rates_b <- merge(rates_b, df_pop_agg_b[, c("region_name", "age_group_id", "age_group_name", "pop_size")], by = c("region_name", "age_group_id", "age_group_name"))
  rates_b[, paste0("V", 1:1000) := lapply(0:999, function(x) get(paste0("draw_", x)) /pop_size)]
  rates_b[, c(paste0("draw_", 0:999))] <- NULL
  
  ## Age-standardized rates by region
  calc_agestd_rates <- function(data_dt){
    draws <- paste0("V", 1:1000)
    dt_tmp <- copy(data_dt)
    age_weights <- get_age_metadata(12, gbd_round_id=6)
    dt_tmp <- merge(dt_tmp, age_weights[, c("age_group_id", "age_group_weight_value")], by = c("age_group_id"), all.x=T)
    dt_tmp[, (draws) := lapply(.SD, function(x) sum(x * age_group_weight_value)), by = "region_name", .SDcols = draws]
    dt_tmp <- unique(dt_tmp, by = c("region_name"))
    dt_tmp[,`:=` (age_group_weight_value = NULL, age_group_id = 27)]
    dt_tmp[, age_group_name := "Age-standardized"] 
    dt_tmp[, mean := rowMeans(.SD), .SDcols = paste0("V", 1:1000)]
    dt_tmp$upper <- apply(dt_tmp[,c(paste0("V", 1:1000))], 1, quantile, probs = 0.975, na.rm = TRUE)
    dt_tmp$lower <- apply(dt_tmp[,c(paste0("V", 1:1000))], 1, quantile, probs = 0.025, na.rm = TRUE)
    dt_tmp[, c(paste0("V", 1:1000))] <- NULL
    return(dt_tmp)
  }
  
  agestdrate_m <- calc_agestd_rates(rates_m)
  agestdrate_m$sex <- "Male"
  agestdrate_f <- calc_agestd_rates(rates_f)
  agestdrate_f$sex <- "Female"
  agestdrate_b <- calc_agestd_rates(rates_b)
  agestdrate_b$sex <- "Both"
  
  agestdrate_df <- plyr::rbind.fill(agestdrate_m, agestdrate_f, agestdrate_b)
  
  
  ## Age-aggregate cases from each region 
  mylist <- split(region_agg_m, region_agg_m$region_name)
  for(i in 1:7){
    #i <- 1 # dev
    df_tmp <- as.data.frame(mylist[[i]])
    name <- as.character(df_tmp[1, 1])
    df_tmp <- as.data.table(df_tmp %>% adorn_totals("row"))
    df_tmp <- df_tmp[region_name=="Total",]
    df_tmp$region_name <- gsub('Total', name, df_tmp$region_name)
    mylist[[i]] <- df_tmp
  }  
  
  region_collapse_m <- do.call(rbind,mylist)
  region_collapse_m <- region_collapse_m[region_collapse_m$age_group_id!=0,]
  
  mylist <- split(region_agg_f, region_agg_f$region_name)
  for(i in 1:7){
    #i <- 1 # dev
    df_tmp <- as.data.frame(mylist[[i]])
    name <- as.character(df_tmp[1, 1])
    df_tmp <- as.data.table(df_tmp %>% adorn_totals("row"))
    df_tmp <- df_tmp[region_name=="Total",]
    df_tmp$region_name <- gsub('Total', name, df_tmp$region_name)
    mylist[[i]] <- df_tmp
  }  
  
  region_collapse_f <- do.call(rbind,mylist)
  region_collapse_f <- region_collapse_f[region_collapse_f$age_group_id!=0,]
  
  
  
  ## Add in row for global total and collapse cases to regional and global mean and UIs
  region_collapse_m <- as.data.table(region_collapse_m %>% adorn_totals("row"))
  region_collapse_m$region_name <- gsub("Total", "Global", region_collapse_m$region_name)
  region_collapse_m$sex <- "Male"
  
  region_collapse_f <- as.data.table(region_collapse_f %>% adorn_totals("row"))
  region_collapse_f$region_name <- gsub("Total", "Global", region_collapse_f$region_name)
  region_collapse_f$sex <- "Female"
  
  region_collapse_b <- plyr::rbind.fill(region_collapse_m, region_collapse_f)
  setDT(region_collapse_b)
  region_collapse_b[,c("age_group_id", "age_group_name", "sex")] <- NULL
  region_collapse_b <- region_collapse_b[, lapply(.SD, sum, na.rm=TRUE), by=c("region_name") ]
  region_collapse_b$sex <- "Both"
  
  #Finalize - cases
  final_dt <- plyr::rbind.fill(region_collapse_m, region_collapse_f, region_collapse_b)
  setDT(final_dt)
  final_dt$age_group_name <- "All age"
  final_dt$age_group_id <- 22
  final_dt[, mean := rowMeans(.SD), .SDcols = paste0("draw_", 0:999)]
  final_dt$upper <- apply(final_dt[,c(paste0("draw_", 0:999))], 1, quantile, probs = 0.975, na.rm = TRUE)
  final_dt$lower <- apply(final_dt[,c(paste0("draw_", 0:999))], 1, quantile, probs = 0.025, na.rm = TRUE)
  final_dt[, c(paste0("draw_", 0:999))] <- NULL
  final_dt$year_id <- year_forecast
  final_dt$cause <- severity
  
  #Finalize - age-standardized
  agestdrate_df$age_group_name <- "Age-standardized"
  agestdrate_df$age_group_id <- 27
  agestdrate_df$year_id <- year_forecast
  agestdrate_df$cause <- severity
  agestdrate_df$pop_size <- NULL
  
  #Combine 
  combo_dt <- plyr::rbind.fill(final_dt, agestdrate_df)
  
  ## Save output
  output_path_reg <- paste0(my_dir, "FILEPATH")
  write.xlsx(combo_dt, paste0(output_path_reg, year_forecast, "FILEPATH"))
}
## Final combined data
df1 <- read.xlsx(paste0(my_dir, "FILEPATH"))
df2 <- read.xlsx(paste0(my_dir, "FILEPATH"))
df3 <- read.xlsx(paste0(my_dir, "FILEPATH"))
df4 <- read.xlsx(paste0(my_dir, "FILEPATH"))
merge_years <- plyr::rbind.fill(df1,df2,df3,df4)
write.xlsx(merge_years, paste0(my_dir, "FILEPATH"))
