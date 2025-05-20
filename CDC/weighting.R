# Do weighting for any NEISS injury types with 2+ GBD n-code matches

rm(list=ls())

library(data.table)
library(dplyr)

invisible(sapply(list.files("FILEPATH", full.names = T), source))

ncodes <- get_rei_metadata(rei_set_id = 7, release_id = 9)
detailed <- ncodes[level == 2]

draws <- paste0("draw_", 0:499)  
mean_ui <- function (draw_dt) {
  dt <- copy(draw_dt)
  dt[, val := rowMeans(.SD), .SDcols = draws]
  dt[, lower := apply(.SD, 1, quantile, probs = 0.025), .SDcols = draws]
  dt[, upper := apply(.SD, 1, quantile, probs = 0.975), .SDcols = draws]
  dt[, (draws) := NULL]
  return(dt)
} 

# dt <- data.table()
# 
# # Get outputs by ncode for years 2015-2020 
# # Incidence, all-ages, both-sex, rates per 1 of n-code injuries
# for (i in detailed$rei_id){
#   print(i)
#   print(paste0("start: ", Sys.time()))
#   outputs <- get_draws(gbd_id_type = c("cause_id","rei_id"),
#                        gbd_id = c(294,i),
#                        location_id = 102,
#                        sex_id = 3,
#                        year_id = c(2015:2020),
#                        age_group_id = 22,
#                        measure_id = 6,
#                        metric_id = 3,
#                        source = "como",
#                        version_id = 1471,
#                        release_id = 9,
#                        num_workers = 4)
# 
#   dt <- rbind(dt, outputs)
#   print(paste0("done: ", Sys.time()))
# }

#fwrite(dt, "FILEPATH")
dt <- fread("FILEPATH")

dt <- merge(dt, ncodes[,.(rei_id,rei,rei_name)], by = c("rei_id"))

# Get corresponding population to use to convert to counts
population <- get_population(age_group_id = 22,
                             location_id = 102,
                             year_id = c(2015:2020),
                             sex_id = 3,
                             release_id = 9)

counts <- merge(dt, population, by = c("year_id","location_id","age_group_id","sex_id"))

# Do draw-level computation to convert to n-code incidence counts summed by year
pattern <- "^draw_\\d+$"
draw_cols <- grep(pattern, names(dt), value = TRUE)
counts <- counts[, (draw_cols) := lapply(.SD, function(x) x * population), .SDcols = draw_cols]
counts_sum <- counts[, lapply(.SD, sum), by = .(rei,rei_id,rei_name), .SDcols = draw_cols]

counts_sum <- merge(counts_sum, ncodes[,.(rei_id,sort_order)], by = "rei_id")
counts_sum_table <- mean_ui(counts_sum)
fwrite(counts_sum_table, "FILEPATH")

# Get map of NEISS-AIP injury types to GBD n-codes (S1 for paper)
# "Yellow" matches need custom weighted DWs
my_map <- read.xlsx("FILEPATH")
setDT(my_map)

yellow <- my_map[color == "Yellow"]
yellow_uniques <- yellow[,.(bdyptg_c,bdypt,diag,`GBD.ncode-1`,`ncode-1.name`,`GBD.ncode-2`,`ncode-2.name`,`GBD.ncode-3`,`ncode-3.name`,`GBD.ncode-4`,`ncode-4.name`,
                            `GBD.ncode-5`,`ncode-5.name`,`GBD.ncode-6`,`ncode-6.name`,`GBD.ncode-7`,`ncode-7.name`,`GBD.ncode-8`,`ncode-8.name`,`GBD.ncode-9`,`ncode-9.name`,`GBD.ncode-10`,`ncode-10.name`)]
colnames(yellow_uniques) <- c("bdyptg_c","bdypt","diag","ncode1","ncode1_name","ncode2","ncode2_name","ncode3","ncode3_name","ncode4","ncode4_name","ncode5","ncode5_name",
                              "ncode6","ncode6_name","ncode7","ncode7_name","ncode8","ncode8_name","ncode9","ncode9_name","ncode10","ncode10_name")
yellow_uniques <- yellow_uniques[order(yellow_uniques$ncode1_name),]
yellow_uniques <- yellow_uniques[diag != "Other"]

dws <- fread("FILEPATH", na.strings = "")

# Drop any "UNTRT" dws
dws <- dws[treatment == "TRT" | is.na(treatment)]

# If both ST & LT present, only keep ST
dws <- dws %>%
  group_by(ncode) %>%
  filter(duration == max(duration) | is.na(duration)) %>%
  ungroup()
setDT(dws)

# Do two matches first
two_matches <- yellow_uniques[is.na(ncode3)]
two_matches <- two_matches[, (draws) := lapply(1:500, function(x) NA)]
two_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

for (i in 1:nrow(two_matches)){
  
  dw1 <- dws[ncode == trimws(two_matches[i]$ncode1)]
  dw1_draws <- dw1[,6:505]
  dw2 <- dws[ncode == trimws(two_matches[i]$ncode2)]
  dw2_draws <- dw2[,6:505]
  
  counts1 <- counts_sum[rei == trimws(two_matches[i]$ncode1)]
  counts1_draws <- counts1[,4:503]
  counts2 <- counts_sum[rei == trimws(two_matches[i]$ncode2)]
  counts2_draws <- counts2[,4:503]
  
  weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws)/(counts1_draws + counts2_draws)
  
  two_matches[i, (draws) := weighted_avg]
  
  print(paste0("done: ",i))
}

# Do three matches 
three_matches <- yellow_uniques[is.na(ncode4) & !is.na(ncode3)]
three_matches <- three_matches[, (draws) := lapply(1:500, function(x) NA)]
three_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

for (i in 1:nrow(three_matches)){
  
  dw1 <- dws[ncode == trimws(three_matches[i]$ncode1)]
  dw1_draws <- dw1[,6:505]
  dw2 <- dws[ncode == trimws(three_matches[i]$ncode2)]
  dw2_draws <- dw2[,6:505]
  dw3 <- dws[ncode == trimws(three_matches[i]$ncode3)]
  dw3_draws <- dw3[,6:505]
  
  counts1 <- counts_sum[rei == trimws(three_matches[i]$ncode1)]
  counts1_draws <- counts1[,4:503]
  counts2 <- counts_sum[rei == trimws(three_matches[i]$ncode2)]
  counts2_draws <- counts2[,4:503]
  counts3 <- counts_sum[rei == trimws(three_matches[i]$ncode3)]
  counts3_draws <- counts3[,4:503]
  
  weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws + dw3_draws*counts3_draws)/(counts1_draws + counts2_draws + counts3_draws)
  
  three_matches[i, (draws) := weighted_avg]
  
  print(paste0("done: ",i))
}

# Do four matches 
four_matches <- yellow_uniques[is.na(ncode5) & !is.na(ncode4)]
four_matches <- four_matches[, (draws) := lapply(1:500, function(x) NA)]
four_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

for (i in 1:nrow(four_matches)){
  
  dw1 <- dws[ncode == trimws(four_matches[i]$ncode1)]
  dw1_draws <- dw1[,6:505]
  dw2 <- dws[ncode == trimws(four_matches[i]$ncode2)]
  dw2_draws <- dw2[,6:505]
  dw3 <- dws[ncode == trimws(four_matches[i]$ncode3)]
  dw3_draws <- dw3[,6:505]
  dw4 <- dws[ncode == trimws(four_matches[i]$ncode4)]
  dw4_draws <- dw4[,6:505]
  
  counts1 <- counts_sum[rei == trimws(four_matches[i]$ncode1)]
  counts1_draws <- counts1[,4:503]
  counts2 <- counts_sum[rei == trimws(four_matches[i]$ncode2)]
  counts2_draws <- counts2[,4:503]
  counts3 <- counts_sum[rei == trimws(four_matches[i]$ncode3)]
  counts3_draws <- counts3[,4:503]
  counts4 <- counts_sum[rei == trimws(four_matches[i]$ncode4)]
  counts4_draws <- counts4[,4:503]
  
  weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws + dw3_draws*counts3_draws + dw4_draws*counts4_draws)/(counts1_draws + counts2_draws + counts3_draws + counts4_draws)
  
  four_matches[i, (draws) := weighted_avg]
  
  print(paste0("done: ",i))
}

# Do five matches
five_matches <- yellow_uniques[is.na(ncode6) & !is.na(ncode5)]
five_matches <- five_matches[, (draws) := lapply(1:500, function(x) NA)]
five_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

dw1 <- dws[ncode == trimws(five_matches$ncode1)]
dw1_draws <- dw1[,6:505]
dw2 <- dws[ncode == trimws(five_matches$ncode2)]
dw2_draws <- dw2[,6:505]
dw3 <- dws[ncode == trimws(five_matches$ncode3)]
dw3_draws <- dw3[,6:505]
dw4 <- dws[ncode == trimws(five_matches$ncode4)]
dw4_draws <- dw4[,6:505]
dw5 <- dws[ncode == trimws(five_matches$ncode5)]
dw5_draws <- dw5[,6:505]

counts1 <- counts_sum[rei == trimws(five_matches$ncode1)]
counts1_draws <- counts1[,4:503]
counts2 <- counts_sum[rei == trimws(five_matches$ncode2)]
counts2_draws <- counts2[,4:503]
counts3 <- counts_sum[rei == trimws(five_matches$ncode3)]
counts3_draws <- counts3[,4:503]
counts4 <- counts_sum[rei == trimws(five_matches$ncode4)]
counts4_draws <- counts4[,4:503]
counts5 <- counts_sum[rei == trimws(five_matches$ncode5)]
counts5_draws <- counts5[,4:503]

weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws + dw3_draws*counts3_draws + dw4_draws*counts4_draws + dw5_draws*counts5_draws)/(counts1_draws + counts2_draws + counts3_draws + counts4_draws + counts5_draws)

five_matches[, (draws) := weighted_avg]

# Do six matches
six_matches <- yellow_uniques[is.na(ncode7) & !is.na(ncode6)]
six_matches <- six_matches[, (draws) := lapply(1:500, function(x) NA)]
six_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

dw1 <- dws[ncode == trimws(six_matches$ncode1)]
dw1_draws <- dw1[,6:505]
dw2 <- dws[ncode == trimws(six_matches$ncode2)]
dw2_draws <- dw2[,6:505]
dw3 <- dws[ncode == trimws(six_matches$ncode3)]
dw3_draws <- dw3[,6:505]
dw4 <- dws[ncode == trimws(six_matches$ncode4)]
dw4_draws <- dw4[,6:505]
dw5 <- dws[ncode == trimws(six_matches$ncode5)]
dw5_draws <- dw5[,6:505]
dw6 <- dws[ncode == trimws(six_matches$ncode6)]
dw6_draws <- dw6[,6:505]

counts1 <- counts_sum[rei == trimws(six_matches$ncode1)]
counts1_draws <- counts1[,4:503]
counts2 <- counts_sum[rei == trimws(six_matches$ncode2)]
counts2_draws <- counts2[,4:503]
counts3 <- counts_sum[rei == trimws(six_matches$ncode3)]
counts3_draws <- counts3[,4:503]
counts4 <- counts_sum[rei == trimws(six_matches$ncode4)]
counts4_draws <- counts4[,4:503]
counts5 <- counts_sum[rei == trimws(six_matches$ncode5)]
counts5_draws <- counts5[,4:503]
counts6 <- counts_sum[rei == trimws(six_matches$ncode6)]
counts6_draws <- counts6[,4:503]

weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws + dw3_draws*counts3_draws + dw4_draws*counts4_draws + dw5_draws*counts5_draws + dw6_draws*counts6_draws)/(counts1_draws + counts2_draws + counts3_draws + counts4_draws + counts5_draws + counts6_draws)

six_matches[, (draws) := weighted_avg]

# Do seven matches
seven_matches <- yellow_uniques[is.na(ncode8) & !is.na(ncode7)]
seven_matches <- seven_matches[, (draws) := lapply(1:500, function(x) NA)]
seven_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

dw1 <- dws[ncode == trimws(seven_matches$ncode1)]
dw1_draws <- dw1[,6:505]
dw2 <- dws[ncode == trimws(seven_matches$ncode2)]
dw2_draws <- dw2[,6:505]
dw3 <- dws[ncode == trimws(seven_matches$ncode3)]
dw3_draws <- dw3[,6:505]
dw4 <- dws[ncode == trimws(seven_matches$ncode4)]
dw4_draws <- dw4[,6:505]
dw5 <- dws[ncode == trimws(seven_matches$ncode5)]
dw5_draws <- dw5[,6:505]
dw6 <- dws[ncode == trimws(seven_matches$ncode6)]
dw6_draws <- dw6[,6:505]
dw7 <- dws[ncode == trimws(seven_matches$ncode7)]
dw7_draws <- dw7[,6:505]

counts1 <- counts_sum[rei == trimws(seven_matches$ncode1)]
counts1_draws <- counts1[,4:503]
counts2 <- counts_sum[rei == trimws(seven_matches$ncode2)]
counts2_draws <- counts2[,4:503]
counts3 <- counts_sum[rei == trimws(seven_matches$ncode3)]
counts3_draws <- counts3[,4:503]
counts4 <- counts_sum[rei == trimws(seven_matches$ncode4)]
counts4_draws <- counts4[,4:503]
counts5 <- counts_sum[rei == trimws(seven_matches$ncode5)]
counts5_draws <- counts5[,4:503]
counts6 <- counts_sum[rei == trimws(seven_matches$ncode6)]
counts6_draws <- counts6[,4:503]
counts7 <- counts_sum[rei == trimws(seven_matches$ncode7)]
counts7_draws <- counts7[,4:503]

weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws + dw3_draws*counts3_draws + dw4_draws*counts4_draws + dw5_draws*counts5_draws + dw6_draws*counts6_draws + dw7_draws*counts7_draws)/(counts1_draws + counts2_draws + counts3_draws + counts4_draws + counts5_draws + counts6_draws + counts7_draws)

seven_matches[, (draws) := weighted_avg]

# Do ten matches
ten_matches <- yellow_uniques[!is.na(ncode10)]
ten_matches <- ten_matches[, (draws) := lapply(1:500, function(x) NA)]
ten_matches[, (draws) := lapply(.SD, as.numeric), .SDcols = draws]

dw1 <- dws[ncode == trimws(ten_matches$ncode1)]
dw1_draws <- dw1[,6:505]
dw2 <- dws[ncode == trimws(ten_matches$ncode2)]
dw2_draws <- dw2[,6:505]
dw3 <- dws[ncode == trimws(ten_matches$ncode3)]
dw3_draws <- dw3[,6:505]
dw4 <- dws[ncode == trimws(ten_matches$ncode4)]
dw4_draws <- dw4[,6:505]
dw5 <- dws[ncode == trimws(ten_matches$ncode5)]
dw5_draws <- dw5[,6:505]
dw6 <- dws[ncode == trimws(ten_matches$ncode6)]
dw6_draws <- dw6[,6:505]
dw7 <- dws[ncode == trimws(ten_matches$ncode7)]
dw7_draws <- dw7[,6:505]
dw8 <- dws[ncode == trimws(ten_matches$ncode8)]
dw8_draws <- dw8[,6:505]
dw9 <- dws[ncode == trimws(ten_matches$ncode9)]
dw9_draws <- dw9[,6:505]
dw10 <- dws[ncode == trimws(ten_matches$ncode10)]
dw10_draws <- dw10[,6:505]

counts1 <- counts_sum[rei == trimws(ten_matches$ncode1)]
counts1_draws <- counts1[,4:503]
counts2 <- counts_sum[rei == trimws(ten_matches$ncode2)]
counts2_draws <- counts2[,4:503]
counts3 <- counts_sum[rei == trimws(ten_matches$ncode3)]
counts3_draws <- counts3[,4:503]
counts4 <- counts_sum[rei == trimws(ten_matches$ncode4)]
counts4_draws <- counts4[,4:503]
counts5 <- counts_sum[rei == trimws(ten_matches$ncode5)]
counts5_draws <- counts5[,4:503]
counts6 <- counts_sum[rei == trimws(ten_matches$ncode6)]
counts6_draws <- counts6[,4:503]
counts7 <- counts_sum[rei == trimws(ten_matches$ncode7)]
counts7_draws <- counts7[,4:503]
counts8 <- counts_sum[rei == trimws(ten_matches$ncode8)]
counts8_draws <- counts8[,4:503]
counts9 <- counts_sum[rei == trimws(ten_matches$ncode9)]
counts9_draws <- counts9[,4:503]
counts10 <- counts_sum[rei == trimws(ten_matches$ncode10)]
counts10_draws <- counts10[,4:503]

weighted_avg <- (dw1_draws*counts1_draws + dw2_draws*counts2_draws + dw3_draws*counts3_draws + dw4_draws*counts4_draws + dw5_draws*counts5_draws + dw6_draws*counts6_draws + dw7_draws*counts7_draws + dw8_draws*counts8_draws + dw9_draws*counts9_draws + dw10_draws*counts10_draws)/(counts1_draws + counts2_draws + counts3_draws + counts4_draws + counts5_draws + counts6_draws + counts7_draws + counts8_draws + counts9_draws + counts10_draws)

ten_matches[, (draws) := weighted_avg]

# Combine all!
all_dws <- rbindlist(list(two_matches,three_matches,four_matches,five_matches,six_matches,seven_matches,ten_matches))

# Summarize!
weighted_dws <- mean_ui(all_dws)
fwrite(weighted_dws, "FILEPATH")