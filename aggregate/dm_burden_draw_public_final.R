################################################################################
# GBD diabetes paper code - pulling YLL, YLD, deaths, prev, dalys
# This code pulls and aggregates the data that made the estimates. 
################################################################################


#path to save docs
path<- filepath

# ---Load functions and arguments-----------------------------------------------

source("address/get_draws.R")
source("address/get_population.R")
source("address/get_age_metadata.R")
source("address/get_demographics.R")
source("address/get_rei_metadata.R")

rei<-get_rei_metadata(release_id = 9, rei_set_id = 2)
rei<-rei[,c('rei_id','parent_id','most_detailed','lancet_label')]
age<-get_age_metadata(release_id = 10)
age1<-age[,c('age_group_name','age_group_weight_value','age_group_id')]

demographics <- get_demographics(gbd_team="epi", gbd_round_id=7)
ages <- unlist(demographics$age_group_id, use.names=F)

#Arguments
argue <- commandArgs(trailingOnly = T)
locations <- as.numeric(argue[1])
loc_launch<-locations

#settings
###demographics
year_id<-c(1990,2010,2021)
cause<-c(587,975,976)
sex_id<-c(1,2)
metric_id<-c(1,3)
measure_id<-c(1,2,3,4,5)
age_group_id <- c(2,3,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,34,235,238,388,389)

#versions
version<-250
como<-1353
codcorrect<-330
compare<-7851
paf<-445

# ---Create dataframe-----------------------------------------------------------

#step1: pulling data

#population #using run_id #348
pop<-get_population(location_id =loc_launch,age_group_id = age_group_id,year_id=year_id,sex_id=sex_id,with_ui=TRUE, release_id=9 )
pop<-pop[,c('age_group_id','location_id','year_id','sex_id','population')]

#YLDs/Prevalence
df_burden1 <- get_draws("cause_id", cause, age_group_id = age_group_id,location_id=loc_launch, year_id=year_id, source="como", 
                version_id=como, release_id = 9, measure_id = measure_id, metric_id = 3,sex_id=sex_id) #sex_specific

#YLLs/Deaths
df_burden2 <- get_draws("cause_id", cause, age_group_id = age_group_id,location_id=loc_launch, year_id=year_id, source="codcorrect", 
                        version_id=codcorrect, release_id = 9, measure_id = measure_id, metric_id = 1,sex_id=sex_id) #sex_specific

#add files to get yll/yld/prev/death together
df_burden<-rbind(df_burden1, df_burden2)

#reshape files to create summary measures
df_burden_l<- melt(setDT(df_burden), id.vars = c("sex_id","year_id","metric_id",'version_id',
                                                     'location_id','measure_id','age_group_id','cause_id'), variable.name = "draw")

#step 2: ensure both numbers and rates are generated
#create the corresponding value
df_burden_long<-merge(df_burden_l,pop,by=c('age_group_id','sex_id','year_id','location_id'),all.x=T)
df_burden_long$updated_num<- ifelse(df_burden_long$metric_id==1, df_burden_long$value/df_burden_long$population, 
                                    df_burden_long$value*df_burden_long$population)

df_fix<-df_burden_long[,c('age_group_id','sex_id','year_id','location_id',"metric_id",'version_id','measure_id','cause_id','draw',
                          'population','updated_num')]
df_fix$metric_id1<-ifelse(df_fix$metric_id==1,3,1) #update the metric version
df_fix<-df_fix[ ,-c("metric_id")]
colnames(df_fix)[colnames(df_fix) == "updated_num"] ="value"
colnames(df_fix)[colnames(df_fix) == "metric_id1"] ="metric_id"
df_burden_long<-df_burden_long[ ,-c("updated_num")]

df<-rbind(df_burden_long, df_fix)

#step 3: Create both sex estimates
sex<-df[df$metric_id==1,]

sex_both<-sex[,j=list(value=sum(value),population=sum(population)), by=c("year_id","metric_id",'draw','version_id',
                                         'location_id','measure_id','age_group_id','cause_id')]
sex_both$sex_id<-3
sex_both_num<-copy(sex_both)
sex_both$values<-sex_both$value/sex_both$population
sex_both$value<-sex_both$values
sex_both<-sex_both[ ,-c("values")]
sex_both$metric_id<-3
sex_both_rate<-copy(sex_both)
df<-rbind(df,sex_both_num,sex_both_rate )

#step 4: create DALYS (sum YLL+YLD)
t<-df[df$measure_id %in% c(3,4) & df$metric_id==1,]

t.2<-t[,j=list(value=sum(value)), by=c("sex_id","year_id","metric_id",'draw','population',
                                         'location_id','age_group_id','cause_id')]
t.2$measure_id<-2
t.2$version_id<-9999
daly_num<-copy(t.2)

#rate (sum YLL+YLD)
t.2$values<-t.2$value/t.2$population
t.2$value<-t.2$values
t.2<-t.2[ ,-c("values")]
t.2$metric_id<-3
daly_rate<-copy(t.2)
df<-rbind(df,daly_num,daly_rate)

#STEP 5: create all age numbers
df_all_age<-df[df$metric_id==1,]
all_age<-df_all_age[,j=list(value=sum(value)),by=c('sex_id','year_id','location_id','metric_id','version_id','measure_id','cause_id',
                                           'draw')]
all_age$age_group_id<-22
all_age$population<-NA
df<-rbind(df,all_age)

#STEP 6: create age-std rates
age_std<-df[df$metric_id==3 & df$age_group_id!=22,]
age_std<-merge(age_std, age1, by=c('age_group_id'),all.x=T)
age_std$rate<-age_std$value*age_std$age_group_weight_value
age_std_rate<-age_std[,j=list(value=sum(rate)), by=c("sex_id","year_id","metric_id",'draw','location_id',
                                                     'cause_id','measure_id','version_id')]
age_std_rate$age_group_id<-27
age_std_rate$population<-9999

#STEP 7: create <5 year groups
age_u5<-df[df$metric_id==3 & df$age_group_id %in% c(2,3,34,238,388,389),]
age_u5<-merge(age_u5, age1, by=c('age_group_id'),all.x=T)
age_u5$rate<-age_u5$value*age_u5$age_group_weight_value
age_u5_rate<-age_u5[,j=list(value=sum(rate)), by=c("sex_id","year_id","metric_id",'draw','location_id',
                                                     'cause_id','measure_id','version_id')]
age_u5_rate$age_group_id<-9999
age_u5_rate$population<-9999

#Bind age std and <5 years to main file
df<-rbind(df,age_std_rate,age_u5_rate)

#Write files
write.csv(df, paste0(path,'draws/disease/burden_',loc_launch,'.csv'))


# ---ATTRIBUTABLE DALYS---------------------------------------------------------

burden <- get_draws("rei_id", 
                    c(85,86,87, 98, 99, 100, 102, 104, 105, 108, 110, 111, 112, 113, 116, 117, 118, 
                      119,125, 169, 202, 203, 331, 337 ,338 ,380), 
                    age_group_id = c( 2,3,6,7,8,9,10,11,12, 13, 14, 15, 16, 17, 18, 19, 20,30 , 31, 
                                      32, 34, 235, 238, 388, 389),
                    location_id=loc_launch, 
                    year_id=year_id, source="burdenator", 
                    version_id=250, release_id = 9, measure_id = 2, metric_id =c(1,2), sex_id =c(1,2))

#reshape files to create summary measures
burden_l<- melt(setDT(burden), id.vars = c("sex_id","year_id","metric_id",'version_id','rei_id',
                                           'location_id','measure_id','age_group_id','cause_id'), variable.name = "draw")

#subset to just diabetes
burden_l<-burden_l[burden_l$cause_id %in% cause,]
burden_l<-merge(burden_l, rei, by=c('rei_id'),all.x=T)

# ---calculate dalys------------------------------------------------------------

#bring in number of dalys
dalys<-df[df$measure_id==2 & df$metric_id==1,]

#Custom PAF: DALYS attributable to risk factors other than FPG
detailed<-burden_l[(burden_l$most_detailed==1|burden_l$rei_id==108) & burden_l$metric_id==2 & burden_l$rei_id!=105,]
paf_agg<-detailed[,j=list(value=(1-value)),  by=c('age_group_id','rei_id','location_id','sex_id','metric_id','draw',
                                                  'cause_id','year_id','lancet_label','measure_id')]
paf_agg1<-setDT(paf_agg)[, .(value = prod(value)), .(age_group_id,location_id,cause_id,year_id,sex_id,measure_id,metric_id,draw)]
paf_agg1$risk_val<-1-paf_agg1$value

#larger groups to do risk specific grouping of paf
groups<-burden_l[burden_l$lancet_label %in% c('Environmental/occupational risks','Dietary risks','Tobacco') 
                 & burden_l$metric_id==2 ,]


#dalys due to RF (total)
daly_rf<-merge(dalys, paf_agg1,by=c('age_group_id','sex_id','year_id','location_id','cause_id','draw','measure_id'),all.y=T)
daly_rf$value<-daly_rf$risk_val*daly_rf$value.x
daly_rf_sum<-daly_rf[,j=list(value=sum(value)),
                      by=c("year_id",'location_id','cause_id','measure_id','draw')]
daly_rf_sum$type<-'num_dalys_d2_rf'
daly_rf_sum1<-daly_rf_sum[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)),
                     by=c("year_id", 'location_id','cause_id','measure_id','type')] 

# % of dalys due to rf
daly_rf_pct<-daly_rf[,j=list(value=sum(value)/sum(value.x)),
                     by=c("year_id",'location_id','cause_id','measure_id','draw')]
daly_rf_pct$type<-'pct_dalys_d2_rf'
daly_rf_pct1<-daly_rf_pct[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)),
                          by=c("year_id", 'location_id','cause_id','measure_id','type')] 

# detailed level rf attribution
dalys<-dalys[dalys$cause_id == 976 & dalys$age_group_id %in% c(8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,235) & 
               dalys$sex_id %in% c(1,2) & dalys$metric_id==1,]
detailed<-detailed[detailed$cause_id == 976,]
rf_detailed<-merge(dalys, detailed,by=c('age_group_id','sex_id','year_id','location_id','cause_id','draw','measure_id'),all.y=T)
rf_groups<-merge(dalys, groups,by=c('age_group_id','sex_id','year_id','location_id','cause_id','draw','measure_id'),all.y=T)
rf_groups<-rf_groups[rf_groups$cause_id == 976 & rf_groups$age_group_id %in% c(8,9,10,11,12,13,14,15,16,17,18,19,20,30,31,32,235)]
rf_detailed<-rbind(rf_detailed, rf_groups)
rf_detailed$value<-rf_detailed$value.y*rf_detailed$value.x
rf_detailed_calc<-rf_detailed[,j=list(value=sum(value),value.daly=sum(value.x)),
                           by=c("year_id", 'location_id','cause_id','measure_id', 'draw','lancet_label','rei_id')]
rf_detailed_calc$paf_pct<-rf_detailed_calc$value/rf_detailed_calc$value.daly
rf_detailed_calc$type<-'pafs'
rf_detailed_calc_pct<-rf_detailed_calc[,j=list(value=mean(paf_pct), lower=quantile(paf_pct, probs=0.025), upper=quantile(paf_pct,probs=0.975)),
                           by=c("year_id",'location_id','cause_id','measure_id','lancet_label','rei_id','type')]

#Percent change
paf_pct_change<-rf_detailed_calc[,j=list(value=((paf_pct[year_id==2021]-paf_pct[year_id==1990])/paf_pct[year_id==1990])), 
                                 by=c("cause_id",'draw','location_id','measure_id','rei_id','lancet_label')]
paf_pct_change$type<-'paf pct change'
paf_pct_change$value[paf_pct_change$value=='Inf']<-0
paf_pct_change$value[is.na(paf_pct_change$value)]<-0

pct_change<-paf_pct_change[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)), 
                  by=c("cause_id",'rei_id','type', 'location_id','measure_id','lancet_label')]
pct_change$year_id<-""
pct_change$type<-'pct_change 1990-2021'

#create files to save
burden_draw<-rbind(daly_rf_sum,daly_rf_pct,rf_detailed_calc,paf_pct_change, fill=T)
burden_summary<-rbind(daly_rf_sum1,daly_rf_pct1,rf_detailed_calc_pct,pct_change, fill=T)
write.csv(burden_draw, paste0(path,'draws/risk_factors/burden_',loc_launch,'.csv'))
write.csv(burden_summary, paste0(path,'summary/risk_factors/burden_',loc_launch,'.csv'))

# ---ANALYZE DATA : collapse draw files into mean & UI--------------------------

# Mean & UI for numbers and rates
mean<-df[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)), by=c("sex_id","year_id","metric_id",#'population',
                                        'location_id','age_group_id','cause_id','measure_id')]
mean$type<-'mean_ui'

#Proportion T2
df_15<-df[!df$age_group_id %in% c(2,3,4,5,6,7,34,238,388,389),]
ratio<-df_15[,j=list(value=(value[cause_id==976]/value[cause_id==587])), by=c("sex_id","year_id","metric_id",'population','draw',
                                                                      'location_id','age_group_id','measure_id')]
                                                                      
ratio_mean<-ratio[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)), 
                  by=c("sex_id","year_id","metric_id",#'population',
                                        'location_id','age_group_id','measure_id')]
ratio_mean$cause_id<-""

ratio_mean$type<-'prop_t2'

#Percent change
ratio<-df[,j=list(value=((value[year_id==2021]-value[year_id==1990])/value[year_id==1990])), by=c("sex_id","cause_id","metric_id",'draw',
                                                                           'location_id','age_group_id','measure_id')]

ratio$value[ratio$value=='Inf']<-0
ratio$value[is.na(ratio$value)]<-0
pct_change<-ratio[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)), 
                  by=c("sex_id","cause_id","metric_id",#'population',
                       'location_id','age_group_id','measure_id')]
pct_change$year_id<-""
pct_change$type<-'pct_change 1990-2021'

#Male/Female ratio
sex<-df[df$sex_id %in% c(1,2) & df$age_group_id==27,]
ratio<-sex[,j=list(value=value[sex_id==1]/value[sex_id==2]), by=c("cause_id","year_id","metric_id",'draw',
                                                                           'location_id','age_group_id','measure_id')]


ratio_sex<-ratio[,j=list(value=mean(value), lower=quantile(value, probs=0.025), upper=quantile(value,probs=0.975)), 
                  by=c("cause_id","year_id","metric_id",
                       'location_id','age_group_id','measure_id')]
ratio_sex$type<-'male_female ratio'
ratio_sex$sex_id<-""

summary<-rbind(mean, ratio_sex, ratio_mean,pct_change)

write.csv(summary, paste0(path,'summary/disease/burden_',loc_launch,'.csv'))


# ---Forecasts------------------------------------------------------------------

#pull in forecast
forecast<-filepath
#forecast rate
forecast_rate<-readRDS(paste0(forecast,'diabetes_rate_per_1.rds'))
forecast_rate<-forecast_rate[forecast_rate$location_id==loc_launch & forecast_rate$draw<100 & forecast_rate$age_group_id==27,]
forecast_rate$draw1<-paste0('draw_',forecast_rate$draw)

#forecast pop
forecast_count<-readRDS(paste0(forecast,'diabetes_count.rds'))
forecast_count<-forecast_count[forecast_count$location_id==loc_launch & forecast_count$draw<100,]
forecast_count$draw1<-paste0('draw_',forecast_count$draw)

#both forecast files
forecast_all<-rbind(forecast_rate, forecast_count)

forecast_all$cause_id[forecast_all$acause=='diabetes']<-587
forecast_all$cause_id[forecast_all$acause=='diabetes_typ1']<-975
forecast_all$cause_id[forecast_all$acause=='diabetes_typ2']<-976

#pull in gbd estimate
gbd<-df[df$measure_id==5  & df$sex_id==3 & ((df$age_group_id %in% c(27) & df$metric_id==3)|(df$age_group_id==22 & df$metric_id==1)),]
gbd$draw1<-gbd$draw

#intercept shift
int_forecast<-forecast_all[forecast_all$year_id==2021,]
int_gbd<-gbd[gbd$year_id==2021,]
int<-merge(int_gbd,int_forecast, by=c('age_group_id','draw1','location_id','sex_id','cause_id','year_id'),all=T)

#calculate difference between gbd 2021 estimate and forecast. gbd - forecast
int$diff<-(int$value.x-int$value.y)

#apply the difference to forecasted value
int<-int[,c('age_group_id','draw1','location_id','sex_id','cause_id','value.x','value.y','diff' )]

int_update<-merge(int,forecast_all,by=c('age_group_id','draw1','location_id','sex_id','cause_id'), all=T)
int_update$updated_forecast<-int_update$value+int_update$diff
int_update<-int_update[,c('age_group_id','draw1','location_id','sex_id','cause_id','year_id','updated_forecast','value')]
int_update<-int_update[int_update$year_id !=2021,]

#summary of forecasted estimate
gbd$updated_forecast<-gbd$value

forcast_plot<-rbind(gbd, int_update, fill=T)
forcast_plot$type<-'mean_ui'

#percent change between 1990-2021 and 2021-2050
forecast_pct<-forcast_plot[,j=list(value= (updated_forecast[year_id==2050] - updated_forecast[year_id==2021])/updated_forecast[year_id==2021]), 
                           by=c("sex_id","cause_id",'location_id','age_group_id','draw1')]
forecast_pct$type<-'percent change 2021-2050'
forecast_pct1<-forecast_pct[,j=list(value=mean( value), lower=quantile( value, probs=0.025), 
                                   upper=quantile( value,probs=0.975)), 
                           by=c("sex_id","cause_id",'location_id','age_group_id','type')]
#forecast mean & Ui
forecast_sum<-forcast_plot[,j=list(value=mean( updated_forecast), lower=quantile( updated_forecast, probs=0.025), 
                                   upper=quantile( updated_forecast,probs=0.975)), 
                                            by=c("sex_id","cause_id",'location_id','age_group_id','year_id','type')]

forcast_plot<-rbind(forcast_plot,forecast_pct,fill=T)
forecast_sum<-rbind(forecast_sum,forecast_pct1,fill=T)

write.csv(forcast_plot,paste0(path,'draws/forecast/burden_',loc_launch,'.csv'))
write.csv(forecast_sum,paste0(path,'summary/forecast/burden_',loc_launch,'.csv'))
