# Analyzing microdata with FPG/HbA1c for sensitivity test

library(ggplot2)

# read in microdata with FPG/HbA1c
df <- as.data.table(read.csv("FILEPATH/microdata.csv"))

# subsetting
df <- df[!is.na(fpg_mmol)&!is.na(a1c), ]

df <- df[fpg_mmol<2.5|fpg_mmol>25, is_outlier := 1]
df <- df[a1c<2|a1c>20, is_outlier := 1]
df <- df[is.na(is_outlier), is_outlier := 0]

df <- df[is_outlier == 0, ]

# linear regression
linreg <- lm(fpg_mmol~a1c, data = df)

summary(linreg)

ggplot(df, aes(x = a1c, y = fpg_mmol)) +
  geom_point() +
  stat_smooth(method = "lm")

# predict FPG equivalent of HbA1c 8%
newdata <- data.frame(a1c = 8)

predict(linreg, newdata)


# other scatter plots
ggplot(df, aes(x = a1c, y = fpg_mmol)) +
  geom_point() + facet_wrap(~location_name) +
  geom_abline(slope = 1, intercept = 0)
