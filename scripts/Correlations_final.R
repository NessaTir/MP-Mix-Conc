# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the statistical analysis of the relationship of all processed data and the concentrations. 
# Relationship will be tested for beeing
# 1) linear
#   - untransformed
#   - log-log transformed
#   - log-log transformed with low concentrations
#   - log-log transformed with high concentrations
# 2) logarithmic
# 3) exponential

# This script is closely connected to 'Plots - 4.5. Correlation Plot'

# ----- 2. Load in needed packages ---------------------------------------------
library(tidyverse)

# read in the data files
library(readxl)

# work with summarized data, e.g., means
library(rstatix)

# test for linear correlation
library(ggpubr)

# using AIC to find the best fitted test
library(AICcmodavg)


# ----- 3. Read in needed data files ------------------------------------------- 
## ---- 3.01. Surface ----------------------------------------------------------
surface <- read.csv2("processed/surface_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.02. Volume -----------------------------------------------------------
volume <- read.csv2("processed/volume_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.03. Calcification ----------------------------------------------------
calcification <- read.csv2("processed/calcification_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.04. Necrosis ---------------------------------------------------------
necrosis <- read.csv2("processed/necrosis_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.05. Polypactivity ----------------------------------------------------
polypactivity <- read.csv2("processed/polypactivity_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.06. YII --------------------------------------------------------------
YII <- read.csv2("processed/YII_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.07. FvFm -------------------------------------------------------------
FvFm <- read.csv2("processed/FvFm_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.08. rETRmax ----------------------------------------------------------
rETR <- read.csv2("processed/rETR_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.09. Ek ---------------------------------------------------------------
Ek <- read.csv2("processed/Ek_all.csv") %>%
  mutate(conc = as.numeric(conc))

## ---- 3.10. Alpha ------------------------------------------------------------
alpha <- read.csv2("processed/alpha_all.csv") %>%
  mutate(conc = as.numeric(conc))


# ----- 4. Looking for Relationships -------------------------------------------
## ---- 4.01. Surface ----------------------------------------------------------
### --- 4.01.1. Pve ------------------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

surface_Pve_wocon <- surface_Pve %>%
  subset(conc != "0")

# untransformed
Surf_Pve1 <- lm((value) ~ (conc), data = surface_Pve)
#view the output of the model
summary(Surf_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = surface_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -54.244 -14.692  -2.047  10.365  60.831 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 35.50787    2.57335    13.8   <2e-16 ***
#   conc        -0.18893    0.05725    -3.3   0.0014 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 21.22 on 88 degrees of freedom
# Multiple R-squared:  0.1101,	Adjusted R-squared:    0.1 
# F-statistic: 10.89 on 1 and 88 DF,  p-value: 0.001398 

# log-log transformed
Surf_Pve2 <- lm(value_log ~ conc_log, data = surface_Pve)
#view the output of the model
summary(Surf_Pve2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = surface_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.99039 -0.42016  0.07274  0.47442  1.17521 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.47227    0.09815  35.378  < 2e-16 ***
#   conc_log    -0.11769    0.04303  -2.735  0.00761 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.6864 on 84 degrees of freedom
# Multiple R-squared:  0.08176,	Adjusted R-squared:  0.07083 
# F-statistic: 7.479 on 1 and 84 DF,  p-value: 0.007612 

# log-log transformed - early linearity
surface_Pve_wocon100 <- surface_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
Surf_Pve3 <- lm(value_log ~ log(conc), data = surface_Pve_wocon100)
#view the output of the model
summary(Surf_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = surface_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.0273 -0.3533  0.1084  0.3741  1.1383 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 3.427629   0.087607  39.125   <2e-16 ***
#   log(conc)   0.006433   0.046820   0.137    0.891    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.6376 on 51 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.00037,	Adjusted R-squared:  -0.01923 
# F-statistic: 0.01888 on 1 and 51 DF,  p-value: 0.8913


# log-log transformed - late linearity
surface_Pve_wocon01 <- surface_Pve_wocon %>%
  subset(conc != "0.1")
Surf_Pve4 <- lm(value_log ~ log(conc), data = surface_Pve_wocon01)
#view the output of the model
summary(Surf_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ log(conc), data = surface_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.0450 -0.4402  0.0411  0.5478  1.1242 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.44529    0.15827  21.768   <2e-16 ***
#   log(conc)   -0.11023    0.05454  -2.021   0.0488 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.7314 on 49 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.07694,	Adjusted R-squared:  0.0581 
# F-statistic: 4.084 on 1 and 49 DF,  p-value: 0.04876

# Logarithmic
#fit the model
Surf_Pve5 <- lm(value ~ conc_log, data = surface_Pve)
#view the output of the model
summary(Surf_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = surface_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -55.775 -14.267  -2.437  13.344  61.656 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   37.038      3.039  12.189  < 2e-16 ***
#   conc_log      -3.671      1.295  -2.836  0.00567 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 21.53 on 88 degrees of freedom
# Multiple R-squared:  0.08372,	Adjusted R-squared:  0.07331 
# F-statistic: 8.041 on 1 and 88 DF,  p-value: 0.005674

# Exponential
#fit the model
Surf_Pve6 <- lm(log(value) ~ conc, data = surface_Pve)
#view the output of the model
summary(Surf_Pve6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = surface_Pve)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.02461 -0.44326  0.09223  0.44034  1.19492 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.431428   0.082492   41.60  < 2e-16 ***
#   conc        -0.006507   0.001902   -3.42 0.000968 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.6711 on 84 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.1222,	Adjusted R-squared:  0.1118 
# F-statistic:  11.7 on 1 and 84 DF,  p-value: 0.0009676


# AIC test
models_Surf_Pve <- list(Surf_Pve2, Surf_Pve3, Surf_Pve4)
model.names <- c('Surf_Pve2', 'Surf_Pve3', 'Surf_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_Surf_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Surf_Pve3 3 107.16       0.00   0.99   0.99 -50.33
# Surf_Pve4 3 117.30      10.15   0.01   1.00 -55.40
# Surf_Pve2 3 183.62      76.46   0.00   1.00 -88.66


### --- 4.01.2. Spi ------------------------------------------------------------
# Surface Spi
surface_Spi <- subset(surface, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

surface_Spi_wocon <- surface_Spi %>%
  subset(conc != "0")

# untransformed
Surf_Spi1 <- lm((value) ~ (conc), data = surface_Spi)
#view the output of the model
summary(Surf_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = surface_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -111.935   -9.132    5.044   19.470   57.807 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 34.79340    3.89116   8.942 5.41e-14 ***
#   conc        -0.22205    0.08657  -2.565    0.012 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 32.09 on 88 degrees of freedom
# Multiple R-squared:  0.06956,	Adjusted R-squared:  0.05899 
# F-statistic: 6.579 on 1 and 88 DF,  p-value: 0.01201 

# log-log transformed
Surf_Spi2 <- lm(value_log ~ conc_log, data = surface_Spi)
#view the output of the model
summary(Surf_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = surface_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.4782 -0.1617  0.2830  0.5601  1.1232 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.45418    0.14645   23.59   <2e-16 ***
#   conc_log    -0.07409    0.06389   -1.16     0.25    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.9904 on 80 degrees of freedom
# (8 observations deleted due to missingness)
# Multiple R-squared:  0.01653,	Adjusted R-squared:  0.004236 
# F-statistic: 1.345 on 1 and 80 DF,  p-value: 0.2497


# log-log transformed - early linearity
surface_Spi_wocon100 <- surface_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
Surf_Spi3 <- lm(value_log ~ log(conc), data = surface_Spi_wocon100)
#view the output of the model
summary(Surf_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ log(conc), data = surface_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.4841 -0.1057  0.1606  0.4838  1.1173 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.48799    0.12396  28.138   <2e-16 ***
#   conc_log    -0.08928    0.06528  -1.368    0.178    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.8755 on 48 degrees of freedom
# Multiple R-squared:  0.0375,	Adjusted R-squared:  0.01745 
# F-statistic:  1.87 on 1 and 48 DF,  p-value: 0.1778


# log-log transformed - late linearity
surface_Spi_wocon01 <- surface_Spi_wocon %>%
  subset(conc != "0.1")
Surf_Spi4 <- lm(value_log ~ log(conc), data = surface_Spi_wocon01)
#view the output of the model
summary(Surf_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = surface_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.5079 -0.2910  0.1891  0.5262  1.0935 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   3.5515     0.2045  17.366   <2e-16 ***
#   conc_log     -0.1065     0.0704  -1.513    0.137    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.9022 on 47 degrees of freedom
# Multiple R-squared:  0.04645,	Adjusted R-squared:  0.02616 
# F-statistic: 2.289 on 1 and 47 DF,  p-value: 0.137


# Logarithmic
#fit the model
Surf_Spi5 <- lm(value ~ conc_log, data = surface_Spi)
#view the output of the model
summary(Surf_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = volume_Pve_wocon)
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.044047 -0.017154 -0.002132  0.015144  0.067152 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.059694   0.004490  13.295  < 2e-16 ***
#   conc_log    -0.005106   0.001711  -2.984  0.00392 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.02549 on 70 degrees of freedom
# Multiple R-squared:  0.1128,	Adjusted R-squared:  0.1001 
# F-statistic: 8.902 on 1 and 70 DF,  p-value: 0.00392

# Exponential
#fit the model
Surf_Spi6 <- lm(log(value) ~ conc, data = surface_Spi)
#view the output of the model
summary(Surf_Spi6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = surface_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.5770 -0.1630  0.2878  0.5857  1.1210 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.407289   0.124862  27.288   <2e-16 ***
#   conc        -0.003190   0.002902  -1.099    0.275    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.9912 on 80 degrees of freedom
# (8 observations deleted due to missingness)
# Multiple R-squared:  0.01489,	Adjusted R-squared:  0.002571 
# F-statistic: 1.209 on 1 and 80 DF,  p-value: 0.2749


# AIC test
models_Surf_Spi <- list(Surf_Spi2, Surf_Spi3, Surf_Spi4)
model.names <- c('Surf_Spi2', 'Surf_Spi3', 'Surf_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_Surf_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Surf_Spi3 3 133.08       0.00   0.55   0.55  -63.28
# Surf_Spi4 3 133.46       0.39   0.45   1.00  -63.47
# Surf_Spi2 3 235.40     102.33   0.00   1.00 -114.55



## ---- 4.02. Volume -----------------------------------------------------------
### --- 4.02.1. Pve ------------------------------------------------------------
volume_Pve <- subset(volume, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

volume_Pve_wocon <- volume_Pve %>%
  filter(conc != "0")

# untransformed
Vol_Pve1 <- lm((value) ~ (conc), data = volume_Pve)
#view the output of the model
summary(Vol_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = volume_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.061735 -0.017387 -0.001495  0.017026  0.061432 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.564e-02  3.059e-03  18.192  < 2e-16 ***
#   conc        -2.472e-04  6.805e-05  -3.632 0.000471 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.02522 on 88 degrees of freedom
# Multiple R-squared:  0.1304,	Adjusted R-squared:  0.1205 
# F-statistic: 13.19 on 1 and 88 DF,  p-value: 0.0004715

# log-log transformed
Vol_Pve2 <- lm(value_log ~ conc_log, data = volume_Pve)
#view the output of the model
summary(Vol_Pve2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = volume_Pve_wocon)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.6331 -0.2702  0.1026  0.3899  1.0134 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.96149    0.07967 -37.170   <2e-16 ***
#   conc_log    -0.09099    0.03491  -2.606   0.0108 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5578 on 85 degrees of freedom
# Multiple R-squared:  0.07399,	Adjusted R-squared:  0.0631 
# F-statistic: 6.792 on 1 and 85 DF,  p-value: 0.01081

# log-log transformed - early linearity
volume_Pve_wocon100 <- volume_Pve_wocon %>%
  subset(conc != "100")
Vol_Pve3 <- lm(value_log ~ log(conc), data = volume_Pve_wocon100)
#view the output of the model
summary(Vol_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = volume_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.41334 -0.28985  0.05222  0.36336  0.80837 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.00060    0.07292 -41.147   <2e-16 ***
#   conc_log     0.01127    0.03879   0.291    0.773    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5359 on 52 degrees of freedom
# Multiple R-squared:  0.001621,	Adjusted R-squared:  -0.01758 
# F-statistic: 0.08441 on 1 and 52 DF,  p-value: 0.7726


# log-log transformed - late linearity
volume_Pve_wocon01 <- volume_Pve_wocon %>%
  subset(conc != "0.1")
Vol_Pve4 <- lm(value_log ~ log(conc), data = volume_Pve_wocon01)
#view the output of the model
summary(Vol_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = volume_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.4695 -0.2628  0.1177  0.3899  1.0001 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.94446    0.12882 -22.858   <2e-16 ***
#   conc_log    -0.09639    0.04455  -2.164   0.0353 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5975 on 50 degrees of freedom
# Multiple R-squared:  0.08561,	Adjusted R-squared:  0.06732 
# F-statistic: 4.681 on 1 and 50 DF,  p-value: 0.0353

# Logarithmic
#fit the model
Vol_Pve5 <- lm(value ~ conc_log, data = volume_Pve_wocon)
#view the output of the model
summary(Vol_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = volume_Pve_wocon)
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.044047 -0.017154 -0.002132  0.015144  0.067152 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.059694   0.004490  13.295  < 2e-16 ***
#   conc_log    -0.005106   0.001711  -2.984  0.00392 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.02549 on 70 degrees of freedom
# Multiple R-squared:  0.1128,	Adjusted R-squared:  0.1001 
# F-statistic: 8.902 on 1 and 70 DF,  p-value: 0.00392

# Exponential
#fit the model
Vol_Pve6 <- lm(log(value) ~ conc, data = volume_Pve)
#view the output of the model
summary(Vol_Pve6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = volume_Pve)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.60083 -0.29337  0.08486  0.37631  0.87822 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.993784   0.066732 -44.863  < 2e-16 ***
#   conc        -0.005071   0.001547  -3.277  0.00152 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5462 on 85 degrees of freedom
# (3 observations deleted due to missingness)
# Multiple R-squared:  0.1122,	Adjusted R-squared:  0.1017 
# F-statistic: 10.74 on 1 and 85 DF,  p-value: 0.001518


# AIC test
models_Vol_Pve <- list(Vol_Pve2, Vol_Pve3, Vol_Pve4)
model.names <- c('Vol_Pve2', 'Vol_Pve3', 'Vol_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_Vol_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Vol_Pve3 3  90.31       0.00   0.98   0.98 -41.92
# Vol_Pve4 3  98.46       8.15   0.02   1.00 -45.98
# Vol_Pve2 3 149.58      59.27   0.00   1.00 -71.65


### --- 4.02.2. Spi ------------------------------------------------------------
# volume Spi
volume_Spi <- subset(volume, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

volume_Spi_wocon <- volume_Spi %>%
  filter(conc != "0")

# untransformed
Vol_Spi1 <- lm((value) ~ (conc), data = volume_Spi)
#view the output of the model
summary(Vol_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = volume_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -0.054205 -0.017964 -0.000673  0.012412  0.094862 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.433e-02  3.867e-03  14.051   <2e-16 ***
#   conc        -2.051e-04  8.603e-05  -2.385   0.0192 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.03189 on 88 degrees of freedom
# Multiple R-squared:  0.06069,	Adjusted R-squared:  0.05002 
# F-statistic: 5.686 on 1 and 88 DF,  p-value: 0.01925

# log-log transformed
volume_Spi_wocon <- volume_Spi %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
Vol_Spi2 <- lm(value_log ~ conc_log, data = volume_Spi)
#view the output of the model
summary(Vol_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = volume_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1873 -0.1662  0.2771  0.5270  1.4402 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.20782    0.14476 -22.160   <2e-16 ***
#   conc_log    -0.06203    0.06238  -0.994    0.323    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.018 on 86 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.01137,	Adjusted R-squared:  -0.0001287 
# F-statistic: 0.9888 on 1 and 86 DF,  p-value: 0.3228

# log-log transformed - early linearity
volume_Spi_wocon100 <- volume_Spi_wocon %>%
  subset(conc != "100")
Vol_Spi3 <- lm(value_log ~ log(conc), data = volume_Spi_wocon100)
#view the output of the model
summary(Vol_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ log(conc), data = volume_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1756 -0.2460  0.1798  0.4536  1.1305 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.12067    0.12180 -25.621   <2e-16 ***
#   conc_log     0.04550    0.06418   0.709    0.482    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.8867 on 51 degrees of freedom
# Multiple R-squared:  0.009757,	Adjusted R-squared:  -0.00966 
# F-statistic: 0.5025 on 1 and 51 DF,  p-value: 0.4816


# log-log transformed - late linearity
volume_Spi_wocon01 <- volume_Spi_wocon %>%
  subset(conc != "0.1")
Vol_Spi4 <- lm(value_log ~ log(conc), data = volume_Spi_wocon01)
#view the output of the model
summary(Vol_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = volume_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.01816 -0.24918  0.07437  0.47604  1.35220 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -2.96661    0.19823 -14.965   <2e-16 ***
#   conc_log    -0.13114    0.06694  -1.959   0.0557 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.8988 on 50 degrees of freedom
# Multiple R-squared:  0.07128,	Adjusted R-squared:  0.0527 
# F-statistic: 3.837 on 1 and 50 DF,  p-value: 0.05571


# Logarithmic
#fit the model
Vol_Spi5 <- lm(value ~ conc_log, data = volume_Spi)
#view the output of the model
summary(Vol_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = volume_Spi)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.05289 -0.01741 -0.00217  0.01105  0.10010 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.054858   0.004571  12.002   <2e-16 ***
#   conc_log    -0.003260   0.001948  -1.674   0.0977 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.03239 on 88 degrees of freedom
# Multiple R-squared:  0.03086,	Adjusted R-squared:  0.01985 
# F-statistic: 2.802 on 1 and 88 DF,  p-value: 0.0977

# Exponential
#fit the model
Vol_Spi6 <- lm(log(value) ~ conc, data = volume_Spi)
#view the output of the model
summary(Vol_Spi6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = volume_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1993 -0.1359  0.2803  0.5460  1.3320 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.201219   0.123055 -26.015   <2e-16 ***
#   conc        -0.004720   0.002785  -1.695   0.0937 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.007 on 86 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.03232,	Adjusted R-squared:  0.02106 
# F-statistic: 2.872 on 1 and 86 DF,  p-value: 0.09375


# AIC test
models_Vol_Spi <- list(Vol_Spi2, Vol_Spi3, Vol_Spi4)
model.names <- c('Vol_Spi2', 'Vol_Spi3', 'Vol_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_Vol_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Vol_Spi3 3 133.08       0.00   0.55   0.55 -63.28
# Vol_Spi4 3 133.46       0.39   0.45   1.00 -63.47
# Vol_Spi2 3 164.82      31.74   0.00   1.00 -79.21



## ---- 4.03. Calcification ----------------------------------------------------
### --- 4.03.1. Pve ------------------------------------------------------------
calcification_Pve <- subset(calcification, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

calcification_Pve_wocon <- calcification_Pve %>%
  filter(conc != "0")

# untransformed
Calc_Pve1 <- lm((value) ~ (conc), data = calcification_Pve)
#view the output of the model
summary(Calc_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = calcification_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -75.608 -26.373  -0.341  23.189  79.265 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  87.5865     4.1396  21.158   <2e-16 ***
#   conc         -0.1886     0.0921  -2.048   0.0436 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 34.14 on 88 degrees of freedom
# Multiple R-squared:  0.04548,	Adjusted R-squared:  0.03464 
# F-statistic: 4.193 on 1 and 88 DF,  p-value: 0.04356

# log-log transformed
Calc_Pve2 <- lm(value_log ~ conc_log, data = calcification_Pve)
#view the output of the model
summary(Calc_Pve2)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ conc_log, data = calcification_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.9078 -0.2754  0.1147  0.3467  0.8417 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.39089    0.07503  58.518   <2e-16 ***
#   conc_log    -0.05291    0.03197  -1.655    0.102    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5317 on 88 degrees of freedom
# Multiple R-squared:  0.03018,	Adjusted R-squared:  0.01916 
# F-statistic: 2.739 on 1 and 88 DF,  p-value: 0.1015

# log-log transformed - early linearity
calcification_Pve_wocon100 <- calcification_Pve_wocon %>%
  subset(conc != "100")
Calc_Pve3 <- lm(value_log ~ log(conc), data = calcification_Pve_wocon100)
#view the output of the model
summary(Calc_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = calcification_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.65768 -0.24633  0.04884  0.31059  0.68515 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.401900   0.061093   72.05   <2e-16 ***
#   conc_log    0.008116   0.032496    0.25    0.804    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4489 on 52 degrees of freedom
# Multiple R-squared:  0.001198,	Adjusted R-squared:  -0.01801 
# F-statistic: 0.06237 on 1 and 52 DF,  p-value: 0.8038


# log-log transformed - late linearity
calcification_Pve_wocon01 <- calcification_Pve_wocon %>%
  subset(conc != "0.1")
Calc_Pve4 <- lm(value_log ~ log(conc), data = calcification_Pve_wocon01)
#view the output of the model
summary(Calc_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = calcification_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.5321 -0.3905  0.1490  0.4168  0.8108 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.46974    0.11887  37.602   <2e-16 ***
#   conc_log    -0.07590    0.03999  -1.898   0.0632 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5525 on 52 degrees of freedom
# Multiple R-squared:  0.0648,	Adjusted R-squared:  0.04681 
# F-statistic: 3.603 on 1 and 52 DF,  p-value: 0.06324

# Logarithmic
#fit the model
Calc_Pve5 <- lm(value ~ conc_log, data = calcification_Pve)
#view the output of the model
summary(Calc_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = calcification_Pve)
# Residuals:
#  Min      1Q  Median      3Q     Max 
#-75.947 -26.648   0.564  22.015  84.002 
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   87.925      4.877  18.029   <2e-16 ***
#  conc_log      -2.903      2.078  -1.397    0.166    
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 34.56 on 88 degrees of freedom
#Multiple R-squared:  0.02169,	Adjusted R-squared:  0.01058 
#F-statistic: 1.951 on 1 and 88 DF,  p-value: 0.1659

# Exponential
#fit the model
Calc_Pve6 <- lm(log(value) ~ conc, data = calcification_Pve)
#view the output of the model
summary(Calc_Pve6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = calcification_Pve)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.89431 -0.26309  0.09206  0.38510  0.77965 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.377428   0.063758  68.657   <2e-16 ***
#   conc        -0.003109   0.001419  -2.192    0.031 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.5258 on 88 degrees of freedom
# Multiple R-squared:  0.05177,	Adjusted R-squared:  0.041 
# F-statistic: 4.805 on 1 and 88 DF,  p-value: 0.03102


# AIC test
models_Calc_Pve <- list(Calc_Pve2, Calc_Pve3, Calc_Pve4)
model.names <- c('Calc_Pve2', 'Calc_Pve3', 'Calc_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_Calc_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Calc_Pve3 3  71.19       0.00      1      1 -32.36
# Calc_Pve4 3  93.60      22.41      0      1 -43.56
# Calc_Pve2 3 145.97      74.78      0      1 -69.85


### --- 4.03.2. Spi ------------------------------------------------------------
# calcification Spi
calcification_Spi <- subset(calcification, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

calcification_Spi_wocon <- calcification_Spi %>%
  filter(conc != "0")

# untransformed
Calc_Spi1 <- lm((value) ~ (conc), data = calcification_Spi)
#view the output of the model
summary(Calc_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = calcification_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -96.485 -21.093  -1.005  20.179 111.680 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  91.5509     5.5799  16.407   <2e-16 ***
#   conc         -0.2961     0.1241  -2.385   0.0192 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 46.01 on 88 degrees of freedom
# Multiple R-squared:  0.06071,	Adjusted R-squared:  0.05003 
# F-statistic: 5.688 on 1 and 88 DF,  p-value: 0.01923

# log-log transformed
Calc_Spi2 <- lm(value_log ~ conc_log, data = calcification_Spi)
#view the output of the model
summary(Calc_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = calcification_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.1082 -0.0506  0.1518  0.4678  1.1528 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.41402    0.12103  36.472   <2e-16 ***
#   conc_log    -0.11143    0.05099  -2.185   0.0316 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.8413 on 86 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.0526,	Adjusted R-squared:  0.04158 
# F-statistic: 4.775 on 1 and 86 DF,  p-value: 0.0316

# log-log transformed - early linearity
calcification_Spi_wocon100 <- calcification_Spi_wocon %>%
  subset(conc != "100")
Calc_Spi3 <- lm(value_log ~ log(conc), data = calcification_Spi_wocon100)
#view the output of the model
summary(Calc_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = calcification_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.1379 -0.0760  0.0860  0.3626  1.0031 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.36649    0.11326   38.55   <2e-16 ***
#   conc_log    -0.03039    0.06083   -0.50     0.62    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.8158 on 50 degrees of freedom
# Multiple R-squared:  0.004966,	Adjusted R-squared:  -0.01493 
# F-statistic: 0.2495 on 1 and 50 DF,  p-value: 0.6196


# log-log transformed - late linearity
calcification_Spi_wocon01 <- calcification_Spi_wocon %>%
  subset(conc != "0.1")
Calc_Spi4 <- lm(value_log ~ log(conc), data = calcification_Spi_wocon01)
#view the output of the model
summary(Calc_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = calcification_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.1010 -0.0676  0.1825  0.5080  1.1733 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.32966    0.20722  20.894   <2e-16 ***
#   conc_log    -0.08828    0.06971  -1.266    0.211    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.9631 on 52 degrees of freedom
# Multiple R-squared:  0.02992,	Adjusted R-squared:  0.01126 
# F-statistic: 1.604 on 1 and 52 DF,  p-value: 0.211


# Logarithmic
#fit the model
Calc_Spi5 <- lm(value ~ conc_log, data = calcification_Spi)
#view the output of the model
summary(Calc_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = calcification_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -97.558 -21.596  -1.921  19.973 119.656 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   93.089      6.572  14.164   <2e-16 ***
#   conc_log      -5.202      2.800  -1.858   0.0666 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 46.57 on 88 degrees of freedom
# Multiple R-squared:  0.03774,	Adjusted R-squared:  0.0268 
# F-statistic: 3.451 on 1 and 88 DF,  p-value: 0.06656

# Exponential
#fit the model
Calc_Spi6 <- lm(log(value) ~ conc, data = calcification_Spi)
#view the output of the model
summary(Calc_Spi6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = calcification_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.1297 -0.0774  0.1671  0.4427  1.2442 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.363893   0.102775   42.46   <2e-16 ***
#   conc        -0.005609   0.002261   -2.48   0.0151 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.835 on 86 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.06677,	Adjusted R-squared:  0.05592 
# F-statistic: 6.153 on 1 and 86 DF,  p-value: 0.01507


# AIC test
models_Calc_Spi <- list(Calc_Spi2, Calc_Spi3, Calc_Spi4)
model.names <- c('Calc_Spi2', 'Calc_Spi3', 'Calc_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_Calc_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Calc_Spi3 3 130.85       0.00      1      1 -62.18
# Calc_Spi4 3 153.63      22.77      0      1 -73.57
# Calc_Spi2 3 223.59      92.73      0      1 -108.65


## ---- 4.04. Necrosis ---------------------------------------------------------
### --- 4.04.1. Pve ------------------------------------------------------------
necrosis_Pve <- subset(necrosis, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

necrosis_Pve_wocon <- necrosis_Pve %>%
  filter(conc != "0")

necrosis_Pve_Inf <- necrosis_Pve %>%
  filter(value_log != "-Inf")

# untransformed
Necro_Pve1 <- lm((value) ~ (conc), data = necrosis_Pve)
#view the output of the model
summary(Necro_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = necrosis_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.698 -1.463 -1.339 -1.326 39.334 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  1.32557    0.77966   1.700   0.0926 .
# conc         0.01372    0.01735   0.791   0.4310  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 6.429 on 88 degrees of freedom
# Multiple R-squared:  0.007062,	Adjusted R-squared:  -0.004221 
# F-statistic: 0.6259 on 1 and 88 DF,  p-value: 0.431


# log-log transformed
Necro_Pve2 <- lm(value_log ~ conc_log, data = necrosis_Pve_Inf)
#view the output of the model
summary(Necro_Pve2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = necrosis_Pve_Inf)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8448 -0.5163 -0.0236  0.6951  1.7973 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  2.28671    0.56070   4.078  0.00354 **
#   conc_log    -0.07487    0.21028  -0.356  0.73100   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.282 on 8 degrees of freedom
# Multiple R-squared:  0.0156,	Adjusted R-squared:  -0.1074 
# F-statistic: 0.1268 on 1 and 8 DF,  p-value: 0.731

# log-log transformed - early linearity
necrosis_Pve_wocon100 <- necrosis_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "-Inf")
Necro_Pve3 <- lm(value_log ~ log(conc), data = necrosis_Pve_wocon100)
#view the output of the model
summary(Necro_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = necrosis_Pve_wocon100)
# Residuals:
#      1       2       3       4       5 
# -0.0461  0.9978  0.9978 -0.3316 -1.6179 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   2.0079     0.5610   3.579   0.0373 *
#   conc_log      0.2173     0.3852   0.564   0.6122  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.254 on 3 degrees of freedom
# Multiple R-squared:  0.09587,	Adjusted R-squared:  -0.2055 
# F-statistic: 0.3181 on 1 and 3 DF,  p-value: 0.6122


# log-log transformed - late linearity
necrosis_Pve_wocon01 <- necrosis_Pve_wocon %>%
  subset(conc != "0.1") %>%
  filter(value_log != "-Inf")
Necro_Pve4 <- lm(value_log ~ log(conc), data = necrosis_Pve_wocon01)
#view the output of the model
summary(Necro_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = necrosis_Pve_wocon01)
# Residuals:
#         1        2        3        4        5        7        8 
# -1.90147 -0.55391  0.34495  1.64482  1.63297  0.05946 -1.22682 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   1.6168     0.8223   1.966    0.106
# conc_log      0.1061     0.2620   0.405    0.702
# Residual standard error: 1.478 on 5 degrees of freedom
# Multiple R-squared:  0.03175,	Adjusted R-squared:  -0.1619 
# F-statistic: 0.1639 on 1 and 5 DF,  p-value: 0.7023

# Logarithmic
#fit the model
Necro_Pve5 <- lm(value ~ conc_log, data = necrosis_Pve)
#view the output of the model
summary(Necro_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = necrosis_Pve)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -2.549 -1.882 -1.370 -1.161 39.483 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   1.1614     0.9074   1.280    0.204
# conc_log      0.3007     0.3866   0.778    0.439
# Residual standard error: 6.43 on 88 degrees of freedom
# Multiple R-squared:  0.006825,	Adjusted R-squared:  -0.004461 
# F-statistic: 0.6048 on 1 and 88 DF,  p-value: 0.4389

# Exponential
#fit the model
Necro_Pve6 <- lm(value_log ~ conc, data = necrosis_Pve_Inf)
#view the output of the model
summary(Necro_Pve6)
# OUTPUT: Call:
# lm(formula = value_log ~ conc, data = necrosis_Pve_Inf)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.88429 -0.52672 -0.04259  0.70333  1.87405 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  2.278451   0.491108   4.639  0.00167 **
#   conc        -0.004141   0.008951  -0.463  0.65598   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.275 on 8 degrees of freedom
# Multiple R-squared:  0.02605,	Adjusted R-squared:  -0.09569 
# F-statistic: 0.214 on 1 and 8 DF,  p-value: 0.656


# AIC test
models_Necro_Pve <- list(Necro_Pve2, Necro_Pve3, Necro_Pve4)
model.names <- c('Necro_Pve2', 'Necro_Pve3', 'Necro_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_Necro_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Necro_Pve4 3 36.98       0.00   0.58   0.58 -11.49
# Necro_Pve2 3 41.12       4.14   0.11   0.97 -15.56
# Necro_Pve3 3 43.90       6.92   0.02   1.00  -6.95


### --- 4.04.2. Spi ------------------------------------------------------------
# necrosis Spi
necrosis_Spi <- subset(necrosis, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

necrosis_Spi_wocon <- necrosis_Spi %>%
  filter(conc != "0")

necrosis_Spi_Inf <- necrosis_Spi %>%
  filter(value_log != "-Inf")

# untransformed
Necro_Spi1 <- lm((value) ~ (conc), data = necrosis_Spi)
#view the output of the model
summary(Necro_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = necrosis_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -8.963 -3.677 -3.096 -3.090 64.589 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  3.09011    1.69815   1.820   0.0722 .
# conc         0.05873    0.03778   1.555   0.1237  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 14 on 88 degrees of freedom
# Multiple R-squared:  0.02673,	Adjusted R-squared:  0.01567 
# F-statistic: 2.416 on 1 and 88 DF,  p-value: 0.1237

# log-log transformed
Necro_Spi2 <- lm(value_log ~ conc_log, data = necrosis_Spi_Inf)
#view the output of the model
summary(Necro_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = necrosis_Spi_Inf)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.5412 -1.0392 -0.1687  1.5558  2.1511 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  2.42257    0.61841   3.917  0.00137 **
#   conc_log    -0.05972    0.20663  -0.289  0.77652   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.598 on 15 degrees of freedom
# Multiple R-squared:  0.005538,	Adjusted R-squared:  -0.06076 
# F-statistic: 0.08354 on 1 and 15 DF,  p-value: 0.7765

# log-log transformed - early linearity
necrosis_Spi_wocon100 <- necrosis_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "-Inf")
Necro_Spi3 <- lm(value_log ~ log(conc), data = necrosis_Spi_wocon100)
#view the output of the model
summary(Necro_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = necrosis_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.8921 -0.9047  0.1045  0.9349  2.0776 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   2.0314     0.4981   4.078  0.00354 **
#   conc_log     -0.3198     0.2586  -1.237  0.25122   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.564 on 8 degrees of freedom
# Multiple R-squared:  0.1605,	Adjusted R-squared:  0.0556 
# F-statistic:  1.53 on 1 and 8 DF,  p-value: 0.2512


# log-log transformed - late linearity
necrosis_Spi_wocon01 <- necrosis_Spi_wocon %>%
  subset(conc != "0.1") %>%
  filter(value_log != "-Inf")
Necro_Spi4 <- lm(value_log ~ log(conc), data = necrosis_Spi_wocon01)
#view the output of the model
summary(Necro_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ log(conc), data = necrosis_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.2388 -1.0424 -0.1753  1.5526  2.1479 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  2.20564    0.78979   2.793   0.0175 *
#   conc_log    -0.01206    0.23372  -0.052   0.9598  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.551 on 11 degrees of freedom
# Multiple R-squared:  0.0002419,	Adjusted R-squared:  -0.09065 
# F-statistic: 0.002662 on 1 and 11 DF,  p-value: 0.9598


# Logarithmic
#fit the model
Necro_Spi5 <- lm(value ~ conc_log, data = necrosis_Spi)
#view the output of the model
summary(Necro_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = necrosis_Spi)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -7.203 -5.165 -3.323 -2.961 66.349 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)   2.9610     1.9898   1.488    0.140
# conc_log      0.9191     0.8479   1.084    0.281
# Residual standard error: 14.1 on 88 degrees of freedom
# Multiple R-squared:  0.01318,	Adjusted R-squared:  0.001964 
# F-statistic: 1.175 on 1 and 88 DF,  p-value: 0.2813

# Exponential
#fit the model
Necro_Spi6 <- lm(value_log ~ conc, data = necrosis_Spi_Inf)
#view the output of the model
summary(Necro_Spi6)
# OUTPUT: Call:
# lm(formula = value_log ~ conc, data = necrosis_Spi_Inf)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.28154 -1.38315 -0.07964  1.21871  2.00965 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 2.156881   0.500336   4.311 0.000618 ***
#   conc        0.003340   0.008394   0.398 0.696273    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.594 on 15 degrees of freedom
# Multiple R-squared:  0.01045,	Adjusted R-squared:  -0.05552 
# F-statistic: 0.1584 on 1 and 15 DF,  p-value: 0.6963


# AIC test
models_Necro_Spi <- list(Necro_Spi2, Necro_Spi3, Necro_Spi4)
model.names <- c('Necro_Spi2', 'Necro_Spi3', 'Necro_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_Necro_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Necro_Spi3 3 45.09       0.00   0.99   0.99 -17.55
# Necro_Spi4 3 54.80       9.71   0.01   1.00 -23.07
# Necro_Spi2 3 69.90      24.81   0.00   1.00 -31.03



## ---- 4.05. Polyp activity ---------------------------------------------------
### --- 4.05.1. Pve ------------------------------------------------------------
polypactivity_Pve <- subset(polypactivity, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

polypactivity_Pve_wocon <- polypactivity_Pve %>%
  filter(conc != "0")

# untransformed
Polyp_Pve1 <- lm((value) ~ (conc), data = polypactivity_Pve)
#view the output of the model
summary(Polyp_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = polypactivity_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.43369 -0.09729  0.03164  0.11600  0.35637 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.8486369  0.0179668  47.233  < 2e-16 ***
#   conc        -0.0026056  0.0003997  -6.518 4.31e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1482 on 88 degrees of freedom
# Multiple R-squared:  0.3256,	Adjusted R-squared:  0.3179 
# F-statistic: 42.49 on 1 and 88 DF,  p-value: 4.313e-09

# log-log transformed
Polyp_Pve2 <- lm(value_log ~ conc_log, data = polypactivity_Pve)
#view the output of the model
summary(Polyp_Pve2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = polypactivity_Pve)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.24188 -0.06172  0.02039  0.13389  0.49272 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.12512    0.03365  -3.718 0.000353 ***
#   conc_log    -0.09204    0.01434  -6.418 6.76e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2385 on 88 degrees of freedom
# Multiple R-squared:  0.3188,	Adjusted R-squared:  0.3111 
# F-statistic: 41.19 on 1 and 88 DF,  p-value: 6.757e-09

# log-log transformed - early linearity
polypactivity_Pve_wocon100 <- polypactivity_Pve_wocon %>%
  subset(conc != "100")
Polyp_Pve3 <- lm(value_log ~ log(conc), data = polypactivity_Pve_wocon100)
#view the output of the model
summary(Polyp_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = polypactivity_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.65878 -0.06273  0.01831  0.11558  0.28569 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.20063    0.02506  -8.006 1.24e-10 ***
#   conc_log    -0.03694    0.01333  -2.771  0.00773 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1841 on 52 degrees of freedom
# Multiple R-squared:  0.1287,	Adjusted R-squared:  0.1119 
# F-statistic:  7.68 on 1 and 52 DF,  p-value: 0.007728


# log-log transformed - late linearity
polypactivity_Pve_wocon01 <- polypactivity_Pve_wocon %>%
  subset(conc != "0.1")
Polyp_Pve4 <- lm(value_log ~ log(conc), data = polypactivity_Pve_wocon01)
#view the output of the model
summary(Polyp_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = polypactivity_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.25297 -0.06857  0.00042  0.20571  0.48163 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.18274    0.06122  -2.985  0.00432 ** 
#   conc_log    -0.07731    0.02060  -3.754  0.00044 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2845 on 52 degrees of freedom
# Multiple R-squared:  0.2132,	Adjusted R-squared:  0.1981 
# F-statistic: 14.09 on 1 and 52 DF,  p-value: 0.0004404

# Logarithmic
#fit the model
Polyp_Pve5 <- lm(value ~ conc_log, data = polypactivity_Pve)
#view the output of the model
summary(Polyp_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = polypactivity_Pve)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.44188 -0.07412  0.00793  0.11620  0.33590 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.88380    0.02063  42.843  < 2e-16 ***
#   conc_log    -0.05964    0.00879  -6.785 1.29e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1462 on 88 degrees of freedom
# Multiple R-squared:  0.3435,	Adjusted R-squared:  0.336 
# F-statistic: 46.04 on 1 and 88 DF,  p-value: 1.29e-09

# Exponential
#fit the model
Polyp_Pve6 <- lm(log(value) ~ conc, data = polypactivity_Pve)
#view the output of the model
summary(Polyp_Pve6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = polypactivity_Pve)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.20309 -0.10635  0.05954  0.16130  0.53151 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.1773198  0.0289732  -6.120 2.53e-08 ***
#   conc        -0.0041135  0.0006446  -6.381 7.96e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2389 on 88 degrees of freedom
# Multiple R-squared:  0.3164,	Adjusted R-squared:  0.3086 
# F-statistic: 40.72 on 1 and 88 DF,  p-value: 7.963e-09


# AIC test
models_Polyp_Pve <- list(Polyp_Pve2, Polyp_Pve3, Polyp_Pve4)
model.names <- c('Polyp_Pve2', 'Polyp_Pve3', 'Polyp_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_Polyp_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Polyp_Pve3 3 -25.05       0.00      1      1 15.77
# Polyp_Pve2 3  12.65      37.70      0      1 -3.15
# Polyp_Pve4 3  21.94      46.99      0      1 -7.73


### --- 4.05.2. Spi ------------------------------------------------------------
# polypactivity Spi
polypactivity_Spi <- subset(polypactivity, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

polypactivity_Spi_wocon <- polypactivity_Spi %>%
  filter(conc != "0")

polypactivity_Spi_Inf <- polypactivity_Spi %>%
  filter(value_log != "-Inf")

# untransformed
Polyp_Spi1 <- lm((value) ~ (conc), data = polypactivity_Spi)
#view the output of the model
summary(Polyp_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = polypactivity_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -0.41886 -0.19799 -0.00455  0.17405  0.50738 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.6425804  0.0275045  23.363   <2e-16 ***
#   conc        -0.0014996  0.0006119  -2.451   0.0162 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2268 on 88 degrees of freedom
# Multiple R-squared:  0.06388,	Adjusted R-squared:  0.05325 
# F-statistic: 6.005 on 1 and 88 DF,  p-value: 0.01624

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
Polyp_Spi2 <- lm(value_log ~ conc_log, data = polypactivity_Spi)
#view the output of the model
summary(Polyp_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = polypactivity_Spi)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.35187 -0.33291  0.06984  0.38083  0.84536 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.45459    0.06228  -7.299 1.22e-10 ***
#   conc_log    -0.08467    0.02654  -3.191  0.00197 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4413 on 88 degrees of freedom
# Multiple R-squared:  0.1037,	Adjusted R-squared:  0.0935 
# F-statistic: 10.18 on 1 and 88 DF,  p-value: 0.001968

# log-log transformed - early linearity
polypactivity_Spi_wocon100 <- polypactivity_Spi_wocon %>%
  subset(conc != "100")
Polyp_Spi3 <- lm(value_log ~ log(conc), data = polypactivity_Spi_wocon100)
#view the output of the model
summary(Polyp_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = polypactivity_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.96980 -0.22303  0.02622  0.28296  0.59843 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.53428    0.05227 -10.222 4.85e-14 ***
#   conc_log    -0.02786    0.02780  -1.002    0.321    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3841 on 52 degrees of freedom
# Multiple R-squared:  0.01895,	Adjusted R-squared:  8.268e-05 
# F-statistic: 1.004 on 1 and 52 DF,  p-value: 0.3209


# log-log transformed - late linearity
polypactivity_Spi_wocon01 <- polypactivity_Spi_wocon %>%
  subset(conc != "0.1")
Polyp_Spi4 <- lm(value_log ~ log(conc), data = polypactivity_Spi_wocon01)
#view the output of the model
summary(Polyp_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = polypactivity_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.3496 -0.3962  0.1465  0.3949  0.8476 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.47190    0.10733  -4.397 5.45e-05 ***
#   conc_log    -0.08159    0.03611  -2.260   0.0281 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4988 on 52 degrees of freedom
# Multiple R-squared:  0.08942,	Adjusted R-squared:  0.0719 
# F-statistic: 5.106 on 1 and 52 DF,  p-value: 0.02806


# Logarithmic
#fit the model
Polyp_Spi5 <- lm(value ~ conc_log, data = polypactivity_Spi)
#view the output of the model
summary(Polyp_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ conc_log, data = polypactivity_Spi)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.41787 -0.21032 -0.00065  0.16601  0.49935 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.66473    0.03186  20.862   <2e-16 ***
#   conc_log    -0.03555    0.01358  -2.619   0.0104 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2258 on 88 degrees of freedom
# Multiple R-squared:  0.07229,	Adjusted R-squared:  0.06175 
# F-statistic: 6.857 on 1 and 88 DF,  p-value: 0.01039

# Exponential
#fit the model
Polyp_Spi6 <- lm(log(value) ~ conc, data = polypactivity_Spi)
#view the output of the model
summary(Vol_Spi6)
# OUTPUT: Call:
# lm(formula = log(value) ~ conc, data = volume_Spi)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1993 -0.1359  0.2803  0.5460  1.3320 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.201219   0.123055 -26.015   <2e-16 ***
#   conc        -0.004720   0.002785  -1.695   0.0937 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 1.007 on 86 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.03232,	Adjusted R-squared:  0.02106 
# F-statistic: 2.872 on 1 and 86 DF,  p-value: 0.09375

# Exponential
ks.test(polypactivity_Spi$value, polypactivity_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  polypactivity_Spi$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_Polyp_Spi <- list(Polyp_Spi2, Polyp_Spi3, Polyp_Spi4)
model.names <- c('Polyp_Spi2', 'Polyp_Spi3', 'Polyp_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_Polyp_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Polyp_Spi3 3 54.35       0.00      1      1 -23.93
# Polyp_Spi4 3 82.58      28.23      0      1 -38.05
# Polyp_Spi2 3 97.35      43.00      0      1 -45.50




## ---- 4.06. YII --------------------------------------------------------------
### --- 4.06.1. Pve ------------------------------------------------------------
YII_Pve <- subset(YII, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
YII_Pve1 <- lm((value) ~ (conc), data = YII_Pve)
#view the output of the model
summary(YII_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = YII_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -12.4270  -3.5702   0.1322   2.9918  21.5346 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 102.14475    0.68912 148.225   <2e-16 ***
#   conc          0.03206    0.01533   2.091   0.0394 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 5.683 on 88 degrees of freedom
# Multiple R-squared:  0.04735,	Adjusted R-squared:  0.03652 
# F-statistic: 4.374 on 1 and 88 DF,  p-value: 0.03938

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
YII_Pve_wocon <- YII_Pve %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
YII_Pve2 <- lm(value_log ~ conc_log, data = YII_Pve_wocon)
#view the output of the model
summary(YII_Pve2)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = YII_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.142236 -0.026691  0.002753  0.022866  0.175445 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.624755   0.007039 657.033  < 2e-16 ***
#   conc_log    0.008721   0.002496   3.494 0.000829 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.05452 on 70 degrees of freedom
# Multiple R-squared:  0.1485,	Adjusted R-squared:  0.1364 
# F-statistic: 12.21 on 1 and 70 DF,  p-value: 0.0008289

# log-log transformed - early linearity
YII_Pve_wocon100 <- YII_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
YII_Pve3 <- lm(value_log ~ conc_log, data = YII_Pve_wocon100)
#view the output of the model
summary(YII_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = YII_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.157088 -0.028741  0.005773  0.030356  0.160593 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.628468   0.007809 592.678  < 2e-16 ***
#   conc_log    0.013559   0.004154   3.264  0.00194 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.05739 on 52 degrees of freedom
# Multiple R-squared:  0.1701,	Adjusted R-squared:  0.1541 
# F-statistic: 10.66 on 1 and 52 DF,  p-value: 0.001944


# log-log transformed - late linearity
YII_Pve_wocon01 <- YII_Pve_wocon %>%
  subset(conc != "0.1")
YII_Pve4 <- lm(value_log ~ conc_log, data = YII_Pve_wocon01)
#view the output of the model
summary(YII_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = YII_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.14684 -0.02434 -0.00117  0.01970  0.17084 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.643170   0.011158 416.116   <2e-16 ***
#   conc_log    0.002723   0.003754   0.726    0.471    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.05186 on 52 degrees of freedom
# Multiple R-squared:  0.01002,	Adjusted R-squared:  -0.009017 
# F-statistic: 0.5264 on 1 and 52 DF,  p-value: 0.4714

# Logarithmic
plot(YII_Pve$conc, YII_Pve$value)
#fit the model
YII_Pve5 <- lm(value ~ conc_log, data = YII_Pve_wocon)
#view the output of the model
summary(YII_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = YII_Pve_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -13.9522  -2.8509   0.1149   2.1959  19.7962 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 102.1607     0.7322 139.520  < 2e-16 ***
#   conc_log      0.8873     0.2596   3.417  0.00106 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 5.672 on 70 degrees of freedom
# Multiple R-squared:  0.143,	Adjusted R-squared:  0.1307 
# F-statistic: 11.68 on 1 and 70 DF,  p-value: 0.001057

# Exponential
ks.test(YII_Pve$value, YII_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  YII_Pve$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_YII_Pve <- list(YII_Pve2, YII_Pve3, YII_Pve4)
model.names <- c('YII_Pve2', 'YII_Pve3', 'YII_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_YII_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# YII_Pve2 3 -210.26       0.00      1      1 108.31
# YII_Pve4 3 -161.91      48.36      0      1  84.19
# YII_Pve3 3 -150.97      59.29      0      1  78.72


### --- 4.06.2. Spi ------------------------------------------------------------
# YII Spi
YII_Spi <- subset(YII, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
YII_Spi1 <- lm((value) ~ (conc), data = YII_Spi)
#view the output of the model
summary(YII_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = YII_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -16.5670  -4.8838  -0.2954   4.2738  23.8142 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 103.93158    0.89176 116.546   <2e-16 ***
#   conc          0.01954    0.01984   0.985    0.327    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 7.354 on 88 degrees of freedom
# Multiple R-squared:  0.0109,	Adjusted R-squared:  -0.0003393 
# F-statistic: 0.9698 on 1 and 88 DF,  p-value: 0.3274

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
YII_Spi_wocon <- YII_Spi %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
YII_Spi2 <- lm(value_log ~ conc_log, data = YII_Spi_wocon)
#view the output of the model
summary(YII_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = YII_Spi_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.156151 -0.039911 -0.001296  0.036971  0.195982 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.640694   0.009153 507.025   <2e-16 ***
#   conc_log    0.006469   0.003246   1.993   0.0501 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.0709 on 70 degrees of freedom
# Multiple R-squared:  0.0537,	Adjusted R-squared:  0.04018 
# F-statistic: 3.973 on 1 and 70 DF,  p-value: 0.05015

# log-log transformed - early linearity
YII_Spi_wocon100 <- YII_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
YII_Spi3 <- lm(value_log ~ conc_log, data = YII_Spi_wocon100)
#view the output of the model
summary(YII_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = YII_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.174222 -0.036416 -0.000538  0.033110  0.177911 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.645212   0.010230 454.068   <2e-16 ***
#   conc_log    0.012355   0.005441   2.271   0.0273 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.07518 on 52 degrees of freedom
# Multiple R-squared:  0.0902,	Adjusted R-squared:  0.0727 
# F-statistic: 5.155 on 1 and 52 DF,  p-value: 0.02735


# log-log transformed - late linearity
YII_Spi_wocon01 <- YII_Spi_wocon %>%
  subset(conc != "0.1")
YII_Spi4 <- lm(value_log ~ conc_log, data = YII_Spi_wocon01)
#view the output of the model
summary(YII_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = YII_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.155641 -0.049500  0.000609  0.032584  0.196491 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.638655   0.015689 295.659   <2e-16 ***
#   conc_log    0.007133   0.005278   1.352    0.182    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.07292 on 52 degrees of freedom
# Multiple R-squared:  0.03393,	Adjusted R-squared:  0.01536 
# F-statistic: 1.827 on 1 and 52 DF,  p-value: 0.1824


# Logarithmic
plot(YII_Spi$conc, YII_Spi$value)
#fit the model
YII_Spi5 <- lm(value ~ log(conc), data = YII_Spi_wocon)
#view the output of the model
summary(YII_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = YII_Spi_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -15.4772  -4.2948  -0.3825   3.6747  22.4974 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 103.8810     0.9690 107.199   <2e-16 ***
#   log(conc)     0.6787     0.3436   1.975   0.0522 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 7.506 on 70 degrees of freedom
# Multiple R-squared:  0.05279,	Adjusted R-squared:  0.03926 
# F-statistic: 3.901 on 1 and 70 DF,  p-value: 0.0522

# Exponential
ks.test(YII_Spi$value, YII_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  YII_Spi$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_YII_Spi <- list(YII_Spi2, YII_Spi3, YII_Spi4)
model.names <- c('YII_Spi2', 'YII_Spi3', 'YII_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_YII_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# YII_Spi2 3 -172.45       0.00      1      1 89.40
# YII_Spi4 3 -125.10      47.34      0      1 65.79
# YII_Spi3 3 -121.81      50.64      0      1 64.14





## ---- 4.07. FvFm -------------------------------------------------------------
### --- 4.07.1. Pve ------------------------------------------------------------
FvFm_Pve <- subset(FvFm, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
FvFm_Pve1 <- lm((value) ~ (conc), data = FvFm_Pve)
#view the output of the model
summary(FvFm_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = FvFm_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -70.682 -36.743  -8.475  17.206 260.117 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 106.5496     6.8482  15.559   <2e-16 ***
#   conc          0.2916     0.1524   1.914   0.0589 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 56.47 on 88 degrees of freedom
# Multiple R-squared:  0.03997,	Adjusted R-squared:  0.02906 
# F-statistic: 3.663 on 1 and 88 DF,  p-value: 0.05887

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
FvFm_Pve_wocon <- FvFm_Pve %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
FvFm_Pve2 <- lm(value_log ~ conc_log, data = FvFm_Pve_wocon)
#view the output of the model
summary(FvFm_Pve2)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = FvFm_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.00274 -0.31856 -0.01877  0.21841  1.29693 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.59066    0.05550   82.72   <2e-16 ***
#   conc_log     0.03700    0.01968    1.88   0.0643 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4299 on 70 degrees of freedom
# Multiple R-squared:  0.04806,	Adjusted R-squared:  0.03446 
# F-statistic: 3.534 on 1 and 70 DF,  p-value: 0.06427

# log-log transformed - early linearity
FvFm_Pve_wocon100 <- FvFm_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
FvFm_Pve3 <- lm(value_log ~ conc_log, data = FvFm_Pve_wocon100)
#view the output of the model
summary(FvFm_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = FvFm_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.96833 -0.32002  0.01725  0.23465  1.33134 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.556253   0.060869  74.853   <2e-16 ***
#   conc_log    -0.007833   0.032376  -0.242     0.81    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4473 on 52 degrees of freedom
# Multiple R-squared:  0.001125,	Adjusted R-squared:  -0.01808 
# F-statistic: 0.05854 on 1 and 52 DF,  p-value: 0.8098


# log-log transformed - late linearity
FvFm_Pve_wocon01 <- FvFm_Pve_wocon %>%
  subset(conc != "0.1")
FvFm_Pve4 <- lm(value_log ~ conc_log, data = FvFm_Pve_wocon01)
#view the output of the model
summary(FvFm_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = FvFm_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.87560 -0.26526 -0.00667  0.23552  1.42407 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.46352    0.09198  48.528   <2e-16 ***
#   conc_log     0.07841    0.03094   2.534   0.0143 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4275 on 52 degrees of freedom
# Multiple R-squared:  0.1099,	Adjusted R-squared:  0.09281 
# F-statistic: 6.422 on 1 and 52 DF,  p-value: 0.01432

# Logarithmic
plot(FvFm_Pve$conc, FvFm_Pve$value)
#fit the model
FvFm_Pve5 <- lm(value ~ conc_log, data = FvFm_Pve_wocon)
#view the output of the model
summary(FvFm_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = FvFm_Pve_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -73.03 -37.03 -12.33  14.81 251.35 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  109.187      6.884  15.860   <2e-16 ***
#   conc_log       3.538      2.441   1.449    0.152    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 53.33 on 70 degrees of freedom
# Multiple R-squared:  0.02913,	Adjusted R-squared:  0.01526 
# F-statistic:   2.1 on 1 and 70 DF,  p-value: 0.1517

# Exponential
ks.test(FvFm_Pve$value, FvFm_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  FvFm_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_FvFm_Pve <- list(FvFm_Pve2, FvFm_Pve3, FvFm_Pve4)
model.names <- c('FvFm_Pve2', 'FvFm_Pve3', 'FvFm_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_FvFm_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# FvFm_Pve4 3 65.90       0.00   0.92   0.92 -29.71
# FvFm_Pve3 3 70.80       4.90   0.08   1.00 -32.16
# FvFm_Pve2 3 87.08      21.18   0.00   1.00 -40.36


### --- 4.07.2. Spi ------------------------------------------------------------
# FvFm Spi
FvFm_Spi <- subset(FvFm, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
FvFm_Spi1 <- lm((value) ~ (conc), data = FvFm_Spi)
#view the output of the model
summary(FvFm_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = FvFm_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -83.74 -38.15 -12.03  30.45 288.31 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 124.6165     6.8254  18.258   <2e-16 ***
#   conc         -0.2634     0.1519  -1.734   0.0863 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 56.29 on 88 degrees of freedom
# Multiple R-squared:  0.03306,	Adjusted R-squared:  0.02207 
# F-statistic: 3.008 on 1 and 88 DF,  p-value: 0.08633

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
FvFm_Spi_wocon <- FvFm_Spi %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
FvFm_Spi2 <- lm(value_log ~ conc_log, data = FvFm_Spi_wocon)
#view the output of the model
summary(FvFm_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = FvFm_Spi_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.97289 -0.23591 -0.00889  0.27210  1.25004 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.71215    0.05142  91.639   <2e-16 ***
#   conc_log    -0.02650    0.01823  -1.453    0.151    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3983 on 70 degrees of freedom
# Multiple R-squared:  0.02929,	Adjusted R-squared:  0.01542 
# F-statistic: 2.112 on 1 and 70 DF,  p-value: 0.1506

# log-log transformed - early linearity
FvFm_Spi_wocon100 <- FvFm_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
FvFm_Spi3 <- lm(value_log ~ conc_log, data = FvFm_Spi_wocon100)
#view the output of the model
summary(FvFm_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = FvFm_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.73327 -0.23881 -0.03938  0.23225  1.29788 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.736072   0.053726  88.152   <2e-16 ***
#   conc_log    0.004664   0.028577   0.163    0.871    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3948 on 52 degrees of freedom
# Multiple R-squared:  0.0005119,	Adjusted R-squared:  -0.01871 
# F-statistic: 0.02663 on 1 and 52 DF,  p-value: 0.871


# log-log transformed - late linearity
FvFm_Spi_wocon01 <- FvFm_Spi_wocon %>%
  subset(conc != "0.1")
FvFm_Spi4 <- lm(value_log ~ conc_log, data = FvFm_Spi_wocon01)
#view the output of the model
summary(FvFm_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = FvFm_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.96409 -0.23485  0.02376  0.27049  0.71689 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.72975    0.07885  59.983   <2e-16 ***
#   conc_log    -0.03223    0.02653  -1.215     0.23    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3665 on 52 degrees of freedom
# Multiple R-squared:  0.02761,	Adjusted R-squared:  0.008911 
# F-statistic: 1.477 on 1 and 52 DF,  p-value: 0.2298


# Logarithmic
plot(FvFm_Spi$conc, FvFm_Spi$value)
#fit the model
FvFm_Spi5 <- lm(value ~ log(conc), data = FvFm_Spi_wocon)
#view the output of the model
summary(FvFm_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = FvFm_Spi_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -76.59 -33.82 -10.21  25.66 282.14 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  121.725      6.847  17.779   <2e-16 ***
#   log(conc)     -3.924      2.428  -1.616    0.111    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 53.03 on 70 degrees of freedom
# Multiple R-squared:  0.03598,	Adjusted R-squared:  0.02221 
# F-statistic: 2.613 on 1 and 70 DF,  p-value: 0.1105

# Exponential
ks.test(FvFm_Spi$value, FvFm_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  FvFm_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_FvFm_Spi <- list(FvFm_Spi2, FvFm_Spi3, FvFm_Spi4)
model.names <- c('FvFm_Spi2', 'FvFm_Spi3', 'FvFm_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_FvFm_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# FvFm_Spi4 3 49.27       0.00   0.98   0.98 -21.40
# FvFm_Spi3 3 57.32       8.04   0.02   1.00 -25.42
# FvFm_Spi2 3 76.09      26.82   0.00   1.00 -34.87



## ---- 4.08. rETRmax ----------------------------------------------------------
### --- 4.08.1. Pve ------------------------------------------------------------
rETR_Pve <- subset(rETR, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
rETR_Pve1 <- lm((value) ~ (conc), data = rETR_Pve)
#view the output of the model
summary(rETR_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = rETR_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -109.345  -33.855    2.065   27.208  124.708 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 179.58987    6.28012  28.597   <2e-16 ***
#   conc         -0.02225    0.13972  -0.159    0.874    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 51.79 on 88 degrees of freedom
# Multiple R-squared:  0.000288,	Adjusted R-squared:  -0.01107 
# F-statistic: 0.02535 on 1 and 88 DF,  p-value: 0.8739

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
rETR_Pve_wocon <- rETR_Pve %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
rETR_Pve2 <- lm(value_log ~ conc_log, data = rETR_Pve_wocon)
#view the output of the model
summary(rETR_Pve2)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = rETR_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.87688 -0.16597  0.05915  0.20400  0.53229 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 5.121795   0.041030 124.831   <2e-16 ***
#   conc_log    0.001698   0.014549   0.117    0.907    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3178 on 70 degrees of freedom
# Multiple R-squared:  0.0001945,	Adjusted R-squared:  -0.01409 
# F-statistic: 0.01362 on 1 and 70 DF,  p-value: 0.9074

# log-log transformed - early linearity
rETR_Pve_wocon100 <- rETR_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
rETR_Pve3 <- lm(value_log ~ conc_log, data = rETR_Pve_wocon100)
#view the output of the model
summary(rETR_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = rETR_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.86680 -0.18634  0.04533  0.19797  0.53481 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.119274   0.043179 118.559   <2e-16 ***
#   conc_log    -0.001587   0.022967  -0.069    0.945    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3173 on 52 degrees of freedom
# Multiple R-squared:  9.179e-05,	Adjusted R-squared:  -0.01914 
# F-statistic: 0.004774 on 1 and 52 DF,  p-value: 0.9452


# log-log transformed - late linearity
rETR_Pve_wocon01 <- rETR_Pve_wocon %>%
  subset(conc != "0.1")
rETR_Pve4 <- lm(value_log ~ conc_log, data = rETR_Pve_wocon01)
#view the output of the model
summary(rETR_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = rETR_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.88329 -0.14544  0.06062  0.20193  0.50664 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.147443   0.071735  71.756   <2e-16 ***
#   conc_log    -0.006656   0.024132  -0.276    0.784    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3334 on 52 degrees of freedom
# Multiple R-squared:  0.001461,	Adjusted R-squared:  -0.01774 
# F-statistic: 0.07608 on 1 and 52 DF,  p-value: 0.7838

# Logarithmic
plot(rETR_Pve$conc, rETR_Pve$value)
#fit the model
rETR_Pve5 <- lm(value ~ conc_log, data = rETR_Pve_wocon)
#view the output of the model
summary(rETR_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = rETR_Pve_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -106.379  -33.431    2.466   30.469  110.385 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 175.0712     6.4759  27.034   <2e-16 ***
#   conc_log      0.5778     2.2964   0.252    0.802    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 50.16 on 70 degrees of freedom
# Multiple R-squared:  0.0009035,	Adjusted R-squared:  -0.01337 
# F-statistic: 0.06331 on 1 and 70 DF,  p-value: 0.8021

# Exponential
ks.test(rETR_Pve$value, rETR_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  rETR_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_rETR_Pve <- list(rETR_Pve2, rETR_Pve3, rETR_Pve4)
model.names <- c('rETR_Pve2', 'rETR_Pve3', 'rETR_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_rETR_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# rETR_Pve3 3 33.71       0.00   0.93   0.93 -13.62
# rETR_Pve4 3 39.06       5.34   0.06   0.99 -16.29
# rETR_Pve2 3 43.59       9.87   0.01   1.00 -18.62


### --- 4.08.2. Spi ------------------------------------------------------------
# rETR Spi
rETR_Spi <- subset(rETR, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
rETR_Spi1 <- lm((value) ~ (conc), data = rETR_Spi)
#view the output of the model
summary(rETR_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = rETR_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -122.732  -34.476   -5.678   40.913   99.453 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.446e+02  6.006e+00  24.074   <2e-16 ***
#   conc        3.696e-03  1.336e-01   0.028    0.978    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 49.52 on 88 degrees of freedom
# Multiple R-squared:  8.695e-06,	Adjusted R-squared:  -0.01135 
# F-statistic: 0.0007651 on 1 and 88 DF,  p-value: 0.978

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
rETR_Spi_wocon <- rETR_Spi %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
rETR_Spi2 <- lm(value_log ~ conc_log, data = rETR_Spi_wocon)
#view the output of the model
summary(rETR_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = rETR_Spi_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.7887 -0.1681  0.0210  0.2904  0.6255 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.870374   0.052939  91.999   <2e-16 ***
#   conc_log    0.004178   0.018772   0.223    0.825    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4101 on 70 degrees of freedom
# Multiple R-squared:  0.0007072,	Adjusted R-squared:  -0.01357 
# F-statistic: 0.04954 on 1 and 70 DF,  p-value: 0.8245

# log-log transformed - early linearity
rETR_Spi_wocon100 <- rETR_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
rETR_Spi3 <- lm(value_log ~ conc_log, data = rETR_Spi_wocon100)
#view the output of the model
summary(rETR_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = rETR_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.19697 -0.18613 -0.02382  0.22540  0.62930 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.872275   0.049568  98.295   <2e-16 ***
#   conc_log    0.006655   0.026365   0.252    0.802    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3642 on 52 degrees of freedom
# Multiple R-squared:  0.001224,	Adjusted R-squared:  -0.01798 
# F-statistic: 0.06371 on 1 and 52 DF,  p-value: 0.8017


# log-log transformed - late linearity
rETR_Spi_wocon01 <- rETR_Spi_wocon %>%
  subset(conc != "0.1")
rETR_Spi4 <- lm(value_log ~ conc_log, data = rETR_Spi_wocon01)
#view the output of the model
summary(rETR_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = rETR_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.79329 -0.15359  0.06252  0.28121  0.61977 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.861165   0.092385  52.619   <2e-16 ***
#   conc_log    0.007178   0.031078   0.231    0.818    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.4294 on 52 degrees of freedom
# Multiple R-squared:  0.001025,	Adjusted R-squared:  -0.01819 
# F-statistic: 0.05334 on 1 and 52 DF,  p-value: 0.8183


# Logarithmic
plot(rETR_Spi$conc, rETR_Spi$value)
#fit the model
rETR_Spi5 <- lm(value ~ log(conc), data = rETR_Spi_wocon)
#view the output of the model
summary(rETR_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = rETR_Spi_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -121.701  -28.846   -6.634   34.754  104.452 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  139.239      6.198  22.467   <2e-16 ***
#   log(conc)      1.016      2.198   0.462    0.645    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 48.01 on 70 degrees of freedom
# Multiple R-squared:  0.003046,	Adjusted R-squared:  -0.0112 
# F-statistic: 0.2138 on 1 and 70 DF,  p-value: 0.6452

# Exponential
ks.test(rETR_Spi$value, rETR_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  rETR_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_rETR_Spi <- list(rETR_Spi2, rETR_Spi3, rETR_Spi4)
model.names <- c('rETR_Spi2', 'rETR_Spi3', 'rETR_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_rETR_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# rETR_Spi3 3 48.62       0.00      1      1 -21.07
# rETR_Spi4 3 66.38      17.76      0      1 -29.95
# rETR_Spi2 3 80.28      31.67      0      1 -36.97




## ---- 4.09. Ek ---------------------------------------------------------------
### --- 4.09.1. Pve ------------------------------------------------------------
Ek_Pve <- subset(Ek, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
Ek_Pve1 <- lm((value) ~ (conc), data = Ek_Pve)
#view the output of the model
summary(Ek_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = Ek_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -187.740  -61.108   -0.411   58.759  255.741 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 319.31588   11.54176  27.666   <2e-16 ***
#   conc         -0.07852    0.25679  -0.306     0.76    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 95.18 on 88 degrees of freedom
# Multiple R-squared:  0.001061,	Adjusted R-squared:  -0.01029 
# F-statistic: 0.0935 on 1 and 88 DF,  p-value: 0.7605

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
Ek_Pve_wocon <- Ek_Pve %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
Ek_Pve2 <- lm(value_log ~ conc_log, data = Ek_Pve_wocon)
#view the output of the model
summary(Ek_Pve2)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = Ek_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.82094 -0.15392  0.04632  0.23632  0.61544 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.6968168  0.0410720 138.703   <2e-16 ***
#   conc_log    -0.0009916  0.0145641  -0.068    0.946    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3181 on 70 degrees of freedom
# Multiple R-squared:  6.622e-05,	Adjusted R-squared:  -0.01422 
# F-statistic: 0.004636 on 1 and 70 DF,  p-value: 0.9459

# log-log transformed - early linearity
Ek_Pve_wocon100 <- Ek_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
Ek_Pve3 <- lm(value_log ~ conc_log, data = Ek_Pve_wocon100)
#view the output of the model
summary(Ek_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = Ek_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.81120 -0.16149  0.04859  0.17861  0.61788 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.694383   0.043228 131.728   <2e-16 ***
#   conc_log    -0.004163   0.022993  -0.181    0.857    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3177 on 52 degrees of freedom
# Multiple R-squared:  0.0006299,	Adjusted R-squared:  -0.01859 
# F-statistic: 0.03278 on 1 and 52 DF,  p-value: 0.857


# log-log transformed - late linearity
Ek_Pve_wocon01 <- Ek_Pve_wocon %>%
  subset(conc != "0.1")
Ek_Pve4 <- lm(value_log ~ conc_log, data = Ek_Pve_wocon01)
#view the output of the model
summary(Ek_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = Ek_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.83018 -0.12148  0.03205  0.26583  0.57847 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.73379    0.07220  79.415   <2e-16 ***
#   conc_log    -0.01304    0.02429  -0.537    0.594    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3356 on 52 degrees of freedom
# Multiple R-squared:  0.005509,	Adjusted R-squared:  -0.01362 
# F-statistic: 0.288 on 1 and 52 DF,  p-value: 0.5938

# Logarithmic
plot(Ek_Pve$conc, Ek_Pve$value)
#fit the model
Ek_Pve5 <- lm(value ~ conc_log, data = Ek_Pve_wocon)
#view the output of the model
summary(Ek_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = Ek_Pve_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -181.308  -54.581   -0.762   63.331  239.949 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 311.3399    11.9502  26.053   <2e-16 ***
#   conc_log      0.3295     4.2376   0.078    0.938    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 92.57 on 70 degrees of freedom
# Multiple R-squared:  8.635e-05,	Adjusted R-squared:  -0.0142 
# F-statistic: 0.006045 on 1 and 70 DF,  p-value: 0.9382

# Exponential
ks.test(Ek_Pve$value, Ek_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  Ek_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_Ek_Pve <- list(Ek_Pve2, Ek_Pve3, Ek_Pve4)
model.names <- c('Ek_Pve2', 'Ek_Pve3', 'Ek_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_Ek_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Ek_Pve3 3 33.84       0.00   0.94   0.94 -13.68
# Ek_Pve4 3 39.76       5.92   0.05   0.99 -16.64
# Ek_Pve2 3 43.74       9.90   0.01   1.00 -18.69


### --- 4.09.2. Spi ------------------------------------------------------------
# Ek Spi
Ek_Spi <- subset(Ek, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
Ek_Spi1 <- lm((value) ~ (conc), data = Ek_Spi)
#view the output of the model
summary(Ek_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = Ek_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -185.15  -61.34    2.11   47.24  166.52 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 276.25397    9.37645  29.463   <2e-16 ***
#   conc          0.05082    0.20861   0.244    0.808    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 77.32 on 88 degrees of freedom
# Multiple R-squared:  0.0006741,	Adjusted R-squared:  -0.01068 
# F-statistic: 0.05936 on 1 and 88 DF,  p-value: 0.8081

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
Ek_Spi_wocon <- Ek_Spi %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
Ek_Spi2 <- lm(value_log ~ conc_log, data = Ek_Spi_wocon)
#view the output of the model
summary(Ek_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = Ek_Spi_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.0048 -0.1878  0.0053  0.1994  0.5313 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 5.554811   0.038252  145.22   <2e-16 ***
#   conc_log    0.003531   0.013564    0.26    0.795    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2963 on 70 degrees of freedom
# Multiple R-squared:  0.0009669,	Adjusted R-squared:  -0.01331 
# F-statistic: 0.06775 on 1 and 70 DF,  p-value: 0.7954

# log-log transformed - early linearity
Ek_Spi_wocon100 <- Ek_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
Ek_Spi3 <- lm(value_log ~ conc_log, data = Ek_Spi_wocon100)
#view the output of the model
summary(Ek_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = Ek_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.62223 -0.18967  0.01774  0.16377  0.55116 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.549836   0.036443 152.286   <2e-16 ***
#   conc_log    -0.002952   0.019384  -0.152     0.88    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.2678 on 52 degrees of freedom
# Multiple R-squared:  0.0004457,	Adjusted R-squared:  -0.01878 
# F-statistic: 0.02319 on 1 and 52 DF,  p-value: 0.8796


# log-log transformed - late linearity
Ek_Spi_wocon01 <- Ek_Spi_wocon %>%
  subset(conc != "0.1")
Ek_Spi4 <- lm(value_log ~ conc_log, data = Ek_Spi_wocon01)
#view the output of the model
summary(Ek_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = Ek_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.00834 -0.18700  0.00979  0.20177  0.53306 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.54763    0.06619  83.816   <2e-16 ***
#   conc_log     0.00587    0.02227   0.264    0.793    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.3076 on 52 degrees of freedom
# Multiple R-squared:  0.001335,	Adjusted R-squared:  -0.01787 
# F-statistic: 0.0695 on 1 and 52 DF,  p-value: 0.7931


# Logarithmic
plot(Ek_Spi$conc, Ek_Spi$value)
#fit the model
Ek_Spi5 <- lm(value ~ log(conc), data = Ek_Spi_wocon)
#view the output of the model
summary(Ek_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = Ek_Spi_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -180.234  -54.928   -9.418   45.653  170.911 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  268.314      9.841  27.266   <2e-16 ***
#   log(conc)      1.761      3.490   0.505    0.615    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 76.23 on 70 degrees of freedom
# Multiple R-squared:  0.003624,	Adjusted R-squared:  -0.01061 
# F-statistic: 0.2546 on 1 and 70 DF,  p-value: 0.6154

# Exponential
ks.test(Ek_Spi$value, Ek_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  Ek_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_Ek_Spi <- list(Ek_Spi2, Ek_Spi3, Ek_Spi4)
model.names <- c('Ek_Spi2', 'Ek_Spi3', 'Ek_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_Ek_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# Ek_Spi3 3 15.40       0.00      1      1  -4.46
# Ek_Spi4 3 30.37      14.97      0      1 -11.94
# Ek_Spi2 3 33.49      18.09      0      1 -13.57





## ---- 4.10. Alpha ------------------------------------------------------------
### --- 4.10.1. Pve ------------------------------------------------------------
alpha_Pve <- subset(alpha, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
alpha_Pve1 <- lm((value) ~ (conc), data = alpha_Pve)
#view the output of the model
summary(alpha_Pve1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = alpha_Pve)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.139955 -0.034859  0.006227  0.035919  0.130337 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 5.654e-01  6.490e-03  87.118   <2e-16 ***
#   conc        7.233e-05  1.444e-04   0.501    0.618    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.05352 on 88 degrees of freedom
# Multiple R-squared:  0.002843,	Adjusted R-squared:  -0.008488 
# F-statistic: 0.2509 on 1 and 88 DF,  p-value: 0.6177

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
alpha_Pve_wocon <- alpha_Pve %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
alpha_Pve2 <- lm(value_log ~ conc_log, data = alpha_Pve_wocon)
#view the output of the model
summary(alpha_Pve2)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = alpha_Pve_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.27724 -0.05646  0.01523  0.06489  0.20708 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.575022   0.011882 -48.394   <2e-16 ***
#   conc_log     0.002690   0.004213   0.638    0.525    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.09204 on 70 degrees of freedom
# Multiple R-squared:  0.005787,	Adjusted R-squared:  -0.008416 
# F-statistic: 0.4075 on 1 and 70 DF,  p-value: 0.5253

# log-log transformed - early linearity
alpha_Pve_wocon100 <- alpha_Pve_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
alpha_Pve3 <- lm(value_log ~ conc_log, data = alpha_Pve_wocon100)
#view the output of the model
summary(alpha_Pve3)
# OUTPUT: 	Call:
# lm(formula = log(value) ~ log(conc), data = alpha_Pve_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.27715 -0.05266  0.01322  0.06318  0.20742 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.575109   0.012677 -45.367   <2e-16 ***
#   conc_log     0.002576   0.006743   0.382    0.704    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.09316 on 52 degrees of freedom
# Multiple R-squared:  0.002799,	Adjusted R-squared:  -0.01638 
# F-statistic: 0.1459 on 1 and 52 DF,  p-value: 0.704


# log-log transformed - late linearity
alpha_Pve_wocon01 <- alpha_Pve_wocon %>%
  subset(conc != "0.1")
alpha_Pve4 <- lm(value_log ~ conc_log, data = alpha_Pve_wocon01)
#view the output of the model
summary(alpha_Pve4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = alpha_Pve_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.26591 -0.05247  0.02130  0.05798  0.20991 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.586349   0.019482 -30.096   <2e-16 ***
#   conc_log     0.006379   0.006554   0.973    0.335    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.09055 on 52 degrees of freedom
# Multiple R-squared:  0.01789,	Adjusted R-squared:  -0.0009942 
# F-statistic: 0.9474 on 1 and 52 DF,  p-value: 0.3349

# Logarithmic
plot(alpha_Pve$conc, alpha_Pve$value)
#fit the model
alpha_Pve5 <- lm(value ~ conc_log, data = alpha_Pve_wocon)
#view the output of the model
summary(alpha_Pve5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = alpha_Pve_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -0.138607 -0.033285  0.006113  0.035589  0.128022 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.565059   0.006600  85.620   <2e-16 ***
#   conc_log    0.001465   0.002340   0.626    0.533    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.05112 on 70 degrees of freedom
# Multiple R-squared:  0.005569,	Adjusted R-squared:  -0.008638 
# F-statistic: 0.392 on 1 and 70 DF,  p-value: 0.5333

# Exponential
ks.test(alpha_Pve$value, alpha_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  alpha_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_alpha_Pve <- list(alpha_Pve2, alpha_Pve3, alpha_Pve4)
model.names <- c('alpha_Pve2', 'alpha_Pve3', 'alpha_Pve4')

# test for best fitted log-log model
aictab(cand.set = models_alpha_Pve, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# alpha_Pve2 3 -134.87       0.00      1      1 70.61
# alpha_Pve4 3 -101.72      33.15      0      1 54.10
# alpha_Pve3 3  -98.65      36.22      0      1 52.56


### --- 4.10.2. Spi ------------------------------------------------------------
# alpha Spi
alpha_Spi <- subset(alpha, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
alpha_Spi1 <- lm((value) ~ (conc), data = alpha_Spi)
#view the output of the model
summary(alpha_Spi1)
# OUTPUT: 	Call:
# lm(formula = (value) ~ (conc), data = alpha_Spi)
# Residuals:
#        Min       1Q   Median       3Q      Max 
# -0.275674 -0.036443 -0.000562  0.060744  0.134579 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.5167637  0.0100269  51.538   <2e-16 ***
#   conc        -0.0001010  0.0002231  -0.453    0.652    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.08269 on 88 degrees of freedom
# Multiple R-squared:  0.002324,	Adjusted R-squared:  -0.009014 
# F-statistic: 0.205 on 1 and 88 DF,  p-value: 0.6519

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
alpha_Spi_wocon <- alpha_Spi %>%
  subset(conc != "0") %>%
  filter(value_log != "NaN")
alpha_Spi2 <- lm(value_log ~ conc_log, data = alpha_Spi_wocon)
#view the output of the model
summary(alpha_Spi2)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = alpha_Spi_wocon)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.78393 -0.08051  0.01636  0.12040  0.25719 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.6844366  0.0234416 -29.198   <2e-16 ***
#   conc_log     0.0006476  0.0083124   0.078    0.938    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1816 on 70 degrees of freedom
# Multiple R-squared:  8.671e-05,	Adjusted R-squared:  -0.0142 
# F-statistic: 0.00607 on 1 and 70 DF,  p-value: 0.9381

# log-log transformed - early linearity
alpha_Spi_wocon100 <- alpha_Spi_wocon %>%
  subset(conc != "100") %>%
  filter(value_log != "NaN")
alpha_Spi3 <- lm(value_log ~ conc_log, data = alpha_Spi_wocon100)
#view the output of the model
summary(alpha_Spi3)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = alpha_Spi_wocon100)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.57474 -0.08234  0.00287  0.11591  0.27095 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.677561   0.022719 -29.824   <2e-16 ***
#   conc_log     0.009606   0.012084   0.795     0.43    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1669 on 52 degrees of freedom
# Multiple R-squared:  0.01201,	Adjusted R-squared:  -0.006993 
# F-statistic: 0.6319 on 1 and 52 DF,  p-value: 0.4303


# log-log transformed - late linearity
alpha_Spi_wocon01 <- alpha_Spi_wocon %>%
  subset(conc != "0.1")
alpha_Spi4 <- lm(value_log ~ conc_log, data = alpha_Spi_wocon01)
#view the output of the model
summary(alpha_Spi4)
# OUTPUT: 	Call:
# lm(formula = value_log ~ conc_log, data = alpha_Spi_wocon01)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.78494 -0.05358  0.01763  0.11810  0.24248 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.686465   0.040124 -17.108   <2e-16 ***
#   conc_log     0.001308   0.013498   0.097    0.923    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1865 on 52 degrees of freedom
# Multiple R-squared:  0.0001806,	Adjusted R-squared:  -0.01905 
# F-statistic: 0.009394 on 1 and 52 DF,  p-value: 0.9232


# Logarithmic
plot(alpha_Spi$conc, alpha_Spi$value)
#fit the model
alpha_Spi5 <- lm(value ~ log(conc), data = alpha_Spi_wocon)
#view the output of the model
summary(alpha_Spi5)
# OUTPUT: Call:
# lm(formula = value ~ log(conc), data = alpha_Spi_wocon)
# Residuals:
#   Min     1Q Median     3Q    Max 
# -0.282665 -0.046321  0.000879  0.057081  0.140588 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.5117148  0.0106483  48.056   <2e-16 ***
#   log(conc)   0.0004212  0.0037759   0.112    0.911    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.08248 on 70 degrees of freedom
# Multiple R-squared:  0.0001777,	Adjusted R-squared:  -0.01411 
# F-statistic: 0.01244 on 1 and 70 DF,  p-value: 0.9115

# Exponential
ks.test(alpha_Spi$value, alpha_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  alpha_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation


# AIC test
models_alpha_Spi <- list(alpha_Spi2, alpha_Spi3, alpha_Spi4)
model.names <- c('alpha_Spi2', 'alpha_Spi3', 'alpha_Spi4')

# test for best fitted log-log model
aictab(cand.set = models_alpha_Spi, modnames = model.names)
# OUTPUT: Model selection based on AICc:
# K   AICc Delta_AICc AICcWt Cum.Wt     LL
# alpha_Spi2 3 -37.02       0.00   0.67   0.67 21.69
# alpha_Spi3 3 -35.64       1.38   0.33   1.00 21.06
# alpha_Spi4 3 -23.69      13.33   0.00   1.00 15.09