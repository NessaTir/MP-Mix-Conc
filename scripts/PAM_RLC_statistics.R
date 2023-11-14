# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the data processing and statistical analyzes of the additional parameters of
# photosynthetic efficiency of the corals photosymbionts, derived from rapid light curves (RLCs).
# This script builts up on the processing of the RLC using the script 'PAM_RLC_data_processing'
# The parameters were measured using pulse amplitude modulated fluorometry (PAM) and include
#   a) relative electron transport rate (rETRmax)
#   b) efficiency of light capture (α)
#   c)  and minimum saturating irradiance (Ek)
#   Statistical analyzes will be conducted using LMER and GLMER 
#   together with a holm adjusted glht summary



# ----- 2. Load in needed packages ---------------------------------------------
# to easily clean data, to read in .rds files 
library(tidyverse)

# for statistical analyses using LMER and GLMER
library(lme4)

# for statistical testing
library(multcomp)

# check model fits statistically
library(performance)
library(rstatix)
# check model fits visually using qqplot, use glht for lmer testing
library(car)



# ---- 3. Read in needed data files --------------------------------------------
## --- 3.1. Coral identity table -----------------------------------------------
# read in list with overview of all corals used and their treatments etc.
corals <- read_csv2("in/coral_treatments.csv") %>%
  mutate(treat = as.factor(treat), # column for categorical model
         conc = as.numeric(conc)) # column for continuous model 

## --- 3.2. RLC data table -----------------------------------------------------
# read in table with processed RLC parameters
rETR_parameter <- read_rds("processed/rETR_parameter.rds")



# ----- 4. Prepare data for statistical analyzes -------------------------------
## ---- 4.1. Reformat data table -----------------------------------------------
# prepare RLC data table for clean merge with coral info table
RLC_parameter <- rETR_parameter %>% 
  separate(ID_time, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep = "_", remove = FALSE) %>%
  # remove columns to avoid doubling when merged
  select(-spec, -col, -tank)

# merge RLC parameter with coral info table
RLC_parameter <-  merge(corals, RLC_parameter, by = 'ID') %>%
  mutate(treat = as.factor(treat),
         conc = as.numeric(conc))

# level treatments
RLC_parameter$treat <- factor(RLC_parameter$treat, 
                              levels = c("control", "0.1", "1",
                                         "10", "100"))

# reaname timepoints for clean statistical analyses
RLC_parameter <- RLC_parameter %>%
  mutate(tp = case_when(tp == "t0"~ "0",
                        tp == "t1"~ "1",
                        tp == "t2"~ "2",
                        tp == "t3"~ "3"),
         tp = as.numeric(tp))



# ----- 5. Statistics for RLC parameter ----------------------------------------
# Create subsets for the statistical analyses
# create a subset with data of Pve, t0 excluded for continuous model --> Overall effect
Pve_overall_effect <- subset(RLC_parameter, spec == "Pve" & tp!= "0")
# create a subset with data of Pve, t1 excluded for continuous model --> Specific effect: t1
Pve_RLC_t1 <- subset(RLC_parameter, spec == "Pve" & tp == "1")
# create a subset with data of Pve, t2 excluded for continuous model --> Specific effect: t2
Pve_RLC_t2 <- subset(RLC_parameter, spec == "Pve" & tp == "2")
# create a subset with data of Pve, t3 excluded for continuous model --> Specific effect: t3
Pve_RLC_t3 <- subset(RLC_parameter, spec == "Pve" & tp == "3")

# create a subset with data of Spi, t0 excluded for continuous model --> Overall effect
Spi_overall_effect <- subset(RLC_parameter, spec == "Spi" & tp!= "0")
# create a subset with data of Spi, t1 excluded for continuous model --> Specific effect: t1
Spi_RLC_t1 <- subset(RLC_parameter, spec == "Spi" & tp == "1")
# create a subset with data of Spi, t2 excluded for continuous model --> Specific effect: t2
Spi_RLC_t2 <- subset(RLC_parameter, spec == "Spi" & tp == "2")
# create a subset with data of Spi, t3 excluded for continuous model --> Specific effect: t3
Spi_RLC_t3 <- subset(RLC_parameter, spec == "Spi" & tp == "3")

## ---- 5.1. rETRmax -----------------------------------------------------------
### --- 5.1.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.1.1.1 Overall effect -------------------------------------------------
# write LMER
model1_Pve <- lmer(scale(rETRmax) ~ conc + (1|col) + (1|tp), data = Pve_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Pve))          # good fit
shapiro_test(residuals(model1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Pve)     0.992   0.176
check_normality(model1_Pve)
# OK: residuals appear as normally distributed (p = 0.165).
check_heteroscedasticity(model1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.861).

# get summary of LMER
cftest(model1_Pve)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(rETRmax) ~ conc + (1 | col) + (1 | tp), 
#           data = Pve_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0  0.0139180  0.2561474   0.054    0.957
# conc == 0        -0.0006264  0.0013274  -0.472    0.637
# (Univariate p values reported)


#### -- 5.1.1.2 Specific effects -----------------------------------------------
# ------------- t1
# write LMER
# use LMER because of good graphical fit!
model_t1_Pve <- lmer(scale(rETRmax) ~ treat + (1|col), data = Pve_RLC_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.971  0.0436
check_normality(model_t1_Pve)
# Warning: Non-normality of residuals detected (p = 0.044).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.904).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(rETRmax) ~ treat + (1 | col), data = Pve_RLC_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.04935    0.28137   0.175    1.000
# 1 - control == 0   -0.20869    0.28137  -0.742    1.000
# 10 - control == 0  -0.14952    0.28137  -0.531    1.000
# 100 - control == 0 -0.47268    0.28137  -1.680    0.837
# 1 - 0.1 == 0       -0.25805    0.28137  -0.917    1.000
# 10 - 0.1 == 0      -0.19887    0.28137  -0.707    1.000
# 100 - 0.1 == 0     -0.52203    0.28137  -1.855    0.636
# 10 - 1 == 0         0.05918    0.28137   0.210    1.000
# 100 - 1 == 0       -0.26398    0.28137  -0.938    1.000
# 100 - 10 == 0      -0.32316    0.28137  -1.149    1.000
# (Adjusted p values reported -- holm method)


# ------------- t2
# write LMER
model_t2_Pve <- lmer(scale(rETRmax) ~ treat + (1|col), data = Pve_RLC_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.988   0.598
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.598).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.626).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(rETRmax) ~ treat + (1 | col), data = Pve_RLC_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.3278274  0.2731202  -1.200        1
# 1 - control == 0   -0.2853632  0.2731202  -1.045        1
# 10 - control == 0  -0.2846304  0.2731202  -1.042        1
# 100 - control == 0  0.0839061  0.2731202   0.307        1
# 1 - 0.1 == 0        0.0424642  0.2731202   0.155        1
# 10 - 0.1 == 0       0.0431970  0.2731202   0.158        1
# 100 - 0.1 == 0      0.4117335  0.2731202   1.508        1
# 10 - 1 == 0         0.0007328  0.2731202   0.003        1
# 100 - 1 == 0        0.3692693  0.2731202   1.352        1
# 100 - 10 == 0       0.3685365  0.2731202   1.349        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# write LMER 
model_t3_Pve <- lmer(scale(rETRmax) ~ treat + (1|col), data = Pve_RLC_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.988   0.569
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.569).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.804).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(rETRmax) ~ treat + (1 | col), data = Pve_RLC_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.45109    0.27975  -1.612        1
# 1 - control == 0   -0.15836    0.27975  -0.566        1
# 10 - control == 0  -0.41575    0.27975  -1.486        1
# 100 - control == 0 -0.27919    0.27975  -0.998        1
# 1 - 0.1 == 0        0.29274    0.27975   1.046        1
# 10 - 0.1 == 0       0.03534    0.27975   0.126        1
# 100 - 0.1 == 0      0.17190    0.27975   0.614        1
# 10 - 1 == 0        -0.25739    0.27975  -0.920        1
# 100 - 1 == 0       -0.12084    0.27975  -0.432        1
# 100 - 10 == 0       0.13655    0.27975   0.488        1
# (Adjusted p values reported -- holm method)



### --- 5.1.2. Stylophora pistillata -------------------------------------------
#### -- 5.1.2.1 Overall effect -------------------------------------------------
# write LMER
model1_Spi <- lmer(scale(rETRmax) ~ conc + (1|col) + (1|tp), data = Spi_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Spi))          # good fit
shapiro_test(residuals(model1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Spi)     0.995   0.449
check_normality(model1_Spi)
# OK: residuals appear as normally distributed (p = 0.447).
check_heteroscedasticity(model1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.728).

# get summary of LMER
cftest(model1_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(rETRmax) ~ conc + (1 | col) + (1 | tp), 
#           data = Spi_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0 -0.022608   0.319191  -0.071    0.944
# conc == 0         0.001017   0.001306   0.779    0.436
# (Univariate p values reported)


#### -- 5.1.2.2 Specific effects -----------------------------------------------
# ------------- t1
# write LMER
model_t1_Spi <- lmer(scale(rETRmax) ~ treat + (1|col), data = Spi_RLC_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT:  A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.979   0.163  
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.163).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.939).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(rETRmax) ~ treat + (1 | col), data = Spi_RLC_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0  0.39689    0.25796   1.539   0.8673  
# 1 - control == 0    0.69621    0.25796   2.699   0.0696 .
# 10 - control == 0   0.64449    0.25796   2.498   0.1123  
# 100 - control == 0  0.61684    0.25796   2.391   0.1343  
# 1 - 0.1 == 0        0.29933    0.25796   1.160   1.0000  
# 10 - 0.1 == 0       0.24761    0.25796   0.960   1.0000  
# 100 - 0.1 == 0      0.21995    0.25796   0.853   1.0000  
# 10 - 1 == 0        -0.05172    0.25796  -0.201   1.0000  
# 100 - 1 == 0       -0.07937    0.25796  -0.308   1.0000  
# 100 - 10 == 0      -0.02765    0.25796  -0.107   1.0000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# write LMER 
model_t2_Spi <- lmer(scale(rETRmax) ~ treat + (1|col), data = Spi_RLC_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.981   0.213
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.213).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.595).

# get summary of LMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(rETRmax) ~ treat + (1 | col), data = Spi_RLC_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.631750   0.306141   2.064    0.391
# 1 - control == 0    0.190897   0.306141   0.624    1.000
# 10 - control == 0   0.194239   0.306141   0.634    1.000
# 100 - control == 0  0.353810   0.306141   1.156    1.000
# 1 - 0.1 == 0       -0.440853   0.306141  -1.440    1.000
# 10 - 0.1 == 0      -0.437511   0.306141  -1.429    1.000
# 100 - 0.1 == 0     -0.277940   0.306141  -0.908    1.000
# 10 - 1 == 0         0.003342   0.306141   0.011    1.000
# 100 - 1 == 0        0.162913   0.306141   0.532    1.000
# 100 - 10 == 0       0.159571   0.306141   0.521    1.000
# (Adjusted p values reported -- holm method)


# ------------- t3
# write LMER 
model_t3_Spi <- lmer(scale(rETRmax) ~ treat + (1|col), data = Spi_RLC_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.981   0.201
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.201).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.959).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(rETRmax) ~ treat + (1 | col), data = Spi_RLC_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.48011    0.30314  -1.584        1
# 1 - control == 0   -0.46898    0.30314  -1.547        1
# 10 - control == 0  -0.45047    0.30314  -1.486        1
# 100 - control == 0 -0.32790    0.30314  -1.082        1
# 1 - 0.1 == 0        0.01113    0.30314   0.037        1
# 10 - 0.1 == 0       0.02965    0.30314   0.098        1
# 100 - 0.1 == 0      0.15222    0.30314   0.502        1
# 10 - 1 == 0         0.01852    0.30314   0.061        1
# 100 - 1 == 0        0.14108    0.30314   0.465        1
# 100 - 10 == 0       0.12257    0.30314   0.404        1
# (Adjusted p values reported -- holm method)



## ---- 5.2. Ek ----------------------------------------------------------------
### --- 5.2.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.2.1.1 Overall effect -------------------------------------------------
# LMER didn't show a good fit, therefore GLMER is used
model1_Pve <- glmer((Ek) ~ conc + (1|col) + (1|tp), family = poisson, data = Pve_overall_effect)

# get summary of GLMER
cftest(model1_Pve)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (Ek) ~ conc + (1 | col) + (1 | tp), data = Pve_overall_effect, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0  5.803e+00  7.071e-01   8.207 2.22e-16 ***
#   conc == 0        -3.460e-04  8.589e-05  -4.029 5.60e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Univariate p values reported)



#### -- 5.2.1.2 Specific effects -----------------------------------------------
# ------------- t1
# write LMER
# use LMER because of good graphical fit!
model_t1_Pve <- lmer(scale(Ek) ~ treat + (1|col), data = Pve_RLC_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # okay fit
shapiro_test(residuals(model_t1_Pve))    # p < 0.05 = NO Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.951 0.00191
check_normality(model_t1_Pve)
# Warning: Non-normality of residuals detected (p = 0.002).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.695).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Ek) ~ treat + (1 | col), data = Pve_RLC_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.02149    0.24546  -0.088    1.000
# 1 - control == 0   -0.04400    0.24546  -0.179    1.000
# 10 - control == 0  -0.23938    0.24546  -0.975    1.000
# 100 - control == 0 -0.50679    0.24546  -2.065    0.390
# 1 - 0.1 == 0       -0.02251    0.24546  -0.092    1.000
# 10 - 0.1 == 0      -0.21789    0.24546  -0.888    1.000
# 100 - 0.1 == 0     -0.48530    0.24546  -1.977    0.432
# 10 - 1 == 0        -0.19538    0.24546  -0.796    1.000
# 100 - 1 == 0       -0.46280    0.24546  -1.885    0.475
# 100 - 10 == 0      -0.26742    0.24546  -1.089    1.000
# (Adjusted p values reported -- holm method)


# ------------- t2
# write LMER
model_t2_Pve <- lmer(scale(Ek) ~ treat + (1|col), data = Pve_RLC_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.979   0.153
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.153).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.547).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Ek) ~ treat + (1 | col), data = Pve_RLC_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.265400   0.239317  -1.109        1
# 1 - control == 0   -0.275560   0.239317  -1.151        1
# 10 - control == 0  -0.258453   0.239317  -1.080        1
# 100 - control == 0  0.043459   0.239317   0.182        1
# 1 - 0.1 == 0       -0.010160   0.239317  -0.042        1
# 10 - 0.1 == 0       0.006947   0.239317   0.029        1
# 100 - 0.1 == 0      0.308859   0.239317   1.291        1
# 10 - 1 == 0         0.017107   0.239317   0.071        1
# 100 - 1 == 0        0.319019   0.239317   1.333        1
# 100 - 10 == 0       0.301912   0.239317   1.262        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# write LMER 
model_t3_Pve <- lmer(scale(Ek) ~ treat + (1|col), data = Pve_RLC_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.981   0.210
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.210).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.956).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Ek) ~ treat + (1 | col), data = Pve_RLC_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.43787    0.24911  -1.758    0.788
# 1 - control == 0   -0.07624    0.24911  -0.306    1.000
# 10 - control == 0  -0.42795    0.24911  -1.718    0.788
# 100 - control == 0 -0.29393    0.24911  -1.180    1.000
# 1 - 0.1 == 0        0.36164    0.24911   1.452    1.000
# 10 - 0.1 == 0       0.00992    0.24911   0.040    1.000
# 100 - 0.1 == 0      0.14394    0.24911   0.578    1.000
# 10 - 1 == 0        -0.35172    0.24911  -1.412    1.000
# 100 - 1 == 0       -0.21769    0.24911  -0.874    1.000
# 100 - 10 == 0       0.13402    0.24911   0.538    1.000
# (Adjusted p values reported -- holm method)


### --- 5.2.2. Stylophora pistillata -------------------------------------------
#### -- 5.2.2.1 Overall effect -------------------------------------------------
# LMER didn't show a good fit, therefore GLMER is used
model1_Spi <- glmer((Ek) ~ conc + (1|col) + (1|tp), family = poisson, data = Spi_overall_effect)

# check GLMER  
cftest(model1_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (Ek) ~ conc + (1 | col) + (1 | tp), data = Spi_overall_effect, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0 5.725e+00  7.071e-01   8.096 6.66e-16 ***
#   conc == 0        2.224e-04  8.752e-05   2.541    0.011 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Univariate p values reported)



#### -- 5.2.2.2 Specific effects -----------------------------------------------
# ------------- t1
# LMER didn't show a good fit, therefore GLMER is used
model_t1_Spi <- glmer((Ek) ~ treat + (1|col), family = poisson, data = Spi_RLC_t1)

# get summary of GLMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (Ek) ~ treat + (1 | col), data = Spi_RLC_t1, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0 0.043845   0.018029   2.432  0.10513   
#   1 - control == 0   0.061208   0.017953   3.409  0.00521 **
#   10 - control == 0  0.064075   0.017940   3.572  0.00319 **
#   100 - control == 0 0.067530   0.017925   3.767  0.00165 **
#   1 - 0.1 == 0       0.017363   0.017753   0.978  1.00000   
#   10 - 0.1 == 0      0.020230   0.017741   1.140  1.00000   
#   100 - 0.1 == 0     0.023685   0.017725   1.336  1.00000   
#   10 - 1 == 0        0.002868   0.017663   0.162  1.00000   
#   100 - 1 == 0       0.006322   0.017648   0.358  1.00000   
#   100 - 10 == 0      0.003455   0.017635   0.196  1.00000   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# write LMER
model_t2_Spi <- lmer(scale(log(Ek)) ~ treat + (1|col), data = Spi_RLC_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.987   0.515
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.515).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.855).

# get summary of LMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(log(Ek)) ~ treat + (1 | col), data = Spi_RLC_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.50910    0.27721   1.836    0.663
# 1 - control == 0    0.12231    0.27721   0.441    1.000
# 10 - control == 0   0.19656    0.27721   0.709    1.000
# 100 - control == 0  0.28164    0.27721   1.016    1.000
# 1 - 0.1 == 0       -0.38679    0.27721  -1.395    1.000
# 10 - 0.1 == 0      -0.31254    0.27721  -1.127    1.000
# 100 - 0.1 == 0     -0.22746    0.27721  -0.821    1.000
# 10 - 1 == 0         0.07426    0.27721   0.268    1.000
# 100 - 1 == 0        0.15933    0.27721   0.575    1.000
# 100 - 10 == 0       0.08508    0.27721   0.307    1.000
# (Adjusted p values reported -- holm method)


# ------------- t3
# write LMER 
model_t3_Spi <- lmer(scale(Ek) ~ treat + (1|col), data = Spi_RLC_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.989   0.628
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.628).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.837).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Ek) ~ treat + (1 | col), data = Spi_RLC_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.50753    0.28268  -1.795    0.653
# 1 - control == 0   -0.45830    0.28268  -1.621    0.840
# 10 - control == 0  -0.57103    0.28268  -2.020    0.434
# 100 - control == 0 -0.29424    0.28268  -1.041    1.000
# 1 - 0.1 == 0        0.04923    0.28268   0.174    1.000
# 10 - 0.1 == 0      -0.06350    0.28268  -0.225    1.000
# 100 - 0.1 == 0      0.21329    0.28268   0.755    1.000
# 10 - 1 == 0        -0.11273    0.28268  -0.399    1.000
# 100 - 1 == 0        0.16406    0.28268   0.580    1.000
# 100 - 10 == 0       0.27679    0.28268   0.979    1.000
# (Adjusted p values reported -- holm method)



## ---- 5.3. Alpha -------------------------------------------------------------
### --- 5.3.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.3.1.1 Overall effect -------------------------------------------------
# write LMER
model1_Pve <- lmer(scale(alpha^3) ~ conc + (1|col) + (1|tp), data = Pve_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Pve))          # good fit
shapiro_test(residuals(model1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Pve)     0.995   0.568
check_normality(model1_Pve)
# OK: residuals appear as normally distributed (p = 0.568).
check_heteroscedasticity(model1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.409).

# get summary of LMER
cftest(model1_Pve)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(alpha^3) ~ conc + (1 | col) + (1 | tp), 
#           data = Pve_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0 -0.023335   0.243290  -0.096    0.924
# conc == 0         0.001050   0.001323   0.794    0.427
# (Univariate p values reported)


#### -- 5.3.1.2 Specific effects -----------------------------------------------
# ------------- t1
# write LMER
model_t1_Pve <- lmer(scale(alpha^2) ~ treat + (1|col), data = Pve_RLC_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.976  0.0880
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.088).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.916).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(alpha^2) ~ treat + (1 | col), data = Pve_RLC_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.16227    0.26961   0.602    1.000
# 1 - control == 0   -0.39334    0.26961  -1.459    1.000
# 10 - control == 0   0.23056    0.26961   0.855    1.000
# 100 - control == 0 -0.08936    0.26961  -0.331    1.000
# 1 - 0.1 == 0       -0.55562    0.26961  -2.061    0.354
# 10 - 0.1 == 0       0.06829    0.26961   0.253    1.000
# 100 - 0.1 == 0     -0.25164    0.26961  -0.933    1.000
# 10 - 1 == 0         0.62390    0.26961   2.314    0.207
# 100 - 1 == 0        0.30398    0.26961   1.127    1.000
# 100 - 10 == 0      -0.31992    0.26961  -1.187    1.000
# (Adjusted p values reported -- holm method)


# ------------- t2
# write LMER 
model_t2_Pve <- lmer(scale(alpha) ~ treat + (1|col), data = Pve_RLC_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.976  0.0934
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.093).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.710).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(alpha) ~ treat + (1 | col), data = Pve_RLC_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.17130    0.28081  -0.610        1
# 1 - control == 0    0.09286    0.28081   0.331        1
# 10 - control == 0  -0.12846    0.28081  -0.457        1
# 100 - control == 0  0.23755    0.28081   0.846        1
# 1 - 0.1 == 0        0.26416    0.28081   0.941        1
# 10 - 0.1 == 0       0.04285    0.28081   0.153        1
# 100 - 0.1 == 0      0.40885    0.28081   1.456        1
# 10 - 1 == 0        -0.22132    0.28081  -0.788        1
# 100 - 1 == 0        0.14469    0.28081   0.515        1
# 100 - 10 == 0       0.36601    0.28081   1.303        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# write LMER
model_t3_Pve <- lmer(scale(alpha) ~ treat + (1|col), data = Pve_RLC_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.988   0.593
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.593).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.966).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(alpha) ~ treat + (1 | col), data = Pve_RLC_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.02854    0.28547  -0.100        1
# 1 - control == 0   -0.23163    0.28547  -0.811        1
# 10 - control == 0   0.08848    0.28547   0.310        1
# 100 - control == 0  0.07577    0.28547   0.265        1
# 1 - 0.1 == 0       -0.20309    0.28547  -0.711        1
# 10 - 0.1 == 0       0.11702    0.28547   0.410        1
# 100 - 0.1 == 0      0.10431    0.28547   0.365        1
# 10 - 1 == 0         0.32010    0.28547   1.121        1
# 100 - 1 == 0        0.30740    0.28547   1.077        1
# 100 - 10 == 0      -0.01271    0.28547  -0.045        1
# (Adjusted p values reported -- holm method)



### --- 5.3.2. Stylophora pistillata -------------------------------------------
#### -- 5.3.2.1 Overall effect -------------------------------------------------
# write LMER
model1_Spi <- lmer(scale(alpha^3) ~ conc + (1|col) + (1|tp),  data = Spi_overall_effect)
# inspect residuals
qqPlot(residuals(model1_Spi))          # good fit
shapiro_test(residuals(model1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Spi)     0.992   0.133
check_normality(model1_Spi)
# OK: residuals appear as normally distributed (p = 0.141).
check_heteroscedasticity(model1_Spi)
#  OK: Error variance appears to be homoscedastic (p = 0.608).

# get summary of LMER
cftest(model1_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(alpha^3) ~ conc + (1 | col) + (1 | tp), 
#           data = Spi_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0 -0.0173650  0.2906441  -0.060    0.952
# conc == 0         0.0007815  0.0012585   0.621    0.535
# (Univariate p values reported)



#### -- 5.3.2.2 Specific effects -----------------------------------------------
# ------------- t1
# write LMER
model_t1_Spi <- lmer(scale(alpha^3) ~ treat + (1|col), data = Spi_RLC_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.987   0.492
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.492).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.955).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(alpha^3) ~ treat + (1 | col), data = Spi_RLC_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0   0.5068     0.2468   2.053 0.240154    
#   1 - control == 0     1.0624     0.2468   4.305 0.000167 ***
#   10 - control == 0    0.9464     0.2468   3.835 0.001129 ** 
#   100 - control == 0   0.8077     0.2468   3.273 0.008517 ** 
#   1 - 0.1 == 0         0.5556     0.2468   2.252 0.170462    
#   10 - 0.1 == 0        0.4397     0.2468   1.782 0.374002    
#   100 - 0.1 == 0       0.3009     0.2468   1.219 0.890790    
#   10 - 1 == 0         -0.1159     0.2468  -0.470 1.000000    
#   100 - 1 == 0        -0.2547     0.2468  -1.032 0.906009    
#   100 - 10 == 0       -0.1388     0.2468  -0.562 1.000000    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# write LMER
model_t2_Spi <- lmer(scale(alpha) ~ treat + (1|col), data = Spi_RLC_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.977   0.106
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.106).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.797).

# get summary of LMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(alpha) ~ treat + (1 | col), data = Spi_RLC_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.204271   0.271427   0.753        1
# 1 - control == 0    0.142555   0.271427   0.525        1
# 10 - control == 0  -0.006849   0.271427  -0.025        1
# 100 - control == 0  0.175762   0.271427   0.648        1
# 1 - 0.1 == 0       -0.061717   0.271427  -0.227        1
# 10 - 0.1 == 0      -0.211121   0.271427  -0.778        1
# 100 - 0.1 == 0     -0.028509   0.271427  -0.105        1
# 10 - 1 == 0        -0.149404   0.271427  -0.550        1
# 100 - 1 == 0        0.033207   0.271427   0.122        1
# 100 - 10 == 0       0.182611   0.271427   0.673        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# write LMER
model_t3_Spi <- lmer(scale(alpha^2) ~ treat + (1|col), data = Spi_RLC_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.977   0.117
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.117).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.971).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(alpha^2) ~ treat + (1 | col), data = Spi_RLC_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.15989    0.29111  -0.549        1
# 1 - control == 0   -0.25293    0.29111  -0.869        1
# 10 - control == 0   0.02522    0.29111   0.087        1
# 100 - control == 0 -0.22900    0.29111  -0.787        1
# 1 - 0.1 == 0       -0.09304    0.29111  -0.320        1
# 10 - 0.1 == 0       0.18511    0.29111   0.636        1
# 100 - 0.1 == 0     -0.06911    0.29111  -0.237        1
# 10 - 1 == 0         0.27815    0.29111   0.955        1
# 100 - 1 == 0        0.02393    0.29111   0.082        1
# 100 - 10 == 0      -0.25422    0.29111  -0.873        1
# (Adjusted p values reported -- holm method)



## ---- 5.4. Beta --------------------------------------------------------------
# not included in the main manuscript as it doesn't give supportive information
### --- 5.4.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.4.1.1 Overall effect -------------------------------------------------
# LMER didn't show a good fit, therefore GLMER is used
model1_Pve <- glmer((beta+100) ~ conc + (1|col) + (1|tp), family = poisson, data = Pve_overall_effect)

# get summary of GLMER
cftest(model1_Pve)
# OUTPUT:  
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (beta + 100) ~ conc + (1 | col) + (1 | tp), data = Pve_overall_effect, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0 4.605e+00  7.071e-01   6.512  7.4e-11 ***
#   conc == 0        5.563e-08  1.558e-04   0.000        1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Univariate p values reported)


#### -- 5.4.1.2 Specific effects -----------------------------------------------
# ------------- t1
# LMER didn't show a good fit, therefore GLMER is used
model_t1_Pve <- glmer((beta+100) ~ treat + (1|col), family = poisson, data = Pve_RLC_t1)

# get summary of GLMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (beta + 100) ~ treat + (1 | col), data = Pve_RLC_t1, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  9.419e-06  3.333e-02   0.000        1
# 1 - control == 0    4.403e-05  3.333e-02   0.001        1
# 10 - control == 0   3.889e-06  3.333e-02   0.000        1
# 100 - control == 0  2.288e-05  3.333e-02   0.001        1
# 1 - 0.1 == 0        3.461e-05  3.333e-02   0.001        1
# 10 - 0.1 == 0      -5.530e-06  3.333e-02   0.000        1
# 100 - 0.1 == 0      1.346e-05  3.333e-02   0.000        1
# 10 - 1 == 0        -4.014e-05  3.333e-02  -0.001        1
# 100 - 1 == 0       -2.115e-05  3.333e-02  -0.001        1
# 100 - 10 == 0       1.899e-05  3.333e-02   0.001        1
# (Adjusted p values reported -- holm method)


# ------------- t2
# LMER didn't show a good fit, therefore GLMER is used
model_t2_Pve <- glmer((beta+100) ~ treat + (1|col), family = poisson, data = Pve_RLC_t2)

# get summary of GLMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (beta + 100) ~ treat + (1 | col), data = Pve_RLC_t2, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  1.936e-05  3.333e-02   0.001        1
# 1 - control == 0   -4.862e-07  3.333e-02   0.000        1
# 10 - control == 0   9.952e-06  3.333e-02   0.000        1
# 100 - control == 0  1.826e-05  3.333e-02   0.001        1
# 1 - 0.1 == 0       -1.984e-05  3.333e-02  -0.001        1
# 10 - 0.1 == 0      -9.405e-06  3.333e-02   0.000        1
# 100 - 0.1 == 0     -1.098e-06  3.333e-02   0.000        1
# 10 - 1 == 0         1.044e-05  3.333e-02   0.000        1
# 100 - 1 == 0        1.875e-05  3.333e-02   0.001        1
# 100 - 10 == 0       8.307e-06  3.333e-02   0.000        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# LMER didn't show a good fit, therefore GLMER is used
model_t3_Pve <- glmer((beta+100) ~ treat + (1|col), family = poisson, data = Pve_RLC_t3)

# get summary of GLMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (beta + 100) ~ treat + (1 | col), data = Pve_RLC_t3, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  2.836e-06  3.333e-02       0        1
# 1 - control == 0   -3.175e-06  3.333e-02       0        1
# 10 - control == 0   1.149e-05  3.333e-02       0        1
# 100 - control == 0 -7.817e-07  3.333e-02       0        1
# 1 - 0.1 == 0       -6.011e-06  3.333e-02       0        1
# 10 - 0.1 == 0       8.651e-06  3.333e-02       0        1
# 100 - 0.1 == 0     -3.618e-06  3.333e-02       0        1
# 10 - 1 == 0         1.466e-05  3.333e-02       0        1
# 100 - 1 == 0        2.393e-06  3.333e-02       0        1
# 100 - 10 == 0      -1.227e-05  3.333e-02       0        1
# (Adjusted p values reported -- holm method)



### --- 5.4.2. Stylophora pistillata -------------------------------------------
#### -- 5.4.2.1 Overall effect -------------------------------------------------
# LMER didn't show a good fit, therefore GLMER is used
model1_Spi <- glmer((beta+100) ~ conc + (1|col) + (1|tp), family = poisson, data = Spi_overall_effect)

# get summary of GMER
cftest(model1_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (beta + 100) ~ conc + (1 | col) + (1 | tp), data = Spi_overall_effect, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) == 0 4.605e+00  7.071e-01   6.512  7.4e-11 ***
#   conc == 0        9.173e-09  1.558e-04   0.000        1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)



#### -- 5.4.2.2 Specific effects -----------------------------------------------
# ------------- t1
# LMER didn't show a good fit, therefore GLMER is used
model_t1_Spi <- glmer((beta+100) ~ treat + (1|col), family = poisson, data = Spi_RLC_t1)

# get summary of GLMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (beta + 100) ~ treat + (1 | col), data = Spi_RLC_t1, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -2.528e-06  3.333e-02       0        1
# 1 - control == 0    6.898e-07  3.333e-02       0        1
# 10 - control == 0  -1.020e-05  3.333e-02       0        1
# 100 - control == 0 -1.670e-06  3.333e-02       0        1
# 1 - 0.1 == 0        3.218e-06  3.333e-02       0        1
# 10 - 0.1 == 0      -7.677e-06  3.333e-02       0        1
# 100 - 0.1 == 0      8.582e-07  3.333e-02       0        1
# 10 - 1 == 0        -1.089e-05  3.333e-02       0        1
# 100 - 1 == 0       -2.360e-06  3.333e-02       0        1
# 100 - 10 == 0       8.535e-06  3.333e-02       0        1
# (Adjusted p values reported -- holm method)


# ------------- t2
# LMER didn't show a good fit, therefore GLMER is used
model_t2_Spi <- glmer((beta+100) ~ treat + (1|col), family = poisson, data = Spi_RLC_t2)

# get summary of GLMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (beta + 100) ~ treat + (1 | col), data = Spi_RLC_t2, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -3.446e-06  3.333e-02       0        1
# 1 - control == 0    7.259e-06  3.333e-02       0        1
# 10 - control == 0   3.704e-06  3.333e-02       0        1
# 100 - control == 0 -4.883e-07  3.333e-02       0        1
# 1 - 0.1 == 0        1.070e-05  3.333e-02       0        1
# 10 - 0.1 == 0       7.150e-06  3.333e-02       0        1
# 100 - 0.1 == 0      2.957e-06  3.333e-02       0        1
# 10 - 1 == 0        -3.555e-06  3.333e-02       0        1
# 100 - 1 == 0       -7.747e-06  3.333e-02       0        1
# 100 - 10 == 0      -4.193e-06  3.333e-02       0        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# LMER didn't show a good fit, therefore GLMER is used
model_t3_Spi <- glmer((beta+100) ~ treat + (1|col), family = poisson, data = Spi_RLC_t3)

# get summary of GLMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (beta + 100) ~ treat + (1 | col), data = Spi_RLC_t3, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  2.039e-06  3.333e-02       0        1
# 1 - control == 0    4.798e-06  3.333e-02       0        1
# 10 - control == 0   8.261e-07  3.333e-02       0        1
# 100 - control == 0  6.326e-06  3.333e-02       0        1
# 1 - 0.1 == 0        2.759e-06  3.333e-02       0        1
# 10 - 0.1 == 0      -1.213e-06  3.333e-02       0        1
# 100 - 0.1 == 0      4.286e-06  3.333e-02       0        1
# 10 - 1 == 0        -3.972e-06  3.333e-02       0        1
# 100 - 1 == 0        1.528e-06  3.333e-02       0        1
# 100 - 10 == 0       5.500e-06  3.333e-02       0        1
# (Adjusted p values reported -- holm method)



# ---- 6. Write tables ---------------------------------------------------------
## --- 6.1. Means and standard deviation ---------------------------------------
mean_RLC <- RLC_parameter %>%
  # separates species, treatment and timepoints
  group_by(spec, treat, tp) %>% 
  # columns of interest: parameter
  get_summary_stats(alpha, beta, rETRmax, Ek) %>% 
  # remove unnecessary columns
  select(-q1, -q3, -iqr, -mad, -se, -ci, -median) 

# write it into .csv
write_csv2(mean_RLC, "out/RLC_parameter_means.csv")
