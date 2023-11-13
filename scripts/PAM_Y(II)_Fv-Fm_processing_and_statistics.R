# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the data processing and statistical analyzes of the effective (Y(II)) and maximum (Fv/Fm) quantum yield
# of the corals photosymbionts. Both parameters were measured using pulse amplitude modulated fluorometry (PAM)
#   a) Y(II) was measured in light adapted corals
#   b) Fv/Fm was measured in dark adapted corals
# Statistical analyzes will be conducted using LMER and GLMER 
#   together with a holm adjusted glht summary

# Analyzes of further parameters of the photosynthetic efficiency (relative electron transport rate (rETRmax), 
# efficiency of light capture (α), and minimum saturating irradiance (Ek)) of the corals photosymbionts 
# can be found in the script 'RLC_processing_and_statistics'

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



#library(readxl)
#library(qqplotr) # für qqplot mit ggplot, pipe bar
#library(broomExtra) # write models like lm as tibble
#library(magrittr) # A Forward-Pipe Operator
#library(dplyr)
#library(ggpubr)



# ---- 3. Read in needed data files --------------------------------------------
## --- 3.1. Coral identity table -----------------------------------------------
# read in list with overview of all corals used and their treatments etc.
corals <- read_csv2("in/coral_treatments.csv") %>%
  mutate(treat = as.factor(treat), # column for categorical model
         conc = as.numeric(conc)) # column for continuous model 

## --- 3.2. Light adapted PAM data ---------------------------------------------
# read in tables for preparing and analysing Y(II)
# tables are devided by timepoint

# ------------- t0
t0_light <- read_csv2("in/PAM Data/t0_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t0 = "Y(II)")  %>% 
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F, -nr) 
# remove doubled rows
t0_light<- t0_light[-c(49, 50, 51),]

# ------------- t1
t1_light <- read_csv2("in/PAM Data/t1_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t1 = "Y(II)") %>%
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F, -nr)
# Rename wrongly named entries
t1_light[472, 1] <- "Spi_B_8"  
t1_light[473, 1] <- "Spi_B_8"  
t1_light[474, 1] <- "Spi_B_8" 

# ------------- t2
t2_light <- read_csv2("in/PAM Data/t2_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t2 = "Y(II)") %>%
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F, -nr) 
# Rename wrongly named entries
t2_light[136, 1] <- "Pve_D_12" 
t2_light[137, 1] <- "Pve_D_12" 
t2_light[138, 1] <- "Pve_D_12" 

# ------------- t3
t3_light <- read_csv2("in/PAM Data/t3_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t3 = "Y(II)") %>%
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)


## --- 3.3. Dark adapted PAM data ----------------------------------------------
# read in tables for preparing and analysing Fv/Fm
# tables are devided by timepoint

# ------------- t0
t0_dark <- read_csv2("in/PAM Data/t0_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t0 = "Y(II)") %>%
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)

# ------------- t1
t1_dark <- read_csv2("in/PAM Data/t1_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t1 = "Y(II)") %>%
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)

# ------------- t2
t2_dark <- read_csv2("in/PAM Data/t2_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t2 = "Y(II)") %>% 
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F) 

# ------------- t3
t3_dark <- read_csv2("in/PAM Data/t3_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t3 = "Y(II)") %>% 
  # remove unnecessary columns
  select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)



# ----- 4. Prepare data for statistical analyzes -------------------------------
## ---- 4.1. Reformat data table -----------------------------------------------
### --- 4.1.1. Light adapted PAM data ------------------------------------------
# put all tables of Y(II) in long format
# ------------- t0
light_0 <-  merge(corals, t0_light, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         YII_t0 = as.numeric(YII_t0))%>% 
  rename("0" = "YII_t0") %>%
  pivot_longer(cols=c('0'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='YII') # assign new name for the values previously in the column under the headers above

# ------------- t1
light_1 <-  merge(corals, t1_light, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         YII_t1 = as.numeric(YII_t1))%>% 
  rename("1" = "YII_t1") %>%
  pivot_longer(cols=c('1'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='YII') # assign new name for the values previously in the column under the headers above

# ------------- t2
light_2 <-  merge(corals, t2_light, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         YII_t2 = as.numeric(YII_t2))%>% 
  rename("2" = "YII_t2") %>%
  pivot_longer(cols=c('2'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='YII') # assign new name for the values previously in the column under the headers above

# ------------- t3
light_3 <-  merge(corals, t3_light, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         YII_t3 = as.numeric(YII_t3))%>% 
  rename("3" = "YII_t3") %>%
  pivot_longer(cols=c('3'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='YII') # assign new name for the values previously in the column under the headers above

# merge all tables into one for all timepoints
light_all <- rbind(light_0, light_1, light_2, light_3) %>%
  # change column character where necessary
  mutate(tp = as.numeric(tp))

# relevel treatments
light_all$treat <- factor(light_all$treat, 
                         levels = c("control", "0.1", "1",
                                    "10", "100"))



### --- 4.1.2. Dark adapted PAM data ------------------------------------------
# put all tables of Fv/Fm in long format
# ------------- t0
dark_0 <-  merge(corals, t0_dark, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         Fv_Fm_t0 = as.numeric(Fv_Fm_t0))%>% 
  rename("0" = "Fv_Fm_t0") %>%
  pivot_longer(cols=c('0'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='Fv_Fm') # assign new name for the values previously in the column under the headers above

# ------------- t1
dark_1 <-  merge(corals, t1_dark, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         Fv_Fm_t1 = as.numeric(Fv_Fm_t1))%>% 
  rename("1" = "Fv_Fm_t1") %>%
  pivot_longer(cols=c('1'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='Fv_Fm') # assign new name for the values previously in the column under the headers above

# ------------- t2
dark_2 <-  merge(corals, t2_dark, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         Fv_Fm_t2 = as.numeric(Fv_Fm_t2))%>% 
  rename("2" = "Fv_Fm_t2") %>%
  pivot_longer(cols=c('2'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='Fv_Fm') # assign new name for the values previously in the column under the headers above

# ------------- t3
dark_3 <-  merge(corals, t3_dark, by = 'ID', all.x = TRUE) %>%
  # change column character where necessary
  mutate(treat = as.factor(treat),
         Fv_Fm_t3 = as.numeric(Fv_Fm_t3))%>% 
  rename("3" = "Fv_Fm_t3") %>%
  pivot_longer(cols=c('3'), # use previous colums as new entries - categories
               names_to='tp',  # assign new name for the column of tp
               values_to='Fv_Fm') # assign new name for the values previously in the column under the headers above

# merge all tables into one for all timepoints
dark_all <- rbind(dark_0, dark_1, dark_2, dark_3) %>%
  # change column character where necessary
  mutate(tp = as.numeric(tp))

# relevel treatments
dark_all$treat <- factor(dark_all$treat, 
                         levels = c("control", "0.1", "1",
                                    "10", "100"))


# ----- 5. Statistics for PAM light measurements--------------------------------
## ---- 5.1. Y(II) -------------------------------------------------------------
# create table with means of PAM measurements for model creation
Light_mean <- light_all %>%
  group_by(spec, tp, ID, conc, treat, col) %>%
  get_summary_stats(YII, type = "mean")



### --- 5.1.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.1.1.1 Overall effect -------------------------------------------------
# create a subset with data of Pve, t0 excluded for continuous model
Pve_overall_effect <- subset(Light_mean, spec == "Pve" & tp!= "0")

# write LMER
model1_Pve <- lmer(scale(mean^4) ~ conc + (1|col) + (1|tp), data = Pve_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Pve))          # okay fit
shapiro_test(residuals(model1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Pve)     0.991  0.0883
check_normality(model1_Pve)
# OK: residuals appear as normally distributed (p = 0.089).
check_heteroscedasticity(model1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.486).

# get summary of LMER
cftest(model1_Pve)
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(mean^4) ~ conc + (1 | col) + (1 | tp), data = Pve_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept) == 0  0.081855   0.297566   0.275  0.78325   
# conc == 0        -0.003684   0.001162  -3.171  0.00152 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


#### -- 5.1.1.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Pve, to test differences at t1
Pve_light_t1 <- subset(Light_mean, spec == "Pve" & tp == "1")

# write LMER 
model_t1_Pve <- lmer(scale(mean) ~ treat + (1|col), data = Pve_light_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.989   0.643
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.643).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.972).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Pve_light_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0  0.45143    0.24399   1.850  0.32143   
#   1 - control == 0    0.51518    0.24399   2.111  0.20839   
#   10 - control == 0  -0.10166    0.24399  -0.417  1.00000   
#   100 - control == 0 -0.40490    0.24399  -1.660  0.38804   
#   1 - 0.1 == 0        0.06375    0.24399   0.261  1.00000   
#   10 - 0.1 == 0      -0.55308    0.24399  -2.267  0.16380   
#   100 - 0.1 == 0     -0.85633    0.24399  -3.510  0.00404 **
#   10 - 1 == 0        -0.61683    0.24399  -2.528  0.09174 . 
#   100 - 1 == 0       -0.92008    0.24399  -3.771  0.00163 **
#   100 - 10 == 0      -0.30325    0.24399  -1.243  0.64174   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Pve, to test differences at t2
Pve_light_t2 <- subset(Light_mean, spec == "Pve" & tp == "2")

# write LMER
model_t2_Pve <- lmer(scale(mean) ~ treat + (1|col), data = Pve_light_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.983   0.300
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.300).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.944).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Pve_light_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0  0.29656    0.24614   1.205   1.0000  
# 1 - control == 0   -0.09554    0.24614  -0.388   1.0000  
# 10 - control == 0  -0.31846    0.24614  -1.294   1.0000  
# 100 - control == 0 -0.35030    0.24614  -1.423   1.0000  
# 1 - 0.1 == 0       -0.39210    0.24614  -1.593   0.8893  
# 10 - 0.1 == 0      -0.61502    0.24614  -2.499   0.1122  
# 100 - 0.1 == 0     -0.64687    0.24614  -2.628   0.0859 .
# 10 - 1 == 0        -0.22292    0.24614  -0.906   1.0000  
# 100 - 1 == 0       -0.25477    0.24614  -1.035   1.0000  
# 100 - 10 == 0      -0.03185    0.24614  -0.129   1.0000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Pve, to test differences at t3
Pve_light_t3 <- subset(Light_mean, spec == "Pve" & tp == "3")

# write LMER
model_t3_Pve <- lmer(scale(mean^12) ~ treat + (1|col), data = Pve_light_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # okay fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT:  A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.973  0.0601
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.060).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.628).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean^12) ~ treat + (1 | col), data = Pve_light_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.065529   0.221652   0.296        1
# 1 - control == 0    0.357732   0.221652   1.614        1
# 10 - control == 0   0.135670   0.221652   0.612        1
# 100 - control == 0  0.125707   0.221652   0.567        1
# 1 - 0.1 == 0        0.292204   0.221652   1.318        1
# 10 - 0.1 == 0       0.070141   0.221652   0.316        1
# 100 - 0.1 == 0      0.060179   0.221652   0.272        1
# 10 - 1 == 0        -0.222063   0.221652  -1.002        1
# 100 - 1 == 0       -0.232025   0.221652  -1.047        1
# 100 - 10 == 0      -0.009962   0.221652  -0.045        1
# (Adjusted p values reported -- holm method)



### --- 5.1.2. Stylophora pistillata -------------------------------------------
#### -- 5.1.2.1 Overall effect -------------------------------------------------
# create a subset with data of Spi, t0 excluded for continuous model
Spi_overall_effect <- subset(Light_mean, spec == "Spi" & tp!= "0")

# LMER didn't show a good fit, therefore GLMER is used
model1_Spi <- glmer((mean) ~ conc + (1|col) + (1|tp), family = poisson, data = Spi_overall_effect)

# get summary of GLMER
cftest(model1_Spi)
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (mean) ~ conc + (1 | col) + (1 | tp), data = Spi_overall_effect, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0 -0.5181142  0.7129059  -0.727    0.467
# conc == 0        -0.0001463  0.0020294  -0.072    0.943
# (Univariate p values reported)



#### -- 5.1.2.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Spi, to test differences at t1
Spi_light_t1 <- subset(Light_mean, spec == "Spi" & tp == "1")

# write LMER
model_t1_Spi <- lmer(scale(mean^6) ~ treat + (1|col), data = Spi_light_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # okay fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.974  0.0671
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.067).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.278).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean^6) ~ treat + (1 | col), data = Spi_light_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
#   0.1 - control == 0 -0.31084    0.27268  -1.140   1.0000  
#   1 - control == 0   -0.33418    0.27268  -1.226   1.0000  
#   10 - control == 0  -0.61172    0.27268  -2.243   0.2239  
#   100 - control == 0 -0.87271    0.27268  -3.201   0.0137 *
#   1 - 0.1 == 0       -0.02333    0.27268  -0.086   1.0000  
#   10 - 0.1 == 0      -0.30087    0.27268  -1.103   1.0000  
#   100 - 0.1 == 0     -0.56187    0.27268  -2.061   0.3148  
#   10 - 1 == 0        -0.27754    0.27268  -1.018   1.0000  
#   100 - 1 == 0       -0.53854    0.27268  -1.975   0.3379  
#   100 - 10 == 0      -0.26100    0.27268  -0.957   1.0000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Spi, to test differences at t2
Spi_light_t2 <- subset(Light_mean, spec == "Spi" & tp == "2")

# write LMER
model_t2_Spi <- lmer(scale(mean) ~ treat + (1|col), data = Spi_light_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.979   0.153
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.153).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.888).

# get summary of LMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Spi_light_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.35316    0.23040  -1.533    0.877
# 1 - control == 0   -0.57510    0.23040  -2.496    0.126
# 10 - control == 0  -0.44966    0.23040  -1.952    0.459
# 100 - control == 0 -0.42650    0.23040  -1.851    0.513
# 1 - 0.1 == 0       -0.22193    0.23040  -0.963    1.000
# 10 - 0.1 == 0      -0.09649    0.23040  -0.419    1.000
# 100 - 0.1 == 0     -0.07333    0.23040  -0.318    1.000
# 10 - 1 == 0         0.12544    0.23040   0.544    1.000
# 100 - 1 == 0        0.14860    0.23040   0.645    1.000
# 100 - 10 == 0       0.02316    0.23040   0.101    1.000
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Spi, to test differences at t3
Spi_light_t3 <- subset(Light_mean, spec == "Spi" & tp == "3")

# write LMER
model_t3_Spi <- lmer(scale(mean^3) ~ treat + (1|col), data = Spi_light_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # okay fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.974  0.0685
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.069).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.946).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean^3) ~ treat + (1 | col), data = Spi_light_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.15640    0.20201  -0.774    1.000
# 1 - control == 0   -0.29544    0.20201  -1.463    1.000
# 10 - control == 0  -0.49370    0.20201  -2.444    0.145
# 100 - control == 0 -0.35866    0.20201  -1.776    0.682
# 1 - 0.1 == 0       -0.13904    0.20201  -0.688    1.000
# 10 - 0.1 == 0      -0.33730    0.20201  -1.670    0.760
# 100 - 0.1 == 0     -0.20226    0.20201  -1.001    1.000
# 10 - 1 == 0        -0.19826    0.20201  -0.981    1.000
# 100 - 1 == 0       -0.06322    0.20201  -0.313    1.000
# 100 - 10 == 0       0.13504    0.20201   0.668    1.000
# (Adjusted p values reported -- holm method)



## ---- 5.2. Fv/Fm -------------------------------------------------------------
### --- 5.2.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.2.1.1 Overall effect -------------------------------------------------
# create a subset with data of Pve, t0 excluded for continuous model
Pve_overall_effect <- subset(dark_all, spec == "Pve" & tp!= "0")

# write LMER
model1_Pve <- lmer(scale(Fv_Fm) ~ conc + (1|col) + (1|tp), data = Pve_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Pve))          # good fit
shapiro_test(residuals(model1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Pve)     0.990  0.0748
check_normality(model1_Pve)
# OK: residuals appear as normally distributed (p = 0.077).
check_heteroscedasticity(model1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.799).

# get summary of LMER
cftest(model1_Pve)
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(Fv_Fm) ~ conc + (1 | col) + (1 | tp), data = Pve_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0 -0.041710   0.243778  -0.171    0.864
# conc == 0         0.001877   0.001338   1.403    0.161
# (Univariate p values reported)



#### -- 5.2.1.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Pve, to test differences at t1
Pve_dark_t1 <- subset(dark_all, spec == "Pve" & tp == "1")

# write LMER
model_t1_Pve <- lmer(log(Fv_Fm) ~ treat + (1|col), data = Pve_dark_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.990   0.749
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.749).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.756).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = log(Fv_Fm) ~ treat + (1 | col), data = Pve_dark_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.055616   0.094472  -0.589        1
# 1 - control == 0   -0.051225   0.094472  -0.542        1
# 10 - control == 0  -0.079695   0.094472  -0.844        1
# 100 - control == 0 -0.123419   0.094472  -1.306        1
# 1 - 0.1 == 0        0.004391   0.094472   0.046        1
# 10 - 0.1 == 0      -0.024079   0.094472  -0.255        1
# 100 - 0.1 == 0     -0.067803   0.094472  -0.718        1
# 10 - 1 == 0        -0.028470   0.094472  -0.301        1
# 100 - 1 == 0       -0.072194   0.094472  -0.764        1
# 100 - 10 == 0      -0.043724   0.094472  -0.463        1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Pve, to test differences at t2
Pve_dark_t2 <- subset(dark_all, spec == "Pve" & tp == "2")

# write LMER
model_t2_Pve <- lmer(scale(Fv_Fm) ~ treat + (1|col), data = Pve_dark_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT:  A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.979   0.142
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.142).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.818).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Fv_Fm) ~ treat + (1 | col), data = Pve_dark_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.26879    0.30089   0.893        1
# 1 - control == 0    0.22940    0.30089   0.762        1
# 10 - control == 0   0.38475    0.30089   1.279        1
# 100 - control == 0  0.35230    0.30089   1.171        1
# 1 - 0.1 == 0       -0.03939    0.30089  -0.131        1
# 10 - 0.1 == 0       0.11596    0.30089   0.385        1
# 100 - 0.1 == 0      0.08350    0.30089   0.278        1
# 10 - 1 == 0         0.15535    0.30089   0.516        1
# 100 - 1 == 0        0.12289    0.30089   0.408        1
# 100 - 10 == 0      -0.03246    0.30089  -0.108        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Pve, to test differences at t3
Pve_dark_t3 <- subset(dark_all, spec == "Pve" & tp == "3")

# write LMER
model_t3_Pve <- lmer(scale(Fv_Fm) ~ treat + (1|col), data = Pve_dark_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.985   0.394
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.394).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.749).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Fv_Fm) ~ treat + (1 | col), data = Pve_dark_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0   0.3585     0.2866   1.251   1.0000  
# 1 - control == 0     0.1613     0.2866   0.563   1.0000  
# 10 - control == 0    0.4647     0.2866   1.621   0.7346  
# 100 - control == 0   0.8628     0.2866   3.010   0.0261 *
#   1 - 0.1 == 0        -0.1973     0.2866  -0.688   1.0000  
# 10 - 0.1 == 0        0.1062     0.2866   0.370   1.0000  
# 100 - 0.1 == 0       0.5043     0.2866   1.759   0.6282  
# 10 - 1 == 0          0.3034     0.2866   1.059   1.0000  
# 100 - 1 == 0         0.7015     0.2866   2.448   0.1294  
# 100 - 10 == 0        0.3981     0.2866   1.389   0.9891  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



### --- 5.2.2. Stylophora pistillata -------------------------------------------
#### -- 5.2.2.1 Overall effect -------------------------------------------------
# create a subset with data of Spi, t0 excluded for continuous model
Spi_overall_effect <- subset(dark_all, spec == "Spi" & tp!= "0")

# write LMER
model1_Spi <- lmer(scale(Fv_Fm) ~ conc + (1|col) + (1|tp), data = Spi_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Spi))          # good fit
shapiro_test(residuals(model1_Spi))    # p > 0.05 Normality of residuals
# OUTPUT:  A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Spi)     0.993   0.244
check_normality(model1_Spi)
# OK: residuals appear as normally distributed (p = 0.249).
check_heteroscedasticity(model1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.889).

# get summary of LMER
cftest(model1_Spi)
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(Fv_Fm) ~ conc + (1 | col) + (1 | tp), data = Spi_overall_effect)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0  0.024652   0.100095   0.246    0.805
# conc == 0        -0.001109   0.001538  -0.721    0.471
# (Univariate p values reported)


#### -- 5.2.2.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Spi, to test differences at t1
Spi_dark_t1 <- subset(dark_all, spec == "Spi" & tp == "1")

# write LMER
model_t1_Spi <- lmer(scale(Fv_Fm) ~ treat + (1|col), data = Spi_dark_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.990   0.728
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.728).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.641).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Fv_Fm) ~ treat + (1 | col), data = Spi_dark_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.29129    0.32933   0.884        1
# 1 - control == 0   -0.01998    0.32933  -0.061        1
# 10 - control == 0   0.18077    0.32933   0.549        1
# 100 - control == 0  0.01155    0.32933   0.035        1
# 1 - 0.1 == 0       -0.31128    0.32933  -0.945        1
# 10 - 0.1 == 0      -0.11052    0.32933  -0.336        1
# 100 - 0.1 == 0     -0.27974    0.32933  -0.849        1
# 10 - 1 == 0         0.20075    0.32933   0.610        1
# 100 - 1 == 0        0.03153    0.32933   0.096        1
# 100 - 10 == 0      -0.16922    0.32933  -0.514        1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Spi, to test differences at t2
Spi_dark_t2 <- subset(dark_all, spec == "Spi" & tp == "2")

# write LMER
model_t2_Spi <- lmer(scale(Fv_Fm) ~ treat + (1|col), data = Spi_dark_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.988   0.552
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.552).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.747).

# get summary of LMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Fv_Fm) ~ treat + (1 | col), data = Spi_dark_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.316147   0.329381  -0.960        1
# 1 - control == 0    0.007781   0.329381   0.024        1
# 10 - control == 0   0.055045   0.329381   0.167        1
# 100 - control == 0 -0.006052   0.329381  -0.018        1
# 1 - 0.1 == 0        0.323928   0.329381   0.983        1
# 10 - 0.1 == 0       0.371192   0.329381   1.127        1
# 100 - 0.1 == 0      0.310095   0.329381   0.941        1
# 10 - 1 == 0         0.047264   0.329381   0.143        1
# 100 - 1 == 0       -0.013833   0.329381  -0.042        1
# 100 - 10 == 0      -0.061097   0.329381  -0.185        1
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Spi, to test differences at t3
Spi_dark_t3 <- subset(dark_all, spec == "Spi" & tp == "3")

# write LMER
model_t3_Spi <- lmer(scale(Fv_Fm) ~ treat + (1|col), data = Spi_dark_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.976  0.0975
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.097).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.951).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(Fv_Fm) ~ treat + (1 | col), data = Spi_dark_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.38431    0.33310   1.154        1
# 1 - control == 0    0.01615    0.33310   0.048        1
# 10 - control == 0   0.42445    0.33310   1.274        1
# 100 - control == 0 -0.09121    0.33310  -0.274        1
# 1 - 0.1 == 0       -0.36816    0.33310  -1.105        1
# 10 - 0.1 == 0       0.04014    0.33310   0.121        1
# 100 - 0.1 == 0     -0.47552    0.33310  -1.428        1
# 10 - 1 == 0         0.40830    0.33310   1.226        1
# 100 - 1 == 0       -0.10736    0.33310  -0.322        1
# 100 - 10 == 0      -0.51566    0.33310  -1.548        1
# (Adjusted p values reported -- holm method)



# ----- 6. Write tables --------------------------------------------------------
## ---- 6.1. Table of all light adapted PAM data -------------------------------
write_rds(light_all, "processed/light_all.rds")

## ---- 6.2. Table of all dark adapted PAM data --------------------------------
write_rds(dark_all, "processed/dark_all")
