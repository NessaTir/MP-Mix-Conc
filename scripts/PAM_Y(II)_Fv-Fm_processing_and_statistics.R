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
library(dplyr)

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
corals <- read_csv("in/coral_treatments.csv") %>%
  dplyr::mutate(treat = as.factor(treat), # column for categorical model
         conc = as.numeric(conc)) # column for continuous model 

## --- 3.2. Light adapted PAM data ---------------------------------------------
# read in tables for preparing and analysing Y(II)
# tables are devided by timepoint

# ------------- t0
t0_light <- read_csv2("in/PAM Data/t0_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t0 = "Y(II)")  %>% 
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F, -nr) 
# remove doubled rows
t0_light<- t0_light[-c(49, 50, 51),]

# ------------- t1
t1_light <- read_csv2("in/PAM Data/t1_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t1 = "Y(II)") %>%
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F, -nr)
# Rename wrongly named entries
t1_light[472, 1] <- "Spi_B_8"  
t1_light[473, 1] <- "Spi_B_8"  
t1_light[474, 1] <- "Spi_B_8" 

# ------------- t2
t2_light <- read_csv2("in/PAM Data/t2_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t2 = "Y(II)") %>%
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F, -nr) 
# Rename wrongly named entries
t2_light[136, 1] <- "Pve_D_12" 
t2_light[137, 1] <- "Pve_D_12" 
t2_light[138, 1] <- "Pve_D_12" 

# ------------- t3
t3_light <- read_csv2("in/PAM Data/t3_Light_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", YII_t3 = "Y(II)") %>%
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)


## --- 3.3. Dark adapted PAM data ----------------------------------------------
# read in tables for preparing and analysing Fv/Fm
# tables are devided by timepoint

# ------------- t0
t0_dark <- read_csv2("in/PAM Data/t0_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t0 = "Y(II)") %>%
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)

# ------------- t1
t1_dark <- read_csv2("in/PAM Data/t1_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t1 = "Y(II)") %>%
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)

# ------------- t2
t2_dark <- read_csv2("in/PAM Data/t2_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t2 = "Y(II)") %>% 
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F) 

# ------------- t3
t3_dark <- read_csv2("in/PAM Data/t3_dark_all.csv") %>%
  # clean column names for better merge with other timepoints
  rename(ID = "No.", Fv_Fm_t3 = "Y(II)") %>% 
  # remove unnecessary columns
  dplyr::select(-ML, -Temp., -PAR, -ETR, -Fm, -"# t", -Date, -Time, -F)



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

# calculate relative YII values
YII_0 <- subset(Light_mean, tp == "0") %>%
  rename(mean_t0 = mean) %>%
  # remove doubled coulms for clean merge
  dplyr::select(-tp)

YII_1_3 <- subset(Light_mean, tp != "0") %>%
  # remove doubled coulms for clean merge
  dplyr::select(-spec, -col, -treat, -conc, -variable, -n)

Light_relative <- full_join(YII_0, YII_1_3, by = "ID")

Light_relative <- Light_relative %>% 
  mutate(relativeYII = 100/mean_t0*mean) 



### --- 5.1.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.1.1.1 Overall effect -------------------------------------------------
# create a subset with data of Pve, t0 excluded for continuous model
Pve_overall_effect <- subset(Light_relative, spec == "Pve")

# LMER didn't show a good fit - write GLMER
model1_Pve <- glmer((relativeYII) ~ conc + (1|col) + (1|tp), family = poisson, data = Pve_overall_effect)

# get summary of GLMER
cftest(model1_Pve)
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (relativeYII) ~ conc + (1 | col) + (1 | tp), family = poisson, data = Pve_overall_effect)
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0 4.6301317  0.7071406   6.548 5.84e-11 ***
#  conc == 0        0.0001469  0.0001529   0.961    0.337    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


#### -- 5.1.1.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Pve, to test differences at t1
Pve_light_t1 <- subset(Light_mean, spec == "Pve" & tp == "1")

# LMER didn't show a good fit - write GLMER
model_t1_Pve <- glmer((relativeYII) ~ treat + (1|col), family = poisson, data = Pve_light_t1)

# get summary of GLMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = relativeYII ~ treat + (1 | col), family = poisson, data = Pve_light_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0 -0.007604   0.033031  -0.230        1
#   1 - control == 0    0.032411   0.032704   0.991        1
#   10 - control == 0   0.020609   0.032799   0.628        1
#   100 - control == 0  0.007310   0.032908   0.222        1
#   1 - 0.1 == 0        0.040015   0.032767   1.221        1
#   10 - 0.1 == 0       0.028212   0.032862   0.858        1
#   100 - 0.1 == 0      0.014914   0.032971   0.452        1
#   10 - 1 == 0        -0.011803   0.032534  -0.363        1
#   100 - 1 == 0       -0.025101   0.032643  -0.769        1
#   100 - 10 == 0      -0.013299   0.032739  -0.406        1 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Pve, to test differences at t2
Pve_light_t2 <- subset(Light_relative, spec == "Pve" & tp == "2")

# write LMER
model_t2_Pve <- lmer(scale(log(relativeYII)) ~ treat + (1|col), data = Pve_light_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.977   0.107
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.107).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.941).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Pve_light_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0 -0.30011    0.32548  -0.922        1
# 1 - control == 0    0.02977    0.32548   0.091        1
# 10 - control == 0   0.17322    0.32548   0.532        1
# 100 - control == 0  0.20008    0.32548   0.615        1
# 1 - 0.1 == 0        0.32989    0.32548   1.014        1
# 10 - 0.1 == 0       0.47334    0.32548   1.454        1
# 100 - 0.1 == 0      0.50020    0.32548   1.537        1
# 10 - 1 == 0         0.14345    0.32548   0.441        1
# 100 - 1 == 0        0.17031    0.32548   0.523        1
# 100 - 10 == 0       0.02686    0.32548   0.083        1  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Pve, to test differences at t3
Pve_light_t3 <- subset(Light_relative, spec == "Pve" & tp == "3")

# write LMER
model_t3_Pve <- lmer(scale(log(relativeYII)) ~ treat + (1|col), data = Pve_light_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # okay fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT:  A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.977  0.105
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.105).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.807).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean^12) ~ treat + (1 | col), data = Pve_light_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.51753    0.28761  -1.799 0.359762    
# 1 - control == 0    0.38411    0.28761   1.336 0.726824    
# 10 - control == 0   0.60006    0.28761   2.086 0.240423    
# 100 - control == 0  0.60858    0.28761   2.116 0.240423    
# 1 - 0.1 == 0        0.90164    0.28761   3.135 0.013751 *  
# 10 - 0.1 == 0       1.11759    0.28761   3.886 0.000918 ***
# 100 - 0.1 == 0      1.12611    0.28761   3.915 0.000903 ***
# 10 - 1 == 0         0.21595    0.28761   0.751 1.000000    
# 100 - 1 == 0        0.22447    0.28761   0.780 1.000000    
# 100 - 10 == 0       0.00852    0.28761   0.030 1.000000   
# (Adjusted p values reported -- holm method)



### --- 5.1.2. Stylophora pistillata -------------------------------------------
#### -- 5.1.2.1 Overall effect -------------------------------------------------
# create a subset with data of Spi, t0 excluded for continuous model
Spi_overall_effect <- subset(Light_relative, spec == "Spi" & tp!= "0")

# LMER didn't show a good fit - GLMER is used
model1_Spi <- glmer(relativeYII ~ conc + (1|col) + (1|tp), family = poisson, data = Spi_overall_effect)

# get summary of GLMER
cftest(model1_Spi)
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (relativeYII) ~ conc + (1 | col) + (1 | tp), data = Spi_overall_effect, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) == 0 4.6417451  0.7071402   6.564 5.23e-11 ***
# conc == 0        0.0001344  0.0001521   0.883    0.377    
# (Univariate p values reported)



#### -- 5.1.2.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Spi, to test differences at t1
Spi_light_t1 <- subset(Light_relative, spec == "Spi" & tp == "1")

# write LMER
model_t1_Spi <- lmer(scale(relativeYII) ~ treat + (1|col), data = Spi_light_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # okay fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.981  0.219
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.219).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.901).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(relativeYII) ~ treat + (1 | col), data = Spi_light_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0 -0.31601    0.31310  -1.009  1.00000   
# 1 - control == 0   -0.18789    0.31310  -0.600  1.00000   
# 10 - control == 0   0.72875    0.31310   2.328  0.13956   
# 100 - control == 0 -0.05311    0.31310  -0.170  1.00000   
# 1 - 0.1 == 0        0.12812    0.31310   0.409  1.00000   
# 10 - 0.1 == 0       1.04476    0.31310   3.337  0.00847 **
# 100 - 0.1 == 0      0.26290    0.31310   0.840  1.00000   
# 10 - 1 == 0         0.91664    0.31310   2.928  0.03074 * 
# 100 - 1 == 0        0.13478    0.31310   0.430  1.00000   
# 100 - 10 == 0      -0.78186    0.31310  -2.497  0.10016  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Spi, to test differences at t2
Spi_light_t2 <- subset(Light_relative, spec == "Spi" & tp == "2")

# write LMER
model_t2_Spi <- lmer(scale(log(relativeYII)) ~ treat + (1|col), data = Spi_light_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.974   0.0635
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.063).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.666).

# get summary of LMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(relativeYII) ~ treat + (1 | col), data = Spi_light_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  -0.2426     0.2935  -0.827 0.867283    
# 1 - control == 0    -0.3759     0.2935  -1.281 0.801279    
# 10 - control == 0    0.8272     0.2935   2.818 0.038610 *  
# 100 - control == 0   0.3111     0.2935   1.060 0.867283    
# 1 - 0.1 == 0        -0.1333     0.2935  -0.454 0.867283    
# 10 - 0.1 == 0        1.0698     0.2935   3.645 0.002407 ** 
# 100 - 0.1 == 0       0.5537     0.2935   1.887 0.355223    
# 10 - 1 == 0          1.2031     0.2935   4.099 0.000415 ***
# 100 - 1 == 0         0.6870     0.2935   2.341 0.134715    
# 100 - 10 == 0       -0.5161     0.2935  -1.758 0.393482    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Spi, to test differences at t3
Spi_light_t3 <- subset(Light_relative, spec == "Spi" & tp == "3")

# write LMER
model_t3_Spi <- lmer(scale(relativeYII) ~ treat + (1|col), data = Spi_light_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # okay fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.975  0.0747
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.075).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.479).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(relativeYII) ~ treat + (1 | col), data = Spi_light_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.09878    0.28816  -0.343   1.0000  
# 1 - control == 0   -0.15296    0.28816  -0.531   1.0000  
# 10 - control == 0   0.74990    0.28816   2.602   0.0741 .
# 100 - control == 0  0.30876    0.28816   1.071   1.0000  
# 1 - 0.1 == 0       -0.05419    0.28816  -0.188   1.0000  
# 10 - 0.1 == 0       0.84867    0.28816   2.945   0.0291 *
# 100 - 0.1 == 0      0.40753    0.28816   1.414   0.7865  
# 10 - 1 == 0         0.90286    0.28816   3.133   0.0173 *
# 100 - 1 == 0        0.46172    0.28816   1.602   0.7637  
# 100 - 10 == 0      -0.44114    0.28816  -1.531   0.7637  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



## ---- 5.2. Fv/Fm -------------------------------------------------------------
# create table with means of PAM measurements for model creation
Dark_mean <- dark_all %>%
  group_by(spec, tp, ID, conc, treat, col) %>%
  get_summary_stats(Fv_Fm, type = "mean")

# calculate relative YII values
Fv_Fm_0 <- subset(Dark_mean, tp == "0") %>%
  rename(mean_t0 = mean) %>%
  # remove doubled coulms for clean merge
  dplyr::select(-tp)

Fv_Fm_1_3 <- subset(Dark_mean, tp != "0") %>%
  # remove doubled coulms for clean merge
  dplyr::select(-spec, -col, -treat, -conc, -variable, -n)

Dark_relative <- full_join(Fv_Fm_0, Fv_Fm_1_3, by = "ID")

Dark_relative <- Dark_relative %>% 
  mutate(relativeFv_Fm = 100/mean_t0*mean) 


### --- 5.2.1. Pocillopora verrucosa -------------------------------------------
#### -- 5.2.1.1 Overall effect -------------------------------------------------
Pve_overall_effect <- subset(Dark_relative, spec == "Pve" & tp!= "0")

# write LMER
model1_Pve <- lmer(log(relativeFv_Fm) ~ conc + (1|col) + (1|tp), data = Pve_overall_effect)

# inspect residuals
qqPlot(residuals(model1_Pve))          # good fit
shapiro_test(residuals(model1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Pve)     0.993  0.196
check_normality(model1_Pve)
# OK: residuals appear as normally distributed (p = 0.202).
check_heteroscedasticity(model1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.668).

# get summary of LMER
cftest(model1_Pve)
# OUTPUT: 
#  Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = log(relativeFv_Fm) ~ conc + (1 | col) + (1 | tp), 
#           data = Pve_overall_effect)
# Linear Hypotheses:
#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0 4.5918317  0.0480226  95.618   <2e-16 ***
#   conc == 0        0.0013425  0.0006784   1.979   0.0478 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  (Univariate p values reported)



#### -- 5.2.1.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Pve, to test differences at t1
Pve_dark_t1 <- subset(Dark_relative, spec == "Pve" & tp == "1")

# write LMER
model_t1_Pve <- lmer(log(relativeFv_Fm) ~ treat + (1|col), data = Pve_dark_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.982   0.267
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.267).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.940).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = log(relativeFv_Fm) ~ treat + (1 | col), data = Pve_dark_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.202415   0.148279  -1.365    1.000
# 1 - control == 0   -0.242063   0.148279  -1.632    0.923
# 10 - control == 0  -0.307794   0.148279  -2.076    0.379
# 100 - control == 0 -0.198319   0.148279  -1.337    1.000
# 1 - 0.1 == 0       -0.039648   0.148279  -0.267    1.000
# 10 - 0.1 == 0      -0.105379   0.148279  -0.711    1.000
# 100 - 0.1 == 0      0.004095   0.148279   0.028    1.000
# 10 - 1 == 0        -0.065731   0.148279  -0.443    1.000
# 100 - 1 == 0        0.043743   0.148279   0.295    1.000
# 100 - 10 == 0       0.109475   0.148279   0.738    1.000
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Pve, to test differences at t2
Pve_dark_t2 <- subset(Dark_relative, spec == "Pve" & tp == "2")

# LMER didn't show good fit - GLMER was used
model_t2_Pve <- glmer((relativeFv_Fm) ~ treat + (1|col), family = poisson, data = Pve_dark_t2)

# get summary of GLMER
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
Pve_dark_t3 <- subset(Dark_relative, spec == "Pve" & tp == "3")

# write LMER
model_t3_Pve <- lmer(log(relativeFv_Fm) ~ treat + (1|col), data = Pve_dark_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.975   0.0817
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.781).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = log(relativeFv_Fm) ~ treat + (1 | col), data = Pve_dark_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.008981   0.145580   0.062    1.000
# 1 - control == 0   -0.088674   0.145580  -0.609    1.000
# 10 - control == 0  -0.027093   0.145580  -0.186    1.000
# 100 - control == 0  0.272418   0.145580   1.871    0.490
# 1 - 0.1 == 0       -0.097656   0.145580  -0.671    1.000
# 10 - 0.1 == 0      -0.036075   0.145580  -0.248    1.000
# 100 - 0.1 == 0      0.263437   0.145580   1.810    0.493
# 10 - 1 == 0         0.061581   0.145580   0.423    1.000
# 100 - 1 == 0        0.361093   0.145580   2.480    0.131
# 100 - 10 == 0       0.299511   0.145580   2.057    0.357
# (Adjusted p values reported -- holm method)



### --- 5.2.2. Stylophora pistillata -------------------------------------------
#### -- 5.2.2.1 Overall effect -------------------------------------------------
# create a subset with data of Spi, t0 excluded for continuous model
Spi_overall_effect <- subset(Dark_relative, spec == "Spi" & tp!= "0")

# LMER didn't show good fit - GLMER was used
model1_Spi <- glmer((relativeFv_Fm) ~ conc + (1|col) + (1|tp), family = "poisson", data = Spi_overall_effect)

# get summary of GLMER
cftest(model1_Spi)
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (relativeFv_Fm) ~ conc + (1 | col) + (1 | tp), 
# data = Spi_overall_effect, family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0  4.8023397  0.7071354   6.791 1.11e-11 ***
#   conc == 0        -0.0012998  0.0001484  -8.759  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Univariate p values reported)


#### -- 5.2.2.2 Specific effects -----------------------------------------------
# ------------- t1
# create a subset with data of Spi, to test differences at t1
Spi_dark_t1 <- subset(Dark_relative, spec == "Spi" & tp == "1")

# write LMER
model_t1_Spi <- lmer(log(relativeFv_Fm) ~ treat + (1|col), data = Spi_dark_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.974   0.0688
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.069).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.822).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(relativeFv_Fm) ~ treat + (1 | col), data = Spi_dark_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.01512    0.14216   0.106        1
# 1 - control == 0   -0.05006    0.14216  -0.352        1
# 10 - control == 0  -0.02789    0.14216  -0.196        1
# 100 - control == 0 -0.13677    0.14216  -0.962        1
# 1 - 0.1 == 0       -0.06519    0.14216  -0.459        1
# 10 - 0.1 == 0      -0.04301    0.14216  -0.303        1
# 100 - 0.1 == 0     -0.15190    0.14216  -1.068        1
# 10 - 1 == 0         0.02217    0.14216   0.156        1
# 100 - 1 == 0       -0.08671    0.14216  -0.610        1
# 100 - 10 == 0      -0.10888    0.14216  -0.766        1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Spi, to test differences at t2
Spi_dark_t2 <- subset(Dark_relative, spec == "Spi" & tp == "2")

# LMER didn't show good fit - GLMER was used
model_t2_Spi <- glmer((relativeFv_Fm) ~ treat + (1|col), family = "poisson", data = Spi_dark_t2)

# get summary of GLMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (relativeFv_Fm) ~ treat + (1 | col), data = Spi_dark_t2, 
#       family = "poisson")
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0 -0.14483    0.03064  -4.727 2.03e-05 ***
# 1 - control == 0   -0.10199    0.03029  -3.367  0.00456 ** 
# 10 - control == 0  -0.03756    0.02979  -1.261  0.51704    
# 100 - control == 0 -0.18529    0.03097  -5.982 2.20e-08 ***
# 1 - 0.1 == 0        0.04284    0.03139   1.365  0.51704    
# 10 - 0.1 == 0       0.10727    0.03091   3.471  0.00363 ** 
# 100 - 0.1 == 0     -0.04046    0.03205  -1.262  0.51704    
# 10 - 1 == 0         0.06443    0.03056   2.108  0.14009    
# 100 - 1 == 0       -0.08330    0.03172  -2.626  0.04319 *  
# 100 - 10 == 0      -0.14773    0.03124  -4.729 2.03e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t3
Spi_dark_t3 <- subset(Dark_relative, spec == "Spi" & tp == "3")

# write LMER
model_t3_Spi <- lmer(log(relativeFv_Fm) ~ treat + (1|col), data = Spi_dark_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.9992  0.848
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.848).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.960).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = log(relativeFv_Fm) ~ treat + (1 | col), data = Spi_dark_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.07836    0.14464   0.542    1.000
# 1 - control == 0   -0.01482    0.14464  -0.102    1.000
# 10 - control == 0   0.09984    0.14464   0.690    1.000
# 100 - control == 0 -0.16326    0.14464  -1.129    1.000
# 1 - 0.1 == 0       -0.09319    0.14464  -0.644    1.000
# 10 - 0.1 == 0       0.02148    0.14464   0.148    1.000
# 100 - 0.1 == 0     -0.24162    0.14464  -1.671    0.853
# 10 - 1 == 0         0.11466    0.14464   0.793    1.000
# 100 - 1 == 0       -0.14844    0.14464  -1.026    1.000
# 100 - 10 == 0      -0.26310    0.14464  -1.819    0.689
# (Adjusted p values reported -- holm method)



# ----- 6. Write tables --------------------------------------------------------
## ---- 6.1. Table of all light adapted PAM data -------------------------------
write_rds(light_all, "processed/light_all.rds")

## ---- 6.2. Table of all dark adapted PAM data --------------------------------
write_rds(dark_all, "processed/dark_all.rds")
