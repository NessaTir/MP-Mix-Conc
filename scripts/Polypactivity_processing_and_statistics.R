# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the data processing and statistical analyzes of the corals polyp activities.
# Polyp activity was assessed at each sampling timepoint (t0, t1, t2, and t3) at three times a day (11 am, 2 pm and 5 pm).
# The polyps were assessed visually and each coral fragment was categorised as active, moderately activ, or inactive.
# Statistical analyzes will be conducted using the exact fisher test.

# ----- 2. Load in needed packages ---------------------------------------------
# to easily clean data, to read in .rds files 
library(tidyverse)

# to easily read in all data files
library(readxl)


# for statistical testing
library(multcomp)

library(rstatix)

# for statistical analyses using LMER and GLMER
library(lme4)

# check model fits statistically
library(performance)
library(rstatix)
# check model fits visually using qqplot
library(car)


# ----- 3. Read in needed data files ------------------------------------------- 
## ---- 3.1. Coral identity table ----------------------------------------------
# read in list with overview of all corals used and their treatments etc.
corals <- read.csv2("in/coral_treatments.csv", sep=",") %>%
  # modify character of some columns, where necessary
  mutate(treat = as.factor(treat), # column for categorical model
         conc = as.numeric(conc)) # column for continuous model 

## ---- 3.2. Polyp activity tables ---------------------------------------------
# read in Data table of polyp activity of Pocillopora verrucosa - wide format
polyp_data_Pve <- read.csv2("in/Polypactivity_Pve.csv") %>%
  # rename column of coral ID to merge with coral info table
  rename(ID = coral) %>%
  # remove unnecessary columns for clear merge
  dplyr::select(-col, -tank, -origin) 

# read in Data table of polyp activity of Stylophora pistillata - wide format
polyp_data_Spi <- read.csv2("in/Polypactivity_Spi.csv") %>%
  # rename column of coral ID to merge with coral info table
  rename(ID = coral) %>%
  # remove unnecessary columns for clear merge
  dplyr::select(-col, -tank, -origin) 

# bring tables of Pve and Spi together
polyp_data <- rbind(polyp_data_Pve, polyp_data_Spi)
# bring tables of polyp activity together coral information table - wide format
polyp_data_wide <-  merge(corals, polyp_data, by = 'ID', all.x = TRUE)


# ----- 4. Prepare data for statistical analyzes -------------------------------
## ---- 4.1. Reformat data table -----------------------------------------------
Polyps <- polyp_data_wide  %>% 
  # rename columns to get continuous timepoints for statistical analyses
  rename("0" = "PA_t0") %>%
  rename("1" = "PA_t1") %>% 
  rename("2" = "PA_t2") %>%
  rename("3" = "PA_t3") %>%
  # bring table into long format
  pivot_longer(cols=c('0', '1', '2', '3'), 
               names_to='tp',  # assign new name for the column of timepoints
               values_to='activity') # assign new name for the values previously in the column under the headers above

# level polyp activity
Polyps$activity <- factor(Polyps$activity, 
                          levels = c("a", "ma", "ia"))

# level treatment, important for visualisation 
Polyps$treat <- factor(Polyps$treat, 
                       levels = c("control", "0.1", "1", "10", "100"))


## ---- 4.2. Convert categories of polyp activity in numbers -------------------
# convert active into            '1'
# convert moderately active into '0.5'
# convert inactive into          '0'
Polyps <- Polyps %>% 
  # create a new column (ranks) with converted categories
  mutate(ranks = case_when(activity == "a"~ "1",
                           activity == "ma"~ "0.5",
                           activity == "ia"~ "0"),
         ranks = as.numeric(ranks),
         tp = as.numeric(tp))



# ----- 5. Statistical analyses ------------------------------------------------
# 1. an overall analysis is conducted per species 
#   to evaluate whether an overall concentration dependent effect was observed.
#   Therefore, the treatment is used as continuous numerical value (conc), set as fixed factor.
#   Here the colony (col) and the time (time) is set as random factor.
# 2. statistical analyses are split by species and time
#   to evaluate specific differences over the course of the experiment, that might be shadowed in the overall analyses.
#   Therefore, the treatment is used as categorical value (treat), set as fixed factor.
#   Here the colony (col) is set as random factor.

# For statistical analyses LMER was used. If LMER didn't fit test assumptions, GLMER was used instead

## ---- 5.1. Pocillopora verrucosa ---------------------------------------------

### --- 5.1.1. Overall effect --------------------------------------------------
# create a subset with data of Pve, t0 excluded for continuous model
Pve_overall_effect <- subset(Polyps, spec == "Pve" & tp!= "0")

# create a table with mean values
Pve_sum_1 <- Pve_overall_effect %>%
  group_by(tp, ID, conc, col) %>%
  get_summary_stats(ranks, type = "mean")

# LMER didn't show a good fit, therefore GLMER is used
model1_Pve <- glmer((mean+100) ~ conc + (1|col) + (1|tp), family = poisson, data = Pve_sum_1)

# get summary of GLMER
cftest(model1_Pve)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (mean + 100) ~ conc + (1 | col) + (1 | tp), data = Pve_sum_1, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) == 0  4.614e+00  7.071e-01   6.524 6.83e-11 ***
#   conc == 0        -2.586e-05  1.553e-04  -0.167    0.868    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


### --- 5.1.2. Specific effects ------------------------------------------------
# ------------- t0
# create a subset with data of Pve, to test differences at t0
Pve_t0 <- subset(Polyps, spec == "Pve" & tp == "0")

# create a table with mean values
Pve_sum_t0 <- Pve_t0 %>%
  group_by(tp, ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# LMER didn't show a good fit, therefore GLMER is used
model_t0_Pve <- glmer((mean+100) ~ treat + (1|col), family = poisson, data = Pve_sum_t0)

# get summary of GLMER
summary(glht(model_t0_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (mean + 100) ~ treat + (1 | col), data = Pve_sum_t0, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  5.520e-04  3.319e-02   0.017        1
# 1 - control == 0    1.010e-03  3.319e-02   0.030        1
# 10 - control == 0   7.349e-04  3.319e-02   0.022        1
# 100 - control == 0  5.509e-04  3.319e-02   0.017        1
# 1 - 0.1 == 0        4.581e-04  3.319e-02   0.014        1
# 10 - 0.1 == 0       1.828e-04  3.319e-02   0.006        1
# 100 - 0.1 == 0     -1.102e-06  3.319e-02   0.000        1
# 10 - 1 == 0        -2.753e-04  3.318e-02  -0.008        1
# 100 - 1 == 0       -4.592e-04  3.319e-02  -0.014        1
# 100 - 10 == 0      -1.839e-04  3.319e-02  -0.006        1
# (Adjusted p values reported -- holm method)


# ------------- t1
# create a subset with data of Pve, to test differences at t1
Pve_t1 <- subset(Polyps, spec == "Pve" & tp == "1")

# create a table with mean values
Pve_sum_t1 <- Pve_t1 %>%
  group_by(ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model_t1_Pve <- lmer(mean^2 ~ treat + (1|col), data = Pve_sum_t1)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.991   0.796
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.796).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.667).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = mean^2 ~ treat + (1 | col), data = Pve_sum_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0  0.03239    0.08726   0.371   0.7105    
#   1 - control == 0   -0.12653    0.08726  -1.450   0.4411    
#   10 - control == 0  -0.24078    0.08726  -2.759   0.0348 *  
#   100 - control == 0 -0.39951    0.08726  -4.578 4.22e-05 ***
#   1 - 0.1 == 0       -0.15892    0.08726  -1.821   0.3429    
#   10 - 0.1 == 0      -0.27317    0.08726  -3.130   0.0140 *  
#   100 - 0.1 == 0     -0.43190    0.08726  -4.950 7.44e-06 ***
#   10 - 1 == 0        -0.11425    0.08726  -1.309   0.4411    
#   100 - 1 == 0       -0.27298    0.08726  -3.128   0.0140 *  
#   100 - 10 == 0      -0.15873    0.08726  -1.819   0.3429    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Pve, to test differences at t2
Pve_t2 <- subset(Polyps, spec == "Pve" & tp == "2")

# create a table with mean values
Pve_sum_t2 <- Pve_t2 %>%
  group_by(ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model_t2_Pve <- lmer(scale(mean) ~ treat + (1|col), data = Pve_sum_t2)

# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.982   0.248
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.248).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.518).

# get summary of LMER
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Pve_sum_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0 -0.03893    0.18966  -0.205   1.0000    
#   1 - control == 0   -0.50799    0.18966  -2.678   0.0444 *  
#   10 - control == 0  -0.15620    0.18966  -0.824   1.0000    
#   100 - control == 0 -1.99327    0.18966 -10.510  < 2e-16 ***
#   1 - 0.1 == 0       -0.46906    0.18966  -2.473   0.0670 .  
#   10 - 0.1 == 0      -0.11727    0.18966  -0.618   1.0000    
#   100 - 0.1 == 0     -1.95434    0.18966 -10.304  < 2e-16 ***
#   10 - 1 == 0         0.35180    0.18966   1.855   0.2545    
#   100 - 1 == 0       -1.48528    0.18966  -7.831 3.42e-14 ***
#   100 - 10 == 0      -1.83708    0.18966  -9.686  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Pve, to test differences at t3
Pve_t3 <- subset(Polyps, spec == "Pve" & tp == "3")

# create a table with mean values
Pve_sum_t3 <- Pve_t3 %>%
  group_by(ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model_t3_Pve <- lmer(mean^2 ~ treat + (1|col), data = Pve_sum_t3)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.981   0.206
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.206).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.592).

# get summary of LMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = mean^2 ~ treat + (1 | col), data = Pve_sum_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
#   0.1 - control == 0  0.155969   0.091407   1.706   0.7036  
#   1 - control == 0    0.009278   0.091407   0.101   1.0000  
#   10 - control == 0  -0.078630   0.091407  -0.860   1.0000  
#   100 - control == 0 -0.111018   0.091407  -1.215   1.0000  
#   1 - 0.1 == 0       -0.146691   0.091407  -1.605   0.7597  
#   10 - 0.1 == 0      -0.234599   0.091407  -2.567   0.0925 .
#   100 - 0.1 == 0     -0.266988   0.091407  -2.921   0.0349 *
#   10 - 1 == 0        -0.087907   0.091407  -0.962   1.0000  
#   100 - 1 == 0       -0.120296   0.091407  -1.316   1.0000  
#   100 - 10 == 0      -0.032389   0.091407  -0.354   1.0000  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



## ---- 5.2. Stylophora pistillata ---------------------------------------------

### --- 5.2.1. Overall effect --------------------------------------------------
# create a subset with data of Spi t0 excluded for continuous model
Spi_overall_effect <- subset(Polyps, spec == "Spi" & tp!= "0")

# create a table with mean values
Spi_sum_1 <- Spi_overall_effect %>%
  group_by(tp, ID, conc, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model1_Spi <- lmer(scale(mean) ~ conc + (1|col) + (1|tp), data = Spi_sum_1)

# inspect residuals
qqPlot(residuals(model1_Spi))          # good fit
shapiro_test(residuals(model1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable              statistic p.value
# <chr>                     <dbl>   <dbl>
#   1 residuals(model1_Spi)     0.990  0.0510
check_normality(model1_Spi)
# OK: residuals appear as normally distributed (p = 0.052).
check_heteroscedasticity(model1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.672).

cftest(model1_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = scale(mean) ~ conc + (1 | col) + (1 | tp), data = Spi_sum_1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0  0.120578   0.316213   0.381    0.703    
# conc == 0        -0.005427   0.001086  -4.995 5.89e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


### --- 5.2.2. Specific effects ------------------------------------------------
# ------------- t0
# create a subset with data of Spi, to test differences at t0
Spi_t0 <- subset(Polyps, spec == "Spi" & tp == "0")

# create a table with mean values
Spi_sum_t0 <- Spi_t0 %>%
  group_by(tp, ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model_t0_Spi <- lmer(scale(mean) ~ treat + (1|col), data = Spi_sum_t0)

# inspect residuals
qqPlot(residuals(model_t0_Spi))          # good fit
shapiro_test(residuals(model_t0_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t0_Spi)     0.992   0.893
check_normality(model_t0_Spi)
# OK: residuals appear as normally distributed (p = 0.893).
check_heteroscedasticity(model_t0_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.639).

# get summary of LMER
summary(glht(model_t0_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Spi_sum_t0)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.29955    0.20369  -1.471        1
# 1 - control == 0   -0.03750    0.20369  -0.184        1
# 10 - control == 0  -0.15000    0.20369  -0.736        1
# 100 - control == 0 -0.11205    0.20369  -0.550        1
# 1 - 0.1 == 0        0.26205    0.20369   1.286        1
# 10 - 0.1 == 0       0.14955    0.20369   0.734        1
# 100 - 0.1 == 0      0.18750    0.20369   0.921        1
# 10 - 1 == 0        -0.11250    0.20369  -0.552        1
# 100 - 1 == 0       -0.07455    0.20369  -0.366        1
# 100 - 10 == 0       0.03795    0.20369   0.186        1
# (Adjusted p values reported -- holm method)


# ------------- t1
# create a subset with data of Spi, to test differences at t1
Spi_t1 <- subset(Polyps, spec == "Spi" & tp == "1")

# create a table with mean values
Spi_sum_t1 <- Spi_t1 %>%
  group_by(ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER 
model_t1_Spi <- lmer(scale(mean) ~ treat + (1|col), data = Spi_sum_t1)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.981   0.209
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.209).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.994).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Spi_sum_t1)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0 -0.03628    0.17168  -0.211   0.8326    
#   1 - control == 0    0.17965    0.17168   1.046   0.8303    
#   10 - control == 0  -0.25264    0.17168  -1.472   0.7056    
#   100 - control == 0 -0.79182    0.17168  -4.612 3.58e-05 ***
#   1 - 0.1 == 0        0.21593    0.17168   1.258   0.8303    
#   10 - 0.1 == 0      -0.21636    0.17168  -1.260   0.8303    
#   100 - 0.1 == 0     -0.75554    0.17168  -4.401 8.62e-05 ***
#   10 - 1 == 0        -0.43229    0.17168  -2.518   0.0708 .  
#   100 - 1 == 0       -0.97148    0.17168  -5.659 1.52e-07 ***
#   100 - 10 == 0      -0.53918    0.17168  -3.141   0.0118 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2
# create a subset with data of Spi, to test differences at t2
Spi_t2 <- subset(Polyps, spec == "Spi" & tp == "2")

# create a table with mean values
Spi_sum_t2 <- Spi_t2 %>%
  group_by(ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model_t2_Spi <- lmer(scale(mean) ~ treat + (1|col), data = Spi_sum_t2)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.989   0.692
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.692).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.680).

# get summary of lmer
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(mean) ~ treat + (1 | col), data = Spi_sum_t2)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0  -0.1005     0.2196  -0.458  0.88945   
#   1 - control == 0    -0.6713     0.2196  -3.057  0.01789 * 
#   10 - control == 0   -0.4362     0.2196  -1.986  0.28200   
#   100 - control == 0  -0.8391     0.2196  -3.821  0.00133 **
#   1 - 0.1 == 0        -0.5708     0.2196  -2.599  0.06542 . 
#   10 - 0.1 == 0       -0.3357     0.2196  -1.528  0.50557   
#   100 - 0.1 == 0      -0.7386     0.2196  -3.363  0.00693 **
#   10 - 1 == 0          0.2351     0.2196   1.071  0.85295   
#   100 - 1 == 0        -0.1678     0.2196  -0.764  0.88945   
#   100 - 10 == 0       -0.4029     0.2196  -1.835  0.33259   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t3
# create a subset with data of Spi, to test differences at t3
Spi_t3 <- subset(Polyps, spec == "Spi" & tp == "3")

# create a table with mean values
Spi_sum_t3 <- Spi_t3 %>%
  group_by(ID, treat, col) %>%
  get_summary_stats(ranks, type = "mean")

# write LMER
model_t3_Spi <- lmer(mean^2 ~ treat + (1|col), data = Spi_sum_t3)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.990   0.728
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.728).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.416).

# get summary of LMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = mean^2 ~ treat + (1 | col), data = Spi_sum_t3)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.057327   0.087021  -0.659    1.000
# 1 - control == 0    0.129426   0.087021   1.487    0.959
# 10 - control == 0  -0.074204   0.087021  -0.853    1.000
# 100 - control == 0 -0.066401   0.087021  -0.763    1.000
# 1 - 0.1 == 0        0.186753   0.087021   2.146    0.255
# 10 - 0.1 == 0      -0.016877   0.087021  -0.194    1.000
# 100 - 0.1 == 0     -0.009074   0.087021  -0.104    1.000
# 10 - 1 == 0        -0.203630   0.087021  -2.340    0.193
# 100 - 1 == 0       -0.195827   0.087021  -2.250    0.220
# 100 - 10 == 0       0.007803   0.087021   0.090    1.000
# (Adjusted p values reported -- holm method)


# ----- 6. Write tables --------------------------------------------------------
## ---- 6.1. Table of polyp activities -----------------------------------------
# for creation of plots for visualization > Script 'Plots'
write_rds(Polyps, "processed/polyp_activity.rds")

## ---- 6.2. Tables of summary -------------------------------------------------
# create a table to summarize the mean polyp activity for P. verrucosa
Pve <- subset(Polyps, spec == "Pve")
Pve_sum <- Pve %>%
  group_by(tp, ID, conc, col) %>%
  get_summary_stats(ranks, type = "mean")

write_csv2(Pve_sum, "out/polyp_mean_Pve.csv")

# create a table to summarize the mean polyp activity for S. pistillata
Spi <- subset(Polyps, spec == "Spi")
Spi_sum <- Spi %>%
  group_by(tp, ID, conc, col) %>%
  get_summary_stats(ranks, type = "mean")

write_csv2(Spi_sum, "out/polyp_mean_Spi.csv")

