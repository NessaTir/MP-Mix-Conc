# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the statistical analyzes of the corals growth rates. Growth was measured in
#   a) Living tissue area growth
#       --> Necrosis was analyzed in 'Necrosis_statistics'
#   b) Volume growth
#   c) Calcification growth
# This script builds up on data tables produced in the script 'Growth_data_processing'
#   here grwoth rates were calculated for each 4 week time block of microplastic exposure
#     t0-t1: week 0 to 4 after the addition of microplastic
#     t1-t2: week 4 to 8 after the addition of microplastic
#     t2-t3: week 8 to 12 after the addition of microplastic
# Statistical analyzes will be conducted using LMER and GLMER 
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
# check model fits visually using qqplot
library(car)


# ---- 3. Read in needed data files --------------------------------------------
##---- 3.1. Surface rates ------------------------------------------------------
surface <- read_rds("processed/surface_growth.rds") %>%
  mutate(treat = as.factor(treat),
         conc = as.numeric(conc))

##---- 3.2. volume rates -------------------------------------------------------
volume <- read_rds("processed/volume_growth.rds") %>%
  mutate(treat = as.factor(treat),
         conc = as.numeric(conc))

##---- 3.3. calcifiaction rates ------------------------------------------------
calcification <- read_rds("processed/weight_growth.rds") %>%
  mutate(treat = as.factor(treat),
         conc = as.numeric(conc))


# ---- 4. Prepare data for statistical analyzes --------------------------------
# relevel treatments in all tables
## --- 4.1. Surface ------------------------------------------------------------
surface$treat <- factor(surface$treat, 
                              levels = c("control", "0.1", "1",
                                         "10", "100"))
levels(surface$treat)

## --- 4.2. Volume -------------------------------------------------------------
volume$treat <- factor(volume$treat, 
                        levels = c("control", "0.1", "1",
                                   "10", "100"))
levels(volume$treat)

## --- 4.3. Calcification ------------------------------------------------------
calcification$treat <- factor(calcification$treat, 
                        levels = c("control", "0.1", "1",
                                   "10", "100"))
levels(calcification$treat)


# ---- 5. Statistical analyses -------------------------------------------------
# 1. an overall analysis is conducted per growth parameter and species 
#   to evaluate whether an overall concentration dependent effect was observed.
#   Therefore, the treatment is used as continous numerical value (conc), set as fixed factor.
#   Here the colony (col) and the time (time) is set as random factor.
# 2. statistical analyses are split by growth parameter, species and time
#   to evaluate specific differences over the course of the experiment, that might be shadowed in the overall analyses.
#   Therefore, the treatment is used as categorical value (treat), set as fixed factor.
#   Here the colony (col) is set as random factor.

# For statistical analyses LMER was used. If LMER did't fit test assumtions, GLMER was used instead

## --- 5.1. Surface growth -----------------------------------------------------
### -- 5.1.1. Pocillopora verrucosa --------------------------------------------

#### - 5.1.1.1. Overall effect -------------------------------------------------
# create a subset with data of Pve
Pve_surf <- subset(surface, spec == "Pve")

# LMER didn't show a good fit, therefore GLMER is used
model_Pve <- glmer(((surface_growth+100)) ~ conc + (1|col) + (1|time), family = "poisson", data = Pve_surf)
# summary of tested with the GLMER differences
cftest(model_Pve)

# OUTPUT:
#         Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = ((surface_growth + 100)) ~ conc + (1 | col) + 
#             (1 | time), data = Pve_surf, family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0  4.7155330  0.7071380   6.668 2.58e-11 ***
#   conc == 0        -0.0005805  0.0001508  -3.851 0.000118 ***
#  ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Univariate p values reported)


#### - 5.1.1.2. Specific effects -----------------------------------------------
# ------------- t0-t1
# create a subset with data of Pve for t0 - t1
Pve_1_surf <- subset(Pve_surf, time == "1")

# LMER didn't show a good fit, therefore GLMER is used
model_t1_Pve <- glmer(((surface_growth+100)) ~ treat + (1|col), family = poisson, data = Pve_1_surf)

# get summary of glmer
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))

# OUTPUT:
#	 Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = ((surface_growth + 100)) ~ treat + (1 | col), 
#            data = Pve_1_surf, family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.023561   0.032513  -0.725        1
# 1 - control == 0   -0.033689   0.032597  -1.033        1
# 10 - control == 0  -0.039761   0.033139  -1.200        1
# 100 - control == 0 -0.043056   0.032675  -1.318        1
# 1 - 0.1 == 0       -0.010127   0.032787  -0.309        1
# 10 - 0.1 == 0      -0.016199   0.033326  -0.486        1
# 100 - 0.1 == 0     -0.019495   0.032865  -0.593        1
# 10 - 1 == 0        -0.006072   0.033408  -0.182        1
# 100 - 1 == 0       -0.009368   0.032947  -0.284        1
# 100 - 10 == 0      -0.003296   0.033484  -0.098        1
#  (Adjusted p values reported -- holm method)


# ------------- t1-t2
# create a subset with data of Pve for t1 - t2
Pve_2_surf <- subset(Pve_surf, time == "2")

# write lmer 
# transform data to fit LMER
model_t2_Pve <- lmer(log(surface_growth+10) ~ treat + (1|col), data = Pve_2_surf)

# check fit of LMER
# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.992   0.849
check_normality(model_t2_Pve)
# OUTPUT: OK: residuals appear as normally distributed (p = 0.849).
check_heteroscedasticity(model_t2_Pve)
# OUTPUT: OK: Error variance appears to be homoscedastic (p = 0.992).

# get summary of lmer
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = log(surface_growth + 10) ~ treat + (1 | col), 
#           data = Pve_2_surf)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0  0.192014   0.135272   1.419   1.0000  
# 1 - control == 0    0.144563   0.135272   1.069   1.0000  
# 10 - control == 0   0.187751   0.135272   1.388   1.0000  
# 100 - control == 0 -0.187242   0.135272  -1.384   1.0000  
# 1 - 0.1 == 0       -0.047451   0.135272  -0.351   1.0000  
# 10 - 0.1 == 0      -0.004263   0.135272  -0.032   1.0000  
# 100 - 0.1 == 0     -0.379256   0.135272  -2.804   0.0505 .
# 10 - 1 == 0         0.043188   0.135272   0.319   1.0000  
# 100 - 1 == 0       -0.331805   0.135272  -2.453   0.1134  
# 100 - 10 == 0      -0.374993   0.135272  -2.772   0.0505 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2-t3
# create a subset with data of Pve for t2 - t3
Pve_3_surf <- subset(Pve_surf, time == "3")

# LMER didn't show a good fit, therefore GLMER is used
model_t3_Pve <- glmer((surface_growth+100) ~ treat + (1|col), family = poisson, data = Pve_3_surf)

# get summary of GLMER
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (surface_growth + 100) ~ treat + (1 | col), data = Pve_3_surf, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0  0.003051   0.030683   0.099  1.00000   
#   1 - control == 0    0.006757   0.030655   0.220  1.00000   
#   10 - control == 0   0.009861   0.030631   0.322  1.00000   
#   100 - control == 0 -0.095538   0.031467  -3.036  0.01677 * 
#   1 - 0.1 == 0        0.003706   0.030632   0.121  1.00000   
#   10 - 0.1 == 0       0.006810   0.030608   0.222  1.00000   
#   100 - 0.1 == 0     -0.098590   0.031444  -3.135  0.01373 * 
#   10 - 1 == 0         0.003104   0.030579   0.102  1.00000   
#   100 - 1 == 0       -0.102295   0.031416  -3.256  0.01017 * 
#   100 - 10 == 0      -0.105399   0.031393  -3.357  0.00787 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



### -- 5.1.2. Stylophora pistillata --------------------------------------------
#### - 5.1.2.1. Overall effect -------------------------------------------------
# create a subset with data of Spi
Spi_surf <- subset(surface, spec == "Spi")

# LMER didn't show a good fit, therefore GLMER is used
model_Spi <- glmer(((surface_growth+100)) ~ conc + (1|col) + (1|time), family = "poisson", data = Spi_surf)
# summary of tested with the GLMER differences
cftest(model_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = ((surface_growth + 100)) ~ conc + (1 | col) + 
#              (1 | time), data = Spi_surf, family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) == 0  4.7116149  0.7071394   6.663 2.68e-11 ***
#   conc == 0        -0.0006939  0.0001522  -4.561 5.10e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


#### - 5.1.1.2. Specific effects -----------------------------------------------
# ------------- t0-t1
# create a subset with data of Spi for t0 - t1
Spi_1_surf <- subset(Spi_surf, time == "1")

# LMER didn't show a good fit, therefore GLMER is used
model_t1_Spi <- glmer((surface_growth+100) ~ treat + (1|col), family = poisson, data = Spi_1_surf)

# get summary of glmer
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
       test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (surface_growth + 100) ~ treat + (1 | col), data = Spi_1_surf, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.024409   0.033668  -0.725    1.000
# 1 - control == 0   -0.012891   0.033159  -0.389    1.000
# 10 - control == 0  -0.021947   0.033220  -0.661    1.000
# 100 - control == 0 -0.061072   0.033046  -1.848    0.646
# 1 - 0.1 == 0        0.011518   0.033891   0.340    1.000
# 10 - 0.1 == 0       0.002462   0.033876   0.073    1.000
# 100 - 0.1 == 0     -0.036663   0.033767  -1.086    1.000
# 10 - 1 == 0        -0.009056   0.033374  -0.271    1.000
# 100 - 1 == 0       -0.048181   0.033200  -1.451    1.000
# 100 - 10 == 0      -0.039125   0.033261  -1.176    1.000
# (Adjusted p values reported -- holm method)


# ------------- t1-t2
# create a subset with data of Spi for t1 - t2
Spi_2_surf <- subset(Spi_surf, time == "2")

# LMER didn't show a good fit, therefore GLMER is used
model_t2_Spi <- glmer((surface_growth+100) ~ treat + (1|col), family = poisson, data = Spi_2_surf)

# get summary of GLMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (surface_growth + 100) ~ treat + (1 | col), data = Spi_2_surf, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0  0.007146   0.032412   0.220   1.0000  
# 1 - control == 0   -0.013474   0.032571  -0.414   1.0000  
# 10 - control == 0   0.021606   0.032303   0.669   1.0000  
# 100 - control == 0 -0.070316   0.033020  -2.129   0.2657  
# 1 - 0.1 == 0       -0.020620   0.031595  -0.653   1.0000  
# 10 - 0.1 == 0       0.014460   0.031319   0.462   1.0000  
# 100 - 0.1 == 0     -0.077462   0.032058  -2.416   0.1411  
# 10 - 1 == 0         0.035080   0.031483   1.114   1.0000  
# 100 - 1 == 0       -0.056842   0.032218  -1.764   0.5438  
# 100 - 10 == 0       -0.091922   0.031948  -2.877   0.0401 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2-t3
# create a subset with data of Spi for t2 - t3
Spi_3_surf <- subset(Spi_surf, time == "3")

# LMER didn't show a good fit, therefore GLMER is used
model_t3_Spi <- glmer((surface_growth+100) ~ treat + (1|col), family = poisson, data = Spi_3_surf)

# get summary of GLMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (surface_growth + 100) ~ treat + (1 | col), data = Spi_3_surf, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
# 0.1 - control == 0 -0.003798   0.031499  -0.121  1.00000   
# 1 - control == 0   -0.007852   0.031530  -0.249  1.00000   
# 10 - control == 0   0.030748   0.031238   0.984  1.00000   
# 100 - control == 0 -0.082715   0.032121  -2.575  0.09019 . 
# 1 - 0.1 == 0       -0.004054   0.031083  -0.130  1.00000   
# 10 - 0.1 == 0       0.034545   0.030786   1.122  1.00000   
# 100 - 0.1 == 0     -0.078917   0.031682  -2.491  0.10194   
# 10 - 1 == 0         0.038600   0.030818   1.253  1.00000   
# 100 - 1 == 0       -0.074862   0.031713  -2.361  0.12771   
# 100 - 10 == 0      -0.113462   0.031423  -3.611  0.00305 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


## --- 5.2. Volume growth ------------------------------------------------------
### -- 5.2.1. Pocillopora verrucosa --------------------------------------------
#### - 5.2.1.1. Overall effect -------------------------------------------------
# create a subset with data of Pve
Pve_vol <- subset(volume, spec == "Pve")

hist(log((Pve_vol$volume_growth+0.02)))

# Example usage:
# Assuming 'your_data' is the vector or column containing your data
Pve_vol$volume_growth_clean <- remove_outliers(Pve_vol$volume_growth)


# LMER shows a good fit, therefore GLMER is used
model_Pve <- lmer(log(((volume_growth_clean+0.02))) ~ conc + (1|col) + (1|time), data = Pve_vol)
# summary of tested with the LMER differences
cftest(model_Pve)

# inspect residuals
qqPlot(residuals(model_Pve))          # good fit
shapiro_test(residuals(model_Pve))    # p > 0.05 = Non-Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
# 1 residuals(model_Pve)     0.991  0.0986
check_normality(model_Pve)
# OK: residuals appear as normally distributed (p = 0.101).
check_heteroscedasticity(model_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.976).

# OUTPUT:
#Better update yourself, I messed it up :D



#### - 5.2.1.2. Specific effects -----------------------------------------------
# ------------- t0-t1
# create a subset with data of Pve for t0 - t1
Pve_1_vol <- subset(Pve_vol, time == "1")

# write LMER 
model_t1_Pve <- lmer(scale(volume_growth) ~ treat + (1|col), data = Pve_1_vol)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.980   0.177
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.180).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.813).

# get summary of lmer
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(volume_growth) ~ treat + (1 | col), data = Pve_1_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0  0.15365    0.25988   0.591    1.000
# 1 - control == 0   -0.28609    0.25988  -1.101    1.000
# 10 - control == 0   0.20520    0.26388   0.778    1.000
# 100 - control == 0 -0.14669    0.25988  -0.564    1.000
# 1 - 0.1 == 0       -0.43973    0.25988  -1.692    0.816
# 10 - 0.1 == 0       0.05156    0.26388   0.195    1.000
# 100 - 0.1 == 0     -0.30034    0.25988  -1.156    1.000
# 10 - 1 == 0         0.49129    0.26388   1.862    0.626
# 100 - 1 == 0        0.13939    0.25988   0.536    1.000
# 100 - 10 == 0      -0.35190    0.26388  -1.334    1.000
# (Adjusted p values reported -- holm method)


# ------------- t1-t2
# create a subset with data of Pve for t1 - t2
Pve_2_vol <- subset(Pve_vol, time == "2")

# write lmer 
model_t2_Pve <- lmer(scale(volume_growth) ~ treat + (1|col), data = Pve_2_vol)
# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.976   0.101
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.101).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.662).

# get summary of lmer
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#  Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(volume_growth) ~ treat + (1 | col), data = Pve_2_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0  0.01391    0.32179   0.043   1.0000  
# 1 - control == 0    0.17001    0.32179   0.528   1.0000  
# 10 - control == 0   0.29140    0.32179   0.906   1.0000  
# 100 - control == 0 -0.66058    0.32179  -2.053   0.2886  
# 1 - 0.1 == 0        0.15610    0.32179   0.485   1.0000  
# 10 - 0.1 == 0       0.27749    0.32179   0.862   1.0000  
# 100 - 0.1 == 0     -0.67449    0.32179  -2.096   0.2886  
# 10 - 1 == 0         0.12139    0.32179   0.377   1.0000  
# 100 - 1 == 0       -0.83059    0.32179  -2.581   0.0886 .
# 100 - 10 == 0      -0.95198    0.32179  -2.958   0.0309 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Adjusted p values reported -- holm method)


# ------------- t2-t3
# create a subset with data of Pve for t2 - t3
Pve_3_vol <- subset(Pve_vol, time == "3")

# write lmer 
model_t3_Pve <- lmer(scale(volume_growth) ~ treat + (1|col), data = Pve_3_vol)

# inspect residuals
qqPlot(residuals(model_t3_Pve))          # very good fit
shapiro_test(residuals(model_t3_Pve))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Pve)     0.993   0.913
check_normality(model_t3_Pve)
# OK: residuals appear as normally distributed (p = 0.913).
check_heteroscedasticity(model_t3_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.741).

# get summary of lmer
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(volume_growth) ~ treat + (1 | col), data = Pve_3_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0  0.08653    0.27300   0.317  1.00000   
#   1 - control == 0    0.22495    0.27300   0.824  1.00000   
#   10 - control == 0   0.14882    0.27300   0.545  1.00000   
#   100 - control == 0 -0.82043    0.27300  -3.005  0.01857 * 
#   1 - 0.1 == 0        0.13841    0.27300   0.507  1.00000   
#   10 - 0.1 == 0       0.06229    0.27300   0.228  1.00000   
#   100 - 0.1 == 0     -0.90696    0.27300  -3.322  0.00714 **
#   10 - 1 == 0        -0.07612    0.27300  -0.279  1.00000   
#   100 - 1 == 0       -1.04537    0.27300  -3.829  0.00129 **
#   100 - 10 == 0      -0.96925    0.27300  -3.550  0.00346 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Adjusted p values reported -- holm method)



### -- 5.2.2. Stylophora pistillata --------------------------------------------
#### - 5.2.2.1. Overall effect -------------------------------------------------
# create a subset with data of Pve
Spi_vol <- subset(volume, spec == "Spi")

# write LMER
model_Spi <- lmer((scale(volume_growth+100)) ~ conc + (1|col) + (1|time), data = Spi_vol)

# # inspect residuals
qqPlot(residuals(model_Spi))          # good fit
shapiro_test(residuals(model_Spi))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable             statistic p.value
# <chr>                    <dbl>   <dbl>
#   1 residuals(model_Spi)     0.993   0.273
check_normality(model_Spi)
# OK: residuals appear as normally distributed (p = 0.273).
check_heteroscedasticity(model_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.478).

# summary of tested with the LMER differences
cftest(model_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: lmer(formula = (scale(volume_growth + 100)) ~ conc + (1 | col) + 
#             (1 | time), data = Spi_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0  0.087749   0.425073   0.206    0.836    
# conc == 0        -0.004844   0.000962  -5.036 4.76e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


#### - 5.2.2.2. Specific effects -----------------------------------------------
# ------------- t0-t1
# create a subset with data of Spi for t0 - t1
Spi_1_vol <- subset(Spi_vol, time == "1")

# write lmer 
model_t1_Spi <- lmer(scale(volume_growth) ~ treat + (1|col), data = Spi_1_vol)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.987   0.580
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.582).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.757).

# get summary of lmer
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(volume_growth) ~ treat + (1 | col), data = Spi_1_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0 -0.23255    0.24331  -0.956  0.97023    
#   1 - control == 0    0.48949    0.23986   2.041  0.24768    
#   10 - control == 0   0.39856    0.23988   1.662  0.38643    
#   100 - control == 0 -0.46969    0.23624  -1.988  0.24768    
#   1 - 0.1 == 0        0.72203    0.24382   2.961  0.02450 *  
#   10 - 0.1 == 0       0.63111    0.24331   2.594  0.06644 .  
#   100 - 0.1 == 0     -0.23715    0.24016  -0.987  0.97023    
#   10 - 1 == 0        -0.09092    0.23986  -0.379  0.97023    
#   100 - 1 == 0       -0.95918    0.23623  -4.060  0.00049 ***
#   100 - 10 == 0      -0.86826    0.23624  -3.675  0.00214 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t1-t2
# create a subset with data of Spi for t1 - t2
Spi_2_vol <- subset(Spi_vol, time == "2")

# write lmer 
model_t2_Spi <- lmer(scale(volume_growth) ~ treat + (1|col), data = Spi_2_vol)

# inspect residuals
qqPlot(residuals(model_t2_Spi))          # good fit
shapiro_test(residuals(model_t2_Spi))    # p > 0.05 = Normality
# OUTOUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Spi)     0.991   0.819
check_normality(model_t2_Spi)
# OK: residuals appear as normally distributed (p = 0.818).
check_heteroscedasticity(model_t2_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.944).

# get summary of lmer
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(volume_growth) ~ treat + (1 | col), data = Spi_2_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0  0.292397   0.192606   1.518 0.487947    
#   1 - control == 0   -0.155664   0.192606  -0.808 0.837952    
#   10 - control == 0   0.297863   0.192606   1.546 0.487947    
#   100 - control == 0 -0.509709   0.192606  -2.646 0.065087 .  
#   1 - 0.1 == 0       -0.448061   0.186575  -2.402 0.105459    
#   10 - 0.1 == 0       0.005466   0.186575   0.029 0.976629    
#   100 - 0.1 == 0     -0.802106   0.186575  -4.299 0.000154 ***
#   10 - 1 == 0         0.453527   0.186575   2.431 0.105459    
#   100 - 1 == 0       -0.354044   0.186575  -1.898 0.288748    
#   100 - 10 == 0      -0.807572   0.186575  -4.328 0.000150 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2-t3
# create a subset with data of Spi for t1 - t2
Spi_3_vol <- subset(Spi_vol, time == "3")

# write lmer 
model_t3_Spi <- lmer(scale(volume_growth) ~ treat + (1|col), data = Spi_3_vol)

# inspect residuals
qqPlot(residuals(model_t3_Spi))          # good fit
shapiro_test(residuals(model_t3_Spi))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t3_Spi)     0.989   0.636
check_normality(model_t3_Spi)
# OK: residuals appear as normally distributed (p = 0.635).
check_heteroscedasticity(model_t3_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.743).

# get summary of lmer
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(volume_growth) ~ treat + (1 | col), data = Spi_3_vol)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
#   0.1 - control == 0 -0.03329    0.18630  -0.179   1.0000  
#   1 - control == 0    0.08314    0.18630   0.446   1.0000  
#   10 - control == 0   0.10168    0.18630   0.546   1.0000  
#   100 - control == 0 -0.48932    0.18630  -2.627   0.0690 .
#   1 - 0.1 == 0        0.11643    0.18346   0.635   1.0000  
#   10 - 0.1 == 0       0.13497    0.18346   0.736   1.0000  
#   100 - 0.1 == 0     -0.45603    0.18346  -2.486   0.0905 .
#   10 - 1 == 0         0.01854    0.18346   0.101   1.0000  
#   100 - 1 == 0       -0.57246    0.18346  -3.120   0.0163 *
#   100 - 10 == 0      -0.59100    0.18346  -3.221   0.0128 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



## --- 5.3. Weight growth ------------------------------------------------------
### -- 5.3.1. Pocillopora verrucosa --------------------------------------------
#### - 5.3.1.1. Overall effect -------------------------------------------------
# create a subset with data of Pve
Pve_calc <- subset(calcification, spec == "Pve")

# LMER didn't show a good fit, therefore GLMER is used
model_ot_Pve <- glmer(((weight_growth+100)) ~ conc + (1|col) + (1|time), family = "poisson", data = Pve_calc)
# summary of tested with the GLMER differences
cftest(model_ot_Pve)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = ((weight_growth + 100)) ~ conc + (1 | col) + 
#              (1 | time), data = Pve_calc, family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0  4.8575085  0.7071339   6.869 6.45e-12 ***
#   conc == 0        -0.0005011  0.0001398  -3.584 0.000338 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)



#### - 5.3.1.2. Specific effects -----------------------------------------------
# ------------- t0-t1
# create a subset with data of Pve for t0 - t1
Pve_1_calc <- subset(Pve_calc, time == "1")
# write LMER
model_t1_Pve <- lmer(log(weight_growth^2) ~ treat + (1|col), data = Pve_1_calc)

# inspect residuals
qqPlot(residuals(model_t1_Pve))          # good fit
shapiro_test(residuals(model_t1_Pve))    # p > 0.05 = Normality
# A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Pve)     0.977   0.111
check_normality(model_t1_Pve)
# OK: residuals appear as normally distributed (p = 0.109).
check_heteroscedasticity(model_t1_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.628).

# get summary of LMER
summary(glht(model_t1_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = log(weight_growth^2) ~ treat + (1 | col), data = Pve_1_calc)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)
# 0.1 - control == 0 -0.46537    0.34545  -1.347    1.000
# 1 - control == 0   -0.18372    0.34545  -0.532    1.000
# 10 - control == 0   0.11884    0.35077   0.339    1.000
# 100 - control == 0  0.06551    0.34545   0.190    1.000
# 1 - 0.1 == 0        0.28164    0.34545   0.815    1.000
# 10 - 0.1 == 0       0.58421    0.35077   1.665    0.958
# 100 - 0.1 == 0      0.53088    0.34545   1.537    1.000
# 10 - 1 == 0         0.30257    0.35077   0.863    1.000
# 100 - 1 == 0        0.24924    0.34545   0.721    1.000
# 100 - 10 == 0      -0.05333    0.35077  -0.152    1.000
# (Adjusted p values reported -- holm method)


# ------------- t1-t2
# create a subset with data of Pve for t1 - t2
Pve_2_calc <- subset(Pve_calc, time == "2")

# write lmer 
model_t2_Pve <- lmer(scale(weight_growth) ~ treat + (1|col), data = Pve_2_calc)
# inspect residuals
qqPlot(residuals(model_t2_Pve))          # good fit
shapiro_test(residuals(model_t2_Pve))    # p > 0.05 = Normality
# OUTPUT: # A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t2_Pve)     0.987   0.541
check_normality(model_t2_Pve)
# OK: residuals appear as normally distributed (p = 0.541).
check_heteroscedasticity(model_t2_Pve)
# OK: Error variance appears to be homoscedastic (p = 0.820).

# get summary of lmer
summary(glht(model_t2_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = scale(weight_growth) ~ treat + (1 | col), data = Pve_2_calc)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0  0.184362   0.262377   0.703 1.000000    
#   1 - control == 0    0.188211   0.262377   0.717 1.000000    
#   10 - control == 0   0.477703   0.262377   1.821 0.411936    
#   100 - control == 0 -0.665507   0.262377  -2.536 0.078387 .  
#   1 - 0.1 == 0        0.003849   0.262377   0.015 1.000000    
#   10 - 0.1 == 0       0.293341   0.262377   1.118 1.000000    
#   100 - 0.1 == 0     -0.849870   0.262377  -3.239 0.010249 *  
#   10 - 1 == 0         0.289492   0.262377   1.103 1.000000    
#   100 - 1 == 0       -0.853719   0.262377  -3.254 0.010249 *  
#   100 - 10 == 0      -1.143211   0.262377  -4.357 0.000132 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2-t3
# create a subset with data of Pve for t2 - t3
Pve_3_calc <- subset(Pve_calc, time == "3")

# LMER didn't show a good fit, therefore GL MER is used
model_t3_Pve <- glmer(((weight_growth+100)) ~ treat + (1|col), family = "poisson", data = Pve_3_calc)

# get summary of lmer
summary(glht(model_t3_Pve, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = ((weight_growth + 100)) ~ treat + (1 | col), 
#            data = Pve_3_calc, family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)   
#   0.1 - control == 0  0.025612   0.028026   0.914  1.00000   
#   1 - control == 0    0.045302   0.027891   1.624  0.62591   
#   10 - control == 0   0.017626   0.028081   0.628  1.00000   
#   100 - control == 0 -0.059748   0.028635  -2.087  0.25853   
#   1 - 0.1 == 0        0.019690   0.027710   0.711  1.00000   
#   10 - 0.1 == 0      -0.007986   0.027901  -0.286  1.00000   
#   100 - 0.1 == 0     -0.085360   0.028459  -2.999  0.02435 * 
#   10 - 1 == 0        -0.027675   0.027766  -0.997  1.00000   
#   100 - 1 == 0       -0.105049   0.028326  -3.709  0.00208 **
#   100 - 10 == 0      -0.077374   0.028514  -2.714  0.05325 . 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



### -- 5.3.2. Stylophora pistillata --------------------------------------------
#### - 5.3.2.1. Overall effect -------------------------------------------------
# create a subset with data of Spi
Spi_calc <- subset(calcification, spec == "Spi") 

# LMER didn't show a good fit, therefore GLMER is used
model_Spi <- glmer(((weight_growth+100)) ~ conc + (1|col) + (1|time), family = "poisson", data = Spi_calc)
# summary of tested with the GLMER differences
cftest(model_Spi)
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = ((weight_growth + 100)) ~ conc + (1 | col) + 
#              (1 | time), data = Spi_calc, family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept) == 0  4.8670609  0.7071345   6.883 5.87e-12 ***
#   conc == 0        -0.0008224  0.0001410  -5.831 5.52e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)



#### - 5.3.2.2. Specific effects -----------------------------------------------
# ------------- t0-t1
# create a subset with data of Spi for t0 - t1
Spi_1_calc <- subset(Spi_calc, time == "1")

# write LMER
model_t1_Spi <- lmer((weight_growth) ~ treat + (1|col), data = Spi_1_calc)

# inspect residuals
qqPlot(residuals(model_t1_Spi))          # good fit
shapiro_test(residuals(model_t1_Spi))    # p > 0.05 = Normality
# OUTPUT: A tibble: 1 x 3
# variable                statistic p.value
# <chr>                       <dbl>   <dbl>
#   1 residuals(model_t1_Spi)     0.984   0.353
check_normality(model_t1_Spi)
# OK: residuals appear as normally distributed (p = 0.341).
check_heteroscedasticity(model_t1_Spi)
# OK: Error variance appears to be homoscedastic (p = 0.945).

# get summary of LMER
summary(glht(model_t1_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 	 
#      Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lmer(formula = (weight_growth) ~ treat + (1 | col), data = Spi_1_calc)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)  
# 0.1 - control == 0   0.3137     2.0787   0.151   1.0000  
# 1 - control == 0     0.9165     2.0447   0.448   1.0000  
# 10 - control == 0    2.7344     2.0448   1.337   1.0000  
# 100 - control == 0  -3.2236     2.0136  -1.601   0.7657  
# 1 - 0.1 == 0         0.6027     2.1105   0.286   1.0000  
# 10 - 0.1 == 0        2.4206     2.1060   1.149   1.0000  
# 100 - 0.1 == 0      -3.5373     2.0787  -1.702   0.7105  
# 10 - 1 == 0          1.8179     2.0763   0.876   1.0000  
# 100 - 1 == 0        -4.1401     2.0447  -2.025   0.3860  
# 100 - 10 == 0       -5.9580     2.0448  -2.914   0.0357 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Adjusted p values reported -- holm method)


# ------------- t1-t2
# create a subset with data of Spi for t1 - t2
Spi_2_calc <- subset(Spi_calc, time == "2")

# LMER didn't show a good fit, therefore GLMER is used
model_t2_Spi <- glmer((weight_growth+100) ~ treat + (1|col), family = poisson, data = Spi_2_calc)

# get summary of GLMER
summary(glht(model_t2_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (weight_growth + 100) ~ treat + (1 | col), data = Spi_2_calc, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0 -0.0196363  0.0291494  -0.674 1.000000    
#   1 - control == 0   -0.0370122  0.0292750  -1.264 1.000000    
#   10 - control == 0  -0.0366020  0.0292720  -1.250 1.000000    
#   100 - control == 0 -0.1247890  0.0299350  -4.169 0.000306 ***
#   1 - 0.1 == 0       -0.0173759  0.0290608  -0.598 1.000000    
#   10 - 0.1 == 0      -0.0169657  0.0290578  -0.584 1.000000    
#   100 - 0.1 == 0     -0.1051527  0.0297255  -3.537 0.003636 ** 
#   10 - 1 == 0         0.0004102  0.0291837   0.014 1.000000    
#   100 - 1 == 0       -0.0877768  0.0298486  -2.941 0.025033 *  
#   100 - 10 == 0      -0.0881870  0.0298457  -2.955 0.025033 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# ------------- t2-t3
# create a subset with data of Spi for t2 - t3
Spi_3_calc <- subset(Spi_calc, time == "3")

# LMER didn't show a good fit, therefore GLMER is used
model_t3_Spi <- glmer((weight_growth+100) ~ treat + (1|col), family = poisson, data = Spi_3_calc)

# get summary of GMER
summary(glht(model_t3_Spi, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#   Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (weight_growth + 100) ~ treat + (1 | col), data = Spi_3_calc, 
#            family = poisson)
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
#   0.1 - control == 0  0.006266   0.028724   0.218  1.00000    
#   1 - control == 0    0.021594   0.028618   0.755  1.00000    
#   10 - control == 0   0.051205   0.028417   1.802  0.42939    
#   100 - control == 0 -0.084254   0.029376  -2.868  0.02890 *  
#   1 - 0.1 == 0        0.015328   0.028154   0.544  1.00000    
#   10 - 0.1 == 0       0.044939   0.027949   1.608  0.53932    
#   100 - 0.1 == 0     -0.090520   0.028923  -3.130  0.01400 *  
#   10 - 1 == 0         0.029611   0.027840   1.064  1.00000    
#   100 - 1 == 0       -0.105848   0.028818  -3.673  0.00216 ** 
#   100 - 10 == 0      -0.135459   0.028618  -4.733 2.21e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)



# ---- 6. Write tables ---------------------------------------------------------
## --- 6.1. Means and standard deviation ---------------------------------------
### -- 6.1.1. Surface ----------------------------------------------------------
mean_surface_growth <- surface %>% # table
  group_by(spec, treat, time) %>% # separates species, treatment and timepoints
  get_summary_stats(surface_growth) %>% # columns of interest: parameter
  select(-q1, -q3, -iqr, -mad, -se, -ci, -median) # because only selecting what wanted did not work: remove unnecessary columns
# write it into .csv
write_csv2(mean_surface_growth, "out/mean_surface_growth.csv")

### -- 6.1.2. Volume -----------------------------------------------------------
mean_volume_growth <- volume %>% # table
  group_by(spec, treat, time) %>% # separates species, treatment and timepoints
  get_summary_stats(volume_growth) %>% # columns of interest: parameter
  select(-q1, -q3, -iqr, -mad, -se, -ci, -median) # because only selecting what wanted did not work: remove unnecessary columns
# write it into .csv
write_csv2(mean_volume_growth, "out/mean_volume_growth.csv")

### -- 6.1.3. Calcification ----------------------------------------------------
mean_weight_growth <- calcification %>% # table
  group_by(spec, treat, time) %>% # separates species, treatment and timepoints
  get_summary_stats(weight_growth) %>% # columns of interest: parameter
  select(-q1, -q3, -iqr, -mad, -se, -ci, -median) # because only selecting what wanted did not work: remove unnecessary columns
# write it into .csv
write_csv2(mean_weight_growth, "out/mean_weight_growth.csv")
