# ----- 1. Explanation of this script ----
# This script focuses on the assessment of coral necrosis.
# Occurence of necrosis was assessed visually as yes/no for every sampling timepoint
# Relative necrosis was measured using 3D Scanning
# This script builds up on data tables produced in the script 'Growth_data_processing'

# ----- 2. Load in needed packages -------
# to easily clean data, to read in .rds files 
library(tidyverse)

# to easily read in all data files
library(readxl)

# to conduct the exact fisher test function fisher_test
library(rstatix)
library(dplyr)

# for statistical testing
library(multcomp)
library(lme4)

# check model fits visually using qqplot
library(car)

# ----- 3. Read in needed data files ----- 
# Occurence of necrosis for each coral and each sampling timepoint
# table originally in wide format 
# needs to be in long format to process with additional information
Necro_occurence <- read.csv2("in/necrosis.csv") %>%
  # rename columns to get continuous timepoints for later analyses
  rename("0" = "necro_t0") %>%
  rename("1" = "necro_t1") %>% 
  rename("2" = "necro_t2") %>% 
  rename("3" = "necro_t3") %>%  
  # bring table into long format
  pivot_longer(cols = c('0', '1', '2', '3'), # use previous columns as new entries - categories
               names_to = 'tp',  # assign new name for the column of tp
               values_to = 'necro_occured') # assign new name for the values previously in the column under the headers above

# Percent of difference between living tissue area and surface area
# to bring together with the observed occurrences
Necro_percent <- read_rds("processed/necrosis_percent.rds")


# ---- 4. Prepare data for statistical analyzes ----
Necro_occurence <- Necro_occurence %>%
  # create a new column to get a clean merge
  unite(ID_tp, c(ID, tp), sep="_", remove=FALSE) %>%
  # leave out some columns to have a cleaner merge with the percent table
  dplyr::select(-ID, -spec, -col, -origin, -tank, -tp, -treat, -conc)
Necro_percent <- Necro_percent %>%
  # create a new column to get a clean merge
  unite(ID_tp, c(ID, tp), sep="_", remove=FALSE) %>%
  # leave out some columns to have a cleaner merge with the percent table
  dplyr::select(-ID, -spec, -col, -tank, -tp)

# merge into one table
Necrosis <-   merge(Necro_percent, Necro_occurence, by = 'ID_tp', all.x = TRUE) %>%
  separate(ID_tp, c('spec', 'col', 'tank', 'tp')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) %>%
  # make timepoint column numeric for further analyses
  mutate(tp = as.numeric(tp)) %>% 
  # make concentration column numeric for further analyses
  mutate(conc = as.numeric(conc))


# create additional column with renamed treatment, in case glmer not possible with treat column
Necrosis <- Necrosis %>% 
  mutate(treat = as.factor(treat),
         treatment = case_when(treat == "control" ~ "control",
                               treat == 0.1 ~ "A",
                               treat == 1 ~ "B",
                               treat == 10 ~ "C",
                               TRUE ~ "D"))

# create a new column (necro_per) for necrosis corrected for occurence
Necrosis <- Necrosis %>%
  mutate(necro_per = necrosis,
         # if no necrosis was observed enter '0'
         necro_per = case_when(necro_occured == 'no' ~ "0"),
         # format new column as numeric
         necro_per = as.numeric(necro_per))

# replace NAs of positive necrotic occurence with the relative necrotic surface area from 'necrosis'
Necrosis$necro_per <- ifelse(is.na(Necrosis$necro_per),
                             Necrosis$necrosis, Necrosis$necro_per)

# delete old, now in "necro_per" corrected column of necrosis "necrosis" to avoid confusion
Necrosis <- Necrosis %>%
  dplyr::select(-necrosis)

# check levels of certain columns to evaluate for releveling 
# use treatment as categories not numbers (if numbers necessary: use "conc" column)
Necrosis$treat <- factor(Necrosis$treat, 
                         levels = c("control", "0.1", "1", "10", "100"))
# levels(Necrosis$treat)



# ---- 5. Statistical analyses -------------------------------------------------
# Both Species will be assessed separately
# For each timepoint statistical analyzes will be conducted separately, 
# comparing the frequency of occurence of the used categories between the 5 diffferent treatments using LMER and GLMER

# --> Direct comparison
# For t0 in both species for each treatment 100 % of the corals showed no necrosis (see 'Summary_necrosis')
# Therefore for t0 no statistical analyses needed

### -- 5.1.1. Pocillopora verrucosa --------------------------------------------

#### - 5.1.1.1. Overall effect -------------------------------------------------
# create a subset with data of Pve for t3
# to see overall effect (concentration-dependency) at the end of the experiment
# for overall effect at 12 weeks of exposure: Pve_necrosis <- subset(Necrosis, spec == "Pve" & tp == "3")
Pve_necrosis <- subset(Necrosis, spec == "Pve")

# LMER didn't show a good fit - GLMER is used
model_Pve <- glmer((necro_per) ~ conc + (1|col) + (1|tp), family = "poisson", data = Pve_necrosis)
# summary of tested with the GLMER differences
cftest(model_Pve)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (necro_per) ~ conc + (1 | col), data = Pve_necrosis, 
# family = "poisson")
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0 -2.467888   0.479712  -5.145 2.68e-07 ***
# conc == 0         0.003524   0.001353   2.604  0.00921 ** 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Univariate p values reported)


#### - 5.1.2.1. Specific effect ------------------------------------------------
# t1 - supplementary data
Pve_t1 <- subset(Necrosis, spec == "Pve" & tp == "1")
# LMER didn't show a good fit - GLMER is used
model_Pve_t1 <- glmer((necro_per) ~ treat + (1|col), family = "poisson", data = Pve_t1)
summary(glht(model_Pve_t1, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (necro_per) ~ treat + (1 | col), data = Pve_t1, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0  5.602e-07  2.832e+03   0.000   1.0000  
# 1 - control == 0    1.805e+01  2.003e+03   0.009   1.0000  
# 10 - control == 0   1.884e+01  2.003e+03   0.009   1.0000  
# 100 - control == 0  1.832e+01  2.003e+03   0.009   1.0000  
# 1 - 0.1 == 0        1.805e+01  2.003e+03   0.009   1.0000  
# 10 - 0.1 == 0       1.884e+01  2.003e+03   0.009   1.0000  
# 100 - 0.1 == 0      1.832e+01  2.003e+03   0.009   1.0000  
# 10 - 1 == 0         7.870e-01  2.906e-01   2.708   0.0676 .
# 100 - 1 == 0        2.694e-01  3.199e-01   0.842   1.0000  
# 100 - 10 == 0      -5.176e-01  2.659e-01  -1.946   0.4647  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# t2 - supplementary data
Pve_t2 <- subset(Necrosis, spec == "Pve" & tp == "2")
# LMER didn't show a good fit - GLMER is used
model_Pve_t2 <- glmer((necro_per) ~ treat + (1|col), family = "poisson", data = Pve_t2)
summary(glht(model_Pve_t2, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (necro_per) ~ treat + (1 | col), data = Pve_t2, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0  -1.1978     0.5147  -2.327 0.079850 .  
# 1 - control == 0    -0.1767     0.3671  -0.481 0.630281    
# 10 - control == 0    0.8145     0.2977   2.736 0.037321 *  
# 100 - control == 0  -2.0561     0.7359  -2.794 0.036445 *  
# 1 - 0.1 == 0         1.0211     0.5261   1.941 0.156847    
# 10 - 0.1 == 0        2.0124     0.4803   4.189 0.000280 ***
# 100 - 0.1 == 0      -0.8583     0.8268  -1.038 0.598481    
# 10 - 1 == 0          0.9912     0.3170   3.127 0.014148 *  
# 100 - 1 == 0        -1.8794     0.7439  -2.526 0.057627 .  
# 100 - 10 == 0       -2.8707     0.7123  -4.030 0.000501 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# t3 - main part
Pve_t3 <- subset(Necrosis, spec == "Pve" & tp == "3")
# LMER didn't show a good fit - GLMER is used
model_Pve_t3 <- glmer((necro_per) ~ treat + (1|col), family = "poisson", data = Pve_t3)
summary(glht(model_Pve_t3, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (necro_per) ~ treatment + (1 | col), data = Pve_t3, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0  -1.1651     0.3273  -3.560 0.002966 ** 
# 1 - control == 0    -1.0361     0.3118  -3.323 0.006237 ** 
# 10 - control == 0   -0.1645     0.2355  -0.698 1.000000    
# 100 - control == 0   0.2002     0.2152   0.930 1.000000    
# 1 - 0.1 == 0         0.1290     0.3917   0.329 1.000000    
# 10 - 0.1 == 0        1.0006     0.3342   2.994 0.016502 *  
# 100 - 0.1 == 0       1.3653     0.3201   4.265 0.000200 ***
# 10 - 1 == 0          0.8716     0.3190   2.732 0.031479 *  
# 100 - 1 == 0         1.2363     0.3043   4.063 0.000437 ***
# 100 - 10 == 0        0.3648     0.2255   1.617 0.423183
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


### -- 5.1.2. Stylophora pistillata --------------------------------------------

#### - 5.1.2.1. Overall effect -------------------------------------------------
# create a subset with data of Spi for t3
# to see overall effect (concentration-dependency) at the end of the experiment
# for overall effect after 12 weeks exposure: Spi_necrosis <- subset(Necrosis, spec == "Spi" & tp == "3")
Spi_necrosis <- subset(Necrosis, spec == "Spi")

# LMER didn't show a good fit - GLMER is used
model_Spi <- glmer((necro_per) ~ conc + (1|col) + (1|time), family = "poisson", data = Spi_necrosis)
# summary of tested with the GLMER differences
cftest(model_Spi)
# OUTPUT:
# Simultaneous Tests for General Linear Hypotheses
# Fit: glmer(formula = (necro_per) ~ conc + (1 | col), data = Spi_necrosis, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) == 0 -1.2744503  0.4367526  -2.918  0.00352 ** 
# conc == 0         0.0097649  0.0007041  13.868  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#   (Univariate p values reported)


#### - 5.1.2.1. Specific effect ------------------------------------------------
# t1 - supplementary data
Spi_t1 <- subset(Necrosis, spec == "Spi" & tp == "1")
# LMER didn't show a good fit - GLMER is used
model_Spi_t1 <- glmer((necro_per) ~ treat + (1|col), family = "poisson", data = Spi_t1)
summary(glht(model_Spi_t1, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (necro_per) ~ treat + (1 | col), data = Spi_t1, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0  -17.42693 1265.94855  -0.014   1.0000    
# 1 - control == 0     -0.01803    0.29559  -0.061   1.0000    
# 10 - control == 0     0.63040    0.25757   2.448   0.1007    
# 100 - control == 0    1.09208    0.24045   4.542 5.02e-05 ***
# 1 - 0.1 == 0         17.40890 1265.94855   0.014   1.0000    
# 10 - 0.1 == 0        18.05733 1265.94854   0.014   1.0000    
# 100 - 0.1 == 0       18.51901 1265.94854   0.015   1.0000    
# 10 - 1 == 0           0.64843    0.25909   2.503   0.0986 .  
# 100 - 1 == 0          1.11011    0.24209   4.586 4.53e-05 ***
# 100 - 10 == 0         0.46168    0.19384   2.382   0.1034    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)

# t2 - supplementary data
Spi_t2 <- subset(Necrosis, spec == "Spi" & tp == "2")
# LMER didn't show a good fit - GLMER is used
model_Spi_t2 <- glmer((necro_per) ~ treat + (1|col), family = "poisson", data = Spi_t2)
summary(glht(model_Spi_t2, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT:
#Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (necro_per) ~ treat + (1 | col), data = Spi_t2, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0  0.3149242  0.1855084   1.698  0.26810    
# 1 - control == 0    0.3151177  0.1855008   1.699  0.26810    
# 10 - control == 0  -0.4681168  0.2272976  -2.059  0.15779    
# 100 - control == 0  0.8765001  0.1678526   5.222 1.59e-06 ***
# 1 - 0.1 == 0        0.0001936  0.1703998   0.001  0.99909    
# 10 - 0.1 == 0      -0.7830409  0.2151505  -3.640  0.00163 ** 
# 100 - 0.1 == 0      0.5615759  0.1509968   3.719  0.00160 ** 
# 10 - 1 == 0        -0.7832345  0.2151439  -3.641  0.00163 ** 
# 100 - 1 == 0        0.5613824  0.1509875   3.718  0.00160 ** 
# 100 - 10 == 0       1.3446169  0.2001271   6.719 1.83e-10 ***   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)


# t3 - main part
# LMER didn't show a good fit - GLMER is used
Spi_t3 <- subset(Necrosis, spec == "Spi" & tp == "3")
# LMER didn't show a good fit - GLMER is used
model_Spi_t3 <- glmer((necro_per) ~ treat + (1|col), family = "poisson", data = Spi_t3)
summary(glht(model_Spi_t3, linfct = mcp(treat = "Tukey")), 
        test = adjusted("holm"))
# OUTPUT: 
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: glmer(formula = (necro_per) ~ treatment + (1 | col), data = Spi_necrosis, 
#            family = "poisson")
# Linear Hypotheses:
#   Estimate Std. Error z value Pr(>|z|)    
# 0.1 - control == 0  0.58302    0.18586   3.137 0.003415 ** 
# 1 - control == 0    0.67631    0.18287   3.698 0.000798 ***
# 10 - control == 0  -1.11418    0.29953  -3.720 0.000798 ***
# 100 - control == 0  1.30373    0.16789   7.765 7.39e-14 ***
# 1 - 0.1 == 0        0.09329    0.15378   0.607 0.544074    
# 10 - 0.1 == 0      -1.69720    0.28271  -6.003 1.35e-08 ***
# 100 - 0.1 == 0      0.72071    0.13563   5.314 6.44e-07 ***
# 10 - 1 == 0        -1.79049    0.28076  -6.377 1.44e-09 ***
# 100 - 1 == 0        0.62742    0.13150   4.771 9.15e-06 ***
# 100 - 10 == 0       2.41791    0.27124   8.914  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- holm method)




# ---- 6. Write tables ---------------------------------------------------------
# ---- 6.1. Necrosis percent table ---------------------------------------------
# Percent necrosis for all
write.csv2(Necrosis, "out/necrosis_percent.csv")

# ---- 6.1. Categorize severity of necrosis ------------------------------------
# imporant for visualization
# according to Marshall and Schuttenberg 2006
# low: 1-10%
# moderate: > 10-50%
# high: > 50%
Necrosis_category <- Necrosis %>% 
  mutate(cat = case_when(necro_per >= 50 ~ "high",
                         necro_per >= 10 ~ "moderate",
                         necro_per >= 1 ~ "low",
                         TRUE ~ "none"))



## --- 6.2.  Summary of necrosis occurences ------------------------------------
# level the categories
Necrosis_category$cat <- factor(Necrosis_category$cat, 
                                levels = c( "none", "low", "moderate", "high"))

# create a table with percentages of the categories (cat)
# per Species (spec), treatment (treat), timepoint (tp)
Summary_necrosis <- Necrosis_category %>%
  freq_table(spec, treat, tp, cat, na.rm = T)

write.csv2(Summary_necrosis, "out/Summary_necrosis.csv")

