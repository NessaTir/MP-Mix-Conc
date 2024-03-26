# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the statistical analysis of the relationship of all processed data and the concentrations. 
# Relationship will be tested for beeing
# 1) linear
# 2) exponential

# This script is closely connected to 'Plots - 4.5. Correlation Plot'

# ----- 2. Load in needed packages ---------------------------------------------
library(tidyverse)

# read in the data files
library(readxl)

# work with summarized data, e.g., means
library(rstatix)

# test for linear correlation
library(ggpubr)


# ----- 3. Read in needed data files ------------------------------------------- 
## ---- 3.01. Surface ----------------------------------------------------------
surface <- read.csv2("processed/surface_all.csv")

## ---- 3.02. Volume -----------------------------------------------------------
volume <- read.csv2("processed/volume_all.csv")

## ---- 3.03. Calcification ----------------------------------------------------
calcification <- read.csv2("processed/calcification_all.csv")

## ---- 3.04. Necrosis ---------------------------------------------------------
necrosis <- read.csv2("processed/necrosis_all.csv")

## ---- 3.05. Polypactivity ----------------------------------------------------
polypactivity <- read.csv2("processed/polypactivity_all.csv")

## ---- 3.06. YII --------------------------------------------------------------
YII <- read.csv2("processed/YII_all.csv")

## ---- 3.07. FvFm -------------------------------------------------------------
FvFm <- read.csv2("processed/FvFm_all.csv")

## ---- 3.08. rETRmax ----------------------------------------------------------
rETR <- read.csv2("processed/rETR_all.csv")

## ---- 3.09. Ek ---------------------------------------------------------------
Ek <- read.csv2("processed/Ek_all.csv")

## ---- 3.10. Alpha ------------------------------------------------------------
alpha <- read.csv2("processed/alpha_all.csv")


# ----- 4. Linear Relationship -------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(surface_Pve$value_log, surface_Pve$conc, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve$value_log and surface_Pve$conc
# t = -3.2999, df = 88, p-value = 0.001398
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5042723 -0.1339487
# sample estimates:
#  cor 
# -0.3318348 

# Surface Spi
surface_Spi <- subset(surface_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(surface_Spi$value_log, surface_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Spi$value_log and surface_Spi$conc
# t = -1.0995, df = 80, p-value = 0.2749
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.33026771  0.09758609
# sample estimates:
#   cor 
# -0.1220047 


# Volume Pve
volume_Pve <- subset(volume_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(volume_Pve$value_log, volume_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Pve$value_log and volume_Pve$conc
# t = -3.2774, df = 85, p-value = 0.001518
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5096395 -0.1337375
# sample estimates:
#   cor 
# -0.3349497 

# Volume Spi
volume_Spi <- subset(volume_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(volume_Spi$value_log, volume_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Spi$value_log and volume_Spi$conc
# t = -1.6947, df = 86, p-value = 0.09375
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.37508857  0.03083436
# sample estimates:
#   cor 
# -0.1797688 

# calcification Pve
calcification_Pve <- subset(calcification_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(calcification_Pve$value_log, calcification_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Pve$value_log and calcification_Pve$conc
# t = -2.1919, df = 88, p-value = 0.03102
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.41506538 -0.02145178
# sample estimates:
#   cor 
# -0.2275321 

# calcification Spi
calcification_Spi <- subset(calcification_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(calcification_Spi$value_log, calcification_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Spi$value_log and calcification_Spi$conc
# t = -2.4805, df = 86, p-value = 0.01507
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.44381813 -0.05175334
# sample estimates:
#   cor 
# -0.2583948  


# necrosis Pve
necrosis_Pve <- subset(necrosis_all, spec == "Pve" &
                         # remove 0 otherwise log will transform values to unusable
                         value != "0") %>%
  mutate(value_log = log(value))
cor.test(necrosis_Pve$value_log, necrosis_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Pve$value and necrosis_Pve$conc
# t = -0.46259, df = 8, p-value = 0.656
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.7180594  0.5211858
# sample estimates:
#   cor 
# -0.1614063  

# necrosis Spi
necrosis_Spi <- subset(necrosis_all, spec == "Spi" &
     # remove 0 otherwise log will transform values to unusable
                         value != "0") %>%
  mutate(value_log = log(value))
cor.test(necrosis_Spi$value_log, necrosis_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# t = 0.39795, df = 15, p-value = 0.6963
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3979845  0.5555635
# sample estimates:
#   cor 
# 0.1022123 


# polypactivity Pve
polypactivity_Pve <- subset(polypactivity_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(polypactivity_Pve$value_log, polypactivity_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Pve$value_log and polypactivity_Pve$conc
# t = -6.3813, df = 88, p-value = 7.963e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6892596 -0.4022107
# sample estimates:
#   cor 
# -0.5624528  

# polypactivity Spi
polypactivity_Spi <- subset(polypactivity_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(polypactivity_Spi$value_log, polypactivity_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Spi$value_log and polypactivity_Spi$conc
# t = -3.1088, df = 88, p-value = 0.002531
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4897595 -0.1149729
# sample estimates:
#   cor 
# -0.314574 


# YII Pve
YII_Pve <- subset(YII_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(YII_Pve$value_log, YII_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Pve$value_log and YII_Pve$conc
# t = 2.1449, df = 88, p-value = 0.03472
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.01656603 0.41101168
# sample estimates:
#   cor 
# 0.2228925  

# YII Spi
YII_Spi <- subset(YII_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(YII_Spi$value_log, YII_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Spi$value_log and YII_Spi$conc
# t = 1.0509, df = 88, p-value = 0.2962
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.09802393  0.31124258
# sample estimates:
#   cor 
# 0.1113271 


# FvFm Pve
FvFm_Pve <- subset(FvFm_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(FvFm_Pve$value_log, FvFm_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Pve$value_log and FvFm_Pve$conc
# t = 2.6152, df = 88, p-value = 0.01049
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.0650695 0.4505754
# sample estimates:
#   cor 
# 0.2685419

# FvFm Spi
FvFm_Spi <- subset(FvFm_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(FvFm_Spi$value_log, FvFm_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Spi$value_log and FvFm_Spi$conc
# t = -1.7347, df = 88, p-value = 0.08629
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.37481445  0.02624158
# sample estimates:
#   cor 
# -0.1818377  


# rETRmax Pve
rETR_Pve <- subset(rETR_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(rETR_Pve$value_log, rETR_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Pve$value_log and rETR_Pve$conc
# t = -0.14923, df = 88, p-value = 0.8817
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2222648  0.1918172
# sample estimates:
#   cor 
# -0.01590576 

# rETRmax Spi
rETR_Spi <- subset(rETR_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(rETR_Spi$value_log, rETR_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Spi$value_log and rETR_Spi$conc
# t = -0.28161, df = 88, p-value = 0.7789
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.2356329  0.1781924
# sample estimates:
# cor 
# -0.03000602 


# Ek Pve
Ek_Pve <- subset(Ek_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(Ek_Pve$value_log, Ek_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Pve$value and Ek_Pve$conc
# t = -0.31022, df = 88, p-value = 0.7571
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2385099  0.1752393
# sample estimates:
#   cor 
# -0.03305129 

# Ek Spi
Ek_Spi <- subset(Ek_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(Ek_Spi$value_log, Ek_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Spi$value and Ek_Spi$conc
# t = -0.066881, df = 88, p-value = 0.9468
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2139046  0.2002574
# sample estimates:
#   cor 
# -0.007129341 


# alpha Pve
alpha_Pve <- subset(alpha_all, spec == "Pve") %>%
  mutate(value_log = log(value))
cor.test(alpha_Pve$value_log, alpha_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Pve$value_log and alpha_Pve$conc
# t = 0.52306, df = 88, p-value = 0.6022
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1531848  0.2597686
# sample estimates:
#   cor 
# 0.05567239 

# alpha Spi
alpha_Spi <- subset(alpha_all, spec == "Spi") %>%
  mutate(value_log = log(value))
cor.test(alpha_Spi$value_log, alpha_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Spi$value_log and alpha_Spi$conc
# t = -0.51993, df = 88, p-value = 0.6044
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2594578  0.1535102
# sample estimates:
#   cor 
# -0.05534019 


# ----- 5. Other Relationships -------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface_all, spec == "Pve")

# Exponential
ks.test(surface_Pve$value, surface_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  surface_Pve$value
# D = 0.93825, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation

# In general?
cor.test(surface_Pve$value_exp, surface_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: Spearman's rank correlation rho
# data:  surface_Pve$value and surface_Pve$conc
# S = 151098, p-value = 0.0206
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.2437544 
# correlation


# Surface Spi
surface_Spi <- subset(surface_all, spec == "Spi")

# Exponential
ks.test(surface_Spi$value, surface_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  surface_Spi$value
# D = 0.96667, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation

# In general?
cor.test(surface_Spi$value, surface_Spi$conc, method = "spearman", exact=FALSE)
# OUTPUT: Spearman's rank correlation rho
# data:  surface_Spi$value and surface_Spi$conc
# S = 145807, p-value = 0.0585
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#  rho 
# -0.2002052 
# no correlation


# Volume Pve
volume_Pve <- subset(volume_all, spec == "Pve")

# Exponential
ks.test(volume_Pve$value, volume_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  volume_Pve$value
# D = 0.96685, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation

# In general?
cor.test(surface_Pve$value, surface_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: Spearman's rank correlation rho
# data:  surface_Pve$value and surface_Pve$conc
# S = 151098, p-value = 0.0206
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.2437544 
# correlation


# Volume Spi
volume_Spi <- subset(volume_all, spec == "Spi")

# Exponential
ks.test(volume_Spi$value, volume_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  volume_Spi$value
# D = 0.96667, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(volume_Spi$value, volume_Spi$conc,  method = "spearman", exact=FALSE)
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Spi$value and volume_Spi$conc
# Spearman's rank correlation rho
# data:  surface_Spi$value and surface_Spi$conc
# S = 145807, p-value = 0.0585
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.2002052 
# no correlation


# calcification Pve
calcification_Pve <- subset(calcification_all, spec == "Pve")

# Exponential
ks.test(calcification_Pve$value, calcification_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  calcification_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(calcification_Pve$value, calcification_Pve$conc, method = "spearman", exact=FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  calcification_Pve$value and calcification_Pve$conc
# S = 131552, p-value = 0.4375
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.0828644 
# no correlation


# calcification Spi
calcification_Spi <- subset(calcification_all, spec == "Spi")

# Exponential
ks.test(calcification_Spi$value, calcification_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  calcification_Spi$value
# D = 0.96667, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(calcification_Spi$value, calcification_Spi$conc, method = "spearman", exact=FALSE)
# OUTPUT: Spearman's rank correlation rho
# data:  calcification_Spi$value and calcification_Spi$conc
# S = 147938, p-value = 0.03924
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.2177459 
# correlation


# necrosis Pve
necrosis_Pve <- subset(necrosis_all, spec == "Pve")

# Exponential
ks.test(necrosis_Pve$value, necrosis_Pve$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  necrosis_Pve$value
# D = 0.96667, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(necrosis_Pve$value, necrosis_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  necrosis_Pve$value and necrosis_Pve$conc
# S = 116098, p-value = 0.6781
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.04434294 
# no correlation


# necrosis Spi
necrosis_Spi <- subset(necrosis_all, spec == "Spi")

# Exponential
ks.test(necrosis_Spi$value, necrosis_Spi$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  necrosis_Spi$value
# D = 0.96667, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(necrosis_Spi$value, necrosis_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  necrosis_Spi$value and necrosis_Spi$conc
# S = 95824, p-value = 0.04566
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2112297 
# correlation


# polypactivity Pve
polypactivity_Pve <- subset(polypactivity_all, spec == "Pve")

# Exponential
ks.test(polypactivity_Pve$value, polypactivity_Pve$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  polypactivity_Spi$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(polypactivity_Pve$value, polypactivity_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  polypactivity_Pve$value and polypactivity_Pve$conc
# S = 185283, p-value = 1.073e-07
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.5251547 
# correlation


# polypactivity Spi
polypactivity_Spi <- subset(polypactivity_all, spec == "Spi")

# Exponential
ks.test(polypactivity_Spi$value, polypactivity_Spi$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  polypactivity_Spi$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(polypactivity_Spi$value, polypactivity_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  polypactivity_Spi$value and polypactivity_Spi$conc
# S = 149155, p-value = 0.03085
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.2277612 
# correlation


# YII Pve
YII_Pve <- subset(YII_all, spec == "Pve")

# Exponential
ks.test(YII_Pve$value, YII_Pve$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  YII_Pve$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(YII_Pve$value, YII_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  YII_Pve$value and YII_Pve$conc
# S = 78572, p-value = 0.0006377
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.3532365 
# correlation


# YII Spi
YII_Spi <- subset(YII_all, spec == "Spi")

# Exponential
ks.test(YII_Spi$value, YII_Spi$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  YII_Spi$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(YII_Spi$value, YII_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Pearson's product-moment correlation
# data:  YIISpearman's rank correlation rho
# data:  YII_Spi$value and YII_Spi$conc
# S = 96979, p-value = 0.05657
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2017207 
# no correlation


# FvFm Pve
FvFm_Pve <- subset(FvFm_all, spec == "Pve")

# Exponential
ks.test(FvFm_Pve$value, FvFm_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  FvFm_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(FvFm_Pve$value, FvFm_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: Spearman's rank correlation rho
# data:  FvFm_Pve$value and FvFm_Pve$conc
# S = 98045, p-value = 0.06845
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.192947 
# no correlation


# FvFm Spi
FvFm_Spi <- subset(FvFm_all, spec == "Spi")

# Exponential
ks.test(FvFm_Spi$value, FvFm_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  FvFm_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(FvFm_Spi$value, FvFm_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  FvFm_Spi$value and FvFm_Spi$conc
# S = 127106, p-value = 0.665
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.046271 
# no correlation


# rETRmax Pve
rETR_Pve <- subset(rETR_all, spec == "Pve")

# Exponential
ks.test(rETR_Pve$value, rETR_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  rETR_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(rETR_Pve$value, rETR_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  rETR_Pve$value and rETR_Pve$conc
# S = 127731, p-value = 0.6303
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.05141222 
# no correlation


# rETRmax Spi
rETR_Spi <- subset(rETR_all, spec == "Spi")

# Exponential
ks.test(rETR_Spi$value, rETR_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  rETR_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(rETR_Spi$value, rETR_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  rETR_Spi$value and rETR_Spi$conc
# S = 127033, p-value = 0.6691
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.04566615 
# no correlation


# Ek Pve
Ek_Pve <- subset(Ek_all, spec == "Pve")

# Exponential
ks.test(Ek_Pve$value, Ek_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  Ek_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(Ek_Pve$value, Ek_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  Ek_Pve$value and Ek_Pve$conc
# S = 127584, p-value = 0.6384
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.05020252 
# no correlation


# Ek Spi
Ek_Spi <- subset(Ek_all, spec == "Spi")

# Exponential
ks.test(Ek_Spi$value, Ek_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  Ek_Spi$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(Ek_Spi$value, Ek_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  Ek_Spi$value and Ek_Spi$conc
# S = 132764, p-value = 0.3841
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.09284442 
# no correlation


# alpha Pve
alpha_Pve <- subset(alpha_all, spec == "Pve")

# Exponential
ks.test(alpha_Pve$value, alpha_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  alpha_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(alpha_Pve$value, alpha_Pve$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  alpha_Pve$value and alpha_Pve$conc
# S = 119281, p-value = 0.8652
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.01814549
# no correlation


# alpha Spi
alpha_Spi <- subset(alpha_all, spec == "Spi")

# Exponential
ks.test(alpha_Spi$value, alpha_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  alpha_Spie$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# In general?
cor.test(alpha_Spi$value, alpha_Spi$conc, method = "spearman", exact = FALSE)
# OUTPUT: 	Spearman's rank correlation rho
# data:  alpha_Spi$value and alpha_Spi$conc
# S = 124902, p-value = 0.7924
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.02812551 
# no correlation
