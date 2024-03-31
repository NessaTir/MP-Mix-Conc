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

# for GAM
library(mgcv)


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
## ---- 4.01. Surface ----------------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(surface_Pve$value, surface_Pve$conc, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve$value_log and surface_Pve$conc
# t = -3.2999, df = 88, p-value = 0.001398
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5042723 -0.1339487
# sample estimates:
#  cor 
# -0.3318348 

# log-log transformed
# exclude conc 0 mg/l, as log(0) not possible
surface_Pve_wocon <- surface_Pve %>%
  subset(conc != "0")
cor.test(surface_Pve_wocon$value_log, surface_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve_wocon$value_log and surface_Pve_wocon$conc_log
# t = -2.5401, df = 67, p-value = 0.01341
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.49811574 -0.06419875
# sample estimates:
#   cor 
# -0.2963761 

# log-log transformed - early linearity
surface_Pve_wocon100 <- surface_Pve_wocon %>%
  subset(conc != "100")
cor.test(surface_Pve_wocon100$value_log, surface_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve_wocon100$value_log and surface_Pve_wocon100$conc_log
# t = 0.1374, df = 51, p-value = 0.8913
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2523701  0.2880321
# sample estimates:
#   cor 
# 0.01923585 

# log-log transformed - late linearity
surface_Pve_wocon01 <- surface_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(surface_Pve_wocon01$value_log, surface_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve_wocon01$value_log and surface_Pve_wocon01$conc_log
# t = -2.021, df = 49, p-value = 0.04876
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.513699556 -0.001949801
# sample estimates:
#   cor 
# -0.2773844


# Surface Spi
surface_Spi <- subset(surface, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(surface_Spi$value, surface_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Spi$value and surface_Spi$conc
# t = -2.5649, df = 88, p-value = 0.01201
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4464486 -0.0599237
# sample estimates:
#   cor 
# -0.2637419  

# log-log transformed
surface_Spi_wocon <- surface_Spi %>%
  subset(conc != "0")
cor.test(surface_Spi_wocon$value_log, surface_Spi_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Spi_wocon$value_log and surface_Spi_wocon$conc_log
# t = -2.1079, df = 63, p-value = 0.03902
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.47108258 -0.01362858
# sample estimates:
#   cor 
# -0.2566745  

# log-log transformed - early linearity
surface_Spi_wocon100 <- surface_Spi_wocon %>%
  subset(conc != "100")
cor.test(surface_Spi_wocon100$value_log, surface_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Spi_wocon100$value_log and surface_Spi_wocon100$conc_log
# t = -1.3676, df = 48, p-value = 0.1778
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4478609  0.0895187
# sample estimates:
#   cor 
# -0.1936544  

# log-log transformed - late linearity
surface_Spi_wocon01 <- surface_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(surface_Spi_wocon01$value_log, surface_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Spi_wocon01$value_log and surface_Spi_wocon01$conc_log
# t = -1.5131, df = 47, p-value = 0.137
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4683343  0.0699131
# sample estimates:
#   cor 
# -0.2155205 


## ---- 4.02. Volume ----------------------------------------------------------
# Volume Pve
volume_Pve <- subset(volume, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(volume_Pve$value, volume_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Pve$value and volume_Pve$conc
# t = -3.6324, df = 88, p-value = 0.0004715
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5286500 -0.1664462
# sample estimates:
#   cor 
# -0.3610907 

# log-log transformed
volume_Pve_wocon <- volume_Pve %>%
  subset(conc != "0")
cor.test(volume_Pve_wocon$value_log, volume_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Pve_wocon$value_log and volume_Pve_wocon$conc_log
# t = -2.2445, df = 68, p-value = 0.02805
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.46868590 -0.02948102
# sample estimates:
#   cor 
# -0.2626357

# log-log transformed - early linearity
volume_Pve_wocon100 <- volume_Pve_wocon %>%
  subset(conc != "100")
cor.test(volume_Pve_wocon100$value_log, volume_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Pve_wocon100$value_log and volume_Pve_wocon100$conc_log
# t = 0.29054, df = 52, p-value = 0.7726
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2299815  0.3047339
# sample estimates:
#   cor 
# 0.0402582   

# log-log transformed - late linearity
volume_Pve_wocon01 <- volume_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(volume_Pve_wocon01$value_log, volume_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Pve_wocon01$value_log and volume_Pve_wocon01$conc_log
# t = -2.1636, df = 50, p-value = 0.0353
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5236761 -0.0213985
# sample estimates:
#   cor 
# -0.2925902  


# Volume Spi
volume_Spi <- subset(volume, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(volume_Spi$value, volume_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Spi$value and volume_Spi$conc
# t = -2.3846, df = 88, p-value = 0.01925
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.43144116 -0.04138206
# sample estimates:
#   cor 
# -0.2463619 

# log-log transformed
volume_Spi_wocon <- volume_Spi %>%
  subset(conc != "0")
cor.test(volume_Spi_wocon$value_log, volume_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Spi_wocon$value_log and volume_Spi_wocon$conc_log
# t = -1.2413, df = 68, p-value = 0.2188
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.37085306  0.08924552
# sample estimates:
#   cor 
# -0.1488499 

# log-log transformed - early linearity
volume_Spi_wocon100 <- volume_Spi_wocon %>%
  subset(conc != "100")
cor.test(volume_Spi_wocon100$value_log, volume_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Spi_wocon100$value_log and volume_Spi_wocon100$conc_log
# t = 0.70887, df = 51, p-value = 0.4816
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1762228  0.3594723
# sample estimates:
#   cor 
# 0.0987759 

# log-log transformed - late linearity
volume_Spi_wocon01 <- volume_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(volume_Spi_wocon01$value_log, volume_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  volume_Spi_wocon01$value_log and volume_Spi_wocon01$conc_log
# t = -1.9589, df = 50, p-value = 0.05571
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.503216676  0.006385334
# sample estimates:
#   cor 
# -0.2669802


## ---- 4.03. Calcification ----------------------------------------------------
# calcification Pve
calcification_Pve <- subset(calcification, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(calcification_Pve$value, calcification_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Pve$value and calcification_Pve$conc
# t = -2.0478, df = 88, p-value = 0.04356
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.402581624 -0.006465604
# sample estimates:
#   cor 
# -0.2132712 

# log-log transformed
calcification_Pve_wocon <- calcification_Pve %>%
  subset(conc != "0")
cor.test(calcification_Pve_wocon$value_log, calcification_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Pve_wocon$value_log and calcification_Pve_wocon$conc_log
# t = -1.7739, df = 70, p-value = 0.08042
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.41895155  0.02547866
# sample estimates:
#   cor 
# -0.2074142

# log-log transformed - early linearity
calcification_Pve_wocon100 <- calcification_Pve_wocon %>%
  subset(conc != "100")
cor.test(calcification_Pve_wocon100$value_log, calcification_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Pve_wocon100$value_log and calcification_Pve_wocon100$conc_log
# t = 0.24975, df = 52, p-value = 0.8038
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2353286  0.2995969
# sample estimates:
#   cor 
# 0.03461298  

# log-log transformed - late linearity
calcification_Pve_wocon01 <- calcification_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(calcification_Pve_wocon01$value_log, calcification_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Pve_wocon01$value_log and calcification_Pve_wocon01$conc_log
# t = -1.8981, df = 52, p-value = 0.06324
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.48898468  0.01417368
# sample estimates:
#   cor 
# -0.2545529  


# calcification Spi
calcification_Spi <- subset(calcification, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(calcification_Spi$value, calcification_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Spi$value and calcification_Spi$conc
# t = -2.3848, df = 88, p-value = 0.01923
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.43146410 -0.04141021
# sample estimates:
#   cor 
# -0.2463884 

# log-log transformed
calcification_Spi_wocon <- calcification_Spi %>%
  subset(conc != "0")
cor.test(calcification_Spi_wocon$value_log, calcification_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Spi_wocon$value_log and calcification_Spi_wocon$conc_log
# t = -2.1373, df = 68, p-value = 0.03618
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.45881613 -0.01691368
# sample estimates:
#   cor 
# -0.2508906 

# log-log transformed - early linearity
calcification_Spi_wocon100 <- calcification_Spi_wocon %>%
  subset(conc != "100")
cor.test(calcification_Spi_wocon100$value_log, calcification_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Spi_wocon100$value_log and calcification_Spi_wocon100$conc_log
# t = -0.49953, df = 50, p-value = 0.6196
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3368907  0.2064005
# sample estimates:
#   cor 
# -0.07046911

# log-log transformed - late linearity
calcification_Spi_wocon01 <- calcification_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(calcification_Spi_wocon01$value_log, calcification_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  calcification_Spi_wocon01$value_log and calcification_Spi_wocon01$conc_log
# t = -1.2663, df = 52, p-value = 0.211
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.42121563  0.09940143
# sample estimates:
#   cor 
# -0.1729627 


## ---- 4.04. Necrosis ---------------------------------------------------------
# necrosis Pve
necrosis_Pve <- subset(necrosis, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(necrosis_Pve$value, necrosis_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Pve$value and necrosis_Pve$conc
# t = 0.79114, df = 88, p-value = 0.431
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1252331  0.2861487
# sample estimates:
#   cor 
# 0.08403752 

# log-log transformed
necrosis_Pve_wocon <- necrosis_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(necrosis_Pve_wocon$value_log, necrosis_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Pve_wocon$value_log and necrosis_Pve_wocon$conc_log
# t = 0.078869, df = 6, p-value = 0.9397
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6880956  0.7205150
# sample estimates:
#   cor 
# 0.0321816 

# log-log transformed - early linearity
necrosis_Pve_wocon100 <- necrosis_Pve_wocon %>%
  subset(conc != "100")
cor.test(necrosis_Pve_wocon100$value_log, necrosis_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Pve_wocon100$value_log and necrosis_Pve_wocon100$conc_log
# t = 0.24975, df = 52, p-value = 0.8038
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2353286  0.2995969
# sample estimates:
#   cor 
# 0.03461298  

# log-log transformed - late linearity
necrosis_Pve_wocon01 <- necrosis_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(necrosis_Pve_wocon01$value_log, necrosis_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Pve_wocon01$value_log and necrosis_Pve_wocon01$conc_log
# t = -1.8981, df = 52, p-value = 0.06324
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.48898468  0.01417368
# sample estimates:
#   cor 
# -0.2545529  


# necrosis Spi
necrosis_Spi <- subset(necrosis, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(necrosis_Spi$value, necrosis_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Spi$value and necrosis_Spi$conc
# t = -2.3848, df = 88, p-value = 0.01923
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.43146410 -0.04141021
# sample estimates:
#   cor 
# -0.2463884 

# log-log transformed
necrosis_Spi_wocon <- necrosis_Spi %>%
  subset(conc != "0")
cor.test(necrosis_Spi_wocon$value_log, necrosis_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Spi_wocon$value_log and necrosis_Spi_wocon$conc_log
# t = -2.1373, df = 68, p-value = 0.03618
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.45881613 -0.01691368
# sample estimates:
#   cor 
# -0.2508906 

# log-log transformed - early linearity
necrosis_Spi_wocon100 <- necrosis_Spi_wocon %>%
  subset(conc != "100")
cor.test(necrosis_Spi_wocon100$value_log, necrosis_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Spi_wocon100$value_log and necrosis_Spi_wocon100$conc_log
# t = -0.49953, df = 50, p-value = 0.6196
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3368907  0.2064005
# sample estimates:
#   cor 
# -0.07046911

# log-log transformed - late linearity
necrosis_Spi_wocon01 <- necrosis_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(necrosis_Spi_wocon01$value_log, necrosis_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  necrosis_Spi_wocon01$value_log and necrosis_Spi_wocon01$conc_log
# t = -1.2663, df = 52, p-value = 0.211
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.42121563  0.09940143
# sample estimates:
#   cor 
# -0.1729627 


## ---- 4.05. Polyp activity ---------------------------------------------------
# polypactivity Pve
polypactivity_Pve <- subset(polypactivity, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(polypactivity_Pve$value, polypactivity_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Pve$value and polypactivity_Pve$conc
# t = -6.5182, df = 88, p-value = 4.313e-09
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6955206 -0.4122429
# sample estimates:
#   cor 
# -0.5706192 

# log-log transformed
polypactivity_Pve_wocon <- polypactivity_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(polypactivity_Pve_wocon$value_log, polypactivity_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Pve_wocon$value_log and polypactivity_Pve_wocon$conc_log
# t = -5.5024, df = 70, p-value = 5.772e-07
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6929418 -0.3641733
# sample estimates:
#   cor 
# -0.5494833 

# log-log transformed - early linearity
polypactivity_Pve_wocon100 <- polypactivity_Pve_wocon %>%
  subset(conc != "100")
cor.test(polypactivity_Pve_wocon100$value_log, polypactivity_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Pve_wocon100$value_log and polypactivity_Pve_wocon100$conc_log
# t = -2.7713, df = 52, p-value = 0.007728
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.5715859 -0.1006337
# sample estimates:
#   cor 
# -0.358728 

# log-log transformed - late linearity
polypactivity_Pve_wocon01 <- polypactivity_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(polypactivity_Pve_wocon01$value_log, polypactivity_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Pve_wocon01$value_log and polypactivity_Pve_wocon01$conc_log
# t = -3.7539, df = 52, p-value = 0.0004404
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.6492409 -0.2213599
# sample estimates:
#   cor 
# -0.4617517  


# polypactivity Spi
polypactivity_Spi <- subset(polypactivity, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(polypactivity_Spi$value, polypactivity_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Spi$value and polypactivity_Spi$conc
# t = -2.4506, df = 88, p-value = 0.01624
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.43697008 -0.04818192
# sample estimates:
#   cor 
# -0.2527511 

# log-log transformed
polypactivity_Spi_wocon <- polypactivity_Spi %>%
  subset(conc != "0")
cor.test(polypactivity_Spi_wocon$value_log, polypactivity_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Spi_wocon$value_log and polypactivity_Spi_wocon$conc_log
# t = -2.583, df = 70, p-value = 0.01189
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4929697 -0.0679671
# sample estimates:
#   cor 
# -0.2949908 

# log-log transformed - early linearity
polypactivity_Spi_wocon100 <- polypactivity_Spi_wocon %>%
  subset(conc != "100")
cor.test(polypactivity_Spi_wocon100$value_log, polypactivity_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Spi_wocon100$value_log and polypactivity_Spi_wocon100$conc_log
# t = -1.0022, df = 52, p-value = 0.3209
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3910042  0.1350841
# sample estimates:
#   cor 
# -0.1376555 

# log-log transformed - late linearity
polypactivity_Spi_wocon01 <- polypactivity_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(polypactivity_Spi_wocon01$value_log, polypactivity_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  polypactivity_Spi_wocon01$value_log and polypactivity_Spi_wocon01$conc_log
# t = -2.2597, df = 52, p-value = 0.02806
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.52476767 -0.03398397
# sample estimates:
#   cor 
# -0.2990235 


## ---- 4.06. YII --------------------------------------------------------------
# YII Pve
YII_Pve <- subset(YII, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(YII_Pve$value, YII_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Pve$value and YII_Pve$conc
# t = 2.0914, df = 88, p-value = 0.03938
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.0110040 0.4063778
# sample estimates:
#   cor 
# 0.2175992  

# log-log transformed
YII_Pve_wocon <- YII_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(YII_Pve_wocon$value_log, YII_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Pve_wocon$value_log and YII_Pve_wocon$conc_log
# t = 3.4942, df = 70, p-value = 0.0008289
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1687800 0.5664734
# sample estimates:
#   cor 
# 0.3853798 

# log-log transformed - early linearity
YII_Pve_wocon100 <- YII_Pve_wocon %>%
  subset(conc != "100")
cor.test(YII_Pve_wocon100$value_log, YII_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Pve_wocon100$value_log and YII_Pve_wocon100$conc_log
# t = 3.2642, df = 52, p-value = 0.001944
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.1625768 0.6125121
# sample estimates:
#   cor 
# 0.4123855 

# log-log transformed - late linearity
YII_Pve_wocon01 <- YII_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(YII_Pve_wocon01$value_log, YII_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Pve_wocon01$value_log and YII_Pve_wocon01$conc_log
# t = 0.72551, df = 52, p-value = 0.4714
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1722729  0.3582629
# sample estimates:
#   cor 
# 0.1001053 


# YII Spi
YII_Spi <- subset(YII, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(YII_Spi$value, YII_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Spi$value and YII_Spi$conc
# t = 0.98479, df = 88, p-value = 0.3274
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1049549  0.3049041
# sample estimates:
#   cor 
# 0.1044055 

# log-log transformed
YII_Spi_wocon <- YII_Spi %>%
  subset(conc != "0")
cor.test(YII_Spi_wocon$value_log, YII_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Spi_wocon$value_log and YII_Spi_wocon$conc_log
# t = 1.9931, df = 70, p-value = 0.05015
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   7.457339e-05 4.397965e-01
# sample estimates:
#   cor 
# 0.2317391 

# log-log transformed - early linearity
YII_Spi_wocon100 <- YII_Spi_wocon %>%
  subset(conc != "100")
cor.test(YII_Spi_wocon100$value_log, YII_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Spi_wocon100$value_log and YII_Spi_wocon100$conc_log
# t = 2.2705, df = 52, p-value = 0.02735
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.03541451 0.52580474
# sample estimates:
#   cor 
# 0.3003272  

# log-log transformed - late linearity
YII_Spi_wocon01 <- YII_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(YII_Spi_wocon01$value_log, YII_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  YII_Spi_wocon01$value_log and YII_Spi_wocon01$conc_log
# t = 1.3515, df = 52, p-value = 0.1824
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.08788263  0.43072768
# sample estimates:
#   cor 
# 0.1842127 


## ---- 4.07. FvFm -------------------------------------------------------------
# FvFm Pve
FvFm_Pve <- subset(FvFm, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(FvFm_Pve$value, FvFm_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Pve$value and FvFm_Pve$conc
# t = 1.914, df = 88, p-value = 0.05887
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.007486955  0.390824915
# sample estimates:
#   cor 
# 0.1999141  

# log-log transformed
FvFm_Pve_wocon <- FvFm_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(FvFm_Pve_wocon$value_log, FvFm_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Pve_wocon$value_log and FvFm_Pve_wocon$conc_log
# t = 1.88, df = 70, p-value = 0.06427
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.01310153  0.42910697
# sample estimates:
#   cor 
# 0.2192324 

# log-log transformed - early linearity
FvFm_Pve_wocon100 <- FvFm_Pve_wocon %>%
  subset(conc != "100")
cor.test(FvFm_Pve_wocon100$value_log, FvFm_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Pve_wocon100$value_log and FvFm_Pve_wocon100$conc_log
# t = -0.24195, df = 52, p-value = 0.8098
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2986131  0.2363490
# sample estimates:
#   cor 
# -0.03353379 

# log-log transformed - late linearity
FvFm_Pve_wocon01 <- FvFm_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(FvFm_Pve_wocon01$value_log, FvFm_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Pve_wocon01$value_log and FvFm_Pve_wocon01$conc_log
# t = 2.5341, df = 52, p-value = 0.01432
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.06999924 0.55044088
# sample estimates:
#   cor 
# 0.3315456  


# FvFm Spi
FvFm_Spi <- subset(FvFm, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(FvFm_Spi$value, FvFm_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Spi$value and FvFm_Spi$conc
# t = -1.7345, df = 88, p-value = 0.08633
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.37479463  0.02626463
# sample estimates:
#   cor 
# -0.1818154 

# log-log transformed
FvFm_Spi_wocon <- FvFm_Spi %>%
  subset(conc != "0")
cor.test(FvFm_Spi_wocon$value_log, FvFm_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Spi_wocon$value_log and FvFm_Spi_wocon$conc_log
# t = -1.4534, df = 70, p-value = 0.1506
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.38745481  0.06301847
# sample estimates:
#   cor 
# -0.1711488  

# log-log transformed - early linearity
FvFm_Spi_wocon100 <- FvFm_Spi_wocon %>%
  subset(conc != "100")
cor.test(FvFm_Spi_wocon100$value_log, FvFm_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Spi_wocon100$value_log and FvFm_Spi_wocon100$conc_log
# t = 0.1632, df = 52, p-value = 0.871
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2466285  0.2886381
# sample estimates:
#   cor 
# 0.02262619  

# log-log transformed - late linearity
FvFm_Spi_wocon01 <- FvFm_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(FvFm_Spi_wocon01$value_log, FvFm_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  FvFm_Spi_wocon01$value_log and FvFm_Spi_wocon01$conc_log
# t = -1.2151, df = 52, p-value = 0.2298
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.4154426  0.1063250
# sample estimates:
#   cor 
# -0.1661663 



## ---- 4.08. rETRmax -------------------------------------------------------------
# rETR Pve
rETR_Pve <- subset(rETR, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(rETR_Pve$value, rETR_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Pve$value and rETR_Pve$conc
# t = -0.15922, df = 88, p-value = 0.8739
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2232773  0.1907908
# sample estimates:
#   cor 
# -0.01697087  

# log-log transformed
rETR_Pve_wocon <- rETR_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(rETR_Pve_wocon$value_log, rETR_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Pve_wocon$value_log and rETR_Pve_wocon$conc_log
# t = 0.1167, df = 70, p-value = 0.9074
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2184274  0.2448244
# sample estimates:
#   cor 
# 0.01394694 

# log-log transformed - early linearity
rETR_Pve_wocon100 <- rETR_Pve_wocon %>%
  subset(conc != "100")
cor.test(rETR_Pve_wocon100$value_log, rETR_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Pve_wocon100$value_log and rETR_Pve_wocon100$conc_log
# t = -0.069092, df = 52, p-value = 0.9452
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2766318  0.2588437
# sample estimates:
#   cor 
# -0.009580939 

# log-log transformed - late linearity
rETR_Pve_wocon01 <- rETR_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(rETR_Pve_wocon01$value_log, rETR_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Pve_wocon01$value_log and rETR_Pve_wocon01$conc_log
# t = -0.27582, df = 52, p-value = 0.7838
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3028826  0.2319122
# sample estimates:
#   cor 
# -0.0382218  


# rETR Spi
rETR_Spi <- subset(rETR, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(rETR_Spi$value, rETR_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Spi$value and rETR_Spi$conc
# t = 0.027661, df = 88, p-value = 0.978
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2042672  0.2099116
# sample estimates:
#   cor 
# 0.002948651 

# log-log transformed
rETR_Spi_wocon <- rETR_Spi %>%
  subset(conc != "0")
cor.test(rETR_Spi_wocon$value_log, rETR_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Spi_wocon$value_log and rETR_Spi_wocon$conc_log
# t = 0.22257, df = 70, p-value = 0.8245
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2063470  0.2566801
# sample estimates:
#   cor 
# 0.02659284  

# log-log transformed - early linearity
rETR_Spi_wocon100 <- rETR_Spi_wocon %>%
  subset(conc != "100")
cor.test(rETR_Spi_wocon100$value_log, rETR_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Spi_wocon100$value_log and rETR_Spi_wocon100$conc_log
# t = 0.2524, df = 52, p-value = 0.8017
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2349812  0.2999315
# sample estimates:
#   cor 
# 0.03498021  

# log-log transformed - late linearity
rETR_Spi_wocon01 <- rETR_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(rETR_Spi_wocon01$value_log, rETR_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  rETR_Spi_wocon01$value_log and rETR_Spi_wocon01$conc_log
# t = 0.23096, df = 52, p-value = 0.8183
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2377866  0.2972250
# sample estimates:
#   cor 
# 0.03201217 



## ---- 4.09. Ek ---------------------------------------------------------------
# Ek Pve
Ek_Pve <- subset(Ek, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(Ek_Pve$value, Ek_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Pve$value and Ek_Pve$conc
# t = -0.30577, df = 88, p-value = 0.7605
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2380632  0.1756983
# sample estimates:
#   cor 
# -0.03257819 

# log-log transformed
Ek_Pve_wocon <- Ek_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(Ek_Pve_wocon$value_log, Ek_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Pve_wocon$value_log and Ek_Pve_wocon$conc_log
# t = -0.068086, df = 70, p-value = 0.9459
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2393549  0.2239532
# sample estimates:
#   cor 
# -0.008137534 

# log-log transformed - early linearity
Ek_Pve_wocon100 <- Ek_Pve_wocon %>%
  subset(conc != "100")
cor.test(Ek_Pve_wocon100$value_log, Ek_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Pve_wocon100$value_log and Ek_Pve_wocon100$conc_log
# t = -0.18104, df = 52, p-value = 0.857
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2909037  0.2443043
# sample estimates:
#   cor 
# -0.02509808 

# log-log transformed - late linearity
Ek_Pve_wocon01 <- Ek_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(Ek_Pve_wocon01$value_log, Ek_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Pve_wocon01$value_log and Ek_Pve_wocon01$conc_log
# t = -0.53669, df = 52, p-value = 0.5938
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.3353164  0.1974651
# sample estimates:
#   cor 
# -0.07421965 


# Ek Spi
Ek_Spi <- subset(Ek, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(Ek_Spi$value, Ek_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Spi$value and Ek_Spi$conc
# t = 0.24363, df = 88, p-value = 0.8081
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1821077  0.2318073
# sample estimates:
#   cor 
# 0.02596256 

# log-log transformed
Ek_Spi_wocon <- Ek_Spi %>%
  subset(conc != "0")
cor.test(Ek_Spi_wocon$value_log, Ek_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Spi_wocon$value_log and Ek_Spi_wocon$conc_log
# t = 0.26028, df = 70, p-value = 0.7954
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2020291  0.2608840
# sample estimates:
#   cor 
# 0.03109481  

# log-log transformed - early linearity
Ek_Spi_wocon100 <- Ek_Spi_wocon %>%
  subset(conc != "100")
cor.test(Ek_Spi_wocon100$value_log, Ek_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Spi_wocon100$value_log and Ek_Spi_wocon100$conc_log
# t = -0.15227, df = 52, p-value = 0.8796
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2872481  0.2480514
# sample estimates:
#   cor 
# -0.0211113  

# log-log transformed - late linearity
Ek_Spi_wocon01 <- Ek_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(Ek_Spi_wocon01$value_log, Ek_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  Ek_Spi_wocon01$value_log and Ek_Spi_wocon01$conc_log
# t = 0.26362, df = 52, p-value = 0.7931
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2335118  0.3013458
# sample estimates:
#   cor 
# 0.03653302 


## ---- 4.10. Alpha ------------------------------------------------------------
# alpha Pve
alpha_Pve <- subset(alpha, spec == "Pve") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))

# untransformed
cor.test(alpha_Pve$value, alpha_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Pve$value and alpha_Pve$conc
# t = 0.50094, df = 88, p-value = 0.6177
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1554841  0.2575708
# sample estimates:
#   cor 
# 0.05332403 

# log-log transformed
alpha_Pve_wocon <- alpha_Pve %>%
  subset(conc != "0" &
           value_log != "-Inf")
cor.test(alpha_Pve_wocon$value_log, alpha_Pve_wocon$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Pve_wocon$value_log and alpha_Pve_wocon$conc_log
# t = 0.63833, df = 70, p-value = 0.5253
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1583865  0.3024124
# sample estimates:
#   cor 
# 0.07607348 

# log-log transformed - early linearity
alpha_Pve_wocon100 <- alpha_Pve_wocon %>%
  subset(conc != "100")
cor.test(alpha_Pve_wocon100$value_log, alpha_Pve_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Pve_wocon100$value_log and alpha_Pve_wocon100$conc_log
# t = 0.38202, df = 52, p-value = 0.704
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2179457  0.3161839
# sample estimates:
#   cor 
# 0.05290212 

# log-log transformed - late linearity
alpha_Pve_wocon01 <- alpha_Pve_wocon %>%
  subset(conc != "0.1")
cor.test(alpha_Pve_wocon01$value_log, alpha_Pve_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Pve_wocon01$value_log and alpha_Pve_wocon01$conc_log
# t = 0.97332, df = 52, p-value = 0.3349
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1389754  0.3876394
# sample estimates:
#   cor 
# 0.1337627 


# alpha Spi
alpha_Spi <- subset(alpha, spec == "Spi") %>%
  mutate(value_log = log(value),
         conc_log = log(conc))
# untransformed
cor.test(alpha_Spi$value, alpha_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Spi$value and alpha_Spi$conc
# t = -0.45272, df = 88, p-value = 0.6519
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2527721  0.1604889
# sample estimates:
#   cor 
# -0.04820436 

# log-log transformed
alpha_Spi_wocon <- alpha_Spi %>%
  subset(conc != "0")
cor.test(alpha_Spi_wocon$value_log, alpha_Spi_wocon$conc_log,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Spi_wocon$value_log and alpha_Spi_wocon$conc_log
# t = 0.077911, df = 70, p-value = 0.9381
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2228376  0.2404615
# sample estimates:
#   cor 
# 0.009311693 

# log-log transformed - early linearity
alpha_Spi_wocon100 <- alpha_Spi_wocon %>%
  subset(conc != "100")
cor.test(alpha_Spi_wocon100$value_log, alpha_Spi_wocon100$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Spi_wocon100$value_log and alpha_Spi_wocon100$conc_log
# t = 0.79494, df = 52, p-value = 0.4303
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.1629676  0.3665796
# sample estimates:
#   cor 
# 0.1095744  

# log-log transformed - late linearity
alpha_Spi_wocon01 <- alpha_Spi_wocon %>%
  subset(conc != "0.1")
cor.test(alpha_Spi_wocon01$value_log, alpha_Spi_wocon01$conc_log, conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  alpha_Spi_wocon01$value_log and alpha_Spi_wocon01$conc_log
# t = 0.09692, df = 52, p-value = 0.9232
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2552398  0.2801915
# sample estimates:
#   cor 
# 0.0134392 




# ----- 5. Other Relationships -------------------------------------------------
## ---- 5.01. Surface ----------------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface, spec == "Pve")

# Exponential
ks.test(surface_Pve$value, surface_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  surface_Pve$value
# D = 0.93825, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation

# Logarithmic
plot(surface_Pve$conc, surface_Pve$value)
#fit the model
model <- lm(surface_Pve_wocon$value ~ log(surface_Pve_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = surface_Pve_wocon$value ~ log(surface_Pve_wocon$conc))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -44.915 -11.806  -3.226  11.903  62.605 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  33.5447     2.6524  12.647  < 2e-16 ***
#   log(surface_Pve_wocon$conc)  -2.5386     0.9406  -2.699  0.00871 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 20.55 on 70 degrees of freedom
# Multiple R-squared:  0.09426,	Adjusted R-squared:  0.08132 
# F-statistic: 7.285 on 1 and 70 DF,  p-value: 0.008709


# Surface Spi
surface_Spi <- subset(surface_all, spec == "Spi")

# Exponential
ks.test(surface_Spi$value, surface_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  surface_Spi$value
# D = 0.96667, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no exp correlation

# Logarithmic
plot(surface_Spi$conc, surface_Spi$value)
#fit the model
model <- lm(surface_Spi_wocon$value ~ log(surface_Spi_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = surface_Spi_wocon$value ~ log(surface_Spi_wocon$conc))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -119.736   -6.760    5.168   18.011   55.468 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   31.529      4.277   7.371 2.62e-10 ***
#   log(surface_Spi_wocon$conc)   -2.419      1.517  -1.595    0.115    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 33.13 on 70 degrees of freedom
# Multiple R-squared:  0.03506,	Adjusted R-squared:  0.02127 
# F-statistic: 2.543 on 1 and 70 DF,  p-value: 0.1153


## ---- 5.02. Volume ----------------------------------------------------------
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

# Logarithmic
plot(volume_Pve$conc, volume_Pve$value)
#fit the model
model <- lm(volume_Pve_wocon$value ~ log(volume_Pve_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = volume_Pve_wocon$value ~ log(volume_Pve_wocon$conc))
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.041033 -0.017318 -0.001417  0.016130  0.068272 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.053140   0.003348  15.872   <2e-16 ***
#   log(volume_Pve_wocon$conc) -0.002957   0.001187  -2.491   0.0151 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.02593 on 70 degrees of freedom
# Multiple R-squared:  0.08142,	Adjusted R-squared:  0.0683 
# F-statistic: 6.205 on 1 and 70 DF,  p-value: 0.01512


# Volume Spi
volume_Spi <- subset(volume_all, spec == "Spi")

# Exponential
ks.test(volume_Spi$value, volume_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  volume_Spi$value
# D = 0.96667, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# Logarithmic
plot(volume_Spi$conc, volume_Spi$value)
#fit the model
model <- lm(volume_Spi_wocon$value ~ log(volume_Spi_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = volume_Spi_wocon$value ~ log(volume_Spi_wocon$conc))
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.056556 -0.018098 -0.005853  0.011293  0.099932 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.052712   0.004239   12.43   <2e-16 ***
#   log(volume_Spi_wocon$conc) -0.002390   0.001503   -1.59    0.116    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.03284 on 70 degrees of freedom
# Multiple R-squared:  0.03486,	Adjusted R-squared:  0.02107 
# F-statistic: 2.528 on 1 and 70 DF,  p-value: 0.1163


## ---- 5.03. Calcification ----------------------------------------------------
# calcification Pve
calcification_Pve <- subset(calcification_all, spec == "Pve")

# Exponential
ks.test(calcification_Pve$value, calcification_Pve$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  calcification_Pve$value
# D = 0.98889, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# Logarithmic
plot(calcification_Pve$conc, calcification_Pve$value)
#fit the model
model <- lm(calcification_Pve_wocon$value ~ log(calcification_Pve_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = calcification_Pve_wocon$value ~ log(calcification_Pve_wocon$conc))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -65.488 -29.306   1.332  22.433  83.632 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         85.970      4.504  19.086   <2e-16 ***
#   log(calcification_Pve_wocon$conc)   -2.013      1.597  -1.261    0.212    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 34.89 on 70 degrees of freedom
# Multiple R-squared:  0.0222,	Adjusted R-squared:  0.008227 
# F-statistic: 1.589 on 1 and 70 DF,  p-value: 0.2117


# calcification Spi
calcification_Spi <- subset(calcification_all, spec == "Spi")

# Exponential
ks.test(calcification_Spi$value, calcification_Spi$conc, y = "pexp")
# OUTPUT: 	Exact one-sample Kolmogorov-Smirnov test
# data:  calcification_Spi$value
# D = 0.96667, p-value = 6.661e-16
# alternative hypothesis: two-sided
# no correlation

# Logarithmic
plot(calcification_Spi$conc, calcification_Spi$value)
#fit the model
model <- lm(calcification_Spi_wocon$value ~ log(calcification_Spi_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = calcification_Spi_wocon$value ~ log(calcification_Spi_wocon$conc))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -100.441  -22.182   -4.081   20.893  120.214 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         87.766      6.199  14.157   <2e-16 ***
#   log(calcification_Spi_wocon$conc)   -3.349      2.198  -1.523    0.132    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 48.02 on 70 degrees of freedom
# Multiple R-squared:  0.03208,	Adjusted R-squared:  0.01826 
# F-statistic:  2.32 on 1 and 70 DF,  p-value: 0.1322



## ---- 5.04. Necrosis ---------------------------------------------------------
# necrosis Pve
necrosis_Pve <- subset(necrosis_all, spec == "Pve")

# Exponential
ks.test(necrosis_Pve$value, necrosis_Pve$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  necrosis_Pve$value
# D = 0.96667, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# Logarithmic
plot(necrosis_Pve$conc, necrosis_Pve$value)
#fit the model
model <- lm(necrosis_Pve_wocon$value ~ log(necrosis_Pve_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = necrosis_Pve_wocon$value ~ log(necrosis_Pve_wocon$conc))
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -16.999 -10.188  -4.330   8.872  23.807 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                    10.559      6.894   1.532    0.177
# log(necrosis_Pve_wocon$conc)    1.665      2.263   0.735    0.490
# Residual standard error: 16.06 on 6 degrees of freedom
# Multiple R-squared:  0.0827,	Adjusted R-squared:  -0.07018 
# F-statistic: 0.5409 on 1 and 6 DF,  p-value: 0.4898


# necrosis Spi
necrosis_Spi <- subset(necrosis, spec == "Spi")

# Exponential
ks.test(necrosis_Spi$value, necrosis_Spi$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  necrosis_Spi$value
# D = 0.96667, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# Logarithmic
plot(necrosis_Spi$conc, necrosis_Spi$value)
#fit the model
model <- lm(necrosis_Spi_wocon$value ~ log(necrosis_Spi_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = necrosis_Spi_wocon$value ~ log(necrosis_Spi_wocon$conc))
# Residuals:
#   Min     1Q Median     3Q    Max 
# -6.385 -5.373 -4.362 -3.350 67.167 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)                    4.3617     1.9295   2.260   0.0269 *
#   log(necrosis_Spi_wocon$conc)   0.4393     0.6842   0.642   0.5230  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 14.95 on 70 degrees of freedom
# Multiple R-squared:  0.005854,	Adjusted R-squared:  -0.008348 
# F-statistic: 0.4122 on 1 and 70 DF,  p-value: 0.523


# polypactivity Pve
polypactivity_Pve <- subset(polypactivity, spec == "Pve")

# Exponential
ks.test(polypactivity_Pve$value, polypactivity_Pve$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  polypactivity_Spi$value
# D = 0.98889, p-value < 2.2e-16
# alternative hypothesis: two-sided
# no correlation

# Logarithmic
plot(polypactivity_Pve$conc, polypactivity_Pve$value)
#fit the model
model <- lm(polypactivity_Pve_wocon$value ~ log(polypactivity_Pve_wocon$conc))
#view the output of the model
summary(model)
# OUTPUT: Call:
# lm(formula = polypactivity_Pve_wocon$value ~ log(polypactivity_Pve_wocon$conc))
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.45895 -0.07006  0.01204  0.08480  0.31883 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                        0.821296   0.019306  42.542  < 2e-16 ***
#   log(polypactivity_Pve_wocon$conc) -0.042491   0.006846  -6.207 3.37e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1495 on 70 degrees of freedom
# Multiple R-squared:  0.355,	Adjusted R-squared:  0.3458 
# F-statistic: 38.53 on 1 and 70 DF,  p-value: 3.367e-08


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
