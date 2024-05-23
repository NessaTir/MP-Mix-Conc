# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the statistical analysis of the relationship of all processed data and the concentrations. 
# Relationship will be tested using Generalized Additive Models (GAM)

# This scripts statistics is closely used in 'Plots - 4.5. Correlation Plot'

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

# GAM
library(mgcv)


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


# ----- 4. GAM: Relationships -------------------------------------------
## ---- 4.01. Surface ----------------------------------------------------------
### --- 4.01.1. Pve ------------------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface, spec == "Pve")

# fit GAM
# function s within the formula to denote the smooth terms
# bs = cr -> cubic regression splines
mod_gam1 = gam(value~ s(conc, bs = "cr", k=5), data = surface_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   31.310      2.237      14   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc)   1      1 10.89  0.0014 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =    0.1   Deviance explained =   11%
# GCV = 460.57  Scale est. = 450.34    n = 90

plot(mod_gam1)


### --- 4.01.2. Spi ------------------------------------------------------------
# Surface Spi
surface_Spi <- subset(surface, spec == "Spi")

# fit GAM
# function s within the formula to denote the smooth terms
# bs = cr -> cubic regression splines
mod_gam1 = gam(value ~ s(conc, bs = "cr", k=5), data = surface_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   29.859      3.382   8.828 9.28e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc)   1      1 6.579   0.012 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.059   Deviance explained = 6.96%
# GCV = 1053.1  Scale est. = 1029.7    n = 90

plot(mod_gam1)



## ---- 4.02. Volume -----------------------------------------------------------
### --- 4.02.1. Pve ------------------------------------------------------------
volume_Pve <- subset(volume, spec == "Pve") 

# fit GAM
#function s within the formula to denote the smooth terms
#bs = cr -> cubic regression splines
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = volume_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.050151   0.002654    18.9   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc) 1.246  1.432 8.332 0.00145 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.124   Deviance explained = 13.6%
# GCV = 0.0006501  Scale est. = 0.00063388  n = 90

plot(mod_gam1)


### --- 4.02.2. Spi ------------------------------------------------------------
# volume Spi
volume_Spi <- subset(volume, spec == "Spi")

# fit GAM
#function s within the formula to denote the smooth terms
#bs = cr -> cubic regression splines
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = volume_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.050151   0.002654    18.9   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc) 1.246  1.432 8.332 0.00145 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.124   Deviance explained = 13.6%
# GCV = 0.0006501  Scale est. = 0.00063388  n = 90

plot(mod_gam1)



## ---- 4.03. Calcification ----------------------------------------------------
### --- 4.03.1. Pve ------------------------------------------------------------
calcification_Pve <- subset(calcification, spec == "Pve") 

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = calcification_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   83.396      3.597   23.18   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc) 1.046   1.09 3.687  0.0494 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0352   Deviance explained = 4.65%
# GCV = 1191.8  Scale est. = 1164.7    n = 90

plot(mod_gam1)


### --- 4.03.2. Spi ------------------------------------------------------------
# calcification Spi
calcification_Spi <- subset(calcification, spec == "Spi") 

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = calcification_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    84.97       4.85   17.52   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc)   1      1 5.688  0.0192 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =   0.05   Deviance explained = 6.07%
# GCV = 2165.4  Scale est. = 2117.3    n = 90

plot(mod_gam1)



## ---- 4.04. Necrosis ---------------------------------------------------------
### --- 4.04.1. Pve ------------------------------------------------------------
necrosis_Pve <- subset(necrosis, spec == "Pve") %>%
  filter(value != 0)

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = necrosis_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   83.396      3.597   23.18   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc) 1.046   1.09 3.687  0.0494 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0352   Deviance explained = 4.65%
# GCV = 1191.8  Scale est. = 1164.7    n = 90

plot(mod_gam1)


### --- 4.04.2. Spi ------------------------------------------------------------
# necrosis Spi
necrosis_Spi <- subset(necrosis, spec == "Spi") %>%
  filter(value != 0)

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = necrosis_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   23.268      5.878   3.958  0.00137 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc) 1.697   1.91 1.055   0.333
# R-sq.(adj) =  0.0845   Deviance explained = 18.2%
# GCV = 698.17  Scale est. = 587.39    n = 17

plot(mod_gam1)



## ---- 4.05. Polyp activity ---------------------------------------------------
### --- 4.05.1. Pve ------------------------------------------------------------
polypactivity_Pve <- subset(polypactivity, spec == "Pve") 

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = polypactivity_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.79074    0.01533   51.59   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(conc) 3.207  3.377 14.77  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.343   Deviance explained = 36.7%
# GCV = 0.022179  Scale est. = 0.021142  n = 90

plot(mod_gam1)


### --- 4.05.2. Spi ------------------------------------------------------------
# polypactivity Spi
polypactivity_Spi <- subset(polypactivity, spec == "Spi")

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = polypactivity_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.60926    0.02391   25.48   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc)   1      1 6.005  0.0162 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0532   Deviance explained = 6.39%
# GCV = 0.052614  Scale est. = 0.051445  n = 90

plot(mod_gam1)



## ---- 4.06. YII --------------------------------------------------------------
### --- 4.06.1. Pve ------------------------------------------------------------
YII_Pve <- subset(YII, spec == "Pve") 

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = YII_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#  (Intercept) 102.8572     0.5656   181.8   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc) 3.687  3.902 4.756 0.00318 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.141   Deviance explained = 17.7%
# GCV = 30.378  Scale est. = 28.795    n = 90

plot(mod_gam1)


### --- 4.06.2. Spi ------------------------------------------------------------
# YII Spi
YII_Spi <- subset(YII, spec == "Spi")

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = YII_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 104.3657     0.7406   140.9   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc) 1.897  2.002 4.707  0.0114 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0868   Deviance explained = 10.6%
# GCV =  51.01  Scale est. = 49.368    n = 90

plot(mod_gam1)



## ---- 4.07. FvFm -------------------------------------------------------------
### --- 4.07.1. Pve ------------------------------------------------------------
FvFm_Pve <- subset(FvFm, spec == "Pve")

# fit GAM
mod_gam1 = gam(value ~ s(conc, bs = "cr", k = 5), data = FvFm_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  113.029      5.953   18.99   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc)   1      1 3.663  0.0589 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0291   Deviance explained =    4%
# GCV = 3261.7  Scale est. = 3189.3    n = 90

plot(mod_gam1)


### --- 4.07.2. Spi ------------------------------------------------------------
# FvFm Spi
FvFm_Spi <- subset(FvFm, spec == "Spi")

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = FvFm_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  118.764      5.933   20.02   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc)   1      1 3.008  0.0863 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0221   Deviance explained = 3.31%
# GCV = 3240.1  Scale est. = 3168.1    n = 90

plot(mod_gam1)



## ---- 4.08. rETRmax ----------------------------------------------------------
### --- 4.08.1. Pve ------------------------------------------------------------
rETR_Pve <- subset(rETR, spec == "Pve")

# fit GAM
mod_gam1 = gam(value ~ s(conc, bs = "cr", k = 5), data = rETR_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  179.096      5.477    32.7   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc) 3.454  3.703 0.354   0.851
# R-sq.(adj) =  -0.0176   Deviance explained = 2.19%
# GCV =   2840  Scale est. = 2699.5    n = 90

plot(mod_gam1)


### --- 4.08.2. Spi ------------------------------------------------------------
# rETR Spi
rETR_Spi <- subset(rETR, spec == "Spi")

# fit GAM
mod_gam1 = gam(value~ s(conc, bs = "cr", k = 5), data = rETR_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  144.663      5.223    27.7   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc) 3.448  3.696 0.312   0.761
# R-sq.(adj) =  -0.0123   Deviance explained = 2.69%
# GCV = 2582.7  Scale est. = 2455.1    n = 90

plot(mod_gam1)



## ---- 4.09. Ek ---------------------------------------------------------------
### --- 4.09.1. Pve ------------------------------------------------------------
Ek_Pve <- subset(Ek, spec == "Pve")

# fit GAM
mod_gam1 = gam(value ~ s(conc, bs = "cr", k = 5), data = Ek_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   317.57      10.05   31.61   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df    F p-value
# s(conc) 3.445  3.693 0.51    0.77
# R-sq.(adj) =  -0.0129   Deviance explained = 2.63%
# GCV =   9554  Scale est. = 9082.1    n = 90

plot(mod_gam1)


### --- 4.09.2. Spi ------------------------------------------------------------
# Ek Spi
Ek_Spi <- subset(Ek, spec == "Spi")

# fit GAM
mod_gam1 = gam(value ~ s(conc, bs = "cr", k = 5), data = Ek_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  277.383      8.118   34.17   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc) 3.52   3.77 0.541   0.612
# R-sq.(adj) =  -0.00266   Deviance explained =  3.7%
# GCV =   6245  Scale est. = 5931.4    n = 90

plot(mod_gam1)



## ---- 4.10. Alpha ------------------------------------------------------------
### --- 4.10.1. Pve ------------------------------------------------------------
alpha_Pve <- subset(alpha, spec == "Pve") 

# fit GAM
mod_gam1 = gam(value ~ s(conc, bs = "cr", k = 5), data = alpha_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|) 
# (Intercept) 0.567001   0.005641   100.5   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc)   1      1 0.251   0.618
# R-sq.(adj) =  -0.00849   Deviance explained = 0.284%
# GCV = 0.0029295  Scale est. = 0.0028644  n = 90


### --- 4.10.2. Spi ------------------------------------------------------------
# alpha Spi
alpha_Spi <- subset(alpha, spec == "Spi")

# fit GAM
mod_gam1 = gam(value ~ s(conc, bs = "cr", k = 5), data = alpha_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|) 
# (Intercept) 0.514520   0.008716   59.03   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc)   1      1 0.205   0.652
# R-sq.(adj) =  -0.00901   Deviance explained = 0.232%
# GCV = 0.0069924  Scale est. = 0.006837  n = 90