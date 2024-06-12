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
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
          x_axis = as.numeric(x_axis),
          Color = case_when(conc == "0" ~  "#4A8696",
                            conc == "0.1" ~ "#FFED85",
                            conc == "1" ~ "#E09F3E",
                            conc == "10" ~ "#9E2A2B",
                            conc == "100" ~ "#540B0E"))

## ---- 3.02. Volume -----------------------------------------------------------
volume <- read.csv2("processed/volume_all.csv") %>%
  mutate(conc = as.numeric(conc),
         , 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.03. Calcification ----------------------------------------------------
calcification <- read.csv2("processed/calcification_all.csv") %>%
  mutate(conc = as.numeric(conc),
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.04. Necrosis ---------------------------------------------------------
necrosis <- read.csv2("processed/necrosis_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.05. Polypactivity ----------------------------------------------------
polypactivity <- read.csv2("processed/polypactivity_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.06. YII --------------------------------------------------------------
YII <- read.csv2("processed/YII_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.07. FvFm -------------------------------------------------------------
FvFm <- read.csv2("processed/FvFm_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.08. rETRmax ----------------------------------------------------------
rETR <- read.csv2("processed/rETR_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.09. Ek ---------------------------------------------------------------
Ek <- read.csv2("processed/Ek_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))

## ---- 3.10. Alpha ------------------------------------------------------------
alpha <- read.csv2("processed/alpha_all.csv") %>%
  mutate(conc = as.numeric(conc), 
         x_axis = case_when(conc == "0" ~ "1",
                            conc == "0.1" ~ "2",
                            conc == "1" ~ "3",
                            conc == "10" ~ "4",
                            conc == "100" ~ "5"),
         x_axis = as.numeric(x_axis),
         Color = case_when(conc == "0" ~  "#4A8696",
                           conc == "0.1" ~ "#FFED85",
                           conc == "1" ~ "#E09F3E",
                           conc == "10" ~ "#9E2A2B",
                           conc == "100" ~ "#540B0E"))


# ----- 4. GAM: Relationships -------------------------------------------
## ---- 4.01. Surface ----------------------------------------------------------
### --- 4.01.1. Pve ------------------------------------------------------------
# Surface Pve
surface_Pve <- subset(surface, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value+100))

# fit GAM
# function s within the formula to denote the smooth terms
# bs = cr -> cubic regression splines
mod_gam1 = gam(value_log~ s(conc_log, bs = "cr", k=5), data = surface_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.86311    0.01727   281.5   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc_log) 1.601  1.916 4.468  0.0101 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.093   Deviance explained = 10.9%
# GCV = 0.027651  Scale est. = 0.026852  n = 90

plot(mod_gam1)

# for statistical analyses using LMER and GLMER
library(lme4)
# for statistical testing
library(multcomp)

model_Pve <- glmer(value_log ~ conc_log + (1|col), 
                   family = poisson, data = surface_Pve)


# create graph
plot_surf_Pve <- ggplot(surface_Pve, aes(x = x_axis, y = value_log, 
                                     color = Color)) +
  scale_color_identity() +
  geom_point() +
  geom_smooth(aes(x = x_axis, y = value_log), method = "gam",
              formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5),
              color = "black") +
  scale_x_continuous(labels= c("control", "0.1", "1", "10", "100"), 
                     breaks = c(1, 2, 3, 4, 5)) +
  scale_y_continuous(breaks = c(3.5, 4, 4.5, 5, 5.5),
                     limits=c(3.25, 5.5)) +
  labs(title = "P. verrucosa",
       x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "log(tissue growth)") +
  theme_bw() +
  theme(plot.title = element_text(size=12, face="italic"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_text(size = 10),
    legend.position = "none")

# show plot
plot_surf_Pve

# add GAM statistics 
plot_surf_Pve <- plot_surf_Pve + annotate("text", x = 4.25, y = 5.45, 
                         label = "p = 0.0101, edf = 1.601", size = 3)
plot_surf_Pve

### --- 4.01.2. Spi ------------------------------------------------------------
# Surface Spi
surface_Spi <- subset(surface, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value+100))

# fit GAM
# function s within the formula to denote the smooth terms
# bs = cr -> cubic regression splines
mod_gam2 = gam(value_log ~ s(conc_log, bs = "cr", k=5), data = surface_Spi)
summary(mod_gam2)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.78010    0.06587   72.57   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc_log) 1.529  1.831 1.998  0.0962 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0424   Deviance explained = 5.89%
# GCV = 0.40178  Scale est. = 0.39049   n = 90

plot(mod_gam2)

# create graph
plot_surf_Spi <- ggplot(surface_Spi, aes(x = x_axis, y = value_log, 
                                         color = Color)) +
  scale_color_identity() +
  geom_point() +
  geom_smooth(aes(x = x_axis, y = value_log), method = "gam",
              formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5),
              color = "black") +
  scale_x_continuous(labels= c("control", "0.1", "1", "10", "100"), 
                     breaks = c(1, 2, 3, 4, 5)) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(3.5, 4, 4.5, 5, 5.5),
                     limits=c(3.25, 5.5)) +
  labs(title = "S. pistillata",
       x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "log(tissue growth)") +
  theme_bw() +
  theme(plot.title = element_text(size=12, face="italic"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "black"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

# show plot
plot_surf_Spi

# add GAM statistics 
plot_surf_Spi<- plot_surf_Spi + annotate("text", x = 4.25, y = 5.45, 
                                          label = "p < 0.0001, edf = 1", size = 3)

plot_surf_Spi

# bring both surface plots together
plot_surface <- ggarrange(plot_surf_Pve, plot_surf_Spi,
                    #labels = c("A", "B", "C"),
                    ncol = 2, nrow = 1,
                    align = "v")
                    #  heights= c(1, 0.85, 1.2))
plot_surface


## ---- 4.02. Volume -----------------------------------------------------------
### --- 4.02.1. Pve ------------------------------------------------------------
volume_Pve <- subset(volume, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value+1))

# fit GAM
#function s within the formula to denote the smooth terms
#bs = cr -> cubic regression splines
mod_gam3 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = volume_Pve)
summary(mod_gam3)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.048611   0.002521   19.28   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc_log) 1.903  2.231 5.756 0.00336 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.122   Deviance explained =   14%
# GCV = 0.00059105  Scale est. = 0.00057198  n = 90

plot(mod_gam3)

# create graph
plot_vol_Pve <- ggplot(volume_Pve, aes(x = x_axis, y = value_log, 
                                         color = Color)) +
  scale_color_identity() +
  geom_point() +
  geom_smooth(aes(x = x_axis, y = fitted(mod_gam3)), method = "gam",
              formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5),
              color = "black") +
  scale_x_continuous(labels= c("control", "0.1", "1", "10", "100"), 
                     breaks = c(1, 2, 3, 4, 5)) +
  scale_y_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15),
                     limits=c(-0.05, 0.15))+
  labs(title = "P. verrucosa",
       x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "log(volume growth)") +
  theme_bw() +
  theme(plot.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "black"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.position = "none")

plot_vol_Pve

# add GAM statistics 
plot_vol_Pve <- plot_vol_Pve + annotate("text", x = 4.25, y = 0.15, 
                                          label = "p = 0.0034, edf = 1.903", size = 3)
plot_vol_Pve

### --- 4.02.2. Spi ------------------------------------------------------------
# volume Spi
volume_Spi <- subset(volume, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value+1))

# fit GAM
#function s within the formula to denote the smooth terms
#bs = cr -> cubic regression splines
mod_gam4 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = volume_Spi)
summary(mod_gam4)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.04810    0.00317   15.17   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc_log) 1.763  2.089 2.654  0.0717 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0509   Deviance explained = 6.97%
# GCV = 0.00093315  Scale est. = 0.0009045  n = 90

plot(mod_gam4)

# create graph
plot_vol_Spi <- ggplot(volume_Spi, aes(x = x_axis, y = value_log, 
                                       color = Color)) +
  scale_color_identity() +
  geom_point() +
  geom_smooth(aes(x = x_axis, y = fitted(mod_gam4)), method = "gam",
              formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5),
              color = "black") +
  scale_x_continuous(labels= c("control", "0.1", "1", "10", "100"), 
                     breaks = c(1, 2, 3, 4, 5)) +
  scale_y_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15),
    limits=c(-0.05, 0.15)) +
  labs(title = "S. pistillata",
       x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "log(volume growth)") +
  theme_bw() +
  theme(plot.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "black"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none")

plot_vol_Spi

# add GAM statistics 
plot_vol_Spi <- plot_vol_Spi + annotate("text", x = 4.25, y = 0.15, 
                                        label = "p = 0.0717, edf = 1.763", size = 3)
plot_vol_Spi


# bring both surface plots together
plot_volume <- ggarrange(plot_vol_Pve, plot_vol_Spi,
                          #labels = c("A", "B", "C"),
                          ncol = 2, nrow = 1,
                          align = "v")
#  heights= c(1, 0.85, 1.2))
plot_volume

## ---- 4.03. Calcification ----------------------------------------------------
### --- 4.03.1. Pve ------------------------------------------------------------
calcification_Pve <- subset(calcification, spec == "Pve")  %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam5 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = calcification_Pve)
summary(mod_gam5)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value_log ~ s(conc_log, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.30834    0.05556   77.55   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log) 1.615  1.931 1.783    0.13
# R-sq.(adj) =  0.0363   Deviance explained = 5.38%
# GCV = 0.28611  Scale est. = 0.2778    n = 90

plot(mod_gam5)

plot_calc_Pve <- ggplot(volume_Pve, aes(x = x_axis, y = value_log, 
                                       color = Color)) +
  scale_color_identity() +
  geom_point() +
  geom_smooth(aes(x = x_axis, y = value_log), method = "gam",
              formula = y ~ s(x, bs = "cs", fx = TRUE, k = 5),
              color = "black") +
  scale_x_continuous(labels= c("control", "0.1", "1", "10", "100"), 
                     breaks = c(1, 2, 3, 4, 5)) +
#  scale_y_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15),
 #                    limits=c(-0.05, 0.15))+
  labs(title = "P. verrucosa",
       x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "log(volume growth)") +
  theme_bw() +
  theme(plot.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(colour = "black"),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.position = "none")

plot_calc_Pve

# add GAM statistics 
plot_vol_Pve <- plot_vol_Pve + annotate("text", x = 4.25, y = 0.15, 
                                        label = "p = 0.0034, edf = 1.903", size = 3)
plot_vol_Pve

### --- 4.03.2. Spi ------------------------------------------------------------
# calcification Spi
calcification_Spi <- subset(calcification, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value+10))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = calcification_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.36902    0.07511   58.17   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log) 1.454  1.738 0.947   0.273
# R-sq.(adj) =  0.0181   Deviance explained = 3.41%
# GCV = 0.52194  Scale est. = 0.50771   n = 90

plot(mod_gam1)



## ---- 4.04. Necrosis ---------------------------------------------------------
### --- 4.04.1. Pve ------------------------------------------------------------
necrosis_Pve <- subset(necrosis, spec == "Pve") %>%
  filter(value != 0) %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = necrosis_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   2.1488     0.3616   5.942 0.000831 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log) 2.66  2.926 0.928    0.47
# R-sq.(adj) =  0.119   Deviance explained = 37.9%
# GCV =  2.063  Scale est. = 1.3079    n = 10

plot(mod_gam1)


### --- 4.04.2. Spi ------------------------------------------------------------
# necrosis Spi
necrosis_Spi <- subset(necrosis, spec == "Spi") %>%
  filter(value != 0)  %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = necrosis_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   2.2833     0.3568   6.399 1.55e-05 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log) 1.808  2.109 1.164   0.313
# R-sq.(adj) =  0.101   Deviance explained = 20.2%
# GCV = 2.5926  Scale est. = 2.1644    n = 17

plot(mod_gam1)



## ---- 4.05. Polyp activity ---------------------------------------------------
### --- 4.05.1. Pve ------------------------------------------------------------
polypactivity_Pve <- subset(polypactivity, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = polypactivity_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.26872    0.02507  -10.72   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F  p-value    
# s(conc_log) 1.335  1.575 24.55 1.68e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.315   Deviance explained = 32.5%
# GCV = 0.058079  Scale est. = 0.056572  n = 90

plot(mod_gam1)


### --- 4.05.2. Spi ------------------------------------------------------------
# polypactivity Spi
polypactivity_Spi <- subset(polypactivity, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = polypactivity_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.58671    0.04652  -12.61   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc_log)   1      1 10.18 0.00197 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0935   Deviance explained = 10.4%
# GCV = 0.19921  Scale est. = 0.19479   n = 90

plot(mod_gam1)



## ---- 4.06. YII --------------------------------------------------------------
### --- 4.06.1. Pve ------------------------------------------------------------
YII_Pve <- subset(YII, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log~ s(conc_log, bs = "cr", k = 5), data = YII_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#  (Intercept) 4.631792   0.005441   851.3   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value   
# s(conc_log) 3.743  3.936 4.754 0.00234 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.147   Deviance explained = 18.2%
# GCV = 0.0028123  Scale est. = 0.002664  n = 90

plot(mod_gam1)


### --- 4.06.2. Spi ------------------------------------------------------------
# YII Spi
YII_Spi <- subset(YII, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = YII_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.645481   0.007084   655.7   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc_log) 2.411  2.731 2.605  0.0425 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0725   Deviance explained = 9.77%
# GCV = 0.004695  Scale est. = 0.004517  n = 90

plot(mod_gam1)



## ---- 4.07. FvFm -------------------------------------------------------------
### --- 4.07.1. Pve ------------------------------------------------------------
FvFm_Pve <- subset(FvFm, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = FvFm_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.6250     0.0456   101.4   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value  
# s(conc_log) 1.572  1.882 2.525  0.0592 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# R-sq.(adj) =  0.0538   Deviance explained = 7.05%
# GCV = 0.19269  Scale est. = 0.18718   n = 90

plot(mod_gam1)


### --- 4.07.2. Spi ------------------------------------------------------------
# FvFm Spi
FvFm_Spi <- subset(FvFm, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = FvFm_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.68164    0.04545     103   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log) 1.359   1.61 0.842   0.301
# R-sq.(adj) =  0.0136   Deviance explained = 2.87%
# GCV = 0.19095  Scale est. = 0.18595   n = 90

plot(mod_gam1)



## ---- 4.08. rETRmax ----------------------------------------------------------
### --- 4.08.1. Pve ------------------------------------------------------------
rETR_Pve <- subset(rETR, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = rETR_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.14208    0.03346   153.7   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log)   1      1 0.198   0.657
# R-sq.(adj) =  -0.00909   Deviance explained = 0.225%
# GCV = 0.10308  Scale est. = 0.10079   n = 90

plot(mod_gam1)


### --- 4.08.2. Spi ------------------------------------------------------------
# rETR Spi
rETR_Spi <- subset(rETR, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = rETR_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.90593    0.04244   115.6   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log)   1      1 0.189   0.665
# R-sq.(adj) =  -0.0092   Deviance explained = 0.214%
# GCV = 0.16581  Scale est. = 0.16212   n = 90

plot(mod_gam1)



## ---- 4.09. Ek ---------------------------------------------------------------
### --- 4.09.1. Pve ------------------------------------------------------------
Ek_Pve <- subset(Ek, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = Ek_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# ((Intercept)  5.71400    0.03331   171.6   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log)   1      1 0.412   0.523
# R-sq.(adj) =  -0.00665   Deviance explained = 0.466%
# GCV = 0.1021  Scale est. = 0.099833  n = 90

plot(mod_gam1)


### --- 4.09.2. Spi ------------------------------------------------------------
# Ek Spi
Ek_Spi <- subset(Ek, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = Ek_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  5.58509    0.03095   180.5   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log) 1.048  1.094 0.246   0.689
# R-sq.(adj) =  -0.00837   Deviance explained = 0.35%
# GCV = 0.088212  Scale est. = 0.086204  n = 90

plot(mod_gam1)



## ---- 4.10. Alpha ------------------------------------------------------------
### --- 4.10.1. Pve ------------------------------------------------------------
alpha_Pve <- subset(alpha, spec == "Pve") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = alpha_Pve)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|) 
# (Intercept) -0.57192    0.01022  -55.94   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log)   1      1 0.401   0.528
# R-sq.(adj) =  -0.00677   Deviance explained = 0.454%
# GCV = 0.0096215  Scale est. = 0.0094077  n = 90

plot(mod_gam1)


### --- 4.10.2. Spi ------------------------------------------------------------
# alpha Spi
alpha_Spi <- subset(alpha, spec == "Spi") %>%
  mutate(conc_log = log(conc+1),
         value_log = log(value))

# fit GAM
mod_gam1 = gam(value_log ~ s(conc_log, bs = "cr", k = 5), data = alpha_Spi)
summary(mod_gam1)
# OUTPUT: Family: gaussian 
# Link function: identity 
# Formula:
#   value ~ s(conc, bs = "cr", k = 5)
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|) 
# (Intercept) -0.67916    0.01904  -35.67   <2e-16 ***
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value
# s(conc_log)   1      1 0.039   0.845
# R-sq.(adj) =  -0.0109   Deviance explained = 0.0439%
# GCV = 0.03337  Scale est. = 0.032628  n = 90

plot(mod_gam1)

