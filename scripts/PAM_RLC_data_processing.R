# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the data processing and statistical analyzes of the additional parameters of
# photosynthetic efficiency of the corals photosymbionts, derived from rapid light curves (RLCs).
# The parameters were measured using pulse amplitude modulated fluorometry (PAM) and include
#   a) relative electron transport rate (rETRmax)
#   b) efficiency of light capture (α)
#   c)  and minimum saturating irradiance (Ek)
# Statistical analyzes will be conducted in 'PAM_RLC_statistics



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

# for RLC parameter calculation
install.packages("phytotools")
library(phytotools)
# Warning in install.packages: package ‘phytotools’ is not available (for R version 3.6.3)

# ---- 3. Read in needed data files --------------------------------------------
## --- 3.1. Coral identity table -----------------------------------------------
corals <- read_csv2("in/coral_treatments.csv") 
treatments <- corals %>%      
  # to avoid doubling spec, col, and tank columns
  select(-spec, -col, -tank)


## --- 3.2. RLC measurement tables ---------------------------------------------
# Tables are devided by timepoint and species

# ------------- Spi
RLC_t0_Spi <- read_csv2("in/PAM Data/t0_RLC_Spi_all.csv") %>%
  rename(Filename = "# Filename") %>%                       #clean column name
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE)

RLC_t1_Spi <- read_csv2("in/PAM Data/t1_RLC_Spi_all.csv")%>%
  rename(Filename = "# Filename") %>%               
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) 

RLC_t2_Spi <- read_csv2("in/PAM Data/t2_RLC_Spi_all.csv")%>%
  rename(Filename = "# Filename") %>%           
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) 

RLC_t3_Spi <- read_csv2("in/PAM Data/t3_RLC_Spi_all.csv")%>%
  rename(Filename = "# Filename") %>%                
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) 

# ------------- Pve
RLC_t0_Pve <- read_csv2("in/PAM Data/t0_RLC_Pve_all.csv")%>%
  rename(Filename = "# Filename") %>%                  
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) 

RLC_t1_Pve <- read_csv2("in/PAM Data/t1_RLC_Pve_all.csv")%>%
  rename(Filename = "# Filename") %>%                 
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) 

RLC_t2_Pve <- read_csv2("in/PAM Data/t2_RLC_Pve_all.csv")%>%
  rename(Filename = "# Filename") %>%                
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) 

RLC_t3_Pve <- read_csv2("in/PAM Data/t3_RLC_Pve_all.csv")%>%
  rename(Filename = "# Filename") %>%                
  # clean column names for better merge with other timepoints
  separate(Filename, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE)



# ----- 4. Prepare data for statistical analyzes -------------------------------
## ---- 4.1. Unite data tables -------------------------------------------------
# ------------- Pve + Spi data
# merge all timepoints of Spi into long format:
RLC_all_raw <-   rbind(RLC_t0_Pve, RLC_t1_Pve, RLC_t2_Pve, RLC_t3_Pve,
                       RLC_t0_Spi, RLC_t1_Spi, RLC_t2_Spi, RLC_t3_Spi) %>%
  # create column with ID_time
  unite(ID_time, c(ID, tp), sep="_", remove=FALSE) %>%
  select(-spec, -col, -tank) # remove do avoid doubling when joining with info
# ------------- RLC data table + coral ID table
RLC_all_raw <- merge(corals, RLC_all_raw, by = 'ID')

# level treatments
RLC_all_raw$treat <- factor(RLC_all_raw$treat, levels = c("control", "0.1", "1", "10", "100"))

# arrange data by ID_time and then with increasing PAR,
# necessary for calculation of RLCs parameter
RLC_all_raw <- RLC_all_raw %>%
  group_by(ID_time) %>%
  arrange(PAR, .by_group = TRUE)


## ---- 4.2. Calculate RLC parameter -------------------------------------------
# calculation according to the Platt, Gallegos & Harrison 1980 model (PGH)
RLC_PGH <- RLC_all_raw 

# for more details see GitHub script/blog: starting at Step 4.
# https://gfricardo.com/2021/04/07/automated-rapid-light-curve-analyses-in-r/
# it says: For corals we use rETR instead of ETR which is simply the PAR x Yield
# e.g., how to check for anomalies

# use phytotools package 
# for more options and explanation see package phytotools manual
# https://cran.r-project.org/web/packages/phytotools/phytotools.pdf
# written by Antti Takolander 09/2019
# Example script to fit rapid light curves using R package "phytotools" and the model by Platt et al. 1980,
# to extract the coefficients/parameters alpha, beta, ETRmax, Ek and ps from the fitted curve,
# and plotting the curves into individual tiff files.
# The script below relies heavily on code examples provided by the authors of the
# phytotools package. See the code examples of the package for equations for fitting/plotting other light curve models.

# process rapid light curve data
# important: PAR, ETR, unique ID for each curve
# rETR = YII*PAR used

# number of unique IDs in the data
ncurves <- length(unique(RLC_PGH$ID_time)) 
# store the unique IDs
ids <- unique(RLC_PGH$ID_time) 

# create a data frame to store the extracted curve parameters after model fitting
rlc.parameters <- tibble(
  ID_time = ids,
  alpha = 0,
  beta = 0,
  rETRmax = 0,
  Ek = 0,
  ps = 0
)


# the loop below loops runs times specified in ncurves, fits the PGH model to the data,
# plots a fitted curve and data into a tiff file (into current working directory), and extracts the model parameters into rlc.parameters data frame.
# for explanation of the parameters and physiological interpretation see the original publication by Platt et al. 1980.

# Potential problems with the loop are missing PAR or ETR values in the rlc data.
# It is a good practice to check the plots of all fitted curves before relying on / further analyzing any of the fitted values.

for (i in 1:ncurves) {
  # extract the ID of the curve to be fitted
  temp.id <- ids[i] 
  
  # to keep track what's happening if the data has many curves
  print(paste("Now fitting curve ", as.character(temp.id))) 
  
  # extract the data of a single curve into a temporary variable
  temp.rlc.data <- RLC_PGH[RLC_PGH$ID_time == temp.id, ] 
  
  PAR <- temp.rlc.data$PAR
  rETR <- temp.rlc.data$rETR
  
  fit <- fitPGH(PAR, rETR, fitmethod = "Port") 
  
  # store the fitted RLC values into temporary variables
  alpha.rlc <- fit$alpha[1]
  beta.rlc <- fit$beta[1]
  ps.rlc <- fit$ps[1]
  
  # store the parameters
  rlc.parameters$id[i] <- temp.id
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g., Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237).
  rETRmax <- ps.rlc * (alpha.rlc / (alpha.rlc + beta.rlc)) * (beta.rlc / (alpha.rlc + beta.rlc))^(beta.rlc / alpha.rlc)
  Ek <- rETRmax / alpha.rlc
  
  # store the variables
  rlc.parameters$rETRmax[i] <- rETRmax
  rlc.parameters$Ek[i] <- Ek
  
  # plot the data,
  plot(x = PAR, y = rETR, main = temp.id)
  
  # plot the model fit
  with(fit, {
    # the PGH model equation
    P <- ps.rlc * (1 - exp(-1 * alpha.rlc * PAR / ps.rlc)) * exp(-1 * beta.rlc * PAR / ps.rlc) 
    lines(PAR, P)
  }) 
  # close the plotting devide. if this is not done, the next run of the loop will override the plot.
  dev.off() }



# now the data frame rlc.parameters contains the fitted values for each curve. 
rETR_parameter <- rlc.parameters

rETR_parameter <- rETR_parameter %>%
  # remove doubled info
  select(-id)



# ----- 5. Write tables --------------------------------------------------------
## ---- 5.1. Table of all processed RLC data -----------------------------------
write_rds(rETR_parameter, "processed/rETR_parameter.rds")