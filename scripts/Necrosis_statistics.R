# ----- 1. Explanation of this script ----
  # This script focuses on the assessment of coral necrosis.
  # Occurence of necrosis was assessed visually as yes/no for every sampling timepoint
  # Relative necrosis was measured using 3D Scanning
  
  
# ----- 2. load in needed packages -------
# to easily clean data
library(tidyverse)

# to easily read in all data files
library(readxl)

library(dplyr)



library(lme4)
library(multcomp)
library(viridis)
library(car) # for glht
library(tidyverse) # to read in table with read_rds
library(readxl)
library(rstatix) # statistik helfer
library(broomExtra) # write models like lm as tibble
library(magrittr) # A Forward-Pipe Operator

library(ggpubr)
library(skimr) # to skim() - check for missing data
library(performance) # check_normality


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
  pivot_longer(cols = c('0', '1', '2', '3'), # use previous colums as new entries - categories
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
  select(-ID, -spec, -col, -origin, -tank, -tp, -treat, -conc)
Necro_percent <- Necro_percent %>%
  # create a new column to get a clean merge
  unite(ID_tp, c(ID, tp), sep="_", remove=FALSE) %>%
  # leave out some columns to have a cleaner merge with the percent table
  select(-ID, -spec, -col, -tank, -tp, -treat, -conc)

# merge into one table
Necrosis <-   merge(Necro_percent, Necro_occurence, by = 'ID_tp', all.x = TRUE) %>%
  separate(ID_tp, c('spec', 'col', 'tank', 'tp')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) %>%
  # make timepoint column numeric for further analyses
  mutate(tp = as.numeric(tp))

# add treatment info to table
Necrosis <- Necrosis %>% 
  mutate(tank = as.factor(tank),
         treat = case_when(tank == 1 ~ "control",
                           tank == 6 ~ "control",
                           tank == 11 ~ "control",
                           tank == 2 ~ "0.1",
                           tank == 7 ~ "0.1",
                           tank == 12 ~ "0.1",
                           tank == 3 ~ "1",
                           tank == 8 ~ "1",
                           tank == 13 ~ "1",
                           tank == 4 ~ "10",
                           tank == 9 ~ "10",
                           tank == 14 ~ "10",
                           TRUE ~ "100"))

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


# check levels of certain columns to evaluate for releveling 
# treatment
# use treatment as categories not numbers (if numbers necessary: use "conc" column)
Necrosis$treat <- factor(Necrosis$treat, 
                            levels = c("control", "0.1", "1", "10", "100"))
levels(Necrosis$treat)


# ---- 4.1. Categorize severity of necrosis ----
# according to Marshall and Schuttenberg 2006
# low: 1-10%
# moderate: 0-50%
# high: > 50%
Necrosis_category <- Necrosis %>% 
  mutate(cat = case_when(necro_per >= 50 ~ "high",
                           necro_per >= 10 ~ "moderate",
                           necro_per >= 1 ~ "low",
                           TRUE ~ "none"))



# ---- 5. Statistical analyses ----
# For t0 in both species for each treatment 100 % of the corals showed no necrosis (see 'Summary_necrosis')
# Therefore for t0 no statistical analyses needed

# Both Species will be assessed separately
# For each timepoint statistical analyzes will be conducted separately, 
# comparing the frequency of occurence of the used categories between the 5 diffferent treatments using the fisher_test
# --> Direct comparison




## ----5.2. Pocillopora verrucosa ----
### -- 5.1.1. t1  ----
# Therefore, create a subset tables

# ---- compared treatments: control - 0.1 mg/l
Pve_t1_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "control" | treat == "0.1")
Pve_t1_1$treat <-  factor(Pve_t1_1$treat,
                          levels = c("control", "0.1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 1 mg/l
Pve_t1_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "control" | treat == "1")
Pve_t1_1$treat <-  factor(Pve_t1_1$treat,
                          levels = c("control", "1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36  0.486 ns 


# ---- compared treatments: control - 10 mg/l
Pve_t1_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "control" | treat == "10") 
Pve_t1_1$treat <-  factor(Pve_t1_1$treat,
                          levels = c("control", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 100 mg/l
Pve_t1_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "control" | treat == "100") 
Pve_t1_1$treat <-  factor(Pve_t1_1$treat,
                          levels = c("control", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 1 mg/l
Pve_t1_2 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "0.1" | treat == "1") 
Pve_t1_2$treat <-  factor(Pve_t1_2$treat,
                          levels = c("0.1", "1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.486 ns 


# ---- compared treatments: 0.1 mg/l - 10 mg/l
Pve_t1_2 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "0.1" | treat == "10") 
Pve_t1_2$treat <-  factor(Pve_t1_2$treat,
                          levels = c("0.1", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 100 mg/l
Pve_t1_2 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "0.1" | treat == "100") 
Pve_t1_2$treat <-  factor(Pve_t1_2$treat,
                          levels = c("0.1", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns


# ---- compared treatments: 1 mg/l - 10 mg/l
Pve_t1_3 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "1" | treat == "10") 
Pve_t1_3$treat <-  factor(Pve_t1_3$treat,
                          levels = c("1", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.486 ns 


# ---- compared treatments: 1 mg/l - 100 mg/l
# 36     0.486 ns
Pve_t1_3 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "1" | treat == "100") 
Pve_t1_3$treat <-  factor(Pve_t1_3$treat,
                          levels = c("1", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.486 ns 


# ---- compared treatments: 10 mg/l - 100 mg/l
Pve_t1_4 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "1") %>%
  filter(treat == "10" | treat == "100") 
Pve_t1_4$treat <-  factor(Pve_t1_4$treat,
                          levels = c("10", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t1_4))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


### -- 5.1.2. t2  ----
# Therefore, create a subset tables

# ---- compared treatments: control - 0.1 mg/l
Pve_t2_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "control" | treat == "0.1") 
Pve_t2_1$treat <-  factor(Pve_t2_1$treat,
                          levels = c("control", "0.1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 1 mg/l
Pve_t2_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "control" | treat == "1")
Pve_t2_1$treat <-  factor(Pve_t2_1$treat,
                          levels = c("control", "1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.486 ns 


# ---- compared treatments: control - 10 mg/l
Pve_t2_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "control" | treat == "10") 
Pve_t2_1$treat <-  factor(Pve_t2_1$treat,
                          levels = c("control", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 100 mg/l
Pve_t2_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "control" | treat == "100") 
Pve_t2_1$treat <-  factor(Pve_t2_1$treat,
                          levels = c("control", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 1 mg/l
Pve_t2_2 <- Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "0.1" | treat == "1") 
Pve_t2_2$treat <-  factor(Pve_t2_2$treat,
                          levels = c("0.1", "1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 10 mg/l
Pve_t2_2 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "0.1" | treat == "10") 
Pve_t2_2$treat <-  factor(Pve_t2_2$treat,
                          levels = c("0.1", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 100 mg/l
Pve_t2_2 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "0.1" | treat == "100") 
Pve_t2_2$treat <-  factor(Pve_t2_2$treat,
                          levels = c("0.1", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 1 mg/l - 10 mg/l
Pve_t2_3 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "1" | treat == "10") 
Pve_t2_3$treat <-  factor(Pve_t2_3$treat,
                          levels = c("1", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 1 mg/l - 100 mg/l
Pve_t2_3 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "1" | treat == "100") 
Pve_t2_3$treat <-  factor(Pve_t2_3$treat,
                          levels = c("1", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 10 mg/l - 100 mg/l
Pve_t2_4 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "2") %>%
  filter(treat == "10" | treat == "100") 
Pve_t2_4$treat <-  factor(Pve_t2_4$treat,
                          levels = c("10", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t2_4))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


### -- 5.1.3. t3  ----
# Therefore, create a subset tables

# ---- compared treatments: control - 0.1 mg/l
Pve_t3_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "control" | treat == "0.1") 
Pve_t3_1$treat <-  factor(Pve_t3_1$treat,
                          levels = c("control", "0.1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 1 mg/l 
Pve_t3_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "control" | treat == "1")
Pve_t3_1$treat <-  factor(Pve_t3_1$treat,
                          levels = c("control", "1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.112 ns 


# ---- compared treatments: control - 10 mg/l 
Pve_t3_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "control" | treat == "10") 
Pve_t3_1$treat <-  factor(Pve_t3_1$treat,
                          levels = c("control", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 100 mg/l 
Pve_t3_1 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "control" | treat == "100") 
Pve_t3_1$treat <-  factor(Pve_t3_1$treat,
                          levels = c("control", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.603 ns 


# ---- compared treatments: 0.1 mg/l - 1 mg/l
Pve_t3_2 <- Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "0.1" | treat == "1") 
Pve_t3_2$treat <-  factor(Pve_t3_2$treat,
                          levels = c("0.1", "1"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.229 ns 


# ---- compared treatments: 0.1 mg/l - 10 mg/l
Pve_t3_2 <- Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "0.1" | treat == "10") 
Pve_t3_2$treat <-  factor(Pve_t3_2$treat,
                          levels = c("0.1", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 100 mg/l
Pve_t3_2 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "0.1" | treat == "100") 
Pve_t3_2$treat <-  factor(Pve_t3_2$treat,
                          levels = c("0.1", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.735 ns 


# ---- compared treatments: 1 mg/l - 10 mg/l    
Pve_t3_3 <- Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "1" | treat == "10") 
Pve_t3_3$treat <-  factor(Pve_t3_3$treat,
                          levels = c("1", "10"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.229 ns 


# ---- compared treatments: 1 mg/l - 100 mg/l
Pve_t3_3 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "1" | treat == "100") 
Pve_t3_3$treat <-  factor(Pve_t3_3$treat,
                          levels = c("1", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 10 mg/l - 100 mg/l
Pve_t3_4 <-  Necrosis_category %>%
  filter(spec == "Pve" &
           tp == "3") %>%
  filter(treat == "10" | treat == "100") 
Pve_t3_4$treat <-  factor(Pve_t3_4$treat,
                          levels = c("10", "100"))
fisher_test(xtabs (~ treat + cat, data = Pve_t3_4))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.735 ns 



## ---- 5.2. Stylophora pistillata ----
### --- 5.2.1. t1  ----
# Therefore, create a subset tables

# ---- compared treatments: control - 0.1 mg/l 
Spi_t1_1 <-   Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "control" | treat == "0.1") 
Spi_t1_1$treat <-  factor(Spi_t1_1$treat,
                          levels = c("control", "0.1"))

fisher_test(xtabs (~ treat + cat, data = Spi_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36     1 ns 


# ---- compared treatments: control - 1 mg/l
Spi_t1_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "control" | treat == "1")
Spi_t1_1$treat <-  factor(Spi_t1_1$treat,
                          levels = c("control", "1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.735 ns 


# ---- compared treatments: control - 10 mg/l
Spi_t1_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "control" | treat == "10") 
Spi_t1_1$treat <-  factor(Spi_t1_1$treat,
                          levels = c("control", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.603 ns 


# ---- compared treatments: control - 100 mg/l
Spi_t1_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "control" | treat == "100") 
Spi_t1_1$treat <-  factor(Spi_t1_1$treat,
                          levels = c("control", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.338 ns 


# ---- compared treatments: 0.1 mg/l - 1 mg/l
Spi_t1_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "0.1" | treat == "1") 
Spi_t1_2$treat <-  factor(Spi_t1_2$treat,
                          levels = c("0.1", "1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.229 ns 


# ---- compared treatments: 0.1 mg/l - 10 mg/l
Spi_t1_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "0.1" | treat == "10") 
Spi_t1_2$treat <-  factor(Spi_t1_2$treat,
                          levels = c("0.1", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.229 ns 


# ---- compared treatments: 0.1 mg/l - 100 mg/l
Spi_t1_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "0.1" | treat == "100") 
Spi_t1_2$treat <-  factor(Spi_t1_2$treat,
                          levels = c("0.1", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.104 ns 


# ---- compared treatments: 1 mg/l - 10 mg/l
Spi_t1_3 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "1" | treat == "10") 
Spi_t1_3$treat <-  factor(Spi_t1_3$treat,
                          levels = c("1", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36  1 ns 


# ---- compared treatments: 1 mg/l - 100 mg/l
Spi_t1_3 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "1" | treat == "100") 
Spi_t1_3$treat <-  factor(Spi_t1_3$treat,
                          levels = c("1", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.692 ns 


# ---- compared treatments: 10 mg/l - 100 mg/l
Spi_t1_4 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "1") %>%
  filter(treat == "10" | treat == "100") 
Spi_t1_4$treat <-  factor(Spi_t1_4$treat,
                          levels = c("10", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t1_4))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


### --- 5.2.2. t2  ----
# Therefore, create a subset tables
# ---- compared treatments: control - 0.1 mg/l 
Spi_t2_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "control" | treat == "0.1") 
Spi_t2_1$treat <-  factor(Spi_t2_1$treat,
                          levels = c("control", "0.1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 1 mg/l  
Spi_t2_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "control" | treat == "1")
Spi_t2_1$treat <-  factor(Spi_t2_1$treat,
                          levels = c("control", "1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  36 0.692 ns  


# ---- compared treatments: control - 10 mg/l 
Spi_t2_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "control" | treat == "10") 
Spi_t2_1$treat <-  factor(Spi_t2_1$treat,
                          levels = c("control", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.795 ns 


# ---- compared treatments: control - 100 mg/l 
Spi_t2_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "control" | treat == "100") 
Spi_t2_1$treat <-  factor(Spi_t2_1$treat,
                          levels = c("control", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    0.539 ns 


# ---- compared treatments: 0.1 mg/l - 1 mg/l
Spi_t2_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "0.1" | treat == "1") 
Spi_t2_2$treat <-  factor(Spi_t2_2$treat,
                          levels = c("0.1", "1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36  0.539 ns 


# ---- compared treatments: 0.1 mg/l - 10 mg/l
Spi_t2_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "0.1" | treat == "10") 
Spi_t2_2$treat <-  factor(Spi_t2_2$treat,
                          levels = c("0.1", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36  0.404 ns 


# ---- compared treatments: 0.1 mg/l - 100 mg/l
Spi_t2_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "0.1" | treat == "100") 
Spi_t2_2$treat <-  factor(Spi_t2_2$treat,
                          levels = c("0.1", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.436 ns 


# ---- compared treatments: 1 mg/l - 10 mg/l
Spi_t2_3 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "1" | treat == "10") 
Spi_t2_3$treat <-  factor(Spi_t2_3$treat,
                          levels = c("1", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 1 mg/l - 100 mg/l
Spi_t2_3 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "1" | treat == "100") 
Spi_t2_3$treat <-  factor(Spi_t2_3$treat,
                          levels = c("1", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 10 mg/l - 100 mg/l
Spi_t2_4 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "2") %>%
  filter(treat == "10" | treat == "100") 
Spi_t2_4$treat <-  factor(Spi_t2_4$treat,
                          levels = c("10", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t2_4))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 



### -- 5.2.3. t3  ----
# Therefore, create a subset tables
# ---- compared treatments: control - 0.1 mg/l
Spi_t3_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "control" | treat == "0.1") 
# level only treatments that are there, otherwise they will get a 0 entry
Spi_t3_1$treat <-  factor(Spi_t3_1$treat,
                          levels = c("control", "0.1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: control - 1 mg/l
Spi_t3_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "control" | treat == "1")
Spi_t3_1$treat <-  factor(Spi_t3_1$treat,
                          levels = c("control", "1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36  0.735 ns 


# ---- compared treatments: control - 10 mg/l
Spi_t3_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "control" | treat == "10") 
Spi_t3_1$treat <-  factor(Spi_t3_1$treat,
                          levels = c("control", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.229 ns 


# ---- compared treatments: control - 100 mg/l
Spi_t3_1 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "control" | treat == "100") 
Spi_t3_1$treat <-  factor(Spi_t3_1$treat,
                          levels = c("control", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_1))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.122 ns 


# ---- compared treatments: 0.1 mg/l - 1 mg/l
Spi_t3_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "0.1" | treat == "1") 
Spi_t3_2$treat <-  factor(Spi_t3_2$treat,
                          levels = c("0.1", "1"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36    1 ns 


# ---- compared treatments: 0.1 mg/l - 10 mg/l
Spi_t3_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "0.1" | treat == "10") 
Spi_t3_2$treat <-  factor(Spi_t3_2$treat,
                          levels = c("0.1", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.229 ns 


# ---- compared treatments: 0.1 mg/l - 100 mg/l
Spi_t3_2 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "0.1" | treat == "100") 
Spi_t3_2$treat <-  factor(Spi_t3_2$treat,
                          levels = c("0.1", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_2))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.305 ns 


# ---- compared treatments: 1 mg/l - 10 mg/l
Spi_t3_3 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "1" | treat == "10") 
Spi_t3_3$treat <-  factor(Spi_t3_3$treat,
                          levels = c("1", "10"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.603 ns 


# ---- compared treatments: 1 mg/l - 100 mg/l
Spi_t3_3 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "1" | treat == "100") 
Spi_t3_3$treat <-  factor(Spi_t3_3$treat,
                          levels = c("1", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_3))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.735 ns 


# ---- compared treatments: 10 mg/l - 100 mg/l
# 36 0.487 ns 
Spi_t3_4 <-  Necrosis_category %>%
  filter(spec == "Spi" &
           tp == "3") %>%
  filter(treat == "10" | treat == "100") 
Spi_t3_4$treat <-  factor(Spi_t3_4$treat,
                          levels = c("10", "100"))
fisher_test(xtabs (~ treat + cat, data = Spi_t3_4))
# Output: A tibble: 1 x 3
#  n     p p.signif
# * <int> <dbl> <chr>   
#  1    36 0.487 ns 



# ---- 6. Write tables ----
## --- 6.1.  Summary of necrosis occurences ----
# level the categories
Necrosis_category$cat <- factor(Necrosis_category$cat, 
                              levels = c( "none", "low", "moderate", "high"))

# create a table with percentages of the categories (cat)
# per Species (spec), treatment (treat), timepoint (tp)
Summary_necrosis <- Necrosis_category %>%
  freq_table(spec, treat, tp, cat, na.rm = T)

write.csv2(Summary_necrosis, "out/Summary_necrosis.csv")
