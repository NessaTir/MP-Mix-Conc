# ----- 1. Explanation of this script ----
# This script focuses on the calculation of the corals growth rates. Growth was measured in
#   a) Living tissue area growth
#       --> Differences between gross surface area and net tissue area might be caused by necrosis (analyzed in 'Necrosis_statistics')
#   b) Volume growth
#   c) Calcification growth


# ----- 2. Load in needed packages -------
# to easily clean data
library(tidyverse)

# to easily read in all data files
library(readxl)

# to calculate buoyant weight of corals into dry weight of corals
library(seacarb)


# ----- 3. Read in needed data files ----- 
## ---- 3.1. coral identity table --------
# read in list with overview of all corals used and their treatments etc.
corals <- read_csv2("in/coral_treatments.csv")

## ---- 3.2. 3D Scan Day Info ----
# read in data needed for correction of growth rates
# days between each 3D Samling are given
correct_days <- read_csv2("in/3D_Dates.csv") %>%
  #remove unnecessary columns
  select(-t0, -t1, -t2, -t3, -Spec, -Colony, -tank, -treatment, -"t0-t3") %>%
  rename(ID = Fragment,
         "days_1" = "t0-t1",
         "days_2" = "t1-t2",
         "days_3" = "t2-t3")

## ---- 3.2. Tissue and Volune -----------------
# Raw data of 3D Scans are provided as .txt files devided by sampling timepoint
# raw tables include necrosis areas, that have to be excluded with aid of further data tables
# for t0 no data available of:
# Spi_D_2
# Spi_D_7
# Spi_B_14 

### --- 3.2.1. 3D Scan raw data ----
# 3D Scanning table for t0
raw_t0 <- read_tsv("in/t0_basic.txt") %>%
  rename(surf_t0 = Surface,
         vol_t0 = Volume) %>%
  select(-Vertices, -Edges, -Faces) %>%
  # split Filename column into different colums to create ID colums. 
    # Get timepoint (tp), species (spec), colony (col) and tank
  separate(Filename, c('tp', 'spec', 'col', 'tank')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep = "_", remove = FALSE) %>%
  # remove spec, col, and tank to get cleaner merge with coral identity table
  select(-tp,-spec, -col, -tank) 

# 3D Scanning table for t1
raw_t1 <- read_tsv("in/t1_basic.txt") %>%
  rename(surf_t1 = Surface,
         vol_t1 = Volume) %>%
  select(-Vertices, -Edges, -Faces) %>%
  # split Filename column into different colums to create ID colums. 
  # Get timepoint (tp), species (spec), colony (col) and tank
  separate(Filename, c('tp', 'spec', 'col', 'tank')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) %>%
  # remove spec, col, and tank to get cleaner merge with coral identity table
  select(-tp,-spec, -col, -tank) 

# 3D Scanning table for t2
raw_t2 <- read_tsv("in/t2_basic.txt") %>%
  rename(surf_t2 = Surface,
         vol_t2 = Volume) %>%
  select(-Vertices, -Edges, -Faces) %>%
  # split Filename column into different colums to create ID colums. 
  # Get timepoint (tp), species (spec), colony (col) and tank
  separate(Filename, c('tp', 'spec', 'col', 'tank')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) %>%
  # remove spec, col, and tank to get cleaner merge with coral identity table
  select(-tp,-spec, -col, -tank) 

# 3D Scanning table for t3
raw_t3 <- read_tsv("in/t3_basic.txt") %>%
  rename(surf_t3 = Surface,
         vol_t3 = Volume) %>%
  select(-Vertices, -Edges, -Faces) %>%
  # split Filename column into different colums to create ID colums. 
  # Get timepoint (tp), species (spec), colony (col) and tank
  separate(Filename, c('tp', 'spec', 'col', 'tank')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) %>%
  # remove spec, col, and tank to get cleaner merge with coral identity table
  select(-tp,-spec, -col, -tank) 


### -- 3.2.2. table with corrected surface data ----
# to exclude necrotic areas from living tissue area
correct_tissue <- read_tsv("in/correct_tissue.txt") %>%
  rename(cor_surf = Surface) %>%
  select(-Vertices, -Edges, -Faces, -Volume) %>%
  # split Filename column into different colums to create ID colums. 
  # Get timepoint (tp), species (spec), colony (col), tank and additional name (corsurf)
  separate(Filename, c('tp', 'spec', 'col', 'tank', 'corsurf')) %>%
  # unite spec, col and tank to get ID
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE) %>%
  # remove unneccessary columns
  select(-spec, -col, -tank, -corsurf)


### -- 3.2.3. Additional corrections -----
# table with corrected and priviosly missing 3D data
corrected_models <- read.csv2("in/rechecked_3Ddata.csv") %>%
  select(-spec, -col, -tank, -treat, -conc, -origin) %>%
  mutate(surf_t0 = as.numeric(as.character(surf_t0)),
         surf_t1 = as.numeric(as.character(surf_t1)),
         surf_t2 = as.numeric(as.character(surf_t2)),
         surf_t3 = as.numeric(as.character(surf_t3)),
         vol_t0 = as.numeric(as.character(vol_t0)),
         vol_t1 = as.numeric(as.character(vol_t1)),
         vol_t2 = as.numeric(as.character(vol_t2)),
         vol_t3 = as.numeric(as.character(vol_t3)),
         living_t0 = as.numeric(as.character(living_t0)),
         living_t1 = as.numeric(as.character(living_t1)),
         living_t2 = as.numeric(as.character(living_t2)),
         living_t3 = as.numeric(as.character(living_t3))) %>%
  rename(surf_t0_new = surf_t0,
         surf_t1_new = surf_t1,
         surf_t2_new = surf_t2,
         surf_t3_new = surf_t3,
         vol_t0_new = vol_t0,
         vol_t1_new = vol_t1,
         vol_t2_new = vol_t2,
         vol_t3_new = vol_t3,
         living_t0_new = living_t0,
         living_t1_new = living_t1,
         living_t2_new = living_t2,
         living_t3_new = living_t3,
         ID = ï..ID) 


## ---- 3.3 Buoyant weight raw Data ----
# Weight table for Pocillopora verrucosa
weight_Pve <- read.csv2("in/Weight_Pve.csv") %>%
  # remove unnecessary columns
  select(-treatment, -col, -origin)

#Weight table fot Stylophora pistillata
weight_Spi <- read.csv2("in/Weight_Spi.csv") %>%
  # rename coral to ID for better merge
  rename(ID = coral)               




# ----- 4. Data processing ---------
## ---- 4.1. Tissue and Volume -----
# merge tables by ID
# create a list of all tables, that need to be merged to merge them together in one
all_tp <- list(corals, raw_t0, raw_t1, raw_t2, raw_t3)
# create the function for merging
create_merge <- function(corals, raw_t0){
  merge(corals, raw_t0, by = 'ID', all.x = TRUE) # all.x = TRUE for left outer join to fill missing data with NA instead of the whole coral
}

# apply created function on list of tables to merge with "Reduce"
all_scans <- Reduce(create_merge, all_tp) 



## ---- 4.1.1 Surface correction ----
# correct living tissue surface area for necrosis
# bring corrected surface data into wide format
correct_tissue <- correct_tissue %>%
  pivot_wider(names_from = tp,
              values_from = cor_surf)

# merge corrected surface data into all scan table
all_scans <- merge(all_scans, correct_tissue, by = 'ID', all.x = TRUE) %>%
  rename(living_t0 = t0,
         living_t1 = t1,
         living_t2 = t2,
         living_t3 = t3) 

# fill in NA entries of corrected columns with already correct surface data
# --> corrected surface data = living tissue area --> living_tx
all_scans$living_t0 <-  ifelse(is.na(all_scans$living_t0), 
                               all_scans$surf_t0, all_scans$living_t0)
all_scans$living_t1 <-  ifelse(is.na(all_scans$living_t1), 
                               all_scans$surf_t1, all_scans$living_t1)
all_scans$living_t2 <-  ifelse(is.na(all_scans$living_t2), 
                               all_scans$surf_t2, all_scans$living_t2)
all_scans$living_t3 <-  ifelse(is.na(all_scans$living_t3), 
                               all_scans$surf_t3, all_scans$living_t3)


# merge corrected model data into all scan table
all_scans <- merge(all_scans, corrected_models, by = 'ID', all.x = TRUE)

# fill in entries of corrected columns with already correct data
# --> corrected data = _new
all_scans$living_t0 <-  ifelse(is.na(all_scans$living_t0_new), 
                               all_scans$living_t0, all_scans$living_t0_new)
all_scans$living_t1 <-  ifelse(is.na(all_scans$living_t1_new), 
                               all_scans$living_t1, all_scans$living_t1_new)
all_scans$living_t2 <-  ifelse(is.na(all_scans$living_t0_new), 
                               all_scans$living_t2, all_scans$living_t2_new)
all_scans$living_t3 <-  ifelse(is.na(all_scans$living_t3_new), 
                               all_scans$living_t3, all_scans$living_t3_new)

all_scans$surf_t0 <-  ifelse(is.na(all_scans$surf_t0_new), 
                             all_scans$surf_t0, all_scans$surf_t0_new)
all_scans$surf_t1 <-  ifelse(is.na(all_scans$surf_t1_new), 
                             all_scans$surf_t1, all_scans$surf_t1_new)
all_scans$surf_t2 <-  ifelse(is.na(all_scans$surf_t0_new), 
                             all_scans$surf_t2, all_scans$surf_t2_new)
all_scans$surf_t3 <-  ifelse(is.na(all_scans$surf_t3_new), 
                             all_scans$surf_t3, all_scans$surf_t3_new)

all_scans$vol_t0 <-  ifelse(is.na(all_scans$vol_t0_new), 
                            all_scans$vol_t0, all_scans$vol_t0_new)
all_scans$vol_t1 <-  ifelse(is.na(all_scans$vol_t1_new), 
                            all_scans$vol_t1, all_scans$vol_t1_new)
all_scans$vol_t2 <-  ifelse(is.na(all_scans$vol_t0_new), 
                            all_scans$vol_t2, all_scans$vol_t2_new)
all_scans$vol_t3 <-  ifelse(is.na(all_scans$vol_t3_new), 
                            all_scans$vol_t3, all_scans$vol_t3_new)

# remove unnecessary columns
all_scans <- all_scans %>%
  select(-living_t0_new, -living_t1_new, -living_t2_new, -living_t3_new,
         -surf_t0_new, -surf_t1_new, -surf_t2_new, -surf_t3_new,
         -vol_t0_new, -vol_t1_new, -vol_t2_new, -vol_t3_new)


# removed outliers/ very bad models from the data,
# to trick R entired with 0 --> change it to NA
all_scans$surf_t0[all_scans$surf_t0 == 0] = NA
all_scans$vol_t0[all_scans$vol_t0 == 0] = NA
all_scans$living_t0[all_scans$living_t0 == 0] = NA
all_scans$surf_t1[all_scans$surf_t1 == 0] = NA
all_scans$vol_t1[all_scans$vol_t1 == 0] = NA
all_scans$living_t1[all_scans$living_t1== 0] = NA
all_scans$surf_t2[all_scans$surf_t2 == 0] = NA
all_scans$vol_t2[all_scans$vol_t2 == 0] = NA
all_scans$living_t2[all_scans$living_t2 == 0] = NA


# the correction of the 3D Scan data (Volume and linving Tissue area) is now completed


## 4.2. Data of weight ----
# make one table including raw weight data of Pve and Spi and coral identity table
weight_raw <-   rbind(weight_Pve, weight_Spi) %>% 
  #change entires to numeric to use it for further calculations
  mutate(t0 = as.numeric(as.character(t0)),
         t1 = as.numeric(as.character(t1)),
         t2 = as.numeric(as.character(t2)),
         t3 = as.numeric(as.character(t3)),
         hook = as.numeric(as.character(hook)))

# subtract the weight of the hooks used for weighting to be able to calculate the dry weight of the coral only
# create a new table for simplicity
weight <- weight_raw %>% 
  mutate(t0_sub = (t0 - hook),
         t1_sub = (t1 - hook),
         t2_sub = (t2 - hook),
         t3_sub = (t3 - hook)) %>%
  select(-t0, -t1, -t2, -t3, -hook)


# convert buoyant weight to weight using the seacarb package 
# according to the method of Jokiel et al. 1978 and Davies 1989
real_weight <- function(buoyant_weight, S, T, P = 0, rho_aragonite = 2930){
  
  x <- seacarb::rho(S = S, T = T, P = P)
  y <- buoyant_weight / (1 - (x / rho_aragonite))
  attributes(y) <- NULL
  y
}

# replace buoyant weight with dry weight 
# Characteristics of waterbody used for buoyant weighing:
# Salinity: 35
# Temperature: 26.5°C
weight$weight_t0 <- real_weight(weight$t0_sub, S=35, T=26.5, P = 0, rho_aragonite = 2930)
weight$weight_t1 <- real_weight(weight$t1_sub, S=35, T=26.5, P = 0, rho_aragonite = 2930)
weight$weight_t2 <- real_weight(weight$t2_sub, S=35, T=26.5, P = 0, rho_aragonite = 2930)
weight$weight_t3 <- real_weight(weight$t3_sub, S=35, T=26.5, P = 0, rho_aragonite = 2930)





## ---- 4.3. Calculate growth rates ----
# for further calculations write a new table using all of the scan data
growth_3D <- all_scans

### --- 4.3.1.1. Living tissue ----
# growth in surface area (%)
growth_3D$growth_surf_1 <- ((100/growth_3D$living_t0)*growth_3D$living_t1)-100
growth_3D$growth_surf_2 <- ((100/growth_3D$living_t1)*growth_3D$living_t2)-100
growth_3D$growth_surf_3 <- ((100/growth_3D$living_t2)*growth_3D$living_t3)-100

### --- 4.3.1.2. Necrosis -----
# necrotic surface area (%)
growth_3D$necrosis_0 <- (100/growth_3D$surf_t0)*(growth_3D$surf_t0-growth_3D$living_t0)
growth_3D$necrosis_1 <- (100/growth_3D$surf_t1)*(growth_3D$surf_t1-growth_3D$living_t1)
growth_3D$necrosis_2 <- (100/growth_3D$surf_t2)*(growth_3D$surf_t2-growth_3D$living_t2)
growth_3D$necrosis_3 <- (100/growth_3D$surf_t3)*(growth_3D$surf_t3-growth_3D$living_t3)
# It is necessary to check whether all negative values are related to necrosis or whether other sources are responsible
# For further assessment of necrosis see Scripts > Necrosis_statistics

### --- 4.3.2. Volume ----
# growth in volume (cm³ per cm²)
# before volume and surface are given in mm³ and cm²
growth_3D$growth_vol_1 <- (growth_3D$vol_t1 - growth_3D$vol_t0)/1000/(growth_3D$living_t0/100)
growth_3D$growth_vol_2 <- (growth_3D$vol_t2 - growth_3D$vol_t1)/1000/(growth_3D$living_t1/100)
growth_3D$growth_vol_3 <- (growth_3D$vol_t3 - growth_3D$vol_t2)/1000/(growth_3D$living_t2/100)


### --- 4.3.3. Correct 3D data for days between sampling ----
# get growth happening in 4 weeks
growth_rates <- merge(growth_3D, correct_days, by = 'ID', all.x = TRUE)

growth_rates <- growth_rates %>%
  mutate(growth_surf_1 = (growth_surf_1/days_1) * 28,
         growth_surf_2 = (growth_surf_2/days_2) * 28,
         growth_surf_3 = (growth_surf_3/days_3) * 28,
         growth_vol_1 = (growth_vol_1/days_1) * 28,
         growth_vol_2 = (growth_vol_2/days_2) * 28,
         growth_vol_3 = (growth_vol_3/days_3) * 28)


### --- 4.3.4 include weight data ----
# growth in dry weight mg per cm² living tissue area
all_growth <- merge(growth_rates, weight, by = 'ID', all.x = TRUE)

all_growth <- all_growth %>%
  # remove unnecessary columns
  select(-days_1, -days_2, -days_3, -t0_sub, -t1_sub, -t2_sub, -t3_sub) %>%
  # growth in weight mg per cm² surface area
  mutate(growth_weight_1 = ((weight_t1 - weight_t0)*1000) / (surf_t0/100),
         growth_weight_2 = ((weight_t2 - weight_t1)*1000) / (surf_t1/100),
         growth_weight_3 = ((weight_t3 - weight_t2)*1000) / (surf_t2/100))


## ---- 4.4. Reformat tables ----
# write separate tables for each of the growth parameters

# surface
surface_rates <- all_growth %>% 
  select(c(ID, spec, col, tank, treat, conc, 
           growth_surf_1, growth_surf_2, growth_surf_3)) %>% 
  # bring tables into long format
  pivot_longer(cols=c(growth_surf_1, growth_surf_2, growth_surf_3), names_to = "time",
               values_to = "surface_growth", values_drop_na = FALSE)  %>% 
  # rename entires of new column
  mutate(time = case_when(time == "growth_surf_1"~ "1",
                        time == "growth_surf_2" ~ "2",
                        time == "growth_surf_3"~ "3"),
         time = as.numeric(as.character(time))) # change to numeric, to adjust x-axis

# volume
volume_rates <- all_growth %>%
  select(c(ID, spec, col, tank, treat, conc,
           growth_vol_1, growth_vol_2, growth_vol_3)) %>% 
  # bring tables into long format
  pivot_longer(cols=c(growth_vol_1, growth_vol_2, growth_vol_3), names_to = "time",
               values_to = "volume_growth", values_drop_na = FALSE)  %>%
  # rename entires of new column
  mutate  (time = case_when(time == "growth_vol_1"~ "1",
                            time == "growth_vol_2" ~ "2",
                            time == "growth_vol_3"~ "3"),
           time = as.numeric(as.character(time)))


# weight
weight_rates <- all_growth %>%
  select(c(ID, spec, col, tank, treat, conc,
           growth_weight_1, growth_weight_2, growth_weight_3)) %>% 
  # bring tables into long format
  pivot_longer(cols=c(growth_weight_1, growth_weight_2, growth_weight_3), names_to = "time",
               values_to = "weight_growth", values_drop_na = FALSE)  %>%
  # rename entires of new column
  mutate  (time = case_when(time == "growth_weight_1"~ "1",
                            time == "growth_weight_2" ~ "2",
                            time == "growth_weight_3"~ "3"),
           time = as.numeric(as.character(time)))


# necrosis
necrosis_percent <- all_growth %>%
  select(c(ID, spec, col, tank, treat, conc,
           necrosis_0, necrosis_1, necrosis_2, necrosis_3)) %>%
  # bring tables into long format
  pivot_longer(cols = c(necrosis_0, necrosis_1, necrosis_2, necrosis_3),
               names_to = 'tp',
               values_to = 'necrosis', values_drop_na = FALSE) %>%
  # rename entires of new column
  mutate  (tp = case_when(tp == "necrosis_0"~ "0",
                          tp == "necrosis_1"~ "1",
                          tp == "necrosis_2" ~ "2",
                          tp == "necrosis_3"~ "3"),
           tp = as.numeric(as.character(tp)))



# ---- 5. Save tables ----
# write tables of all growth data in long format as .rds to easily read them in for following steps
write_rds(surface_rates, "processed/surface_growth.rds")
write_rds(volume_rates, "processed/volume_growth.rds")
write_rds(weight_rates, "processed/weight_growth.rds")
write_rds(necrosis_percent, "processed/necrosis_percent.rds")
