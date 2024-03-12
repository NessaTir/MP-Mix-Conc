# ----- 1. Explanation of this script ------------------------------------------
# This script focuses on the visualisation of all processed data measured. 
# Depicted in I) the main manuscript
#   a) Living tissue area growth
#   b) Volume growth
#   c) Calcification growth
#   d) Necrosis
#   e) Polyp activity
#   f) Effective quantum yield
#   g) relative electron transport rate
# and II) the supplements

# This script builds up on following other scripts
# 'Growth_data_processing'
# 'Necrosis_statistics'
# 'Polyactivity_processing_and_statistics'

# ----- 2. Load in needed packages ---------------------------------------------
library(tidyverse)

# read in the data files
library(readxl)

# create graphs
library(ggplot2)
library(ggpmisc)

# work with summarized data, e.g., means
library(rstatix)

# bring different plots together
library("patchwork")

# for annotations in heatmap
library("ggtext")


# ----- 3. Read in needed data files ------------------------------------------- 
## ---- 3.2. Coral Information -------------------------------------------------
corals <- read.csv2("in/coral_treatments.csv", sep=",") %>%
  # modify character of some columns, where necessary
  mutate(treat = as.factor(treat), # column for categorical model
         conc = as.numeric(conc)) %>% # column for continuous model 
  dplyr::select(ID, treat, conc)

## ---- 3.3. Processed growth tables -------------------------------------------
# Surface
surface <- read_rds("processed/surface_growth.rds")

# Volume
volume <- read_rds("processed/volume_growth.rds")

# Calcification
calcification <- read_rds("processed/weight_growth.rds") 

# Necrosis
# summary
necrosis <- read_csv2("out/Summary_necrosis.csv") %>%
  dplyr::select(-c(1))

# percent
necrosis_per <- read_csv2("out/necrosis_percent.csv") %>%
  dplyr::select(-c(1))


## ---- 3.4. Photosynthetic efficiency -----------------------------------------
# effective quantum yield
YII <- read_rds("processed/light_all.rds")
# relative YII to t0
YII_relative <- read_rds("processed/light_relative.rds")

# maximum quantum yield
FvFm <- read_rds("processed/dark_all.rds")
# relative YII to t0
FvFm_relative <- read_rds("processed/dark_relative.rds")

# relative electron transport rate
rETR <- read_rds("processed/rETR_parameter.rds") %>%
  # clean column names for better merge with other timepoints
  separate(ID_time, c('spec', 'col', 'tank', 'tp')) %>%    
  unite(ID, c(spec, col, tank), sep="_", remove=FALSE)%>% 
  # reaname timepoints for clean statistical analyses
  mutate(tp = case_when(tp == "t0"~ "0",
                        tp == "t1"~ "1",
                        tp == "t2"~ "2",
                        tp == "t3"~ "3"),
         tp = as.numeric(tp))

# merge rETR Table with coral infos
rETR_i <-  merge(corals, rETR, by = 'ID', all.x = TRUE)


## ---- 3.5. Polyp activity ----------------------------------------------------
Polyps <- read_rds("processed/polyp_activity.rds")


## ---- 3.6. z-Values ----------------------------------------------------------
#z_values <- read_csv2("in/z-values_random_t.csv") %>%
# mutate(z_value = as.numeric(z_value),
#       Parameter = as.factor(Parameter),
#      stars = as.factor(stars))

## ---- 3.7. Supplements -------------------------------------------------------
# MP Concentration measurements
Concentrations <- read_csv2("in/Concentration_newformat.csv") %>%
  mutate(per_L = ((count/volume)*1000),
         tank = as.character(tank))

# MP Size measurements
MP_size <- read_csv2("in/MP_sizes.csv")



# ----- 4. Create Graphs -------------------------------------------------------
## ---- 4.1. Build up theme ----------------------------------------------------
# Create labels for species
spec_labs <- c("P. verrucosa", "S. pistillata")
names(spec_labs) <- c("Pve", "Spi")

# Create labels for treatments

treat_labs <-  c("0", "0.1", "1", "10", "100")
names(treat_labs) <- c("control", "0.1", "1", "10", "100")

# Create a color scheme
color_scheme <- c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B", "#540B0E")

alphavalues <- c(0.005, 0.3, 0.6, 1)


## ---- 4.2. Growth parameters -------------------------------------------------
### --- 4.2.1. Surface ---------------------------------------------------------
surface_plot <- surface %>%
  # separate by species, concentration and timepoint
  group_by(spec, conc, time) %>%
  # ignore NAs
  na.omit() %>%
  # use mean
  summarise(mean = mean(surface_growth), 
            # include sd 
            sd = sd(surface_growth)) %>% 
  # use cumulative sum of growth parameter
  mutate(mean = cumsum(mean)) %>%
  ggplot() + 
  facet_grid( ~ spec, 
              labeller = labeller(spec = spec_labs)) +
  geom_line(aes(x = time, y = mean, group = conc, color = conc), 
            size = 0.8, alpha = 0.8, linetype = 'longdash', position=position_dodge(width=0.2)) +
  scale_x_continuous(labels= c("4", "8", "12"), breaks = c(1, 2, 3)) +
  scale_color_manual(values = color_scheme) +
  geom_errorbar(aes(x = time, ymin = mean-sd, ymax = mean+sd, 
                    color = conc, linetype = NULL, width = 0.2), size = 1, alpha = 0.95,
                position = position_dodge(width = 0.2)) +
  geom_point(aes(x = time, y = mean, color = conc), 
             position = position_dodge(width=0.2), size = 1.5) +
  labs(color = expression(paste("Treatment ", mg, "·", L^-1)), fill = "Treatment (mg/l)") +
  ylab(expression(atop(Cumulative~tissue, growth~("%")))) +
  theme_classic() +
  # design theme
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        # strip.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        #axis.title.x = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.position = "top")

# show plot
surface_plot


### --- 4.2.2. Volume ----------------------------------------------------------
volume_plot <- volume %>%
  # separate by species, concentration and timepoint
  group_by(spec, conc, time) %>%
  # ignore NAs
  na.omit() %>% 
  # use mean
  summarise(mean = mean(volume_growth),
            # include sd
            sd = sd(volume_growth)) %>%
  # use cumulative sum of growth parameter
  mutate(mean = cumsum(mean)) %>%
  ggplot() + 
  facet_grid( ~ spec) +
  geom_line(aes(x = time, y = mean, group = conc, color = conc), 
            size = 0.8, alpha = 0.8, linetype = 'longdash', position=position_dodge(width=0.2)) +
  scale_x_continuous(labels= c("4", "8", "12"), breaks = c(1, 2, 3)) +
  scale_color_manual(values = color_scheme) +
  geom_errorbar(aes(x = time, ymin = mean-sd, ymax = mean+sd, 
                    color = conc, linetype = NULL, width = 0.2), size = 1, alpha = 0.95,
                position = position_dodge(width = 0.2)) +
  geom_point(aes(x =time, y = mean, color = conc), 
             position = position_dodge(width=0.2), size = 1.5) +
  labs(color = "Treatment [mg/l]", fill = "Treatment [mg/l]") +
  ylab(expression(atop(Cumulative~volume, growth~(cm^3~cm^-2)))) +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), # remove species labels
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# show plot
volume_plot


### --- 4.2.3. Calcification ---------------------------------------------------
weight_plot <- calcification %>%
  # separate by species, concentration and timepoint
  group_by(spec, conc, time) %>%
  # ignore NAs
  na.omit() %>% 
  # use mean
  summarise(mean = mean(weight_growth),
            # include sd
            sd = sd(weight_growth)) %>%
  # use cumulative sum of growth parameter
  mutate(mean = cumsum(mean)) %>%
  ggplot() + 
  facet_grid( ~ spec) +
  geom_line(aes(x = time, y = mean, group = conc, color = conc), 
            size = 0.8, alpha = 0.8, linetype = 'longdash', position=position_dodge(width=0.2)) +
  scale_x_continuous(labels= c("4", "8", "12"), breaks = c(1, 2, 3)) + 
  scale_color_manual(values = color_scheme) +
  geom_errorbar(aes(x = time, ymin = mean-sd, ymax = mean+sd, 
                    color = conc, linetype = NULL, width = 0.2), size = 1, alpha = 0.95,
                position = position_dodge(width = 0.2)) +
  geom_point(aes(x =time, y = mean, color = conc), 
             position = position_dodge(width=0.2), size = 1.5) +
  labs(color = "Treatment [mg/l]", fill = "Treatment [mg/l]") +
  ylab(expression(atop(Cumulative~calcification, rate~(mg~cm^-2)))) +
  #ylab("Cumulative calcification [mg/cm²]") +
  xlab("Weeks") +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), # remove species labels
        axis.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "none") # remove legend

# show plot
weight_plot



### --- 4.2.4. Necrosis --------------------------------------------------------
# level treatment, important for visualisation 
necrosis$treat <- factor(necrosis$treat,
                         levels = c("control", "0.1", "1", "10", "100"))
necrosis$cat <- factor(necrosis$cat,
                       levels = c("none", "low", "moderate", "high"))
# create legend names for category names
cat_labs <-  c("", "low", "moderate", "high")
names(cat_labs) <- c("none", "low", "moderate", "high")

# Create a color scheme with switched order, so it's correct in the graph
color_scheme_2 <- c("#540B0E", "#9E2A2B", "#E09F3E", "#FFED85", "#4A8696")

necro_plot <- necrosis %>% 
  filter(tp =="3") %>% 
  ggplot(aes(x=treat, y=prop, fill= fct_rev(treat)))+
  geom_bar(stat = "identity", aes(alpha=cat)) +
  facet_grid( ~ spec, 
              labeller = labeller(spec = spec_labs)) +
  scale_alpha_manual("cat", values=alphavalues,
                     labels = cat_labs) +
  scale_fill_manual(guide = 'none', values = color_scheme_2)+
  scale_x_discrete(labels = treat_labs) +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)),
       y = "Necrosis (%)",
       title = "") +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "top")

# show plot
necro_plot



# bring all growth plots together
growth_plot <- surface_plot / volume_plot / weight_plot / necro_plot

# safe all growth plots as one
ggsave("out/growth.png", plot = growth_plot,
       scale = 1, width = 18, height = 26, units = c("cm"),
       dpi = 600, limitsize = TRUE)  



## ---- 4.3. Polyp activity ----------------------------------------------------
polyps <- ggplot(Polyps, aes(x = tp, y = ranks)) +
  facet_grid( ~ spec, 
              labeller = labeller(spec = spec_labs)) +
  geom_smooth(aes(x = tp, y = ranks, group = treat, color = treat, fill = treat)) + 
  scale_x_continuous(labels= c("0", "4", "8", "12")) +
  scale_color_manual(values = color_scheme,
                     labels = treat_labs) +
  scale_fill_manual(values = color_scheme,
                    labels = treat_labs) +
  ylab("Mean polyp activity") + 
  xlab("Weeks of exposure") +
  labs(color = "Treatment (mg/l)", fill = "Treatment (mg/l)") +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "top")


# show plot
polyps

# save plot
ggsave("out/polyps.png", plot = polyps,
       scale = 1, width = 18, height = 14, units = c("cm"),
       dpi = 600, limitsize = TRUE)  



## ---- 4.4. Photosynthetic efficiency -----------------------------------------
# ----- 4.4.1. Effective quantum yield -----------------------------------------
# use standardized (each fragment to it's mean value at t0) values

YII_plot <- ggplot(YII_relative, aes(x = tp, y = relativeYII)) +
  facet_grid( ~ spec, 
              labeller = labeller(spec = spec_labs)) +
  geom_smooth(aes(x = tp, y = relativeYII, group = treat, color = treat, fill = treat)) + 
  # brings lines under the boxplots to see changes over time
  scale_x_continuous(labels= c("4", "8", "12"), breaks = c(1, 2, 3)) + 
  scale_color_manual(values = color_scheme, labels=c('0', '0.1', '1', '10', '100')) +
  scale_fill_manual(values = color_scheme, labels=c('0', '0.1', '1', '10', '100')) +
  ylab("Relative Y(II)") +
  xlab("Weeks of exposure") +
  labs(color = "Treatment (mg/l)", fill = "Treatment (mg/l)") +
  # create design
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "top")

# show plot
YII_plot


# ----- 4.4.2. rETR ------------------------------------------------------------
# ------------ relative electron transport rate
# exclude t0
rETRmax <- subset(rETR_i, tp!= "0")

rETRmax_plot <- ggplot(rETRmax, aes(x = tp, y = rETRmax)) +
  facet_grid( ~ spec, 
              labeller = labeller(spec = spec_labs)) +
  geom_smooth(aes(x = tp, y = rETRmax, group = treat, color = treat, fill = treat)) + 
  # brings lines under the boxplots to see changes over time
  scale_x_continuous(labels= c("4", "8", "12"), breaks = c(1, 2, 3)) +
  scale_color_manual(values = color_scheme, labels=c('0', '0.1', '1', '10', '100')) +
  scale_fill_manual(values = color_scheme, labels=c('0', '0.1', '1', '10', '100')) +
  ylab("rETRmax") + 
  xlab("Weeks of exposure") +
  labs(color = "Treatment (mg/l)", fill = "Treatment (mg/l)") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "none")

# show plot
rETRmax_plot

PAM_plot <-  YII_plot / rETRmax_plot  

# show plot
PAM_plot

# save plot
ggsave("out/PAM_plot_relativeYII.png", plot = PAM_plot,
       scale = 1, width = 18, height = 22, units = c("cm"),
       dpi = 600, limitsize = TRUE)  




## ---- 4.5. Heatmap of z-values -----------------------------------------------
# level Parameter in right order
z_values$Parameter <- factor(z_values$Parameter, 
                             levels = c("alpha", "Ek", "rETRmax", "Fv_Fm", "YII",
                                        "polypactivity", 
                                        "weight", "volume", "surface"))

# add asterids 
z_values_end <- z_values %>%
  mutate(
    label = case_when(
      p > 0.05 ~ "",
      p > 0.01 ~ "*",
      p > 0.001 ~ "**",
      !is.na(p) ~ "***",
      TRUE ~ NA_character_
    )
  )



heatmap <- ggplot(z_values, aes(Species, Parameter, fill = z_value))+
  geom_tile(color= c("white"), size=0.1) +
  # create color gradient
  scale_fill_gradient2(low = "#652966",
                       mid = "#7998cc",
                       high = '#4AC424') +
  scale_y_discrete(labels = c("surface" = "Tissue growth", "volume" = "Volume growth", "weight" = "Calcification",
                              "polypactivity" = "Polypactivity", 
                              "YII" = "Y(II)", "Fv_Fm" = "Fv/Fm", "rETRmax" =" rETRmax", "Ek" = "Ek", "alpha" = "Alpha")) +
  scale_x_discrete(label = spec_labs) +
  geom_text(aes(label = stars), color = c("white"), size = 5,
            fill = NA, label.color = NA) +
  ylab("Parameter") + 
  xlab("") +
  theme_minimal(base_size = 8) +
  theme(strip.background = element_blank(),
        # strip.text.x = element_text(face = "italic", size = 12),
        axis.text.x = element_text(size = 10, face = "italic", angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10))

# show graph 
heatmap

# safe graph
ggsave("out/heatmap.png", plot = heatmap,
       scale = 1, width = 8, height = 12, units = c("cm"),
       dpi = 600, limitsize = TRUE)  

## ---- 4.6. Advanced heatmap as correlation plots -------------------------
# create summary table with cumulative surface growth week 0-12
surface_all <-  surface %>%
  # separate by species and concentration
  group_by(ID, spec, conc) %>%
  # add column with parameter
  # ignore NAs
  na.omit() %>%
  # use mean
  summarise(value = sum(surface_growth)) %>% 
  mutate(parameter = "surface",
         conc = as.numeric(conc)) 

# create summary table with cumulative volume growth week 0-12
volume_all <-  volume %>%
  # separate by species and concentration
  group_by(ID, spec, conc) %>%
  # ignore NAs
  na.omit() %>%
  # use mean
  summarise(value = sum(volume_growth))%>% 
  mutate(parameter = "volume",
         conc = as.numeric(conc)) 

# create summary table with cumulative weight growth week 0-12
calcification_all <-  calcification %>%
  # separate by species and concentration
  group_by(ID, spec, conc) %>%
  # ignore NAs
  na.omit() %>%
  # use mean
  summarise(value = sum(weight_growth))%>% 
  mutate(parameter = "calcification",
         conc = as.numeric(conc)) 

# create summary table with relative necrosis after 12 weeks
necrosis_all <-  necrosis_per %>%
  # select only the last timepoint
  filter(tp=="3") %>%
  # ignore NAs
  na.omit() %>%
  #rename column enrty
  rename(value = necro_per) %>% 
  mutate(parameter = "necrosis",
         conc = as.numeric(conc)) 
necrosis_all <- necrosis_all %>%
  # remove unnecessary colums
  dplyr::select(-necro_occured, -treatment, -col, -tank, -tp, -treat)

# create summary table with relative YII after 12 weeks
YII_all <-  YII_relative %>%
  # select only the last timepoint
  filter(tp=="3") %>%
  # ignore NAs
  na.omit() %>%
  #rename column enrty
  rename(value = relativeYII) %>% 
  # use mean
  mutate(parameter = "YII")   %>% 
  #keep only relevant columns
  dplyr::select(ID, spec, conc, value, parameter)

# create summary table with relative FvFm after 12 weeks
FvFm_all <-  FvFm_relative %>%
  # select only the last timepoint
  filter(tp=="3") %>%
  # ignore NAs
  na.omit() %>%
  #rename column enrty
  rename(value = relativeFv_Fm) %>% 
  # use mean
  mutate(parameter = "FvFm")   %>% 
  #keep only relevant columns
  dplyr::select(ID, spec, conc, value, parameter)

# create summary table with rETR after 12 weeks
rETR_all <-  rETR_i %>%
  # select only the last timepoint
  filter(tp=="3") %>%
  # ignore NAs
  na.omit() %>%
  #rename column enrty
  rename(value = rETRmax) %>% 
  # use mean
  mutate(parameter = "rETRmax")   %>% 
  #keep only relevant columns
  dplyr::select(ID, spec, conc, value, parameter)


# create summary table with alpha after 12 weeks
alpha_all <-  rETR_i %>%
  # select only the last timepoint
  filter(tp=="3") %>%
  # ignore NAs
  na.omit() %>%
  #rename column enrty
  rename(value = alpha) %>% 
  # use mean
  mutate(parameter = "alpha")   %>% 
  #keep only relevant columns
  dplyr::select(ID, spec, conc, value, parameter)


# create summary table with Ek after 12 weeks
Ek_all <-  rETR_i %>%
  # select only the last timepoint
  filter(tp=="3") %>%
  # ignore NAs
  na.omit() %>%
  #rename column enrty
  rename(value = Ek) %>% 
  # use mean
  mutate(parameter = "Ek")   %>% 
  #keep only relevant columns
  dplyr::select(ID, spec, conc, value, parameter)

# create summary table with mean polyp acrtivity over time of exposure
polypactivity_all <-  Polyps %>%
  # select only the last timepoint
  filter(tp!="0") %>%
  # separate by species and concentration
  group_by(ID, spec, conc) %>%
  # ignore NAs
  na.omit() %>%
  # use mean
  summarise(value = mean(ranks))%>% 
  mutate(parameter = "polypactivity",
         conc = as.numeric(conc)) 

# bring all tables together
all_data <- rbind(surface_all, volume_all, calcification_all, necrosis_all,
                  polypactivity_all,
                  YII_all, FvFm_all, 
                  rETR_all, alpha_all, Ek_all)

all_data$parameter <- factor(all_data$parameter, 
                             levels = c("surface", "volume", "calcification", "necrosis", "polypactivity",
                                        "YII", "FvFm", "rETRmax", "Ek", "alpha"))


ggplot(all_data, aes(conc, value)) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  stat_poly_line(color = "black") +
  #scale_color_continuous(treat = color_scheme) +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10')+
  stat_correlation(p.digits = 4,
                   aes(label = paste(after_stat(r.label),
                                     after_stat(p.value.label),
                                     sep = "*\"; \"*")),
                   method = "pearson", conf.level = 0.95) +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Value") +
  # stat_correlation(use_label(c("R", "P")))+
  theme_bw()+
  theme(#legend.position= "none", #c(0.1, 0.9),
    #legend.direction = "horizontal",
    # legend.text = element_text(size = 10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_text(face = "italic", size = 12),
    strip.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12))


# safe graph
ggsave("out/correlations_log.png", plot = last_plot(),
       scale = 1, width = 18, height = 24, units = c("cm"),
       dpi = 600, limitsize = TRUE)  



# ----- 5. Supplements  --------------------------------------------------------
## ---- 5.1. Advanced Hearmap --------------------------------------------------
# Statistics for supplement table
# Surface Pve
surface_Pve <- subset(surface_all, spec == "Pve")
cor.test(surface_Pve$value, surface_Pve$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve$value and surface_Pve$conc
# t = -3.2999, df = 88, p-value = 0.001398
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5042723 -0.1339487
# sample estimates:
#  cor 
# -0.3318348 

# Surface Spi
surface_Spi <- subset(surface_all, spec == "Spi")
cor.test(surface_Spi$value, surface_Spi$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: 	Pearson's product-moment correlation
# data:  surface_Pve$value and surface_Pve$conc
# t = -3.2999, df = 88, p-value = 0.001398
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.5042723 -0.1339487
# sample estimates:
#  cor 
# -0.3318348 






# CONTINUE HERE ---------------------
# standardize each fragment to it's mean value at t0
FvFm_0 <- subset(FvFm, tp== "0") %>% 
  # remove unnecessary columns for clear merge
  dplyr::select(ID, Fv_Fm) %>% 
  rename(FvFm_t0 = Fv_Fm) 

FvFm_relative <- subset(FvFm, tp!= "0")

FvFm_relative <- full_join(FvFm_relative, FvFm_0, by="ID")

FvFm_relative <- FvFm_relative %>% 
  mutate(relativeFvFm = 100/FvFm_t0*Fv_Fm) 
#-------------------------------------


## ---- 5.1. MP concentrations -------------------------------------------------
Conc_curve <- ggplot(Concentrations, aes(x = days, y = per_L)) +
  geom_smooth(aes(x = days, y = per_L, group = conc, color = conc, fill = conc)) + 
  # brings lines under the boxplots to see changes over time
  scale_color_manual(values = color_scheme, labels=c('0', '0.1', '1', '10', '100')) +
  scale_fill_manual(values = color_scheme, labels=c('0', '0.1', '1', '10', '100')) +
  ylab("Particles per litre") +
  xlab("Days of exposure") +
  labs(color = "Treatment (mg/l)", fill = "Treatment (mg/l)") +
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "top")


# safe graph
ggsave("out/Concentrations.png", plot = Conc_curve,
       scale = 1, width = 16, height = 18, units = c("cm"),
       dpi = 600, limitsize = TRUE)



## ---- 5.1.1. Table of MP concentrations --------------------------------------
# get mean and sd of Concentrations
mean_concentration <- Concentrations %>% # table
  group_by(conc) %>% # separates species, treatment and timepoints
  get_summary_stats(per_L) %>% # columns of interest: parameter
  select(-q1, -q3, -iqr, -mad, -se, -ci, -median) # because only selecting what wanted did not work: remove unnecessary columns
# write it into .csv
write_csv2(mean_concentration, "out/mean_concentration.csv")


## ---- 5.2. MP Sizes ----------------------------------------------------------
# get mean and sd of MP Sizes
mean_Size <- MP_size %>% # table
  group_by(Form, Polymer) %>% # separates species, treatment and timepoints
  get_summary_stats(Length, Width) %>% # columns of interest: parameter
  select(-q1, -q3, -iqr, -mad, -se, -ci, -median) # because only selecting what wanted did not work: remove unnecessary columns

# write it into .csv
write_csv2(mean_Size, "out/mean_MP_size.csv")


# ----- 6. Write tables --------------------------------------------------------
## ---- 6.1. MP Size measurements ----------------------------------------------

