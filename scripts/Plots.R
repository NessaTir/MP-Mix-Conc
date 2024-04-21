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


## ---- 3.6. Supplements -------------------------------------------------------
# MP Concentration measurements
Concentrations <- read_csv2("in/Concentration_newformat.csv") %>%
  mutate(per_L = ((count/volume)*1000),
         tank = as.character(tank),
         conc = as.numeric(conc))

# MP Size measurements
# MP_size <- read_csv2("in/MP_sizes.csv")



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
        legend.position = "top",
        plot.margin = margin(t = 2,  # Top margin
                             r = 2,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2)) # Left margin))

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
        legend.position = "none",
        plot.margin = margin(t = 2,  # Top margin
                             r = 2,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2)) # Left margin))

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
  xlab("Weeks of exposure") +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(), # remove species labels
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.position = "none", # remove legend
        plot.margin = margin(t = 2,  # Top margin
                             r = 2,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2)) # Left margin))
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
        legend.position = "top",
        plot.margin = margin(t = 2,  # Top margin
                             r = 2,  # Right margin
                             b = 2,  # Bottom margin
                             l = 2)) # Left margin))

# show plot
necro_plot



# bring all growth plots together
growth_plot <- surface_plot / volume_plot / weight_plot / necro_plot

# safe all growth plots as one
# ggsave("out/growth.png", plot = growth_plot,
#       scale = 1, width = 18, height = 26, units = c("cm"),
#       dpi = 600, limitsize = TRUE)  




## ---- 4.3. Polyp activity ----------------------------------------------------
Polyps_wo0 <-  Polyps %>%
  # select everything but the timepoint before the treatment
  filter(tp!="0")
polyps <- ggplot(Polyps_wo0, aes(x = tp, y = ranks)) +
  facet_grid( ~ spec, 
              labeller = labeller(spec = spec_labs)) +
  geom_smooth(aes(x = tp, y = ranks, group = treat, color = treat, fill = treat)) + 
  scale_x_continuous(labels= c("4", "8", "12"), breaks = c(1, 2, 3)) +
  scale_color_manual(values = color_scheme,
                     labels = treat_labs) +
  scale_fill_manual(values = color_scheme,
                    labels = treat_labs) +
  ylab("Mean polyp activity") + 
  xlab("Weeks of exposure") +
  labs(color = expression(paste("Treatment ", mg, "·", L^-1)), 
       fill = expression(paste("Treatment ", mg, "·", L^-1))) +
  # design theme
  theme_classic() +
  theme(panel.background = element_rect(color = "black"),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_blank(), # remove species labels
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.position = "top",
        plot.margin = margin(t = 2,  # Top margin
                             r = 2,  # Right margin
                             l = 2,  # Left margin
                             b = 2)) # Bottom margin)


# show plot
polyps

# save plot
# ggsave("out/polyps.png", plot = polyps,
#        scale = 1, width = 18, height = 6, units = c("cm"),
#        dpi = 600, limitsize = TRUE)  

QI_I <- surface_plot / volume_plot / weight_plot
# save plot
ggsave("out/QI_I.png", plot = QI_I,
       scale = 1, width = 16, height = 20, units = c("cm"),
       dpi = 600, limitsize = TRUE) 

QI_II <- necro_plot / polyps
# save plot
ggsave("out/QI_II.png", plot = QI_II,
       scale = 1, width = 16, height = 20, units = c("cm"),
       dpi = 600, limitsize = TRUE) 

 

## ---- 4.4. Photosynthetic efficiency -----------------------------------------
### --- 4.4.1. Effective quantum yield -----------------------------------------
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
  labs(color = expression(paste("Treatment ", mg, "·", L^-1)), fill = expression(paste("Treatment ", mg, "·", L^-1))) +
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


### --- 4.4.2. rETR ------------------------------------------------------------
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
       scale = 1, width = 12, height = 20, units = c("cm"),
       dpi = 600, limitsize = TRUE)  



## ---- 4.5. Correlation graph -------------------------------------------------
# ----- SURFACE
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

# write it into .csv
write_csv2(surface_all, "processed/surface_all.csv")

# create the base graph
surf_cor <- ggplot(surface_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Tissue growth (%)") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_text(face = "italic", size = 12),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

# create table with statistical results
annotation_surf_cor <- data.frame(
  parameter = factor(x = c("surface"), 
                     levels = c("surface")),
  value = c(75),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.0014, r = -0.3318 ", 
            "p = 0.2749, r = -0.122"))


# add statistics to graph
SURF <- surf_cor + geom_text(data = annotation_surf_cor,
                              mapping = aes(x = conc, y = value,
                                            label = label),
                              color = "black",
                              size = 3)
SURF

# ----- VOLUME
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

# write it into .csv
write_csv2(volume_all, "processed/volume_all.csv")

# create the base graph
vol_cor <- ggplot(volume_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = expression(paste("Volume growth ", cm^3, "·", cm^-2))) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

# create table with statistical results
annotation_vol_cor <- data.frame(
  parameter = factor(x = c("volume"), 
                     levels = c("volume")),
  value = c(0.15),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.0015, r = -0.3349", 
            "p = 0.0938, r = -0.1798"))

# add statistics to graph
VOL <- vol_cor + geom_text(data = annotation_vol_cor,
                     mapping = aes(x = conc, y = value,
                                   label = label),
                     color = "black",
                     size = 3)
VOL




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

# write it into .csv
write_csv2(calcification_all, "processed/calcification_all.csv")

calc_cor <- ggplot(calcification_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = expression(paste("Calcification ", mg, "·", cm^-2))) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_calc_cor <- data.frame(
  parameter = factor(x = c("calcification"), 
                     levels = c("calcification")),
  value = c(200),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.031, r = -0.2275", 
            "p = 0.0151, r = -0.2584"))


CALC <- calc_cor + geom_text(data = annotation_calc_cor,
                             mapping = aes(x = conc, y = value,
                                           label = label),
                             color = "black",
                             size = 3)
CALC


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

# write it into .csv
write_csv2(necrosis_all, "processed/necrosis_all.csv")


necro_cor <- ggplot(necrosis_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = expression(paste("Necrosis (%) "))) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_necro_cor <- data.frame(
  parameter = factor(x = c("necrosis"), 
                     levels = c("necrosis")),
  value = c(50),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.656, r = -0.1614",
            "p = 0.6963, r = 0.1022"))

NECRO <- necro_cor + geom_text(data = annotation_necro_cor,
                             mapping = aes(x = conc, y = value,
                                           label = label),
                             color = "black",
                             size = 3)
NECRO


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

# write it into .csv
write_csv2(polypactivity_all, "processed/polypactivity_all.csv")

poly_cor <- ggplot(polypactivity_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Polypactivity (mean)") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_poly_cor <- data.frame(
  parameter = factor(x = c("polypactivity"), 
                     levels = c("polypactivity")),
  value = c(1.1),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = <0.0001, r = -0.5625", 
            "p = 0.0025, r = -0.3146"))

POLY <- poly_cor + geom_text(data = annotation_poly_cor,
                           mapping = aes(x = conc, y = value,
                                         label = label),
                           color = "black",
                           size = 3)

POLY


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

# write it into .csv
write_csv2(YII_all, "processed/YII_all.csv")

YII_cor <- ggplot(YII_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Relative Y(II)") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_YII_cor <- data.frame(
  parameter = factor(x = c("YII"), 
                     levels = c("YII")),
  value = c(125),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.0347, r = 0.2229", 
            "p = 0.2962, r = 0.1113"))

YII_cor_plot <- YII_cor + geom_text(data = annotation_YII_cor,
                               mapping = aes(x = conc, y = value,
                                             label = label),
                               color = "black",
                               size = 3)
YII_cor_plot


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

# write it into .csv
write_csv2(FvFm_all, "processed/FvFm_all.csv")

FvFm_cor <- ggplot(FvFm_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Relative FvFm") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_FvFm_cor <- data.frame(
  parameter = factor(x = c("FvFm"), 
                     levels = c("FvFm")),
  value = c(350),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.0105, r = 0.2685", 
            "p = 0.0823, r = -0.1818"))

FvFm_cor_plot <- FvFm_cor + geom_text(data = annotation_FvFm_cor,
                                    mapping = aes(x = conc, y = value,
                                                  label = label),
                                    color = "black",
                                    size = 3)
FvFm_cor_plot


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

# write it into .csv
write_csv2(rETR_all, "processed/rETR_all.csv")

rETR_cor <- ggplot(rETR_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "rETRmax") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_rETR_cor <- data.frame(
  parameter = factor(x = c("rETRmax"), 
                     levels = c("rETRmax")),
  value = c(300),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.8817, r = -0.0159",
            "p = 0.7789, r = -0.03"))

rETR_cor_plot <- rETR_cor + geom_text(data = annotation_rETR_cor,
                                    mapping = aes(x = conc, y = value,
                                                  label = label),
                                    color = "black",
                                    size = 3)
rETR_cor_plot


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

# write it into .csv
write_csv2(Ek_all, "processed/Ek_all.csv")

Ek_cor <- ggplot(Ek_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Ek") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_blank(),
    legend.position = "none")

annotation_Ek_cor <- data.frame(
  parameter = factor(x = c("Ek"), 
                     levels = c("Ek")),
  value = c(550),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.7571, r = -0.0331", 
            "p = 0.9468, r = -0.0071"))

Ek_cor_plot <- Ek_cor + geom_text(data = annotation_Ek_cor,
                                    mapping = aes(x = conc, y = value,
                                                  label = label),
                                    color = "black",
                                    size = 3)
Ek_cor_plot


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

# write it into .csv
write_csv2(alpha_all, "processed/alpha_all.csv")

alpha_cor <- ggplot(alpha_all, aes(conc, value, color = as.factor(conc))) +
  facet_grid(parameter~spec, scales="free", 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  # stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Alpha") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 10),
    legend.position = "none")

ggplot(alpha_all, aes(x = conc, y = value)) +
  facet_grid(~ spec, 
             labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  geom_smooth(aes(x = conc, y = value)) +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log10') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = "Alpha") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 10),
    legend.position = "none")


annotation_alpha_cor <- data.frame(
  parameter = factor(x = c("alpha"), 
                     levels = c("alpha")),
  value = c(0.725),
  conc = c(15),
  spec = c("Pve", "Spi"),
  label = c("p = 0.6022, r = 0.0557", 
            "p = 0.6044, r = -0.0553"))

ALPHA <- alpha_cor + geom_text(data = annotation_alpha_cor,
                           mapping = aes(x = conc, y = value,
                                         label = label),
                           color = "black",
                           size = 3)

ALPHA





Correlation <- # bring all growth plots together
  SURF / VOL / CALC / NECRO / POLY / YII_cor_plot / FvFm_cor_plot / rETR_cor_plot / Ek_cor_plot / ALPHA

Correlation


# safe graph
ggsave("out/Correlation.png", plot = Correlation,
       scale = 1, width = 24, height = 42, units = c("cm"),
       dpi = 600, limitsize = TRUE)  



# ----- 5. Supplements  --------------------------------------------------------
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

Conc_curve

# safe graph
ggsave("out/Concentrations.png", plot = Conc_curve,
       scale = 1, width = 16, height = 18, units = c("cm"),
       dpi = 600, limitsize = TRUE)

# display relationship between added mg/L and measured ppl
Conc <- ggplot(Concentrations, aes(conc, per_L, color = as.factor(conc))) +
#  facet_grid(parameter~spec, scales="free", 
 #            labeller = labeller(spec = spec_labs)) +
  geom_point() +
  scale_color_manual(breaks = c("0", "0.1", "1", "10", "100"),
                     values = c("#4A8696", "#FFED85", "#E09F3E", "#9E2A2B","#540B0E"))  +
  stat_poly_line(color = "black") +
  scale_x_continuous(trans='log', labels= c("0", "0.1", "1", "10", "100"), breaks = c(0, 0.1, 1, 10, 100)) +
  scale_y_continuous(trans='log') +
  labs(x = expression(paste("Treatment ", mg, "·", L^-1)), 
       y = expression(paste("Particles ","·", L^-1))) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.background = element_rect(colour = "black"),
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    legend.position = "none")

cor.test(Concentrations$per_L, Concentrations$conc,  conf.level = 0.95, method = "pearson")
# OUTPUT: Pearson's product-moment correlation
# data:  test$value_log and test$conc
# t = 29.72, df = 7018, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.3134020 0.3549612
# sample estimates:
#   cor 
# 0.3343441 
# --> linear correlation

# Exponential
ks.test(Concentrations$per_L, Concentrations$conc, y = "pexp")
# OUTPUT: 	Asymptotic one-sample Kolmogorov-Smirnov test
# data:  Concentrations$per_L
# D = 0.96923, p-value < 2.2e-16
# alternative hypothesis: two-sided

annotation_conc <- data.frame(
  parameter = factor(x = c("conc"), 
                     levels = c("conc")),
  per_L = c(10000),
  conc = c(15),
  # spec = c("Pve", "Spi"),
  label = c("p < 0.0001, r = 0.3343"))

Conc_graph <- Conc + geom_text(data = annotation_conc,
                               mapping = aes(x = conc, y = per_L,
                                             label = label),
                               color = "black",
                               size = 3)

Conc_graph

# safe graph
ggsave("out/Concentrations_cor.png", plot = Conc_graph,
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

