# Lake Winnipeg walleye length-mass relationships over time
# This corresponds to the manuscript section "Lake Winnipeg Walleye Length-Mass Relationships"

# Note, not all of these libraries may have been used in this script.
library(tidyverse)
library(WRS2)
library(car)
library(emmeans)
library(broom)
library(sjstats)
library(cowplot)
library(jtools)

# Looking at the gill net index data. A much bigger file
# Downloaded from: https://www.gov.mb.ca/sd/fish_and_wildlife/fish/commercial_fishing/netting_data.html
# Sauger and dwarf walleye were removed prior to analysis in R
gill_index <- read_csv("lake_wgp_index_gillnetting2009_2018_KMJ.csv", col_names = TRUE, col_types = "ncffccffcnnffnn")

# Take the log of each weight and length
gill_index <- gill_index %>% 
  mutate(log_ForkLength_mm = log10(Length)) %>% 
  mutate(log_Mass_kg = log10(Weight/1000))

# Reordering gill index data by year
gill_index <- gill_index[order(gill_index$Year),]

# Turning year into a factor for running models
gill_index$Year <- as.factor(gill_index$Year)

# Controlling for the effect of Mesh
Mesh_cov_year <- lm(log_Mass_kg ~ log_ForkLength_mm * (factor(Year) + factor(Site) + factor(Sex)) + factor(Mesh), data = dplyr::filter(gill_index, Length >= 375))

View(dplyr::filter(gill_index, Length >= 375) %>% 
       group_by(Site) %>% 
       tally())

# This is Supplementary Table S1 in the manuscript
Mesh_cov_year
summary(Mesh_cov_year)
summ(Mesh_cov_year, confint = TRUE)
eta_sq(Mesh_cov_year)
sjstats::anova_stats(car::Anova(Mesh_cov_year, type = 3)) %>% dplyr::select(1:7)
emmeans::pmmeans(Mesh_cov_year, "Mesh")
emmeans::pmmeans(Mesh_cov_year, "log_ForkLength_mm", by = "Year")

# This is Supplementary Figure S1 in the manuscript
emmeans(Mesh_cov_year, "Site", by = "Year")
emmeans(Mesh_cov_year, "Site" ~ "Year")
emmip(Mesh_cov_year, Site ~ Year, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Year") + 
  ylab(bquote('Log'[10]~'mass')) +
  theme(text = element_text(size=20))

