#  This script corresponds to the section of the manuscript titled "Spatial Differences Among Walleye in 2017"

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

# Now breaking down the data to before and after the rainbow smelt crash
# Using 2009 and 2010 to represent years prior to the smelt crash, and 2017 and 2018 to represent years following the smelt crash
smelt_crash_index <- dplyr::filter(gill_index, 
                                   Year == "2009" |
                                     Year == "2010" |
                                     Year == "2017" |
                                     Year == "2018")

smelt_crash_index$Smelt_Crash <- ifelse(smelt_crash_index$Year == "2009", "pre",
                                        ifelse(smelt_crash_index$Year == "2010", "pre",
                                               ifelse(smelt_crash_index$Year == "2017", "post",
                                                      ifelse(smelt_crash_index$Year == "2018", "post", "NA"))))

post_smelt_crash_index <- dplyr::filter(smelt_crash_index, Smelt_Crash == "post")

# An initial basin-level model with large fish only >= 375 mm in fork length
basin_post_smelt_lm_2017 <- lm(log_Mass_kg ~ (log_ForkLength_mm * factor(Basin)) + factor(Sex) + factor(Mesh) + factor(Age), data = dplyr::filter(post_smelt_crash_index, Year == "2017", Length >= 375))
summary(basin_post_smelt_lm_2017)
pairs(emtrends(basin_post_smelt_lm_2017, "Basin", var = "log_ForkLength_mm"))
pairs(emmeans(basin_post_smelt_lm_2017, "Basin", var = "log_ForkLength_mm"))
sjstats::anova_stats(car::Anova(basin_post_smelt_lm_2017, type = 3)) %>% dplyr::select(1:7)

basin_length_mass_gill_net <- ggplot(basin_post_smelt_lm_2017$model, aes(x = log_ForkLength_mm, y = log_Mass_kg, color = basin_post_smelt_lm_2017$model$`factor(Basin)`, fill = basin_post_smelt_lm_2017$model$`factor(Basin)`, Group = basin_post_smelt_lm_2017$model$`factor(Basin)`)) +
  scale_color_manual(values = c("#4575b4", "#d73027", "#fdae61"), breaks = c("North", "Central", "South"), labels = c("North Basin", "Channel", "South Basin"), name = "Basin") +
  geom_point(size = 1.2, alpha = 0.5) +
  labs(x = bquote('Log'[10]~'fork length'), y =  bquote('Log'[10]~'mass'), fill = "Basin") +
  geom_smooth(linetype = 0, method = "lm") +
  scale_fill_manual(values = c("#4575b4", "#d73027", "#fdae61"), breaks = c("North", "Central", "South"), labels = c("North Basin", "Channel", "South Basin"), name = "Basin") +
  theme_bw() +
  theme(text = element_text(size=15), legend.justification = c("right", "bottom"), legend.position = c(0.99, 0.01), legend.background = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)), shape = guide_legend(override.aes = list(size = 2))) +
  xlim(2.57, 2.89) +
  ylim(-0.46, 0.72)
# Panel B of figure 5 in the manuscript
basin_length_mass_gill_net

# Repeating the length-mass model with the same terms for the year 2017, but including the Grand Rapids site and all walleye regardless of length. This gives a more overall look at how length and mass relate.
all_site_post_smelt_lm_2017 <- lm(log_Mass_kg ~ (log_ForkLength_mm * factor(Site)) + factor(Sex) + factor(Mesh) + factor(Age), data = dplyr::filter(post_smelt_crash_index, Year == "2017"))
summary(all_site_post_smelt_lm_2017)
sjstats::anova_stats(car::Anova(all_site_post_smelt_lm_2017, type = 3)) %>% dplyr::select(1:7)

View(dplyr::filter(post_smelt_crash_index, Year == "2017") %>% 
       group_by(Site) %>% 
       tally())


pairs(emtrends(all_site_post_smelt_lm_2017, "Site", var = "log_ForkLength_mm"))
pairs(emmeans(all_site_post_smelt_lm_2017, "Site", var = "log_ForkLength_mm"))

post_smelt_crash_index %>% 
  dplyr::filter(Year == "2017") %>%
  group_by(Basin) %>% 
  dplyr::tally()

# Plotting the linear model for all sites and walleye sizes from 2017. This is panel A of Figure 5 in the manuscript.
all_sites_lengthmass_2017 <- ggplot(all_site_post_smelt_lm_2017$model, aes(x = log_ForkLength_mm, y = log_Mass_kg, color = all_site_post_smelt_lm_2017$model$`factor(Site)`, fill = all_site_post_smelt_lm_2017$model$`factor(Site)`, Group = all_site_post_smelt_lm_2017$model$`factor(Site)`)) +
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach"), name = "Site") +
  scale_fill_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach"), name = "Site") +
  geom_point(size = 1.2, alpha = 0.5) +
  labs(x = bquote('Log'[10]~'fork length'), y =  bquote('Log'[10]~'mass'), fill = "Site") +
  geom_smooth(linetype = 0, method = "lm") +
  theme_bw() +
  theme(text = element_text(size=20), legend.justification = c("right", "bottom"), legend.position = c(0.99, 0.01), legend.background = element_blank()) 
all_sites_lengthmass_2017
# Use this line for the legend title if you want a subtitle for Panels A & B on a second line and smaller. Use Inkscape to left justify the subtitle
# expression('Site\n'*scriptstyle(italic('Panels A & B')))
nrow(dplyr::filter(post_smelt_crash_index, Length >= 375, Year == "2017"))


#############################################################
# Running a between-site linear model on the data used for metabolites to explore length-mass relationships
# Selecting only variables of interest from the metabolite sheet
metabolite_combined <- dplyr::select(metabolite_condition, -TotalLength_mm, -log_TotalLength_mm)
# Rename the year column to be consistent with the post smelt sheet
metabolite_combined <- dplyr::rename(metabolite_combined,
                                     Year = Year_tagged)
# Add in a column to show that these fish are from the sequencing study
metabolite_combined <- mutate(metabolite_combined,
                              Study = "sequenced")
# Add a column for Basin information, for combination with the gill net data
metabolite_combined <- mutate(metabolite_combined,
                              Basin = ifelse(Site == "Red_River", "South",
                                             ifelse(Site == "Matheson", "Central",
                                                    ifelse(Site == "Dauphin_River", "North", NA))))
# Reorder the columns to match the gill net data
metabolite_combined <- metabolite_combined[c( "ID","Year","Site","Basin","Sex","ForkLength_mm","Weight_kg","log_ForkLength_mm","log_Mass_kg","Study")]

# Run a linear model looking at length-mass relationships among walleye for whom metabolites were measured
metabolite_lm <- lm(log_Mass_kg ~ (log_ForkLength_mm * Site) + Sex + (log_ForkLength_mm * Sex), data = metabolite_combined)
summary(metabolite_lm)

emtrends(metabolite_lm, list(pairwise ~ Site), var = "log_ForkLength_mm", adjust = "tukey")

emmeans(metabolite_lm, list(pairwise ~ Site), adjust = "tukey")


sjstats::anova_stats(car::Anova(metabolite_lm, type = 3))  %>% dplyr::select(1:7)
# Looking at the linear prediction of log mass based on site. No real differences
emmip(metabolite_lm, log_ForkLength_mm ~ Site, CIs = TRUE)


# Panel C of figure 5 in the manuscript
length_mass_metabolite <- ggplot(metabolite_lm$model, aes(x = log_ForkLength_mm, y = log_Mass_kg, color = metabolite_lm$model$Site, fill = metabolite_lm$model$Site, Group = metabolite_lm$model$Site)) +
  scale_color_manual(values = c("#d73027", "#fdae61", "#4575b4"), breaks = c("Dauphin_River", "Matheson_Island", "Red_River"), labels = c("Dauphin River", "Matheson Island", "Red River"), name = "Site") +
  geom_point(size = 1.2, alpha = 0.5) +
  labs(x = bquote('Log'[10]~'fork length'), y = bquote('Log'[10]~'mass'), fill = "Site") +
  geom_smooth(linetype = 0, method = "lm") +
  scale_fill_manual(values = c("#d73027", "#fdae61", "#4575b4"), breaks = c("Dauphin_River", "Matheson_Island", "Red_River"), labels = c("Dauphin River", "Matheson Island", "Red River"), name = "Site") +
  theme_bw() +
  theme(text = element_text(size=15), legend.justification = c("right", "bottom"), legend.position = c(0.99, 0.01), legend.background = element_blank()) +
  xlim(2.57, 2.89) +
  ylim(-0.46, 0.72)
length_mass_metabolite

###############################################################
# Using a linear model to explore differences between the metabolite and gill net data at the basin-level
# Combine the data
combined_data <- rbind(post_smelt_combined, metabolite_combined)
basin_combined_lm <- lm(log_Mass_kg ~ log_ForkLength_mm + Basin + Study + (log_ForkLength_mm * Sex) + (log_ForkLength_mm * Basin) + (log_ForkLength_mm * Study), data = dplyr::filter(combined_data, Year == 2017, ForkLength_mm >= 375, Site != "Grand_Rapids"))
summary(basin_combined_lm)
display(basin_combined_lm)
eta_sq(basin_combined_lm)
# Table 5 in the manuscript
sjstats::anova_stats(car::Anova(basin_combined_lm, type = 3))  %>% dplyr::select(1:7)
# See if the p value is significantly different between the gill net and metabolite data
pairs(emmeans(basin_combined_lm, list(pairwise ~ Study), adjust = "tukey"))
pairs(emtrends(basin_combined_lm, list(pairwise ~ Study), var = "log_ForkLength_mm", adjust = "tukey"))

# Check how different the basins are 
emmeans(basin_combined_lm, list(pairwise ~ Basin), adjust = "tukey")

# Panel D of Figure 5 in the manuscript
length_mass_combined <- ggplot(basin_combined_lm$model, aes(x = log_ForkLength_mm, y = log_Mass_kg, color = Study, group = Study, fill = Study)) +
  scale_color_manual(name = "Data Set", labels = c("Gill Net Index", "Metabolite"), values = c("#4575b4", "#d73027")) +
  geom_point(size = 1.2, alpha = 0.5) +
  labs(x = bquote('Log'[10]~'fork length'), y =  bquote('Log'[10]~'mass')) + 
  geom_smooth(linetype = 0, method = "lm", se = TRUE) +
  scale_fill_manual(name = "Data Set", labels = c("Gill Net Index", "Metabolite"), values = c("#4575b4", "#d73027")) +
  theme_bw() +
  theme(text = element_text(size=15), legend.justification = c("right", "bottom"), legend.position = c(0.99, 0.01), legend.background = element_blank()) +
  xlim(2.57, 2.89) +
  ylim(-0.46, 0.72)
length_mass_combined


# Combining length-mass plots for gill net index and metabolite data
right_col <- plot_grid(basin_length_mass_gill_net, length_mass_metabolite, length_mass_combined, labels = c('B', 'C', 'D'), label_size = 12, nrow = 3, rel_widths = c(1, 1, 1))
# Figure 5 in the manuscript with all panels put together
plot_grid(all_sites_lengthmass_2017, right_col, label_size = 12, ncol = 2, labels = c('A'))
