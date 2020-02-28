# This script corresponds to the section of the manuscript titled "Lake Winnipeg Walleye Length-at-Age Over Time"
# It first looks at a comparison in the years 2015, 2016, 2017, and 2018
# Then looks at each of 6 sites individually, between the years 20012 and 2018
# Figures 3 and 4 are the outputs of this script. 

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


# Running linear models with the gill net index data to see how size at age differs between the north, central, and south basins
# Using length as the dependent variable because while mass may shrink from a fish starving, length wouldn't and we are interested at size at age
# Turning the gill index age into a numeric and not a factor so it's a covariate and we lose fewer degrees of freedom
gill_index$Age <- as.numeric(gill_index$Age)

# This line shows that ages 9 and above have 14 or fewer fish represented, while ages 8 and below have more 30 in all cases except age 1. Re-running models on just ages where we have enough data, ages 2 - 6
dplyr::filter(gill_index, Year == c("2018")) %>% 
  group_by(Age) %>% 
  tally()

# Only taking the ages 2 through 6 years old, since all other ages have < 20 observations and are likely unreliable, although two Dauphin River age 2 year classes still have low sample sizes with even this limited filtering. 
gill_sufficient_age_data <- gill_index[gill_index$Age == "2" | gill_index$Age == "3" | gill_index$Age == "4" | gill_index$Age == "5" | gill_index$Age == "6",]


# Repeating the linear model but in a site-wise comparison, instead

site_age_index_2017 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2017")))
nrow(dplyr::filter(gill_sufficient_age_data, Year == c("2017")))
summary(site_age_index_2017)
age_2017_label <- c("italic('F')==48.50","italic('p')<2.2%*%10^-16","Adjusted","italic(R)^{2}==0.67")
eta_sq(site_age_index_2017)
pmmeans(site_age_index_2017, "Age", by = "Site")



pairs(emmeans(site_age_index_2017, ~ Site | Age))
pairs(emtrends(site_age_index_2017, ~ Site, var = "Age"))

pairs(emmeans(site_age_index_2017, ~ Site | Age))

emmip(site_age_index_2017, Site ~ Age, CIs = TRUE)
length_2017 <- emmip(site_age_index_2017, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  ylim(2.20, 2.753) 
length_2017
# For these colors, the order of how they appear in scale color manual is 1=Grand Rapids, 2=Riverton, 3=Frog Bay, 4=Grand Beach, 5=Matheson Island,6=Dauphin River


site_age_index_2009 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2009")))
summary(site_age_index_2009)
eta_sq(site_age_index_2009)
pmmeans(site_age_index_2009, "Age", by = "Site")


pairs(emmeans(site_age_index_2009, ~ Site | Age))
# This will be figure 4 in the MS
emmip(site_age_index_2009, Site ~ Age, CIs = TRUE)
emmip(site_age_index_2009, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14))


site_age_index_2018 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2018")))
nrow(dplyr::filter(gill_sufficient_age_data, Year == c("2018")))
summary(site_age_index_2018)
eta_sq(site_age_index_2018)
pairs(pmmeans(site_age_index_2018, "Site", by = "Age"))


pairs(emmeans(site_age_index_2018, ~ Site | Age))
pairs(emtrends(site_age_index_2018, ~ Site, var = "Age"))


pairs(emmeans(site_age_index_2018, ~ Site | Age))
# This will be figure 4 in the MS
emmip(site_age_index_2018, Site ~ Age, CIs = TRUE)
length_2018 <- emmip(site_age_index_2018, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  ylim(2.20, 2.753)
length_2018


site_age_index_2016 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2016")))
nrow(dplyr::filter(gill_sufficient_age_data, Year == c("2016")))
summary(site_age_index_2016)
eta_sq(site_age_index_2016)
pmmeans(site_age_index_2016, "Age", by = "Site")


pairs(emmeans(site_age_index_2016, ~ Site | Age))
# This will be figure 4 in the MS
emmip(site_age_index_2016, Site ~ Age, CIs = TRUE)
length_2016 <- emmip(site_age_index_2016, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  ylim(2.20, 2.753)



site_age_index_2015 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2015")))
nrow(dplyr::filter(gill_sufficient_age_data, Year == c("2015")))
summary(site_age_index_2015)
eta_sq(site_age_index_2015)
pmmeans(site_age_index_2015, "Age", by = "Site")


pairs(emmeans(site_age_index_2015, ~ Site | Age))
pairs(emtrends(site_age_index_2015, ~ Site, var = "Age"))

# This will be figure 4 in the MS
emmip(site_age_index_2015, Site ~ Age, CIs = TRUE)
length_2015 <- emmip(site_age_index_2015, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  ylim(2.20, 2.753)
length_2015


site_age_index_2014 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2014")))
nrow(dplyr::filter(gill_sufficient_age_data, Year == c("2014")))
summary(site_age_index_2014)
eta_sq(site_age_index_2014)
pmmeans(site_age_index_2014, "Age", by = "Site")


pairs(emmeans(site_age_index_2014, ~ Site | Age))
pairs(emtrends(site_age_index_2014, ~ Site, var = "Age"))

# This will be figure 4 in the MS
emmip(site_age_index_2014, Site ~ Age, CIs = TRUE)
length_2014 <- emmip(site_age_index_2014, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  ylim(2.20, 2.753)
length_2014


site_age_index_2013 <- lm(log_ForkLength_mm ~ (factor(Site) * factor(Age)) + factor(Sex) + factor(Mesh), data = dplyr::filter(gill_sufficient_age_data, Year == c("2013")))
nrow(dplyr::filter(gill_sufficient_age_data, Year == c("2013")))
summary(site_age_index_2013)
eta_sq(site_age_index_2013)
pmmeans(site_age_index_2013, "Age", by = "Site")


pairs(emmeans(site_age_index_2013, ~ Site | Age))
pairs(emtrends(site_age_index_2013, ~ Site, var = "Age"))

# This will be figure 4 in the MS
emmip(site_age_index_2013, Site ~ Age, CIs = TRUE)
length_2013 <- emmip(site_age_index_2013, Site ~ Age, CIs = TRUE) + 
  scale_color_manual(values = c("#91bfdb", "#a50f15", "#fdae61", "#de2d26", "#fc8d59", "#4575b4"), breaks = c("Grand_Rapids", "Dauphin_River", "Matheson_Island", "Frog_Bay", "Riverton", "Grand_Beach"), labels = c("Grand Rapids", "Dauphin River", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach")) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  ylim(2.20, 2.753)
length_2013


legend <- get_legend(length_2015)

# Figure 4 in the manuscript
length_plots <- plot_grid(length_2015 + theme(legend.position="none"), length_2016 + theme(legend.position="none"), length_2017 + theme(legend.position="none"), length_2018 + theme(legend.position="none"), labels = c('2015', '2016', '2017', '2018'), label_size = 17, label_x = 0.17, label_y = 0.97, hjust = -0.1, vjust = 1.7, ncol = 2)
plot_grid(length_plots, legend, rel_widths = c(1, 0.2))

#########################################################################################################################################################


# Running a similar analysis of length over time, but with just the Dauphin River fish
dauphin_river <- dplyr::filter(gill_index, Site == "Dauphin_River")

# Checking how much data there is to work with
View(dauphin_river %>% 
       group_by(Year, Age) %>% 
       tally())
nrow(dauphin_river)
dauphin_river <- dplyr::filter(dauphin_river, Year != "2009", Year != "2010", Year != "2011", Age > 1, Age <= 6)

# Running a single linear model with age and year interacting
dauphin_age_year <- lm(log_ForkLength_mm ~ (factor(Age) * factor(Year)) + factor(Sex) + factor(Mesh), data = dauphin_river)
summary(dauphin_age_year)
pairs(emmeans(dauphin_age_year, ~ Year | Age))
emmip(dauphin_age_year, Year ~ Age, CIs = TRUE)

dauphin_plot <- emmip(dauphin_age_year, Year ~ Age, CIs = TRUE) + 
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  scale_color_brewer() +
  ylim(2.20, 2.753) 

dauphin_plot

# Following the same process for the Riverton site
riverton <- dplyr::filter(gill_index, Site == "Riverton")
riverton <- dplyr::filter(riverton, Year != "2009", Year != "2010", Year != "2011", Age > 1, Age <= 6)
nrow(riverton)
riverton_age_year <- lm(log_ForkLength_mm ~ (factor(Age) * factor(Year)) + factor(Sex) + factor(Mesh), data = riverton)
summary(riverton_age_year)
emmip(riverton_age_year, Year ~ Age, CIs = TRUE)


riverton_plot <- emmip(riverton_age_year, Year ~ Age, CIs = TRUE) + 
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  scale_color_brewer(palette = "Reds") +
  ylim(2.20, 2.753) 
riverton_plot

dauphin_plot2 <- emmip(dauphin_age_year, Year ~ Age, CIs = TRUE, plotit = FALSE)
dauphin_plot2

riverton_plot2 <- emmip(riverton_age_year, Year ~ Age, CIs = TRUE, plotit = FALSE)
riverton_plot2

ggplot(data = dauphin_plot2, aes(x=xvar,
                                 y=yvar,
                                 group=Year,
                                 color=Year)) + 
  geom_point(data = dauphin_plot2, size = 0.5) +
  geom_line(data = dauphin_plot2, size = 0.75) +
  geom_linerange(data = dauphin_plot2, aes(ymin = LCL, ymax = UCL), size = 1.25, alpha = 0.5) +
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) + 
  scale_color_brewer() +
  theme(text = element_text(size=20)) +
  geom_line(data = riverton_plot2, aes(x=xvar, y=yvar, group=Year, color = Year), size = 0.75, color = rep(c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d"), 5), inherit.aes = FALSE, alpha = 0.4) +
  geom_point(data = riverton_plot2, size = 0.5, alpha = 0.4) +
  geom_linerange(data = riverton_plot2, aes(ymin = LCL, ymax = UCL), size = 1.25, alpha = 0.4) 


# Trying Matheson Island with the same process
matheson <- dplyr::filter(gill_index, Site == "Matheson_Island")
matheson <- dplyr::filter(matheson, Year != "2009", Year != "2010", Year != "2011", Age > 1, Age <= 6)
nrow(matheson)
matheson_age_year <- lm(log_ForkLength_mm ~ (factor(Age) * factor(Year)) + factor(Sex) + factor(Mesh), data = matheson)
summary(matheson_age_year)

emmip(matheson_age_year, Year ~ Age, CIs = TRUE)
matheson_plot <- emmip(matheson_age_year, Year ~ Age, CIs = TRUE) + 
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  scale_color_brewer(palette = "Oranges") +
  ylim(2.20, 2.753)
matheson_plot



# Grand Rapids
grand_rapids <- dplyr::filter(gill_index, Site == "Grand_Rapids")
grand_rapids <- dplyr::filter(grand_rapids, Year != "2009", Year != "2010", Year != "2011", Age > 1, Age <= 6)
nrow(grand_rapids)
grand_rapids_age_year <- lm(log_ForkLength_mm ~ (factor(Age) * factor(Year)) + factor(Sex) + factor(Mesh), data = grand_rapids)
summary(grand_rapids_age_year)

emmip(grand_rapids_age_year, Year ~ Age, CIs = TRUE)
grand_rapids_plot <- emmip(grand_rapids_age_year, Year ~ Age, CIs = TRUE) + 
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  scale_color_brewer(palette = "Greens") +
  ylim(2.20, 2.753)
grand_rapids_plot

# Frog Bay
frog_bay <- dplyr::filter(gill_index, Site == "Frog_Bay")
frog_bay <- dplyr::filter(frog_bay, Year != "2009", Year != "2010", Year != "2011", Age > 1, Age <= 6)
nrow(frog_bay)
frog_bay_age_year <- lm(log_ForkLength_mm ~ (factor(Age) * factor(Year)) + factor(Sex) + factor(Mesh), data = frog_bay)
summary(frog_bay_age_year)

emmip(frog_bay_age_year, Year ~ Age, CIs = TRUE)
frog_bay_plot <- emmip(frog_bay_age_year, Year ~ Age, CIs = TRUE) + 
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  scale_color_brewer(palette = "YlGnBu") +
  ylim(2.20, 2.753)
frog_bay_plot


# Grand Beach
grand_beach <- dplyr::filter(gill_index, Site == "Grand_Beach")
grand_beach <- dplyr::filter(grand_beach, Year != "2009", Year != "2010", Year != "2011", Age > 1, Age <= 6)
nrow(grand_beach)
grand_beach_age_year <- lm(log_ForkLength_mm ~ (factor(Age) * factor(Year)) + factor(Sex) + factor(Mesh), data = grand_beach)
summary(grand_beach_age_year)

emmip(grand_beach_age_year, Year ~ Age, CIs = TRUE)
grand_beach_plot <- emmip(grand_beach_age_year, Year ~ Age, CIs = TRUE) + 
  theme_bw() +
  xlab("Age") + 
  ylab(bquote('Log'[10]~'fork length')) +
  theme(text = element_text(size=14)) +
  scale_color_brewer(palette = "Purples") +
  ylim(2.20, 2.753)
grand_beach_plot

plot_grid(dauphin_plot, grand_rapids_plot, matheson_plot, frog_bay_plot, riverton_plot, grand_beach_plot, labels = c("Dauphin River", "Grand Rapids", "Matheson Island", "Frog Bay", "Riverton", "Grand Beach"), label_x = 0.12, label_y = 0.95, hjust = -0.1, vjust = 1.7, ncol = 2, greedy = TRUE)
# A is Dauphin, B is Grand Rapids, C is Matheson, D is Frog Bay, E is Riverton, F is Grand Beach

# Figure 3 in the manuscript
plot_grid(dauphin_plot, grand_rapids_plot, matheson_plot, frog_bay_plot, riverton_plot, grand_beach_plot, labels = NA, ncol = 2, greedy = FALSE)