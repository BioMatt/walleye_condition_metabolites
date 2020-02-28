# Spatial differences in 9 metabolites among Lake Winnipeg walleye from 3 sites in 2017
# This script corresponds to the section of the manuscript titled "Modeling Metabolite Differences". At the end is the method used to create Figure 6 with the predicted values for each metabolite by site. 

# Note, not all of these libraries may have been used in this script.
library(tidyverse)
library(WRS2)
library(car)
library(emmeans)
library(broom)
library(sjstats)
library(cowplot)
library(jtools)

# Read the metabolite data
select_metabolites <- read_delim("select_metabolites.txt", delim = "\t")


# An ANCOVA is not appropriate because both the metabolite (methionine) and the independent variable, fork length, are continuous. Plus, the interaction term is significant. So trying regular linear models instead, taking sex and site into account
methionine_avg_lm <- lm(Methionine_Average ~ log_ForkLength_mm + Site, data = select_metabolites)


summary(methionine_avg_lm)
eta_sq(methionine_avg_lm)
sjstats::anova_stats(car::Anova(methionine_avg_lm, type = 3))  %>% dplyr::select(1:7)

methionine_emmeans <- emmeans(methionine_avg_lm, pairwise ~ Site, adjust = "tukey")

methionine_emmeans
plot(methionine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c("blue", "green", "red"), ylab = "Average Methionine Estimated Marginal Mean", PIs = FALSE)


# The linear model for tryptophan
tryptophan_lm <- lm(Tryptophan ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(tryptophan_lm)
eta_sq(tryptophan_lm)

sjstats::anova_stats(car::Anova(tryptophan_lm, type = 3))  %>% dplyr::select(1:7)



tryptophan_emmeans <- emmeans(tryptophan_lm, pairwise ~ Site, adjust = "tukey")
tryptophan_emmeans
plot(tryptophan_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Average Tryptophan Estimated Marginal Mean", PIs = TRUE)


# Plotting points for Tryptophan and predicted levels for each site
ggplot(select_metabolites, aes(x=Site, y=Tryptophan, colour=Site)) +
  geom_point(size=1) +
  geom_line(aes(y=predict(tryptophan_lmer), group=Site, size="Site"), size = 1) +
  scale_size_manual(name="Predictions", values=c("Site"=0.5, "Site"=3)) +
  theme_bw(base_size=22) 


# The linear model for lysine
lysine_avg_lm <- lm(Lysine_Average ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(lysine_avg_lm)
eta_sq(lysine_avg_lm)

sjstats::anova_stats(car::Anova(lysine_avg_lm, type = 3))  %>% dplyr::select(1:7)


summary(lysine_avg_lm)
eta_sq(lysine_avg_lm)

lysine_emmeans <- emmeans(lysine_avg_lm, pairwise ~ Site, adjust = "tukey")
lysine_emmeans
plot(lysine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Average Lysine Estimated Marginal Mean", PIs = TRUE)


# Plotting points for Tryptophan and predicted levels for each site
ggplot(select_metabolites, aes(x=Site, y=Lysine_Average, colour=Site)) +
  geom_point(size=1) +
  geom_line(aes(y=predict(lysine_avg_lmer), group=Site, size="Site")) +
  scale_size_manual(name="Predictions", values=c("Site"=0.5, "Site"=3)) +
  theme_bw(base_size=22) 


# The linear model for trans hydroxyproline
trans_Hydroxyproline_lm <- lm(trans_Hydroxyproline ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(trans_Hydroxyproline_lm)
eta_sq(trans_Hydroxyproline_lm)

sjstats::anova_stats(car::Anova(trans_Hydroxyproline_lm, type = 3))  %>% dplyr::select(1:7)


trans_Hydroxyproline_emmeans <- emmeans(trans_Hydroxyproline_lm, pairwise ~ Site, adjust = "tukey")
trans_Hydroxyproline_emmeans
plot(trans_Hydroxyproline_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Average trans-Hydroxyproline Estimated Marginal Mean", PIs = TRUE)

# Plotting points for Tryptophan and predicted levels for each site
ggplot(select_metabolites, aes(x=Site, y=trans_Hydroxyproline, colour=Site)) +
  geom_point(size=1) +
  geom_line(aes(y=predict(trans_hydroxyproline_lmer), group=Site, size="Site")) +
  scale_size_manual(name="Predictions", values=c("Site"=0.5, "Site"=3)) +
  theme_bw(base_size=22) 


# The linear model for dimethylarginine
total_dimethylarginine_lm <- lm(total_Dimethylarginine ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(total_dimethylarginine_lm)
eta_sq(total_dimethylarginine_lm)

sjstats::anova_stats(car::Anova(total_dimethylarginine_lm, type = 3))  %>% dplyr::select(1:7)



total_dimethylarginine_emmeans <- emmeans(total_dimethylarginine_lm, pairwise ~ Site, adjust = "tukey")
total_dimethylarginine_emmeans
plot(total_dimethylarginine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Average Dimethylarginine Estimated Marginal Mean", PIs = TRUE)

# Plotting points for dmethylarginine and predicted levels for each site
ggplot(select_metabolites, aes(x=Site, y=total_Dimethylarginine, colour=Site)) +
  geom_point(size=1) +
  geom_line(aes(y=predict(total_dimethylarginine_lmer), group=Site, size="Site")) +
  scale_size_manual(name="Predictions", values=c("Site"=0.5, "Site"=3)) +
  theme_bw(base_size=22) 

# The linear model for methionine sulfoxide
methionine_sulfoxide_lm <- lm(Methionine_Sulfoxide ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(methionine_sulfoxide_lm)
eta_sq(methionine_sulfoxide_lm)

sjstats::anova_stats(car::Anova(methionine_sulfoxide_lm, type = 3))  %>% dplyr::select(1:7)

methionine_sulfoxide_emmeans <- emmeans(methionine_sulfoxide_lm, pairwise ~ Site, adjust = "tukey")
methionine_sulfoxide_emmeans
plot(methionine_sulfoxide_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Methionine Sulfoxide Estimated Marginal Mean", PIs = TRUE) + theme_bw()

# The linear model for dimethyl sulfone
dimethyl_sulfone_lm <- lm(dimethyl_sulfone ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(dimethyl_sulfone_lm)
eta_sq(dimethyl_sulfone_lm)

sjstats::anova_stats(car::Anova(dimethyl_sulfone_lm, type = 3))  %>% dplyr::select(1:7)

test <- lm(log_ForkLength_mm ~ Site, data = select_metabolites)
summary(test)

emtrends(dimethyl_sulfone_lm, pairwise ~ Site, var = "log_ForkLength_mm",adjust = "tukey")
dimthyl_sulfone_emmeans <- emmeans(dimethyl_sulfone_lm, pairwise ~ Site, adjust = "tukey")

dimthyl_sulfone_emmeans
plot(dimthyl_sulfone_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Methionine Sulfoxide Estimated Marginal Mean", PIs = TRUE) + theme_bw()

# The linear model for kynurenine
Kynurenine_lm <- lm(Kynurenine ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(Kynurenine_lm)
eta_sq(Kynurenine_lm)


sjstats::anova_stats(car::Anova(Kynurenine_lm, type = 3))  %>% dplyr::select(1:7)

emmeans(Kynurenine_lm, pairwise ~ Site, adjust = "tukey")
emtrends(Kynurenine_lm, pairwise ~ Site, var = "log_ForkLength_mm",adjust = "tukey")

kynurenine_emmeans <- emmeans(Kynurenine_lm, pairwise ~ Site, adjust = "tukey")

kynurenine_emmeans
plot(kynurenine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "Kynurenine Estimated Marginal Mean", PIs = TRUE) + theme_bw()

# The linear model for alpha aminoadipic acid
alpha_aminoadipic_acid_lm <- lm(alpha_aminoadipic_acid ~ log_ForkLength_mm + Site, data = select_metabolites)
summary(alpha_aminoadipic_acid_lm)
eta_sq(alpha_aminoadipic_acid_lm)

sjstats::anova_stats(car::Anova(alpha_aminoadipic_acid_lm, type = 3))  %>% dplyr::select(1:7)

emmeans(alpha_aminoadipic_acid_lm, pairwise ~ Sex, adjust = "tukey")
emtrends(alpha_aminoadipic_acid_lm, pairwise ~ Site, var = "log_ForkLength_mm",adjust = "tukey")

alpha_aminoadipic_acid_lm_emmeans <- emmeans(alpha_aminoadipic_acid_lm, pairwise ~ Site, adjust = "tukey")

alpha_aminoadipic_acid_lm_emmeans
plot(alpha_aminoadipic_acid_lm_emmeans, horizontal = FALSE, comparisons = TRUE, colors = "blue", ylab = "alpha aminoadipic acid Estimated Marginal Mean", PIs = TRUE) + theme_bw()


# Collecting all 6 plots in the same place
avg_methionine_plot <- plot(methionine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Methionine", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))
tryptophan_plot <- plot(tryptophan_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Tryptophan", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))
lysine_plot <- plot(lysine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Lysine", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))

trans_hydroxyproline_plot <- plot(trans_Hydroxyproline_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Hydroxyproline", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))
dimethylarginine_plot <- plot(total_dimethylarginine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Dimethylarginine", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))
methionine_sulfoxide_plot <- plot(methionine_sulfoxide_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Methionine Sulfoxide", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))

dimethyl_sulfone_plot <- plot(dimthyl_sulfone_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Dimethyl Sulfone", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))
kynurenine_plot <- plot(kynurenine_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = "Kynurenine", PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))
kynurenine_plot
alpha_aminoadipic_acid_plot <- plot(alpha_aminoadipic_acid_lm_emmeans, horizontal = FALSE, comparisons = TRUE, colors = c(NA, "#4575b4", "#d73027", "black"), ylab = bquote(alpha*'-aminoadipic acid'), PIs = FALSE) + theme_bw() + theme(text = element_text(size=14))

# Figure 6 in the manuscript.
plot_grid(avg_methionine_plot, tryptophan_plot, lysine_plot, dimethyl_sulfone_plot, kynurenine_plot, alpha_aminoadipic_acid_plot, trans_hydroxyproline_plot, dimethylarginine_plot, methionine_sulfoxide_plot, labels = "AUTO", label_size = 12)



