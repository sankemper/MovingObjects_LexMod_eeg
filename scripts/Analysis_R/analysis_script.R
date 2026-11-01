################################################################################
## MASTER'S THESIS - Statistical testing: Deconvolution data 
## Santina Kemper 
## Last Update: 19.10.2026
## V1
################################################################################

## packages
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(mgcv)
library(gam)
library(readr)
library(itsadug)
library(purrr)
library(plotfunctions)
library(ggplot2)
library(patchwork)
library(performance)
library(insight)
library(see)
library(nlme)
library(tidyr)
library(tidyverse)
library(knitr)
library(kableExtra)
library(tidymodels)

################################################################################
## ============================= read in data ==================================
################################################################################
path <- "Analysis/Unfold_preprocessed/unfold_gamm_long.csv"
data <- read_csv(path)

noun <- data %>%
  filter(event == "N1_crit")

det <- data %>%
  filter(event == "DET1_crit")


adj <- data %>% 
  filter(event == "ADJ1_crit")

## ============================ GLOBAL FUNCTION ================================

theme_house <- function(base_size = 20) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold", colour = "#3a86d4"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

red_dark <- "#B2182B"
red_light <- "#EF8A8A"
blue_dark <- "#2166AC"
blue_light <- "#67A9CF"
## ============================ GRAND AVERAGE ==================================
# windows.N400 = [100 500];
# windows.Nref = [300 800];
# windows.P600 = [500 800];
# windows.ELAN = [100 300];
# roi.N400 = {'Cz','C3','C4','CP1','CP2','Pz','P3','P4'}; 
# roi.Nref = {'F3','FC1','FC5','C3'};
# roi.P600 = {'Pz','P3','P4','CP1','CP2'};
# roi.ELAN = {'Fp1','F3','FC5','FC1'};

roi_N400 = c('Cz','C3','C4','CP1','CP2','Pz','P3','P4') 
roi_Nref = c('F3','FC1','FC5','C3')
roi_P600 = c('Pz','P3','P4','CP1','CP2')
roi_ELAN = c('Fp1','F3','FC5','FC1')

det_N400 <- det %>%
  filter(channel %in% roi_N400)
det_Nref <- det %>%
  filter(channel %in% roi_Nref)
det_P600 <- det %>%
  filter(channel %in% roi_P600)
det_ELAN <- det %>%
  filter(channel %in% roi_ELAN)
n_N400 <- noun %>%
  filter(channel %in% roi_N400)
n_Nref <- noun %>%
  filter(channel %in% roi_Nref)
n_P600 <- noun %>%
  filter(channel %in% roi_P600)
n_ELAN <- noun %>%
  filter(channel %in% roi_ELAN)
adj_N400 <- adj %>%
  filter(channel %in% roi_N400)
adj_Nref <- adj %>%
  filter(channel %in% roi_Nref)
adj_P600 <- adj %>%
  filter(channel %in% roi_P600)
adj_ELAN <- adj %>%
  filter(channel %in% roi_ELAN)

all_N400 <- rbind(det_N400, adj_N400, n_N400)
all_P600 <- rbind(det_P600, adj_P600, n_P600)
all_Nref <- rbind(det_Nref, adj_Nref, n_Nref)
all_ELAN <- rbind(det_ELAN, adj_ELAN, n_ELAN)

event <- c("DET1_crit", "N1_crit", "ADJ1_crit")

grand_roi_N400 <- all_N400 %>%
  group_by(event, informative, mention, modified, subject, time) %>%
  summarise(
    subj_mean = mean(amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(event, informative, mention, modified, time) %>%
  summarise(
    mean_roi = mean(subj_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    event = factor(
      event,
      levels = c("DET1_crit", "ADJ1_crit", "N1_crit"),
      labels = c("DET", "ADJ", "N")
    ),
    cond = interaction(
      factor(mention, levels = c(1, 2)),
      factor(modified, levels = c(0, 1)),
      sep = "_mod"
    )
  )

cond_colors <- c(
  "1_mod0" = "#2166AC", # Mention 1, unmodified:  blue
  "2_mod0" = "#67A9CF",  # Mention 2, unmodified: light blue
  "1_mod1" = "#B2182B",  # Mention 1, modified: red 
  "2_mod1" = "#EF8A8A"   # Mention 2, modified: lighter red
)

p_grand_N400 <- ggplot(grand_roi_N400, aes(x = time, y = mean_roi, color = cond)) +
  annotate(
    "rect",
    xmin = 100, xmax = 500,
    ymin = -Inf, ymax = Inf,
    fill = "lightgrey", alpha = 0.25
  ) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_y_reverse() +
  scale_color_manual(
    name = "Mention x Modification",
    values = cond_colors,
    breaks = c("1_mod0", "2_mod0", "1_mod1", "2_mod1"),
    labels = c(
      "1_mod0" = "1mention_unmodified",
      "2_mod0" = "2mention_unmodified",
      "1_mod1" = "1mention_modified",
      "2_mod1" = "2mention_modified"
    ),
    drop = TRUE
  ) +
  facet_grid(
    informative ~ event,
    labeller = label_both,
    drop = FALSE
  ) +
  labs(
    x = "time (ms)",
    y = "amplitude (\u00b5V)",
    color = "Mention_Modification"
  ) +
  theme_house()

ggsave("Grand_Average/grand_average_N400.png", p_grand_N400, width = 12, height = 6, dpi = 300)

grand_roi_Nref <- all_Nref %>%
  group_by(event, informative, mention, modified, subject, time) %>%
  summarise(
    subj_mean = mean(amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(event, informative, mention, modified, time) %>%
  summarise(
    mean_roi = mean(subj_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    event = factor(
      event,
      levels = c("DET1_crit", "ADJ1_crit", "N1_crit"),
      labels = c("DET", "ADJ", "N")
    ),
    cond = interaction(
      factor(mention, levels = c(1, 2)),
      factor(modified, levels = c(0, 1)),
      sep = "_mod"
    )
  )

p_grand_Nref <- ggplot(grand_roi_Nref, aes(x = time, y = mean_roi, color = cond)) +
  annotate(
    "rect",
    xmin = 300, xmax = 800,
    ymin = -Inf, ymax = Inf,
    fill = "lightgrey", alpha = 0.25
  ) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_y_reverse() +
  scale_color_manual(
    name = "Mention x Modification",
    values = cond_colors,
    breaks = c("1_mod0", "2_mod0", "1_mod1", "2_mod1"),
    labels = c(
      "1_mod0" = "1mention_unmodified",
      "2_mod0" = "2mention_unmodified",
      "1_mod1" = "1mention_modified",
      "2_mod1" = "2mention_modified"
    ),
    drop = TRUE
  ) +
  facet_grid(
    informative ~ event,
    labeller = label_both,
    drop = FALSE
  ) +
  labs(
    x = "time (ms)",
    y = "amplitude (\u00b5V)",
    color = "Mention_Modification"
  ) +
  theme_house()

ggsave("Grand_Average/grand_average_Nref.png", p_grand_Nref, width = 12, height = 6, dpi = 300)


grand_roi_P600 <- all_P600 %>%
  group_by(event, informative, mention, modified, subject, time) %>%
  summarise(
    subj_mean = mean(amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(event, informative, mention, modified, time) %>%
  summarise(
    mean_roi = mean(subj_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    event = factor(
      event,
      levels = c("DET1_crit", "ADJ1_crit", "N1_crit"),
      labels = c("DET", "ADJ", "N")
    ),
    cond = interaction(
      factor(mention, levels = c(1, 2)),
      factor(modified, levels = c(0, 1)),
      sep = "_mod"
    )
  )

p_grand_P600 <- ggplot(grand_roi_P600, aes(x = time, y = mean_roi, color = cond)) +
  annotate(
    "rect",
    xmin = 500, xmax = 800,
    ymin = -Inf, ymax = Inf,
    fill = "lightgrey", alpha = 0.25
  ) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_y_reverse() +
  scale_color_manual(
    name = "Mention x Modification",
    values = cond_colors,
    breaks = c("1_mod0", "2_mod0", "1_mod1", "2_mod1"),
    labels = c(
      "1_mod0" = "1mention_unmodified",
      "2_mod0" = "2mention_unmodified",
      "1_mod1" = "1mention_modified",
      "2_mod1" = "2mention_modified"
    ),
    drop = TRUE
  ) +
  facet_grid(
    informative ~ event,
    labeller = label_both,
    drop = FALSE
  ) +
  labs(
    x = "time (ms)",
    y = "amplitude (\u00b5V)",
    color = "Mention_Modification"
  ) +
  theme_house()

ggsave("Grand_Average/grand_average_P600.png", p_grand_P600, width = 12, height = 6, dpi = 300)

grand_roi_ELAN <- all_ELAN %>%
  group_by(event, informative, mention, modified, subject, time) %>%
  summarise(
    subj_mean = mean(amplitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(event, informative, mention, modified, time) %>%
  summarise(
    mean_roi = mean(subj_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    event = factor(
      event,
      levels = c("DET1_crit", "ADJ1_crit", "N1_crit"),
      labels = c("DET", "ADJ", "N")
    ),
    cond = interaction(
      factor(mention, levels = c(1, 2)),
      factor(modified, levels = c(0, 1)),
      sep = "_mod"
    )
  )

p_grand_ELAN <- ggplot(grand_roi_ELAN, aes(x = time, y = mean_roi, color = cond)) +
  annotate(
    "rect",
    xmin = 100, xmax = 300,
    ymin = -Inf, ymax = Inf,
    fill = "lightgrey", alpha = 0.25
  ) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  scale_y_reverse() +
  scale_color_manual(
    name = "Mention x Modification",
    values = cond_colors,
    breaks = c("1_mod0", "2_mod0", "1_mod1", "2_mod1"),
    labels = c(
      "1_mod0" = "1mention_unmodified",
      "2_mod0" = "2mention_unmodified",
      "1_mod1" = "1mention_modified",
      "2_mod1" = "2mention_modified"
    ),
    drop = TRUE
  ) +
  facet_grid(
    informative ~ event,
    labeller = label_both,
    drop = FALSE
  ) +
  labs(
    x = "time (ms)",
    y = "amplitude (\u00b5V)",
    color = "Mention_Modification"
  ) +
  theme_house()

ggsave("Grand_Average/grand_average_ELAN.png", p_grand_ELAN, width = 12, height = 6, dpi = 300)

## ============================ LINEAR MIXED MODEL =============================
##### === Data sets for each ERP === ######
# windows.N400 = [100 500];
# windows.Nref = [300 800];
# windows.P600 = [500 800];
# windows.ELAN = [100 300];

det_N400_means <- det_N400 %>%
  filter(time >= 100, time <=500) %>%
  group_by(subject, mention, modified, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

det_Nref_means <- det_Nref %>%
  filter(time >= 300, time <=800) %>%
  group_by(subject, mention, modified, informative, ) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

det_P600_means <- det_P600 %>%
  filter(time >= 500, time <=800) %>%
  group_by(subject, mention, modified, informative, ) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

det_ELAN_means <- det_ELAN %>%
  filter(time >= 100, time <=300) %>%
  group_by(subject, mention, modified, informative, ) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

adj_N400_means <- adj_N400 %>%
  filter(time >= 100, time <= 500) %>%       
  group_by(subject, mention, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

adj_Nref_means <- adj_Nref %>%
  filter(time >= 300, time <=800) %>%
  group_by(subject, mention, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

adj_P600_means <- adj_P600 %>%
  filter(time >= 500, time <=800) %>%
  group_by(subject, mention, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

adj_ELAN_means <- adj_ELAN %>%
  filter(time >= 100, time <=300) %>%
  group_by(subject, mention, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

n_N400_means <- n_N400 %>%
  filter(time >= 100, time <= 500) %>%       
  group_by(subject, mention, modified, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

n_Nref_means <- n_Nref %>%
  filter(time >= 300, time <=800) %>%
  group_by(subject, mention, modified, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

n_P600_means <- n_P600 %>%
  filter(time >= 500, time <=800) %>%
  group_by(subject, mention, modified, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")

n_ELAN_means <- n_ELAN %>%
  filter(time >= 100, time <=300) %>%
  group_by(subject, mention, modified, informative) %>%
  summarise(amplitude = mean(amplitude), .groups = "drop")


### ==== DETERMINER ==== ###
##### === Sum coding === #####
det_N400_means$informative <- factor(
  det_N400_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

det_N400_means$mention <- factor(
  det_N400_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

det_N400_means$prev_mod <- factor(
  det_N400_means$prev_mod,
  levels = c(0, 1),
  labels = c("No_prev", "Prev")
)

det_N400_means$modified <- factor(
  det_N400_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

det_N400_means$informative <- relevel(det_N400_means$informative, ref = "Informative")
det_N400_means$mention <- relevel(det_N400_means$mention, ref = "First")
det_N400_means$prev_mod <- relevel(det_N400_means$prev_mod, ref = "Prev")
det_N400_means$modified <- relevel(det_N400_means$modified, ref = "Modified")

contrasts(det_N400_means$informative) <- contr.sum(2) / 2
contrasts(det_N400_means$mention)     <- contr.sum(2) / 2
contrasts(det_N400_means$prev_mod)     <- contr.sum(2) / 2
contrasts(det_N400_means$modified)    <- contr.sum(2) / 2

det_Nref_means$informative <- factor(
  det_Nref_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

det_Nref_means$mention <- factor(
  det_Nref_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

det_Nref_means$modified <- factor(
  det_Nref_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

det_Nref_means$informative <- relevel(det_Nref_means$informative, ref = "Informative")
det_Nref_means$mention <- relevel(det_Nref_means$mention, ref = "First")
det_Nref_means$modified <- relevel(det_Nref_means$modified, ref = "Modified")

contrasts(det_Nref_means$informative) <- contr.sum(2) / 2
contrasts(det_Nref_means$mention)     <- contr.sum(2) / 2
contrasts(det_Nref_means$modified)    <- contr.sum(2) / 2

det_P600_means$informative <- factor(
  det_P600_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

det_P600_means$mention <- factor(
  det_P600_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

det_P600_means$modified <- factor(
  det_P600_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

det_P600_means$informative <- relevel(det_P600_means$informative, ref = "Informative")
det_P600_means$mention <- relevel(det_P600_means$mention, ref = "First")
det_P600_means$modified <- relevel(det_P600_means$modified, ref = "Modified")

contrasts(det_P600_means$informative) <- contr.sum(2) / 2
contrasts(det_P600_means$mention)     <- contr.sum(2) / 2
contrasts(det_P600_means$modified)    <- contr.sum(2) / 2

det_ELAN_means$informative <- factor(
  det_ELAN_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

det_ELAN_means$mention <- factor(
  det_ELAN_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

det_ELAN_means$modified <- factor(
  det_ELAN_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

det_ELAN_means$informative <- relevel(det_ELAN_means$informative, ref = "Informative")
det_ELAN_means$mention <- relevel(det_ELAN_means$mention, ref = "First")
det_ELAN_means$modified <- relevel(det_ELAN_means$modified, ref = "Modified")

contrasts(det_ELAN_means$informative) <- contr.sum(2) / 2
contrasts(det_ELAN_means$mention)     <- contr.sum(2) / 2
contrasts(det_ELAN_means$modified)    <- contr.sum(2) / 2

##### === run linear model === #####
## === N400 ===
## setting up models 
m1_detN400 <- lmer(amplitude ~ modified * informative * mention + (1 + informative | subject),
                          data = det_N400_means)
summary(m1_detN400) 
isSingular(m1_detN400)

# plot + model check 
theme_set(theme_classic(base_size = 6))
diagnostic_plots <- check_model(m1_detN400)
diagnostic_plots <- plot(check_model(m1_detN400, panel = FALSE))

diagnostic_plots <- plot(check_model(m1_detN400, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))
panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_detN400.png", combined,
                width = 12, height = 10, dpi = 300)

combined <- wrap_plots(diagnostic_plots[c(1, 2, 3, 5)], ncol = 2)
ggplot2::ggsave("diagnostics_detN400.png", combined, width = 12, height = 10)

aggregate(amplitude ~ mention, data = det_N400_means, FUN = mean)
aggregate(amplitude ~ modified, data = det_N400_means, FUN = mean)

coefs <- as.data.frame(summary(m1_detN400)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error`
)

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on N400 amplitude.",
    label     = "tab:n400") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 4.74, informativity slope SD = 1.62, correlation = .63; residual SD = 2.50.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("det_N400.tex")


m2_detN400 <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                   data = det_N400_means)
summary(m2_detN400) 

m3_detN400 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * mention | subject),
                   data = det_N400_means)
summary(m3_detN400) 

m4_detN400 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified | subject),
                   data = det_N400_means)
summary(m4_detN400) 

m5_detN400 <- lmer(amplitude ~ modified * informative * mention + (1 + mention * mention | subject),
                   data = det_N400_means)
summary(m5_detN400) 

## model comparison
result1_N400_det <- compare_performance(m5_detN400, m4_detN400, m3_detN400)
print_md(result1_N400_det)
plot(result1_N400_det)

result_test1_N400_det <- test_performance(m5_detN400, m4_detN400, m3_detN400)
print_md(result_test1_N400_det)

result2_N400_det <- compare_performance(m1_detN400, m2_detN400, m3_detN400)
print_md(result2_N400_det)
plot(result2_N400_det)

result_test2_N400_det <- test_performance(m2_detN400, m1_detN400, m3_detN400)
print_md(result_test2_N400_det)

# go with model 1 

# == contrast testing ==
emmeans(m1_detN400, pairwise ~ informative, lmer.df = "satterthwaite") # for contrast testing 
emmeans(m1_detN400,  pairwise ~  modified, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ mention, lmer.df = "satterthwaite")
emmeans(m1_detN400, pairwise ~ informative |mention, lmer.df = "satterthwaite") 
emmeans(m1_detN400,  pairwise ~  modified|mention, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ mention|informative, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ mention|modified, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ informative|modified, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ modified|informative, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ informative|modified*mention, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ modified|informative*mention, lmer.df = "satterthwaite")
emmeans(m1_detN400,  pairwise ~ mention|modified*informative, lmer.df = "satterthwaite")

## === Nref ===
m1_detNref <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                   data = det_Nref_means)
summary(m1_detNref) 
isSingular(m1_detNref)

# plot + model check 
diagnostic_plots <- plot(check_model(m1_detNref, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_detNref.png", combined,
                width = 12, height = 10, dpi = 300)


b <- summary(m1_detNref)$coefficients["modified1:mention1", ]
b["Estimate"] + c(-1, 1) * qt(.975, b["df"]) * b["Std. Error"]

coefs <- as.data.frame(summary(m1_detNref)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error`
)

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on Nref amplitude at the determiner",
    label     = "tab_nref") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 3.08; residual SD = 3.19.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("det_Nref.tex")

# other models 
m2_detNref <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                   data = det_Nref_means)
summary(m2_detNref) 

m3_detNref <- lmer(amplitude ~ modified * informative * mention + (1 + informative * mention | subject),
                   data = det_Nref_means)
summary(m3_detNref) 

m4_detNref <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified | subject),
                   data = det_Nref_means)
summary(m4_detNref) 

m5_detNref <- lmer(amplitude ~ modified * informative * mention + (1 + modified * mention | subject),
                   data = det_Nref_means)
summary(m5_detNref) 

## model comparison
result1_Nref_det <- compare_performance(m5_detNref, m4_detNref, m3_detNref)
print_md(result1_Nref_det)
plot(result1_Nref_det)

result_test1_Nref_det <- test_performance(m5_detNref, m4_detNref, m3_detNref)
print_md(result_test1_Nref_det)

result2_Nref_det <- compare_performance(m2_detNref, m1_detNref, m4_detNref)
print_md(result2_Nref_det)
plot(result2_Nref_det)

result_test2_Nref_det <- test_performance(m2_detNref, m1_detNref, m4_detNref)
print_md(result_test2_Nref_det)

# go with model 1 (simple model)

# == contrast testing ==
emmeans(m1_detNref, pairwise ~ informative) # for contrast testing 
emmeans(m1_detNref,  pairwise ~  modified)
emmeans(m1_detNref,  pairwise ~ mention)
emmeans(m1_detNref, pairwise ~ informative |mention) 
emmeans(m1_detNref,  pairwise ~  modified|mention,  lmer.df = "satterthwaite")
emmeans(m1_detNref,  pairwise ~ mention|informative)
emmeans(m1_detNref,  pairwise ~ mention|modified)
emmeans(m1_detNref,  pairwise ~ informative|modified)
emmeans(m1_detNref,  pairwise ~ modified|informative)

## === P600 ===
m1_detP600 <- lmer(amplitude ~ modified * informative * mention + (1 + informative | subject),
                   data = det_P600_means)
summary(m1_detP600) 

m2_detP600 <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                   data = det_P600_means)
summary(m2_detP600) 

# plot + model check 
diagnostic_plots <- plot(check_model(m2_detP600, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_detP600.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m2_detP600)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
     lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
     upper = coefs$Estimate + tcrit * coefs$`Std. Error` )

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on P600 amplitude at the determiner",
    label     = "tab_P600") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 5.25; residual SD = 5.70.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("det_P600.tex")

# other models 
m3_detP600 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * mention | subject),
                   data = det_P600_means)
summary(m3_detP600) 

m4_detP600 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified  | subject),
                   data = det_P600_means)
summary(m4_detP600) 

m5_detP600 <- lmer(amplitude ~ modified * informative * mention + (1 + mention * modified| subject),
                   data = det_P600_means)
summary(m5_detP600) 

## model comparison
result1_P600_det <- compare_performance(m5_detP600, m4_detP600, m3_detP600)
print_md(result1_P600_det)
plot(result1_P600_det)

result_test1_P600_det <- test_performance(m5_detP600, m4_detP600, m3_detP600)
print_md(result_test1_P600_det)

result2_P600_det <- compare_performance(m2_detNref, m1_detNref, m3_detNref, m4_detNref)
print_md(result2_P600_det)
plot(result2_P600_det)

result_test2_P600_det <- test_performance(m2_detNref, m1_detNref, m4_detNref, m3_detNref)
print_md(result_test2_P600_det)

# go with model 1 is not singular 

# == contrast testing ==
emmeans(m3_detP600, pairwise ~ informative) # for contrast testing 
emmeans(m3_detP600,  pairwise ~  modified)
emmeans(m3_detP600,  pairwise ~ mention)
emmeans(m3_detP600, pairwise ~ informative |mention) 
emmeans(m3_detP600,  pairwise ~  modified|mention)
emmeans(m3_detP600,  pairwise ~ mention|informative)
emmeans(m3_detP600,  pairwise ~ mention|modified)
emmeans(m3_detP600,  pairwise ~ informative|modified)
emmeans(m3_detP600,  pairwise ~ modified|informative)

## === ELAN ===
m1_detELAN <- lmer(amplitude ~ modified * informative * mention + (1 + informative | subject),
                   data = det_ELAN_means)
summary(m1_detELAN) 
aggregate(amplitude ~ mention, data = det_ELAN_means, FUN = mean)

# plot + model check 
diagnostic_plots <- plot(check_model(m1_detELAN, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_detELAN.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m1_detELAN)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error` )

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on ELAN amplitude at the determiner",
    label     = "tab_ELAN") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 5.25, informativity slope SD = 0.52, correlation = .4; residual SD = 2.53",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("det_ELAN.tex")

# other models 
m2_detELAN <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                   data = det_ELAN_means)
summary(m2_detELAN) 

m3_detELAN <- lmer(amplitude ~ modified * informative * mention + (1 + informative * mention | subject),
                   data = det_ELAN_means)
summary(m3_detELAN) 
isSingular(m3_detELAN)

m4_detELAN <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified  | subject),
                   data = det_ELAN_means)
summary(m4_detELAN) 

m5_detELAN <- lmer(amplitude ~ modified * informative * mention + (1 + mention * modified| subject),
                   data = det_ELAN_means)
summary(m5_detELAN) 

## model comparison
result1_ELAN_det <- compare_performance(m5_detELAN, m4_detELAN, m3_detELAN)
print_md(result1_ELAN_det)
plot(result1_ELAN_det)

result_test1_ELAN_det <- test_performance(m5_detELAN, m4_detELAN, m3_detELAN)
print_md(result_test1_ELAN_det)

result2_ELAN_det <- compare_performance(m2_detELAN, m1_detELAN, m3_detELAN, m4_detELAN)
print_md(result2_ELAN_det)
plot(result2_ELAN_det)

result_test2_ELAN_det <- test_performance(m2_detELAN, m1_detELAN, m3_detELAN, m4_detELAN)
print_md(result_test2_ELAN_det )

# go with model 1

# == contrast testing ==
emmeans(m4_detELAN, pairwise ~ informative) # for contrast testing 
emmeans(m4_detELAN,  pairwise ~  modified)
emmeans(m4_detELAN,  pairwise ~ mention)
emmeans(m4_detELAN, pairwise ~ informative |mention) 
emmeans(m4_detELAN,  pairwise ~  modified|mention)
emmeans(m4_detELAN,  pairwise ~ mention|informative)
emmeans(m4_detELAN,  pairwise ~ mention|modified)
emmeans(m4_detELAN,  pairwise ~ informative|modified)
emmeans(m4_detELAN,  pairwise ~ modified|informative)

### ==== ADJECTIVE ==== ###
##### === Sum coding === #####
adj_N400_means$informative <- factor(
  adj_N400_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

adj_N400_means$mention <- factor(
  adj_N400_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

adj_N400_means$informative <- relevel(adj_N400_means$informative, ref = "Informative")
adj_N400_means$mention <- relevel(adj_N400_means$mention, ref = "First")

contrasts(adj_N400_means$informative) <- contr.sum(2) / 2
contrasts(adj_N400_means$mention)     <- contr.sum(2) / 2

adj_Nref_means$informative <- factor(
  adj_Nref_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

adj_Nref_means$mention <- factor(
  adj_Nref_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

adj_Nref_means$informative <- relevel(adj_Nref_means$informative, ref = "Informative")
adj_Nref_means$mention <- relevel(adj_Nref_means$mention, ref = "First")

contrasts(adj_Nref_means$informative) <- contr.sum(2) / 2
contrasts(adj_Nref_means$mention)     <- contr.sum(2) / 2

adj_P600_means$informative <- factor(
  adj_P600_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

adj_P600_means$mention <- factor(
  adj_P600_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

adj_P600_means$informative <- relevel(adj_P600_means$informative, ref = "Informative")
adj_P600_means$mention <- relevel(adj_P600_means$mention, ref = "First")

contrasts(adj_P600_means$informative) <- contr.sum(2) / 2
contrasts(adj_P600_means$mention)     <- contr.sum(2) / 2

adj_ELAN_means$informative <- factor(
  adj_ELAN_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

adj_ELAN_means$mention <- factor(
  adj_ELAN_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

adj_ELAN_means$informative <- relevel(adj_ELAN_means$informative, ref = "Informative")
adj_ELAN_means$mention <- relevel(adj_ELAN_means$mention, ref = "First")

contrasts(adj_ELAN_means$informative) <- contr.sum(2) / 2
contrasts(adj_ELAN_means$mention)     <- contr.sum(2) / 2

##### === run linear model === #####
## === N400 ===
m1_adjN400 <- lmer(amplitude ~ informative * mention + (1 + informative | subject),
                   data = adj_N400_means)
summary(m1_adjN400) 

m2_adjN400 <- lmer(amplitude ~ informative * mention + (1 | subject),
                   data = adj_N400_means)
summary(m2_adjN400) 
aggregate(amplitude ~ mention, data = adj_N400_means, FUN = mean)
aggregate(amplitude ~ modified, data = det_N400_means, FUN = mean)

# plot + model check 
diagnostic_plots <- plot(check_model(m2_adjN400, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)

combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_adjN400.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m2_adjN400)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error`
)

tab <- data.frame(
  Predictor = c("Intercept", "Informative", "Mention",
                "Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on N400 amplitude at the adjective",
    label     = "tab_n400_adj") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 5.41; residual SD = 4.18.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("adj_N400.tex")


model_lm_tidy_adj_N400 <- linear_reg() %>%
  set_engine("lm") %>%
  fit(amplitude ~ informative * mention,
      data = adj_N400_means)
check_model(model_lm_tidy_adj_N400)

# other model 
m3_adjN400 <- lmer(amplitude ~ informative * mention + (1 + mention | subject),
                   data = adj_N400_means)
summary(m3_adjN400) 


## model comparison
result1_N400_adj <- compare_performance(m1_adjN400, m2_adjN400, m3_adjN400)
print_md(result1_N400_adj)
plot(result1_N400_adj)

result_test1_N400_adj <- test_performance(m1_adjN400, m2_adjN400, m3_adjN400)
print_md(result_test1_N400_adj)

# go with model 2 

# == contrast testing ==
emmeans(m2_adjN400, pairwise ~ informative) # for contrast testing 
emmeans(m2_adjN400,  pairwise ~ mention)
emmeans(m2_adjN400, pairwise ~ informative|mention) 
emmeans(m2_adjN400,  pairwise ~ mention|informative)

## === Nref ===
m1_adjNref <- lmer(amplitude ~ informative * mention + (1 + informative | subject),
                   data = adj_Nref_means)
summary(m1_adjNref) 

m2_adjNref <- lmer(amplitude ~ informative * mention + (1  | subject),
                   data = adj_Nref_means)
summary(m2_adjNref) 

m3_adjNref <- lmer(amplitude ~ informative * mention + (1 + mention | subject),
                   data = adj_Nref_means)
summary(m3_adjNref) 

aggregate(amplitude ~ mention, data = m2_adjNref, FUN = mean)

# plot + model check 
diagnostic_plots <- plot(check_model(m3_adjNref, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)

combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_adjNref.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m2_adjNref)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error`
)

tab <- data.frame(
  Predictor = c("Intercept", "Informative", "Mention",
                "Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on Nref amplitude at the adjective",
    label     = "tab_nrerf_adj") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 3.25; residual SD = 2.23.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("adj_Nref.tex")

## model comparison
result1_Nref_adj <- compare_performance(m1_adjNref, m2_adjNref, m3_adjNref)
print_md(result1_Nref_adj)
plot(result1_Nref_adj)

result_test1_Nref_adj <- test_performance(m1_adjNref, m2_adjNref, m3_adjNref)
print_md(result_test1_Nref_adj)

# go with model 2

# == contrast testing ==
emmeans(m3_adjNref, pairwise ~ informative) # for contrast testing 
emmeans(m3_adjNref,  pairwise ~ mention)
emmeans(m3_adjNref, pairwise ~ informative|mention)
emmeans(m3_adjNref,  pairwise ~ mention|informative)

## === P600 ===
m1_adjP600 <- lmer(amplitude ~ informative * mention + (1 + informative | subject),
                   data = adj_P600_means)
summary(m1_adjP600) 

m2_adjP600 <- lmer(amplitude ~ informative * mention + (1 | subject),
                   data = adj_P600_means)
summary(m2_adjP600) 

# plot + model check 
diagnostic_plots <- plot(check_model(m2_adjP600, panel = FALSE))

is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_adjP600.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m2_adjP600)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error`
)

tab <- data.frame(
  Predictor = c("Intercept", "Informative", "Mention",
                "Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on P600 amplitude at the adjective",
    label     = "tab_nrerf_adj") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 4.95; residual SD = 2.73.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("adj_P600.tex")

m3_adjP600 <- lmer(amplitude ~ informative * mention + (1 + mention | subject),
                   data = adj_P600_means)
summary(m3_adjP600) 


## model comparison
result1_P600_adj <- compare_performance(m1_adjP600, m2_adjP600, m3_adjP600)
print_md(result1_P600_adj)
plot(result1_P600_adj)

result_test1_P600_adj <- test_performance(m1_adjP600, m2_adjP600, m3_adjP600)
print_md(result_test1_P600_adj)

# go with model 2

# == contrast testing ==
emmeans(m3_adjP600, pairwise ~ informative) # for contrast testing 
emmeans(m3_adjP600,  pairwise ~ mention)
emmeans(m3_adjP600, pairwise ~ informative|mention) 
emmeans(m3_adjP600,  pairwise ~ mention|informative)

## === ELAN ===
m1_adjELAN <- lmer(amplitude ~ informative * mention + (1 + informative | subject),
                   data = adj_ELAN_means)
summary(m1_adjELAN) 

m2_adjELAN <- lmer(amplitude ~ informative * mention + (1 | subject),
                   data = adj_ELAN_means)
summary(m2_adjELAN) 

m3_adjELAN <- lmer(amplitude ~ informative * mention + (1 + mention | subject),
                   data = adj_ELAN_means)
summary(m3_adjELAN) 
aggregate(amplitude ~ informative, data = adj_ELAN_means, FUN = mean)

# plot + model check 
diagnostic_plots <- plot(check_model(m3_adjELAN, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_adjELAN.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m3_adjELAN)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error`
)

tab <- data.frame(
  Predictor = c("Intercept", "Informative", "Mention",
                "Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on ELAN amplitude at the adjective",
    label     = "tab_ELAN_adj") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 4.95, mention intercept SD = 5.15, correlation = 0.30; residual SD = 2.73.",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("adj_ELAN.tex")


## model comparison
result1_ELAN_adj <- compare_performance(m1_adjELAN, m2_adjELAN, m3_adjELAN)
print_md(result1_ELAN_adj)
plot(result1_ELAN_adj)

result_test1_ELAN_adj <- test_performance(m1_adjELAN, m2_adjELAN, m3_adjELAN)
print_md(result_test1_ELAN_adj)

# go with model 3

# == contrast testing ==
emmeans(m3_adjELAN, pairwise ~ informative) # for contrast testing 
emmeans(m3_adjELAN,  pairwise ~ mention)
emmeans(m3_adjELAN, pairwise ~ informative|mention) 
emmeans(m3_adjELAN,  pairwise ~ mention|informative)

### ==== NOUN ==== ###
##### === Sum coding === #####
n_N400_means$informative <- factor(
  n_N400_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

n_N400_means$mention <- factor(
  n_N400_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

n_N400_means$modified <- factor(
  n_N400_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

n_N400_means$informative <- relevel(n_N400_means$informative, ref = "Informative")
n_N400_means$mention <- relevel(n_N400_means$mention, ref = "First")
n_N400_means$modified <- relevel(n_N400_means$modified, ref = "Modified")

contrasts(n_N400_means$informative) <- contr.sum(2) / 2
contrasts(n_N400_means$mention)     <- contr.sum(2) / 2
contrasts(n_N400_means$modified)    <- contr.sum(2) / 2

n_Nref_means$informative <- factor(
  n_Nref_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

n_Nref_means$mention <- factor(
  n_Nref_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

n_Nref_means$modified <- factor(
  n_Nref_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

n_Nref_means$informative <- relevel(n_Nref_means$informative, ref = "Informative")
n_Nref_means$mention <- relevel(n_Nref_means$mention, ref = "First")
n_Nref_means$modified <- relevel(n_Nref_means$modified, ref = "Modified")

contrasts(n_Nref_means$informative) <- contr.sum(2) / 2
contrasts(n_Nref_means$mention)     <- contr.sum(2) / 2
contrasts(n_Nref_means$modified)    <- contr.sum(2) / 2

n_P600_means$informative <- factor(
  n_P600_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

n_P600_means$mention <- factor(
  n_P600_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

n_P600_means$modified <- factor(
  n_P600_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

n_P600_means$informative <- relevel(n_P600_means$informative, ref = "Informative")
n_P600_means$mention <- relevel(n_P600_means$mention, ref = "First")
n_P600_means$modified <- relevel(n_P600_means$modified, ref = "Modified")

contrasts(n_P600_means$informative) <- contr.sum(2) / 2
contrasts(n_P600_means$mention)     <- contr.sum(2) / 2
contrasts(n_P600_means$modified)    <- contr.sum(2) / 2

n_ELAN_means$informative <- factor(
  n_ELAN_means$informative,
  levels = c(0, 1),
  labels = c("Uninformative", "Informative")
)

n_ELAN_means$mention <- factor(
  n_ELAN_means$mention,
  levels = c(1, 2),
  labels = c("First", "Second")
)

n_ELAN_means$modified <- factor(
  n_ELAN_means$modified,
  levels = c(0, 1),
  labels = c("Unmodified", "Modified")
)

n_ELAN_means$informative <- relevel(n_ELAN_means$informative, ref = "Informative")
n_ELAN_means$mention <- relevel(n_ELAN_means$mention, ref = "First")
n_ELAN_means$modified <- relevel(n_ELAN_means$modified, ref = "Modified")

contrasts(n_ELAN_means$informative) <- contr.sum(2) / 2
contrasts(n_ELAN_means$mention)     <- contr.sum(2) / 2
contrasts(n_ELAN_means$modified)    <- contr.sum(2) / 2

##### === run linear model === #####
## === N400 ===
m1_nN400 <- lmer(amplitude ~ modified * informative * mention + (1 + informative | subject),
                   data = n_N400_means)
summary(m1_nN400) 
isSingular(m1_nN400)
aggregate(amplitude ~ mention, data = det_N400_means, FUN = mean)

m2_nN400 <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                   data = n_N400_means)
summary(m2_nN400) 

m3_nN400 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * mention | subject),
                   data = n_N400_means)
summary(m3_nN400) 

m4_nN400 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified | subject),
                   data = n_N400_means)
summary(m4_nN400) 

m5_nN400 <- lmer(amplitude ~ modified * informative * mention + (1 + mention * mention | subject),
                   data = n_N400_means)
summary(m5_nN400) 

# plot + model check 
diagnostic_plots <- plot(check_model(m1_nN400, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_nN400.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m1_nN400)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error` )

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on N400 amplitude at the noun",
    label     = "tab_N400_n") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 4.67, informativity slope SD = 1.06, correlation = -0.09; residual SD = 3.68",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("n_N400.tex")


## model comparison
result1_N400_n <- compare_performance(m5_nN400, m4_nN400, m3_nN400)
print_md(result1_N400_n)
plot(result1_N400_n)

result_test1_N400_n <- test_performance(m5_nN400, m4_nN400, m3_nN400)
print_md(result_test1_N400_n)

result2_N400_n <- compare_performance(m1_nN400, m2_nN400, m3_nN400)
print_md(result2_N400_n)
plot(result2_N400_n)

result_test2_N400_n <- test_performance(m1_nN400, m2_nN400, m3_nN400)
print_md(result_test2_N400_n)

# go with model 1 

# == contrast testing ==
emmeans(m3_nN400, pairwise ~ informative) # for contrast testing 
emmeans(m3_nN400,  pairwise ~  modified)
emmeans(m3_nN400,  pairwise ~ mention)
emmeans(m3_nN400, pairwise ~ informative |mention) 
emmeans(m3_nN400,  pairwise ~  modified|mention)
emmeans(m3_nN400,  pairwise ~ mention|informative)
emmeans(m3_nN400,  pairwise ~ mention|modified)
emmeans(m3_nN400,  pairwise ~ informative|modified)
emmeans(m3_nN400,  pairwise ~ modified|informative)

## === Nref ===
m1_nNref <- lmer(amplitude ~ modified * informative * mention + (1 + informative | subject),
                   data = n_Nref_means)
summary(m1_nNref) 

aggregate(amplitude ~ mention, data = det_Nref_means, FUN = mean)
# plot + model check 
diagnostic_plots <- plot(check_model(m1_nNref, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_nNref.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m1_nNref)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error` )

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on Nref amplitude at the noun",
    label     = "tab_Nref_n") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 3.22, informativity slope SD = 0.39, correlation = 1.0; residual SD = 2.39",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("n_Nref.tex")

# other models 
m2_nNref <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                 data = n_Nref_means)
summary(m2_nNref) 

m3_nNref <- lmer(amplitude ~ modified * informative * mention + (1 + informative * mention | subject),
                 data = n_Nref_means)
summary(m3_nNref) 

m4_nNref <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified | subject),
                 data = n_Nref_means)
summary(m4_nNref) 

m5_nNref <- lmer(amplitude ~ modified * informative * mention + (1 + mention * mention | subject),
                 data = n_Nref_means)
summary(m5_nNref) 

## model comparison
result1_Nref_n <- compare_performance(m5_nNref, m4_nNref, m3_nNref)
print_md(result1_Nref_n)
plot(result1_Nref_n)

result_test1_Nref_n <- test_performance(m5_nNref, m4_nNref, m3_nNref)
print_md(result_test1_Nref_n)

result2_Nref_n <- compare_performance(m1_nNref, m2_nNref, m3_nNref)
print_md(result2_Nref_n)
plot(result2_Nref_n)

result_test2_Nref_n <- test_performance(m1_nNref, m2_nNref, m3_nNref)
print_md(result_test2_Nref_n)

# go with model 1

# == contrast testing ==
emmeans(m3_nNref, pairwise ~ informative) # for contrast testing 
emmeans(m3_nNref,  pairwise ~  modified)
emmeans(m3_nNref,  pairwise ~ mention)
emmeans(m3_nNref, pairwise ~ informative |mention) 
emmeans(m3_nNref,  pairwise ~  modified|mention)
emmeans(m3_nNref,  pairwise ~ mention|informative)
emmeans(m3_nNref,  pairwise ~ mention|modified)
emmeans(m3_nNref,  pairwise ~ informative|modified)
emmeans(m3_nNref,  pairwise ~ modified|informative)

## === P600 ===
m1_nP600 <- lmer(amplitude ~ modified * informative * mention + (1 + informative| subject),
                   data = n_P600_means)
summary(m1_nP600) 

aggregate(amplitude ~ mention, data = det_P600_means, FUN = mean)

# plot + model check 
diagnostic_plots <- plot(check_model(m1_nP600, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_nP600.png", combined,
                width = 12, height = 10, dpi = 300)

coefs <- as.data.frame(summary(m1_nP600)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error` )

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on P600 amplitude at the noun",
    label     = "tab_P600_n") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 5.23, informativity slope SD = 2.07, correlation = 0.77; residual SD = 3.07",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("n_P600.tex")


# other models 
m2_nP600 <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                 data = n_P600_means)
summary(m2_nP600)

m3_nP600 <- lmer(amplitude ~ modified * informative * mention + (1 + informative*mention| subject),
                 data = n_P600_means)
summary(m3_nP600) 

m4_nP600 <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified | subject),
                 data = n_P600_means)
summary(m4_nP600) 

m5_nP600 <- lmer(amplitude ~ modified * informative * mention + (1 + mention * mention | subject),
                 data = n_P600_means)
summary(m5_nP600) 

## model comparison
result1_P600_n <- compare_performance(m5_nP600, m4_nP600, m3_nP600)
print_md(result1_P600_n)
plot(result1_P600_n)

result_test1_P600_n <- test_performance(m5_nP600, m4_nP600, m3_nP600)
print_md(result_test1_P600_n)

result2_P600_n <- compare_performance(m1_nP600, m2_nP600, m3_nP600)
print_md(result2_P600_n)
plot(result2_P600_n)

result_test2_P600_n <- test_performance(m1_nP600, m2_nP600, m3_nP600)
print_md(result_test2_P600_n)

# go with model 1

# == contrast testing ==
emmeans(m1_nP600, pairwise ~ informative) # for contrast testing 
emmeans(m1_nP600,  pairwise ~  modified)
emmeans(m1_nP600,  pairwise ~ mention)
emmeans(m1_nP600, pairwise ~ informative |mention) 
emmeans(m1_nP600,  pairwise ~  modified|mention)
emmeans(m1_nP600,  pairwise ~ mention|informative)
emmeans(m1_nP600,  pairwise ~ mention|modified)
emmeans(m1_nP600,  pairwise ~ informative|modified)
emmeans(m1_nP600,  pairwise ~ modified|informative)
emmeans(m1_nP600,  pairwise ~  mention|modified + informative)
emmeans(m1_nP600,  pairwise ~  modified|mention + informative)

## === ELAN ===
m1_nELAN <- lmer(amplitude ~ modified * informative * mention + (1 + informative | subject),
                   data = n_ELAN_means)
summary(m1_nELAN) 

m2_nELAN <- lmer(amplitude ~ modified * informative * mention + (1 | subject),
                 data = n_ELAN_means)
summary(m2_nELAN)

aggregate(amplitude ~ mention, data = det_ELAN_means, FUN = mean)

# plot + model check 
diagnostic_plots <- plot(check_model(m2_nELAN, panel = FALSE))
is_gg <- vapply(diagnostic_plots, ggplot2::is_ggplot, logical(1))
diag_gg <- diagnostic_plots[is_gg]
print(names(diag_gg))

panel_names <- c("PP_CHECK", "NCV", "HOMOGENEITY", "QQ")   # edit to match the printout
panel_names <- intersect(panel_names, names(diag_gg))
print(panel_names)
combined <- wrap_plots(diag_gg[panel_names], ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  ggplot2::theme(
    plot.title    = ggplot2::element_blank(),
    plot.subtitle = ggplot2::element_blank(),
    plot.tag      = ggplot2::element_text(size = 18, face = "bold"),
    axis.title    = ggplot2::element_text(size = 20),
    axis.text     = ggplot2::element_text(size = 20),
    legend.text   = ggplot2::element_text(size = 20)
  )

ggplot2::ggsave("diagnostics_nELAN.png", combined,
                width = 12, height = 10, dpi = 300)


coefs <- as.data.frame(summary(m2_nELAN)$coefficients)
tcrit <- qt(.975, coefs$df)
ci <- cbind(
  lower = coefs$Estimate - tcrit * coefs$`Std. Error`,
  upper = coefs$Estimate + tcrit * coefs$`Std. Error` )

tab <- data.frame(
  Predictor = c("Intercept", "Modified", "Informative", "Mention",
                "Modified $\\times$ Informative",
                "Modified $\\times$ Mention",
                "Informative $\\times$ Mention",
                "Modified $\\times$ Informative $\\times$ Mention"),
  b   = sprintf("%.2f", coefs$Estimate),
  SE  = sprintf("%.2f", coefs$`Std. Error`),
  CI  = sprintf("[%.2f, %.2f]", ci[, 1], ci[, 2]),
  df  = sprintf("%.1f", coefs$df),
  t   = sprintf("%.2f", coefs$`t value`),
  p   = ifelse(coefs$`Pr(>|t|)` < .001, "<.001",
               sub("^0", "", sprintf("%.3f", coefs$`Pr(>|t|)`)))
)

kbl(tab,
    format    = "latex",
    booktabs  = TRUE,
    escape    = FALSE,
    align     = c("l", rep("r", 6)),
    caption   = "Fixed effects for the linear mixed model on ELAN amplitude at the noun",
    label     = "tab_ELAN_n") |>
  kable_styling(latex_options = "hold_position") |>
  footnote(general = "Random effects: subject intercept SD = 3.54; residual SD = 4.73",
           general_title = "Note.",
           threeparttable = TRUE) |>
  save_kable("n_ELAN.tex")

# other models
m3_nELAN <- lmer(amplitude ~ modified * informative * mention + (1 + informative*mention| subject),
                 data = n_ELAN_means)
summary(m3_nELAN) 

m4_nELAN <- lmer(amplitude ~ modified * informative * mention + (1 + informative * modified | subject),
                 data = n_ELAN_means)
summary(m4_nELAN) 

m5_nELAN <- lmer(amplitude ~ modified * informative * mention + (1 + mention * mention | subject),
                 data = n_ELAN_means)
summary(m5_nELAN) 

## model comparison
result1_ELAN_n <- compare_performance(m5_nELAN, m4_nELAN, m3_nELAN)
print_md(result1_ELAN_n)
plot(result1_ELAN_n)

result_test1_ELAN_n <- test_performance(m5_nELAN, m4_nELAN, m3_nELAN)
print_md(result_test1_ELAN_n)

result2_ELAN_n <- compare_performance(m1_nELAN, m2_nELAN, m3_nELAN)
print_md(result2_ELAN_n)
plot(result2_ELAN_n)

result_test2_ELAN_n <- test_performance(m1_nELAN, m2_nELAN, m3_nELAN)
print_md(result_test2_ELAN_n)

# go with model 2

# == contrast testing ==
emmeans(m3_nELAN, pairwise ~ informative) # for contrast testing 
emmeans(m3_nELAN,  pairwise ~  modified)
emmeans(m3_nELAN,  pairwise ~ mention)
emmeans(m3_nELAN, pairwise ~ informative |mention) 
emmeans(m3_nELAN,  pairwise ~  modified|mention)
emmeans(m3_nELAN,  pairwise ~ mention|informative)
emmeans(m3_nELAN,  pairwise ~ mention|modified)
emmeans(m3_nELAN,  pairwise ~ informative|modified)
emmeans(m3_nELAN,  pairwise ~ modified|informative)
emmeans(m3_nELAN,  pairwise ~  mention|modified + informative)
emmeans(m3_nELAN,  pairwise ~  modified|mention + informative)

################################################################################
### END OF SCRIPT ###
################################################################################




