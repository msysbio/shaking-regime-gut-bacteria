# Supplementary kinetic figures (Fig S3-S5)

suppressMessages({library(tidyverse); library(patchwork)})
if (!exists("DATA_DIR")) DATA_DIR <- "."

PAL    <- c(Static = "#888888", Pulsed = "#7F77DD", Continuous = "#1D9E75")
NAVY   <- "#1E2761"
ACCENT <- "#B85042"
LEVELS <- c("Static", "Pulsed", "Continuous")

theme_pub <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#F4F6F8", colour = NA),
          strip.text = element_text(face = "bold", colour = NAVY),
          plot.title = element_text(face = "bold", colour = NAVY),
          plot.subtitle = element_text(colour = "#6B7280", size = base - 2),
          plot.caption = element_text(colour = "#6B7280", size = base - 3.5, hjust = 0),
          legend.position = "bottom")
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(DATA_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(DATA_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  cat("wrote", name, "\n")
}

NO_TEST <- paste0(
  "DESCRIPTIVE ONLY - NO STATISTICAL TEST IS REPORTED, AND NONE IS POSSIBLE.\n",
  "The Static plates were read at endpoint only and have no growth curve, so no kinetic metric has a paired same-run control.\n",
  "Static-vs-shaking cannot be computed at all here; Pulsed-vs-Continuous is a pure between-run comparison with no internal\n",
  "control (the two regimes never share a run). Any apparent difference below may be a run effect. See Fig. 3A.")

kin <- read_csv(file.path(DATA_DIR, "MASTER_kinetic_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
ts  <- read_csv(file.path(DATA_DIR, "MASTER_timeseries_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included") |>
  mutate(condition = factor(condition, levels = LEVELS))

base_shift <- ts |>
  group_by(species, condition, experiment, time_h) |>
  summarise(OD = mean(OD_corr), .groups = "drop") |>
  filter(time_h <= 1) |>
  group_by(species, condition, experiment) |>
  summarise(base = median(OD), .groups = "drop")

curves <- ts |>
  group_by(species, condition, experiment, time_h) |>
  summarise(OD = mean(OD_corr), OD_sd = sd(OD_corr), .groups = "drop") |>
  mutate(OD_sd = ifelse(is.na(OD_sd), 0, OD_sd)) |>
  left_join(base_shift, by = c("species", "condition", "experiment")) |>
  mutate(OD_show = OD - base)

fS3 <- ggplot(curves, aes(time_h, OD_show, colour = condition, fill = condition,
                          group = interaction(experiment, condition))) +
  geom_ribbon(data = ~dplyr::filter(.x, time_h >= 1),
              aes(ymin = OD_show - OD_sd, ymax = OD_show + OD_sd), colour = NA, alpha = 0.15) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  facet_wrap(~species, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = PAL[c("Pulsed", "Continuous")], drop = TRUE) +
  scale_fill_manual(values = PAL[c("Pulsed", "Continuous")], drop = TRUE) +
  labs(title = "Growth curves by species and shaking regime",
       subtitle = "One line per run (mean of wells); band = \u00b11 SD across wells; baseline-shifted to start near 0 | Static absent (endpoint-only)",
       x = "Time (h)", y = expression("OD"[600]*" (baseline-shifted)"), colour = NULL, fill = NULL,
       caption = NO_TEST) +
  theme_pub()
save_fig(fS3, "FigS3_growth_curves", 10.0, 6.2)

mu <- ts |>
  group_by(species, condition, experiment, well) |> arrange(time_h, .by_group = TRUE) |>
  summarise(mu_max = {
    od <- pmax(OD_corr, 1e-3); t <- time_h; w <- 5
    if (length(t) > w) max((log(od[(w+1):length(t)]) - log(od[1:(length(t)-w)])) /
                           (t[(w+1):length(t)] - t[1:(length(t)-w)]), na.rm = TRUE) else NA_real_
  }, .groups = "drop")

met <- kin |>
  select(species, condition, experiment, well, AUC_delta, max_OD, time_max_h) |>
  left_join(mu, by = c("species", "condition", "experiment", "well")) |>
  pivot_longer(c(AUC_delta, max_OD, mu_max, time_max_h),
               names_to = "metric", values_to = "value") |>
  mutate(condition = factor(condition, levels = LEVELS),
         metric = factor(metric,
           levels = c("AUC_delta", "max_OD", "mu_max", "time_max_h"),
           labels = c("AUC (baseline-corrected, OD*h)", "Maximum OD",
                      "Max. specific growth rate (1/h)", "Time to maximum OD (h)")))

run_means <- met |>
  group_by(species, condition, metric, experiment) |>
  summarise(m = mean(value, na.rm = TRUE), .groups = "drop")

fS4 <- ggplot(met, aes(species, value, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.65, colour = "grey40",
               position = position_dodge(0.75)) +
  geom_point(data = run_means, aes(y = m, group = condition), shape = 21, size = 1.9,
             fill = "white", colour = "grey20", stroke = 0.5,
             position = position_dodge(0.75)) +
  facet_wrap(~metric, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = PAL[c("Pulsed", "Continuous")], drop = TRUE) +
  labs(title = "Kinetic growth-curve metrics by species and shaking regime",
       subtitle = "Boxes = wells | white points = per-run means | Static absent: no growth curve was recorded",
       x = NULL, y = NULL, fill = NULL,
       caption = paste0(NO_TEST, "\n",
         "AUC is baseline-corrected (AUC_delta) because starting OD differs across runs. FD was run at two starting ODs (0.010, 0.027).")) +
  theme_pub()
save_fig(fS4, "FigS4_kinetic_metrics", 10.5, 7.0)

cov_t <- ts |>
  filter(time_h >= 4) |>
  group_by(species, condition, experiment, time_h) |>
  filter(n() >= 3, mean(OD_corr) > 0) |>
  summarise(cov = sd(OD_corr) / mean(OD_corr) * 100, .groups = "drop")

cov_run <- cov_t |>
  group_by(species, condition, experiment) |>
  summarise(CoV = mean(cov), .groups = "drop") |>
  mutate(flag = CoV > 50)

pA <- ggplot(cov_t, aes(time_h, cov, colour = condition,
                        group = interaction(experiment, condition))) +
  geom_hline(yintercept = 50, linetype = "dotted", colour = ACCENT, linewidth = 0.4) +
  geom_line(linewidth = 0.55, alpha = 0.85) +
  facet_wrap(~species, nrow = 2) +
  scale_y_log10(labels = scales::label_number(accuracy = 1)) +
  scale_colour_manual(values = PAL[c("Pulsed", "Continuous")], drop = TRUE) +
  labs(title = "A. Well-to-well CoV over time", x = "Time (h)", y = "CoV (%), log scale",
       colour = NULL,
       subtitle = "One line per run | t >= 4 h | dotted line = 50% CoV, for reference") +
  theme_pub() + theme(plot.title = element_text(size = 12))

cov_bar <- cov_run |> group_by(species, condition) |>
  summarise(mean_CoV = mean(CoV), .groups = "drop")

pB <- ggplot(cov_bar, aes(condition, mean_CoV, fill = condition)) +
  geom_col(alpha = 0.5, colour = "grey40", width = 0.62) +
  geom_point(data = cov_run, aes(y = CoV), shape = 21, size = 2.2,
             fill = "white", colour = "grey20", stroke = 0.6,
             position = position_jitter(width = 0.07, height = 0)) +
  geom_text(data = cov_run |> filter(flag),
            aes(y = CoV, label = sprintf("%.0f%% (SD > mean)", CoV)),
            vjust = -1.1, size = 2.5, colour = ACCENT, fontface = "bold") +
  facet_wrap(~species, scales = "free_y", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  scale_fill_manual(values = PAL[c("Pulsed", "Continuous")], drop = TRUE) +
  labs(title = "B. Mean CoV over the growth phase", x = NULL, y = "CoV (%)", fill = NULL,
       subtitle = "Bar = mean across runs | points = individual runs") +
  theme_pub() + theme(plot.title = element_text(size = 12), legend.position = "none")

fS5 <- (pA / pB) +
  plot_annotation(
    title = "Well-to-well variability during growth",
    caption = paste0(
      "Pulsed shaking gave the lower well-to-well CoV over the growth phase in 5 of 6 species (E. coli the exception).\n",
      "In M. gnavus, one continuous run exceeds 100% CoV (SD > mean) and inflates its continuous mean, so M. gnavus's ranking rests on that run.\n",
      "That plate should not be read as a measurement, and it alone drives MG's kinetic Continuous mean.\n",
      "Static is absent throughout: it has no growth curve. As above, no test on a kinetic metric is reported here."),
    theme = theme(plot.title = element_text(face = "bold", colour = NAVY, size = 14),
                  plot.caption = element_text(colour = "#6B7280", size = 7, hjust = 0)))
save_fig(fS5, "FigS5_kinetic_CoV", 10.5, 9.5)

cat("\nSupplementary kinetic figures written to", DATA_DIR, "\n")
cat("REMINDER: these are descriptive. Do not add a significance test to any of them.\n")
