# Main figures (Figs 3-5) and Fig S2

suppressMessages({library(tidyverse); library(ggrepel); library(patchwork)})
if (!exists("DATA_DIR")) DATA_DIR <- "."

PAL   <- c(Static = "#888888", Pulsed = "#7F77DD", Continuous = "#1D9E75")
NAVY  <- "#1E2761"
LEVELS <- c("Static", "Pulsed", "Continuous")

theme_pub <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "#F4F6F8", colour = NA),
          strip.text = element_text(face = "bold", colour = NAVY),
          plot.title = element_text(face = "bold", colour = NAVY),
          plot.subtitle = element_text(colour = "#6B7280", size = base - 1),
          plot.caption = NULL,
          legend.position = "bottom")
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(DATA_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(DATA_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  cat("wrote", name, ".png / .pdf\n")
}

kin <- read_csv(file.path(DATA_DIR, "MASTER_kinetic_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
sta <- read_csv(file.path(DATA_DIR, "MASTER_static_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
wells <- bind_rows(kin |> select(species, experiment, condition, endpoint),
                   sta |> select(species, experiment, condition, endpoint)) |>
  mutate(condition = factor(condition, levels = LEVELS))

st <- read_csv(file.path(DATA_DIR, "02_endpoint_stats.csv"), show_col_types=FALSE)
dl <- read_csv(file.path(DATA_DIR, "01_perrun_paired_deltas.csv"), show_col_types=FALSE)
cd <- read_csv(file.path(DATA_DIR, "02_endpoint_celldata.csv"), show_col_types=FALSE)

nlab <- cd |> select(species, condition, n_runs) |>
  pivot_wider(names_from=condition, values_from=n_runs) |>
  mutate(lab = sprintf("%s (%dP, %dC, %dS)", species, Pulsed, Continuous, Static)) |>
  select(species, lab)

est <- st |>
  mutate(regime = str_remove(contrast, " - Static"),
         regime = factor(regime, levels=c("Pulsed","Continuous"),
                         labels=c("Pulsed - Static", "Continuous - Static")),
         lo = estimate - 1.96*SE, hi = estimate + 1.96*SE,
         star = ifelse(is.na(p_bonf), "n/a", sig),
         txt  = ifelse(is.na(estimate), "not estimable",
                       sprintf("%+.3f %s", estimate, sig))) |>
  left_join(nlab, by="species") |>
  mutate(lab = factor(lab, levels=rev(sort(unique(lab)))))

pts <- dl |>
  mutate(regime = factor(condition, levels=c("Pulsed","Continuous"),
                         labels=c("Pulsed - Static", "Continuous - Static"))) |>
  left_join(nlab, by="species") |>
  mutate(lab = factor(lab, levels=levels(est$lab)))

# Background shading: left of zero = static yields more biomass (static grey);
# right of zero = the shaking regime yields more (that regime's colour).
L <- max(abs(c(est$lo, est$hi, pts$delta)), na.rm = TRUE) * 1.08
shade <- tibble(regime = factor(rep(levels(est$regime), each = 2), levels = levels(est$regime)),
                side   = rep(c("Static", "Shaken"), 2),
                xmin   = rep(c(-L, 0), 2), xmax = rep(c(0, L), 2)) |>
  mutate(fill_key = ifelse(side == "Static", "Static",
                           ifelse(regime == "Pulsed - Static", "Pulsed", "Continuous")),
         fill_key = factor(fill_key, levels = c("Static", "Pulsed", "Continuous")))

p <- ggplot(est, aes(estimate, lab, colour=regime)) +
  geom_rect(data = shade, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_key),
            alpha = 0.16) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50", linewidth=0.4) +
  geom_point(data=pts, aes(x=delta, y=lab), shape=21, size=1.9,
             fill="white", stroke=0.5, alpha=0.85,
             position=position_nudge(y=0.22)) +
  geom_errorbar(aes(xmin=lo, xmax=hi), width=0.18, orientation="y", linewidth=0.7, na.rm=TRUE) +
  geom_point(size=3.2, na.rm=TRUE) +
  geom_text(aes(x=Inf, label=txt), hjust=1.05, size=4, colour=NAVY, fontface="bold") +
  facet_wrap(~regime, nrow=1) +
  scale_colour_manual(values=unname(PAL[c("Pulsed","Continuous")]), guide="none") +
  scale_fill_manual(values = PAL, name = "Shaded side: higher endpoint biomass under",
                    guide = guide_legend(override.aes = list(alpha = 0.35))) +
  scale_x_continuous(breaks=c(-0.5, -0.25, 0, 0.25, 0.5), expand=expansion(mult=c(0, 0.42))) +
  labs(title="A",
       subtitle="Effect on endpoint OD (shaken plate vs its matched static plate)",
       x=expression("Difference in endpoint OD"[600]*" (shaken "-" static)"),
       y=NULL,
       caption=NULL) +
  theme_bw(13) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_line(linetype="dotted"),
        strip.background=element_rect(fill="#F4F6F8", colour=NA),
        strip.text=element_text(face="bold", colour=NAVY, size=13),
        plot.subtitle=element_text(colour="#374151", size=12.5),
        axis.text.y=element_text(size=11),
        axis.text.x=element_text(size=11),
        plot.title=element_text(face="bold", colour=NAVY, size=16),
        axis.title.x=element_text(size=13, colour=NAVY, lineheight=1.4),
        legend.position="bottom", legend.title=element_text(size=11.5, colour=NAVY),
        legend.text=element_text(size=11))
f2a_forest <- p

stats2 <- read_csv(file.path(DATA_DIR, "02_endpoint_stats.csv"), show_col_types = FALSE)

run_means <- wells |>
  group_by(species, condition, experiment) |>
  summarise(m = mean(endpoint), .groups = "drop")

ann <- stats2 |>
  mutate(condition = factor(str_remove(contrast, " - Static"), levels = LEVELS),
         label = ifelse(is.na(estimate), "n/a", sprintf("%s\n%+.3f", sig, estimate))) |>
  select(species, condition, label)
ytop <- wells |> group_by(species) |> summarise(y = max(endpoint) * 1.22, .groups = "drop")
ann  <- left_join(ann, ytop, by = "species")

f2a <- ggplot(wells, aes(condition, endpoint, fill = condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55, width = 0.62, colour = "grey35") +
  geom_point(data = run_means, aes(y = m), shape = 21, size = 2.4,
             fill = "white", colour = "grey20", stroke = 0.7,
             position = position_jitter(width = 0.09, height = 0, seed = 1)) +
  geom_text(data = ann, inherit.aes = FALSE,
            aes(x = condition, y = y, label = label),
            size = 3.6, colour = NAVY, lineheight = 0.9) +
  facet_wrap(~species, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = PAL) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  expand_limits(y = 0) +

  geom_text(data = tibble(species = "FD", condition = factor("Static", levels = LEVELS),
                          y = 0.12, lab = "FD: two starting ODs (0.010, 0.027)"),
            inherit.aes = FALSE, aes(x = 2, y = y, label = lab),
            size = 3.4, colour = "#B85042", fontface = "italic") +
  labs(title = "B",
       subtitle = NULL,
       x = NULL, y = expression("Endpoint OD"[600]), fill = NULL,
       caption = NULL) +
  theme_pub(13) + theme(plot.title = element_text(face = "bold", colour = NAVY, size = 16),
                        strip.text = element_text(face = "bold", colour = NAVY, size = 13),
                        axis.text  = element_text(size = 11),
                        axis.title.y = element_text(size = 13))

fig2a <- (f2a_forest / f2a) +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(
    title = "Endpoint OD by shaking regime",
    caption = NULL,
    theme = theme(plot.title = element_text(face = "bold", colour = NAVY, size = 16),
                  plot.caption = NULL))
save_fig(fig2a, "Fig4_endpoint_OD", 10.6, 11.2)

fcp <- read_csv(file.path(DATA_DIR, "03_FC_plate_foldchange.csv"), show_col_types = FALSE)
fct <- read_csv(file.path(DATA_DIR, "03_FC_tests.csv"), show_col_types = FALSE)

fct <- fct |>
  mutate(sig_CI = ifelse(CI_lo > 1 | CI_hi < 1, "*", "not significant"))
lab <- fct |>
  mutate(label = sprintf(
    "geometric mean %.2fx  %s\n95%% CI %.2f-%.2f\n%s\nt-test (log2 FC) p = %.3f",
    geom_mean_FC, sig_CI, CI_lo, CI_hi, concordant, ttest_log2_p))

f2b <- ggplot(fcp, aes(species, fold_change)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_crossbar(data = fct, inherit.aes = FALSE,
                aes(x = species, y = geom_mean_FC, ymin = CI_lo, ymax = CI_hi),
                width = 0.5, fill = "grey92", colour = NAVY, linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(fill = regime), shape = 21, size = 4, stroke = 0.7,
             colour = "grey20", position = position_jitter(width = 0.10, height = 0, seed = 1)) +
  geom_text(data = lab, inherit.aes = FALSE,
            aes(x = species, y = 1.90, label = label),
            size = 3.4, colour = NAVY, lineheight = 1.05) +
  scale_fill_manual(values = PAL[c("Pulsed", "Continuous")]) +
  scale_y_continuous(limits = c(0.75, 2.05), breaks = seq(0.8, 1.6, 0.2)) +
  labs(title = "Cell density (flow cytometry) vs paired static",
       subtitle = "Each point = one plate | bar = geometric mean +/- 95% bootstrap CI",
       x = NULL, y = "Fold change in cells/mL (shaking / static)\nabove 1.0 = shaking yields MORE cells", fill = NULL,
       caption = NULL) +
  theme_pub()
save_fig(f2b, "Fig5_flow_cytometry", 6.8, 6.8)

bat <- read_csv(file.path(DATA_DIR, "04_batch_effects.csv"), show_col_types = FALSE) |>
  mutate(condition = factor(condition, levels = LEVELS),
         species   = factor(species, levels = rev(c("BT", "EC", "MG", "RI", "SC"))),
         txt = ifelse(is.na(batch_pct), "n/a", sprintf("%d%%\n%s", batch_pct, sig)))

f3 <- ggplot(bat, aes(condition, species, fill = batch_pct)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = txt, colour = batch_pct > 55), size = 3.1, lineheight = 0.9,
            fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(low = "#F4F6F8", high = "#B85042", na.value = "grey93",
                      limits = c(0, 100), name = "% variance\nbetween plates") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = NAVY)) +
  labs(title = "A. Between-plate variance",
       subtitle = "Variance from run identity, per regime",
       x = NULL, y = NULL) +
  theme_pub() + theme(panel.grid = element_blank(), legend.position = "right",
                      plot.title = element_text(face = "bold", colour = NAVY, size = 12))

cov_s <- read_csv(file.path(DATA_DIR, "05_CoV_summary.csv"), show_col_types = FALSE) |>
  filter(condition %in% c("Pulsed", "Continuous")) |>
  mutate(condition = factor(condition, levels = c("Continuous", "Pulsed")))
cov_t <- read_csv(file.path(DATA_DIR, "05_CoV_tests.csv"), show_col_types = FALSE) |>
  filter(metric == "endpoint_CoV")

f4 <- ggplot(cov_s, aes(condition, endpoint_CoV, group = species)) +
  geom_line(colour = "grey55", linewidth = 0.6) +
  geom_point(aes(fill = condition), shape = 21, size = 4.2, stroke = 0.7, colour = "grey20") +
  geom_text_repel(data = cov_s |> filter(condition == "Pulsed"),
            aes(label = species), nudge_x = 0.18, direction = "y", hjust = 0,
            segment.colour = "grey75", segment.size = 0.3, min.segment.length = 0,
            size = 3.1, colour = NAVY, fontface = "bold", box.padding = 0.15) +
  scale_fill_manual(values = PAL[c("Continuous", "Pulsed")]) +
  scale_x_discrete(expand = expansion(add = c(0.35, 0.75))) +
  labs(title = "B. Within-plate variability",
       subtitle = sprintf("Lower within-plate CoV under pulsed (%s species)\nWilcoxon p=%.3f | t p=%.3f",
                          cov_t$favours_Pulsed, cov_t$wilcox_p, cov_t$ttest_p),
       x = NULL, y = "Within-plate CoV of endpoint OD (%)", fill = NULL) +
  theme_pub() + theme(legend.position = "none",
                      plot.title = element_text(face = "bold", colour = NAVY, size = 12))

fig3 <- (f3 | f4) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "Between-plate and within-plate variation in endpoint OD",
    caption = NULL,
    theme = theme(plot.title = element_text(face = "bold", colour = NAVY, size = 14),
                  plot.caption = NULL))
save_fig(fig3, "Fig3_variation", 13.4, 6.6)

cat("\nAll figures written to", DATA_DIR, "\n")

b_s  <- read_csv(file.path(DATA_DIR, "04_batch_effects.csv"), show_col_types = FALSE) |>
  filter(!is.na(ICC))
cd_s <- read_csv(file.path(DATA_DIR, "02_endpoint_celldata.csv"), show_col_types = FALSE)
ts_s <- read_csv(file.path(DATA_DIR, "MASTER_timeseries_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")

mu_s <- ts_s |>
  group_by(species, experiment, well) |> arrange(time_h, .by_group = TRUE) |>
  summarise(mu = {
    od <- pmax(OD_corr, 1e-3); t <- time_h; w <- 5
    if (length(t) > w) max((log(od[(w+1):length(t)]) - log(od[1:(length(t)-w)])) /
                           (t[(w+1):length(t)] - t[1:(length(t)-w)]), na.rm = TRUE) else NA
  }, .groups = "drop") |>
  group_by(species) |> summarise(mu_max = mean(mu, na.rm = TRUE), .groups = "drop")

lvl_s <- cd_s |> group_by(species) |>
  summarise(mean_OD = weighted.mean(mean, n_wells), .groups = "drop")

sc_dat <- b_s |> group_by(species) |>
  summarise(sd_between = mean(sqrt(var_plate)), .groups = "drop") |>
  left_join(lvl_s, by = "species") |> left_join(mu_s, by = "species") |>
  mutate(between_CoV = 100 * sd_between / mean_OD) |>
  pivot_longer(c(sd_between, between_CoV), names_to = "metric", values_to = "val") |>
  mutate(metric = factor(metric, levels = c("sd_between", "between_CoV"),
    labels = c("Absolute spread (OD units)",
               "Same spread, as % of that species' own mean OD")))

s2_panel <- function(m, ttl, ylab) {
  ggplot(filter(sc_dat, metric == m), aes(mu_max, val)) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey75",
                linewidth = 0.5, linetype = "dashed", formula = y ~ x) +
    geom_point(size = 3.4, colour = NAVY) +
    ggrepel::geom_text_repel(aes(label = species), size = 3.2, colour = NAVY,
                             fontface = "bold", box.padding = 0.4, seed = 1) +
    scale_y_continuous(expand = expansion(mult = 0.14)) +
    scale_x_continuous(expand = expansion(mult = 0.12)) +
    labs(title = ttl, x = "Mean maximum specific growth rate (1/h)", y = ylab) +
    theme_pub() + theme(plot.title = element_text(size = 11))
}
fS2 <- (s2_panel("Absolute spread (OD units)", "A. Absolute",
                 expression("Between-run SD of endpoint OD"[600])) |
        s2_panel("Same spread, as % of that species' own mean OD", "B. Relative to each species' mean",
                 expression("Between-run CoV of endpoint OD"[600]*" (%)"))) +
  plot_annotation(
    title = "Between-run variation in endpoint OD versus growth rate",
    subtitle = "One point per species | exploratory, n = 5",
    caption = NULL,
    theme = theme(plot.title = element_text(face = "bold", colour = NAVY, size = 13),
                  plot.subtitle = element_text(colour = "#6B7280", size = 10),
                  plot.caption = NULL))
save_fig(fS2, "FigS2_spread_vs_growthrate", 8.4, 5.4)

cat("\nAll figures written to", DATA_DIR, "\n")
