# What the design can and cannot estimate (regime vs batch)

suppressMessages({
  library(tidyverse); library(lme4); library(lmerTest); library(emmeans)
})

emm_options(lmerTest.limit = 100000, pbkrtest.limit = 100000)

if (!exists("DATA_DIR")) DATA_DIR <- "."
SPECIES <- c("BT", "EC", "FD", "MG", "RI", "SC")

kin <- read_csv(file.path(DATA_DIR, "MASTER_kinetic_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
sta <- read_csv(file.path(DATA_DIR, "MASTER_static_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")

dat <- bind_rows(
  kin |> select(experiment, species, condition, well, endpoint),
  sta |> select(experiment, species, condition, well, endpoint)
) |>
  mutate(condition  = factor(condition, levels = c("Static", "Pulsed", "Continuous")),
         experiment = factor(experiment))

pairs_in_run <- dat |>
  distinct(species, experiment, condition) |>
  group_by(species, experiment) |>
  summarise(conds = list(sort(as.character(condition))), .groups = "drop") |>
  mutate(
    has_S_P = map_lgl(conds, ~ all(c("Static", "Pulsed")     %in% .x)),
    has_S_C = map_lgl(conds, ~ all(c("Static", "Continuous") %in% .x)),
    has_P_C = map_lgl(conds, ~ all(c("Pulsed", "Continuous") %in% .x))
  )

cooc <- pairs_in_run |>
  group_by(species) |>
  summarise(
    `runs with Static+Pulsed`     = sum(has_S_P),
    `runs with Static+Continuous` = sum(has_S_C),
    `runs with Pulsed+Continuous` = sum(has_P_C),
    .groups = "drop"
  )

cat("\n=========================================================================\n")
cat("1. CO-OCCURRENCE WITHIN A RUN  (this is the design constraint)\n")
cat("=========================================================================\n")
print(as.data.frame(cooc), row.names = FALSE)

stopifnot(sum(cooc$`runs with Pulsed+Continuous`) == 0)
cat("\n>>> Pulsed and Continuous co-occur on ZERO runs, for every species.\n")
cat(">>> The P-vs-C contrast therefore carries NO within-run information.\n")
cat(">>> No random effect can de-confound it. Do not report it as a regime effect.\n")
write_csv(cooc, file.path(DATA_DIR, "01_cooccurrence.csv"))

k <- kin |> group_by(species, experiment, condition) |>
  summarise(shaken = mean(endpoint), n_shaken = n(), .groups = "drop")
s <- sta |> group_by(species, experiment) |>
  summarise(static = mean(endpoint), n_static = n(), .groups = "drop")

delta <- inner_join(k, s, by = c("species", "experiment")) |>
  mutate(delta = shaken - static) |>
  arrange(species, condition, experiment)

cat("\n\n=========================================================================\n")
cat("2. PER-RUN PAIRED DELTA:  shaken - its own same-day Static\n")
cat("=========================================================================\n")
print(as.data.frame(delta |> mutate(across(c(shaken, static, delta), \(x) round(x, 4)))),
      row.names = FALSE)
write_csv(delta, file.path(DATA_DIR, "01_perrun_paired_deltas.csv"))

repro <- delta |>
  group_by(species, condition) |>
  summarise(
    n_runs     = n(),
    deltas     = paste(sprintf("%+.3f", delta), collapse = "  "),
    mean_delta = round(mean(delta), 4),
    spread     = round(max(delta) - min(delta), 4),
    consistent_sign = n_distinct(sign(delta)) == 1,
    .groups = "drop"
  ) |>
  mutate(verdict = case_when(
    n_runs < 2                 ~ "NOT REPLICATED - cannot estimate",
    !consistent_sign           ~ "NOT REPRODUCIBLE - runs disagree on sign",
    spread > abs(mean_delta)   ~ "WEAK - spread exceeds the effect",
    TRUE                       ~ "reproducible"
  ))

cat("\n\n=========================================================================\n")
cat("3. IS THE REGIME EFFECT REPRODUCIBLE ACROSS RUNS OF THE SAME REGIME?\n")
cat("=========================================================================\n")
print(as.data.frame(repro), row.names = FALSE)
write_csv(repro, file.path(DATA_DIR, "01_regime_reproducibility.csv"))

route_rows <- list()
for (sp in SPECIES) {
  d <- dat |> filter(species == sp) |> droplevels()
  dd <- delta |> filter(species == sp)

  raw_P  <- dd |> filter(condition == "Pulsed")     |> pull(shaken) |> mean()
  raw_C  <- dd |> filter(condition == "Continuous") |> pull(shaken) |> mean()
  pair_P <- dd |> filter(condition == "Pulsed")     |> pull(delta)  |> mean()
  pair_C <- dd |> filter(condition == "Continuous") |> pull(delta)  |> mean()
  base_P <- dd |> filter(condition == "Pulsed")     |> pull(static) |> mean()
  base_C <- dd |> filter(condition == "Continuous") |> pull(static) |> mean()

  mA   <- suppressMessages(suppressWarnings(lmer(endpoint ~ condition + (1 | experiment), data = d)))
  cA   <- as.data.frame(pairs(emmeans(mA, ~condition), adjust = "none")) |>
            filter(contrast == "Pulsed - Continuous")
  bridge_est <- cA$estimate
  bridge_p   <- cA$p.value

  route_rows[[sp]] <- tibble(
    species              = sp,
    static_baseline_gap  = round(base_P - base_C, 4),
    routeA_raw_PminusC   = round(raw_P - raw_C, 4),
    routeB_paired_PminusC= round(pair_P - pair_C, 4),
    routeC_bridge_PminusC= round(bridge_est, 4),
    routeC_p             = round(bridge_p, 4),
    routeC_singular      = isSingular(mA),
    SIGN_DISAGREEMENT    = n_distinct(sign(c(raw_P - raw_C, pair_P - pair_C, bridge_est))) > 1
  )
}
routes <- bind_rows(route_rows)

cat("\n\n=========================================================================\n")
cat("4. P-vs-C: THREE ROUTES. DO THEY AGREE?\n")
cat("=========================================================================\n")
print(as.data.frame(routes), row.names = FALSE)
write_csv(routes, file.path(DATA_DIR, "01_PvsC_route_comparison.csv"))

bad <- routes$species[routes$SIGN_DISAGREEMENT]
cat("\n>>> Species where the routes disagree on the SIGN of P-C: ",
    paste(bad, collapse = ", "), "\n")
cat(">>> `static_baseline_gap` is the difference between the Pulsed runs' and the\n")
cat(">>>  Continuous runs' Static plates - i.e. how much the runs differed BEFORE\n")
cat(">>>  any shaking. Where this is comparable to routeA, the 'regime effect' IS\n")
cat(">>>  the batch effect.\n\n")

cat("CONCLUSION FOR THE METHODS SECTION:\n")
cat("  Pulsed-vs-Continuous is not identifiable in this design. Report Static-vs-\n")
cat("  Pulsed and Static-vs-Continuous (both within-run), and treat Pulsed-vs-\n")
cat("  Continuous as a design limitation, not a null result.\n\n")

sessionInfo()
