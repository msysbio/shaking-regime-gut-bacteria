# Between-plate (run) variance components

suppressMessages({library(tidyverse); library(lme4)})

if (!exists("DATA_DIR")) DATA_DIR <- "."
SPECIES  <- c("BT", "EC", "MG", "RI", "SC")
REGIMES  <- c("Static", "Pulsed", "Continuous")

kin <- read_csv(file.path(DATA_DIR, "MASTER_kinetic_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
sta <- read_csv(file.path(DATA_DIR, "MASTER_static_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")

dat <- bind_rows(
  kin |> select(experiment, species, condition, well, endpoint),
  sta |> select(experiment, species, condition, well, endpoint)
) |> mutate(experiment = factor(experiment))

icc_cell <- function(d) {
  n_plates <- n_distinct(d$experiment)
  if (n_plates < 2) {
    return(tibble(n_plates = n_plates, n_wells = nrow(d), ICC = NA_real_,
                  var_plate = NA_real_, var_resid = NA_real_, p = NA_real_,
                  singular = NA, note = "only 1 plate - no between-plate variance to estimate"))
  }
  m  <- suppressMessages(suppressWarnings(lmer(endpoint ~ 1 + (1 | experiment), data = d,
                                               REML = TRUE)))
  vc <- as.data.frame(VarCorr(m))
  v_plate <- vc$vcov[vc$grp == "experiment"]
  v_resid <- vc$vcov[vc$grp == "Residual"]
  icc <- v_plate / (v_plate + v_resid)

  m1 <- suppressMessages(suppressWarnings(lmer(endpoint ~ 1 + (1 | experiment), data = d, REML = FALSE)))
  m0 <- lm(endpoint ~ 1, data = d)

  lrt <- as.numeric(2 * (logLik(m1) - logLik(m0)))
  p   <- 0.5 * pchisq(max(lrt, 0), df = 1, lower.tail = FALSE)

  tibble(n_plates = n_plates, n_wells = nrow(d),
         ICC = round(icc, 3), var_plate = signif(v_plate, 3), var_resid = signif(v_resid, 3),
         p = round(p, 4), singular = isSingular(m), note = "")
}

res <- expand_grid(species = SPECIES, condition = REGIMES) |>
  mutate(fit = map2(species, condition, \(s, c) {
    d <- dat |> filter(species == s, condition == c)
    if (nrow(d) == 0) return(tibble(n_plates = 0L, n_wells = 0L, ICC = NA_real_,
                                    var_plate = NA_real_, var_resid = NA_real_,
                                    p = NA_real_, singular = NA, note = "no data"))
    icc_cell(d)
  })) |>
  unnest(fit) |>
  mutate(sig = case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
                         p < 0.05 ~ "*", TRUE ~ "n.s."),
         batch_pct = ifelse(is.na(ICC), NA, round(100 * ICC)))

cat("\n=========================================================================\n")
cat("BATCH EFFECT: % of endpoint-OD variance that is plate-to-plate (within regime)\n")
cat("=========================================================================\n\n")
print(as.data.frame(res |> select(species, condition, n_plates, n_wells,
                                  batch_pct, p, sig, singular, note)), row.names = FALSE)

cat("\n--- summary grid (% between-plate variance) ---\n")
grid <- res |>
  mutate(cell = ifelse(is.na(batch_pct), "  -  ", sprintf("%3d%% %s", batch_pct, sig))) |>
  select(species, condition, cell) |>
  pivot_wider(names_from = condition, values_from = cell)
print(as.data.frame(grid), row.names = FALSE)

write_csv(res, file.path(DATA_DIR, "04_batch_effects.csv"))

cat("\n\nHOW TO READ THIS:\n")
cat("  A high % means: which plate a well sat on explains most of its endpoint OD.\n")
cat("  Wherever that is true, ANY comparison made between plates - which is what\n")
cat("  Pulsed-vs-Continuous necessarily is - is unreliable. This table is the\n")
cat("  quantitative justification for withdrawing the P-vs-C biomass comparison.\n")
cat("\nCAVEATS TO STATE IN THE PAPER:\n")
cat("  - 2-3 plates per cell. These ICCs are themselves imprecisely estimated.\n")
cat("  - BT/EC batch is partly confounded with TIME (Jan runs EXP-06/07 vs spring runs).\n")
cat("  - FD excluded: its plates differ by target-OD design, not by batch.\n\n")
sessionInfo()
