# Within-run endpoint-OD contrasts per regime

suppressMessages({
  library(tidyverse); library(lme4); library(lmerTest); library(emmeans)
})
emm_options(lmerTest.limit = 100000, pbkrtest.limit = 100000)

if (!exists("DATA_DIR")) DATA_DIR <- "."
SPECIES <- c("BT", "EC", "FD", "MG", "RI", "SC")
K_BONF  <- 2

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

fit_within_run <- function(d_sp, regime) {

  runs <- d_sp |>
    filter(condition %in% c("Static", regime)) |>
    distinct(experiment, condition) |>
    count(experiment) |>
    filter(n == 2) |>
    pull(experiment)

  d <- d_sp |>
    filter(experiment %in% runs, condition %in% c("Static", regime)) |>
    mutate(condition = factor(condition, levels = c("Static", regime))) |>
    droplevels()

  n_runs <- n_distinct(d$experiment)

  if (n_runs < 2) {
    return(tibble(
      contrast = paste(regime, "- Static"), n_runs = n_runs, n_wells = nrow(d),
      estimate = NA_real_, SE = NA_real_, df = NA_real_, p_raw = NA_real_,
      p_bonf = NA_real_, singular = NA, sens_estimate = NA_real_, sens_p = NA_real_,
      agree = NA, flag = "NOT ESTIMABLE - only 1 run carries this regime"
    ))
  }

  m  <- suppressMessages(suppressWarnings(
          lmer(endpoint ~ condition + (1 | experiment), data = d)))
  em <- emmeans(m, ~condition)
  ct <- as.data.frame(pairs(em, adjust = "none"))
  est <- -ct$estimate
  se  <-  ct$SE
  df  <-  ct$df
  p   <-  ct$p.value

  m2  <- lm(endpoint ~ condition + experiment, data = d)
  ct2 <- as.data.frame(pairs(emmeans(m2, ~condition), adjust = "none"))
  est2 <- -ct2$estimate
  p2   <-  ct2$p.value

  tibble(
    contrast      = paste(regime, "- Static"),
    n_runs        = n_runs,
    n_wells       = nrow(d),
    estimate      = round(est, 4),
    SE            = round(se, 4),
    df            = round(df, 1),
    p_raw         = p,
    p_bonf        = pmin(p * K_BONF, 1),
    singular      = isSingular(m),
    sens_estimate = round(est2, 4),
    sens_p        = p2,

    agree         = {
      sig1 <- pmin(p  * K_BONF, 1) < 0.05
      sig2 <- pmin(p2 * K_BONF, 1) < 0.05
      (sig1 == sig2) && (!(sig1 || sig2) || sign(est) == sign(est2))
    },
    flag          = ""
  )
}

res <- map_dfr(SPECIES, function(sp) {
  d_sp <- dat |> filter(species == sp) |> droplevels()
  bind_rows(
    fit_within_run(d_sp, "Pulsed"),
    fit_within_run(d_sp, "Continuous")
  ) |> mutate(species = sp, .before = 1)
})

stars <- function(p) ifelse(is.na(p), "",
                     ifelse(p < 0.001, "***",
                     ifelse(p < 0.01,  "**",
                     ifelse(p < 0.05,  "*", "n.s."))))
res <- res |> mutate(sig = stars(p_bonf))

cat("\n=========================================================================\n")
cat("ENDPOINT OD - within-run contrasts (Bonferroni k = ", K_BONF, " per species)\n", sep = "")
cat("=========================================================================\n")
print(as.data.frame(res |> select(species, contrast, n_runs, n_wells, estimate, SE, df,
                                  p_bonf, sig, singular)), row.names = FALSE)

cat("\n--- SENSITIVITY: does a FIXED-effect block model give the same answer? ---\n")
print(as.data.frame(res |> select(species, contrast, estimate, sens_estimate,
                                  p_bonf, sens_p, agree)), row.names = FALSE)

disagree <- res |> filter(!is.na(agree), !agree)
if (nrow(disagree) == 0) {
  cat("\n*** Mixed and fixed-effect models AGREE on every estimable contrast. ***\n")
  cat("*** The result does not depend on how `experiment` is modelled. ***\n")
} else {
  cat("\n!!! DISAGREEMENT between mixed and fixed models - inspect before reporting:\n")
  print(as.data.frame(disagree |> select(species, contrast, estimate, sens_estimate, p_bonf, sens_p)))
}

if (any(res$singular %in% TRUE)) {
  cat("\nNOTE: singular fit(s) in:",
      paste(res$species[res$singular %in% TRUE], res$contrast[res$singular %in% TRUE],
            collapse = "; "),
      "\n  -> the between-run variance is estimated as ~0. Lean on the fixed-effect\n",
      "     sensitivity estimate for these; report both.\n")
}
nb <- res |> filter(!is.na(flag), flag != "")
if (nrow(nb) > 0) {
  cat("\nNOT ESTIMABLE:\n"); print(as.data.frame(nb |> select(species, contrast, flag)), row.names = FALSE)
}

write_csv(res, file.path(DATA_DIR, "02_endpoint_stats.csv"))

cells <- dat |>
  group_by(species, condition) |>
  summarise(n_runs = n_distinct(experiment), n_wells = n(),
            mean = round(mean(endpoint), 4), sd = round(sd(endpoint), 4),
            .groups = "drop")
cat("\n--- cell means for Figure 2A ---\n")
print(as.data.frame(cells), row.names = FALSE)
write_csv(cells, file.path(DATA_DIR, "02_endpoint_celldata.csv"))

cat("\nWrote 02_endpoint_stats.csv and 02_endpoint_celldata.csv\n")
cat("REMINDER: do not add a Pulsed-vs-Continuous test to this figure.\n\n")
sessionInfo()
