# Pulsed-vs-continuous paired analysis and flow cytometry

suppressMessages({library(tidyverse); library(readxl)})

if (!exists("DATA_DIR")) DATA_DIR <- "."
FC_FILE <- file.path(DATA_DIR, "SC_RI_flow_cytometry_master.xlsx")
SPECIES <- c("BT", "EC", "FD", "MG", "RI", "SC")

cat("\n#########################################################################\n")
cat("PART A - ENDPOINT OD: Pulsed vs Continuous, difference-in-differences\n")
cat("#########################################################################\n")

d <- read_csv(file.path(DATA_DIR, "01_perrun_paired_deltas.csv"), show_col_types = FALSE)

did <- map_dfr(SPECIES, function(sp) {
  x <- d |> filter(species == sp)
  P <- x |> filter(condition == "Pulsed")     |> pull(delta)
  C <- x |> filter(condition == "Continuous") |> pull(delta)

  if (length(P) < 2 || length(C) < 2) {
    return(tibble(species = sp, n_P_runs = length(P), n_C_runs = length(C),
                  P_delta = NA_real_, C_delta = NA_real_, difference = NA_real_,
                  p_raw = NA_real_, p_bonf = NA_real_, deltas_overlap = NA,
                  verdict = "NOT ESTIMABLE - fewer than 2 runs in one arm"))
  }

  tt <- t.test(P, C)
  tibble(species = sp, n_P_runs = length(P), n_C_runs = length(C),
         P_delta = round(mean(P), 4), C_delta = round(mean(C), 4),
         difference = round(mean(P) - mean(C), 4),
         p_raw = round(tt$p.value, 4),
         p_bonf = NA_real_,
         deltas_overlap = !(max(P) < min(C) || min(P) > max(C)),
         verdict = "")
})

k <- sum(!is.na(did$p_raw))
did <- did |>
  mutate(p_bonf  = ifelse(is.na(p_raw), NA, pmin(p_raw * k, 1)),
         verdict = case_when(
           is.na(p_raw)      ~ verdict,
           p_bonf   < 0.05   ~ "difference survives correction",
           p_raw    < 0.05   ~ "nominal only - does NOT survive correction",
           TRUE              ~ "no difference"))

cat("\nRun-level deltas (shaken - own Static), Bonferroni k =", k, "species tested\n\n")
print(as.data.frame(did), row.names = FALSE)
write_csv(did, file.path(DATA_DIR, "03_OD_PvsC_did.csv"))

cat("\nREAD THIS AS: with 2-3 runs per regime there is almost no power. A null here\n")
cat("means 'this design cannot resolve it', not 'the regimes are equivalent'.\n")

cat("\n\n#########################################################################\n")
cat("PART B - FLOW CYTOMETRY: per-plate fold change (Shaking / Static)\n")
cat("#########################################################################\n")

read_fc <- function(sheet) {
  read_excel(FC_FILE, sheet = sheet, skip = 3) |>
    rename(plate = Plate, regime = Shaking_regime, species = Species,
           condition = Condition, well = Well,
           cells = `Cells_per_mL (formula)`, excluded = Excluded) |>
    filter(!is.na(plate), tolower(excluded) != "yes") |>
    mutate(cells = as.numeric(cells)) |>
    select(plate, regime, species, condition, well, cells)
}
fc <- bind_rows(read_fc("SC_Raw_Wells"), read_fc("RI_Raw_Wells"))

cat("\nwells used (excluded wells already dropped):\n")
print(as.data.frame(fc |> count(species, plate, regime, condition)), row.names = FALSE)

pm <- fc |>
  filter(condition != "Start") |>
  group_by(species, plate, regime, condition) |>
  summarise(cells = mean(cells), n = n(), .groups = "drop") |>
  mutate(arm = if_else(condition == "Static", "static", "shaken")) |>
  select(species, plate, regime, arm, cells) |>
  pivot_wider(names_from = arm, values_from = cells) |>
  mutate(fold_change = shaken / static,
         log2_fc     = log2(fold_change))

cat("\nPer-plate fold change (Shaking / Static):\n")
print(as.data.frame(pm |> mutate(across(c(static, shaken), \(x) signif(x, 3)),
                                 fold_change = round(fold_change, 3),
                                 log2_fc = round(log2_fc, 3))), row.names = FALSE)
write_csv(pm, file.path(DATA_DIR, "03_FC_plate_foldchange.csv"))

set.seed(1)
boot_gm_ci <- function(x, B = 10000) {
  gm <- \(v) exp(mean(log(v)))
  bs <- replicate(B, gm(sample(x, length(x), replace = TRUE)))
  c(gm = gm(x), lo = unname(quantile(bs, 0.025)), hi = unname(quantile(bs, 0.975)))
}

fc_tests <- map_dfr(c("SC", "RI"), function(sp) {
  x  <- pm |> filter(species == sp)
  fcv <- x$fold_change
  ci  <- boot_gm_ci(fcv)

  w  <- suppressWarnings(wilcox.test(x$log2_fc, mu = 0, alternative = "two.sided"))
  tt <- t.test(x$log2_fc, mu = 0)

  P <- x |> filter(regime == "Pulsed")     |> pull(log2_fc)
  C <- x |> filter(regime == "Continuous") |> pull(log2_fc)
  pvc <- if (length(P) >= 2 && length(C) >= 2)
           round(t.test(P, C)$p.value, 3) else NA_real_

  tibble(species = sp, n_plates = nrow(x),
         geom_mean_FC = round(ci["gm"], 4),
         CI_lo = round(ci["lo"], 4), CI_hi = round(ci["hi"], 4),
         concordant = sprintf("%d/%d plates %s 1.0",
                              sum(fcv > 1), length(fcv), ">"),
         wilcox_p_2sided = round(w$p.value, 4),
         wilcox_p_floor  = round(2 / 2^nrow(x), 4),
         ttest_log2_p    = round(tt$p.value, 4),
         PvsC_p          = pvc)
})

cat("\n--- Shaking vs Static (the headline test) + Pulsed vs Continuous ---\n")
print(as.data.frame(fc_tests), row.names = FALSE)
write_csv(fc_tests, file.path(DATA_DIR, "03_FC_tests.csv"))

cat("\nNOTE ON THE WILCOXON: with n plates, the smallest two-sided p obtainable is\n")
cat("2/2^n (see wilcox_p_floor). At n = 5 that is 0.0625; a one-sided test gives\n")
cat("0.031. A p at the floor means ONLY 'every plate went the same way' - the\n")
cat("concordance count and the bootstrap CI are the real evidence. Say which tail\n")
cat("you used in the manuscript, and prefer the CI.\n")

cat("\n\n#########################################################################\n")
cat("WHY THERE IS NO AUC / max-OD / mu_max SECTION IN THIS SCRIPT\n")
cat("#########################################################################\n")
cat("The Static plates were read at ENDPOINT ONLY - they have no growth curve.\n")
cat("Every within-run control in this study is an endpoint value. A kinetic metric\n")
cat("therefore has no same-run reference to be normalised against, so AUC, max OD\n")
cat("and mu_max can only ever be compared BETWEEN runs, where Pulsed and Continuous\n")
cat("are perfectly confounded with batch. No statistical framework recovers this;\n")
cat("the information is not in the data. Report kinetic curves descriptively, and\n")
cat("recommend running the Static control kinetically in future work.\n\n")
sessionInfo()
