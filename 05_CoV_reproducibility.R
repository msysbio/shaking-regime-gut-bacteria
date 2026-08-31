# Within-plate CoV reproducibility

suppressMessages(library(tidyverse))

if (!exists("DATA_DIR")) DATA_DIR <- "."
SPECIES <- c("BT", "EC", "FD", "MG", "RI", "SC")
T_MIN   <- 4

cv <- function(x) sd(x) / mean(x) * 100

kin <- read_csv(file.path(DATA_DIR, "MASTER_kinetic_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
sta <- read_csv(file.path(DATA_DIR, "MASTER_static_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")
ts  <- read_csv(file.path(DATA_DIR, "MASTER_timeseries_annotated.csv"), show_col_types = FALSE) |>
  filter(status == "Included")

endpoint_cov <- bind_rows(
    kin |> select(species, experiment, condition, endpoint),
    sta |> select(species, experiment, condition, endpoint)
  ) |>
  group_by(species, experiment, condition) |>
  filter(n() >= 3, mean(endpoint) > 0) |>
  summarise(n_wells = n(), mean_OD = round(mean(endpoint), 4),
            CoV_endpoint = round(cv(endpoint), 2), .groups = "drop")

kinetic_cov <- ts |>
  filter(time_h >= T_MIN) |>
  group_by(species, experiment, condition, time_h) |>
  filter(n() >= 3, mean(OD_corr) > 0) |>
  summarise(cov_t = cv(OD_corr), .groups = "drop") |>
  group_by(species, experiment, condition) |>
  summarise(n_timepoints = n(), CoV_kinetic = round(mean(cov_t), 2), .groups = "drop")

per_plate <- full_join(endpoint_cov, kinetic_cov,
                       by = c("species", "experiment", "condition")) |>
  arrange(species, condition, experiment)

cat("\n=========================================================================\n")
cat("PER-PLATE CoV (%)   [kinetic is NA for Static - no growth curve]\n")
cat("=========================================================================\n")
print(as.data.frame(per_plate), row.names = FALSE)
write_csv(per_plate, file.path(DATA_DIR, "05_CoV_per_plate.csv"))

summ <- per_plate |>
  group_by(species, condition) |>
  summarise(n_plates = n(),
            endpoint_CoV = round(mean(CoV_endpoint, na.rm = TRUE), 2),
            endpoint_vals = paste(sprintf("%.1f", CoV_endpoint), collapse = ", "),
            kinetic_CoV  = round(mean(CoV_kinetic, na.rm = TRUE), 2),
            .groups = "drop")

cat("\n\n=========================================================================\n")
cat("PER-SPECIES CoV - DESCRIPTIVE ONLY (2-3 plates per cell, no power)\n")
cat("=========================================================================\n")
print(as.data.frame(summ), row.names = FALSE)
write_csv(summ, file.path(DATA_DIR, "05_CoV_summary.csv"))

paired_test <- function(metric) {
  w <- summ |>
    filter(condition %in% c("Pulsed", "Continuous")) |>
    select(species, condition, val = all_of(metric)) |>
    pivot_wider(names_from = condition, values_from = val) |>
    filter(!is.na(Pulsed), !is.na(Continuous)) |>
    mutate(diff = Pulsed - Continuous)

  n <- nrow(w)
  wt <- suppressWarnings(wilcox.test(w$Pulsed, w$Continuous, paired = TRUE))
  tt <- t.test(w$diff)
  favour_pulsed <- sum(w$diff < 0)

  list(table = w,
       row = tibble(metric = metric, n_species = n,
                    favours_Pulsed = sprintf("%d/%d", favour_pulsed, n),
                    mean_diff_pp   = round(mean(w$diff), 2),
                    wilcox_p       = round(wt$p.value, 4),
                    wilcox_p_floor = round(2 / 2^n, 4),
                    ttest_p        = round(tt$p.value, 4),
                    sign_test_p    = round(binom.test(favour_pulsed, n)$p.value, 4)))
}

cat("\n\n=========================================================================\n")
cat("PULSED vs CONTINUOUS, PAIRED ACROSS SPECIES (the only test with any power)\n")
cat("=========================================================================\n")
rows <- list()
for (m in c("endpoint_CoV", "kinetic_CoV")) {
  r <- paired_test(m)
  cat("\n---", m, "--- (negative diff = Pulsed more reproducible)\n")
  print(as.data.frame(r$table |> mutate(across(where(is.numeric), \(x) round(x, 2)))),
        row.names = FALSE)
  rows[[m]] <- r$row
}
tests <- bind_rows(rows)
cat("\n")
print(as.data.frame(tests), row.names = FALSE)
write_csv(tests, file.path(DATA_DIR, "05_CoV_tests.csv"))

cat("\n\nHOW TO REPORT THIS:\n")
cat("  - Per-species CoV: DESCRIPTIVE. 2-3 plates per cell cannot support a test.\n")
cat("  - The across-species paired test is the inferential claim. Note the p-floor:\n")
cat("    with n species, two-sided Wilcoxon cannot return below 2/2^n.\n")
cat("  - CoV cancels plate-level LEVEL shifts, not plate-level NOISE. State that.\n")
cat("  - Do NOT write 'Pulsed lowers CoV in all six species' unless every diff < 0.\n\n")
sessionInfo()
