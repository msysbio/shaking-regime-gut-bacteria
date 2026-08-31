# F. duncaniae growth curves at two starting ODs

library(tidyverse)
library(lme4)
library(lmerTest)

DATA_DIR <- "."

source(file.path(DATA_DIR, "00_load_and_filter.R"))

FD_EXPS <- c("EXP-02", "EXP-03", "EXP-06", "EXP-07")

lg <- load_timeseries() %>%
  filter(species == "FD", experiment %in% FD_EXPS) %>%
  arrange(experiment, well, time_h)

curve_params <- function(d) {
  d <- d %>% arrange(time_h)
  t  <- d$time_h
  od <- d$OD_corr
  t0 <- od[1]
  max_OD <- max(od, na.rm = TRUE)

  auc <- sum(diff(t) * (head(od, -1) + tail(od, -1)) / 2, na.rm = TRUE)

  odf <- pmax(od, 1e-3)
  lod <- log(odf)
  win <- 5
  mu  <- rep(NA_real_, length(t))
  if (length(t) > win) {
    for (i in seq_len(length(t) - win)) {
      j <- i + win
      mu[i] <- (lod[j] - lod[i]) / (t[j] - t[i])
    }
  }
  mu_max <- suppressWarnings(max(mu, na.rm = TRUE))
  t_mumax <- t[which.max(mu)]

  mid <- t0 + 0.5 * (max_OD - t0)
  t_mid <- suppressWarnings(t[which(od >= mid)[1]])
  tibble(t0 = t0, max_OD = max_OD, AUC = auc,
         mu_max = mu_max, t_mumax = t_mumax, t_mid = t_mid)
}

params <- lg %>%
  group_by(experiment, condition, well) %>%
  group_modify(~ curve_params(.x)) %>%
  ungroup() %>%
  mutate(experiment = factor(experiment),
         condition  = factor(condition, levels = c("Pulsed", "Continuous")))

cat("===== per-well curve parameters (FD) =====\n")
print(params, n = 100)
write_csv(params, file.path(DATA_DIR, "FD_curve_params.csv"))

cat("\n===== starting OD (t0) by experiment =====\n")
print(params %>% group_by(experiment, condition) %>%
        summarise(t0_mean = round(mean(t0), 4), t0_sd = round(sd(t0), 4),
                  max_OD = round(mean(max_OD), 3), AUC = round(mean(AUC), 2),
                  mu_max = round(mean(mu_max), 3), .groups = "drop"))

test_outcome <- function(varname) {
  cat("\n==================", varname, "==================\n")
  f_add <- as.formula(paste(varname, "~ t0 + condition + (1 | experiment)"))
  m <- suppressWarnings(lmer(f_add, data = params))
  print(summary(m)$coefficients)
  cat("--- t0 x condition interaction test ---\n")
  f_int <- as.formula(paste(varname, "~ t0 * condition + (1 | experiment)"))
  m2 <- suppressWarnings(lmer(f_int, data = params))
  print(anova(m, m2))
  invisible(NULL)
}

sink(file.path(DATA_DIR, "FD_startingOD_models.txt"), split = TRUE)
for (v in c("max_OD", "AUC", "mu_max")) test_outcome(v)
sink()

cat("\nWrote FD_curve_params.csv and FD_startingOD_models.txt\n")
cat("Read the t0 row in each model: significant => starting OD shapes that outcome.\n\n")
sessionInfo()
