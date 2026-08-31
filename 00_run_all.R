# Run scripts 01-05 in order, with checks and logging

DATA_DIR <- "/Users/u0176884/Desktop/final_metadata_shaking_regime"

setwd(DATA_DIR)

needed <- c("tidyverse", "lme4", "lmerTest", "emmeans", "readxl")
missing <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
invisible(lapply(needed, \(p) suppressMessages(library(p, character.only = TRUE))))

inputs <- c(
  "MASTER_kinetic_annotated.csv",
  "MASTER_static_annotated.csv",
  "MASTER_timeseries_annotated.csv",
  "SC_RI_flow_cytometry_master.xlsx"
)
scripts <- c(
  "01_design_and_identifiability.R",
  "02_endpoint_OD.R",
  "03_PvsC_paired.R",
  "04_batch_effects.R",
  "05_CoV_reproducibility.R"
)

missing_files <- c(inputs, scripts)[!file.exists(c(inputs, scripts))]
if (length(missing_files)) {
  stop("Missing from ", DATA_DIR, ":\n  - ",
       paste(missing_files, collapse = "\n  - "),
       "\n\n(If a data file has a numeric prefix on it, strip the prefix.)")
}
cat("All", length(inputs), "data files and", length(scripts), "scripts found.\n\n")

log_file <- file.path(DATA_DIR, "analysis_log.txt")
con <- file(log_file, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")

cat("SHAKING REGIME ANALYSIS —", format(Sys.time()), "\n")
cat(R.version.string, "\n")
cat(strrep("=", 75), "\n")

results <- data.frame()
for (s in scripts) {
  cat("\n\n", strrep("#", 75), "\n### RUNNING: ", s, "\n", strrep("#", 75), "\n", sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({
    env <- new.env(parent = globalenv())
    source(s, local = env, echo = FALSE)
    TRUE
  }, error = function(e) {
    cat("\n*** FAILED:", conditionMessage(e), "\n")
    FALSE
  })
  results <- rbind(results, data.frame(
    script = s, status = ifelse(ok, "OK", "FAILED"),
    seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)))
  if (!ok) break
}

cat("\n\n", strrep("=", 75), "\nSUMMARY\n", strrep("=", 75), "\n", sep = "")
print(results, row.names = FALSE)

if (all(results$status == "OK")) {
  cat("\nAll scripts completed. Result files written to:\n  ", DATA_DIR, "\n")
  cat("\nKey outputs to look at, in this order:\n")
  cat("  01_PvsC_route_comparison.csv  <- why P-vs-C is not identifiable\n")
  cat("  02_endpoint_stats.csv         <- Fig 2A stats (the solid contrasts)\n")
  cat("  03_FC_tests.csv               <- flow cytometry (SC 1.25x, RI null)\n")
  cat("  04_batch_effects.csv          <- Fig 3 (why P-vs-C had to go)\n")
  cat("  05_CoV_tests.csv              <- the headline: Pulsed 6/6 more reproducible\n")
} else {
  cat("\nSTOPPED at the first failure. Fix it and re-run; nothing downstream ran.\n")
}
cat("\nFull transcript: ", log_file, "\n", sep = "")

sink(type = "message"); sink(); close(con)
