# Load the annotated MASTER_* tables, keep included wells, and check the design table

library(tidyverse)

if (!exists("DATA_DIR")) DATA_DIR <- "."

kin_all <- read_csv(file.path(DATA_DIR, "MASTER_kinetic_annotated.csv"), show_col_types = FALSE)
sta_all <- read_csv(file.path(DATA_DIR, "MASTER_static_annotated.csv"),  show_col_types = FALSE)

kin <- kin_all %>% filter(status == "Included")
sta <- sta_all %>% filter(status == "Included")

dat <- bind_rows(
  kin %>% select(experiment, species, condition, well, replicate, endpoint),
  sta %>% select(experiment, species, condition, well, replicate, endpoint)
) %>%
  mutate(condition  = factor(condition, levels = c("Static", "Pulsed", "Continuous")),
         experiment = factor(experiment),
         species    = factor(species))

load_timeseries <- function() {
  read_csv(file.path(DATA_DIR, "MASTER_timeseries_annotated.csv"), show_col_types = FALSE) %>%
    filter(status == "Included") %>%
    mutate(condition = factor(condition, levels = c("Static", "Pulsed", "Continuous")))
}

excluded_log <- bind_rows(
  kin_all %>% filter(status != "Included") %>% mutate(plate = "Kinetic") %>%
    select(plate, species, experiment, exclusion_category, exclusion_reason),
  sta_all %>% filter(status != "Included") %>% mutate(plate = "Static") %>%
    select(plate, species, experiment, exclusion_category, exclusion_reason)
) %>% distinct()

expected <- tribble(
  ~species, ~Pulsed, ~Continuous, ~Static,
  "BT", 3L, 3L, 5L,   "EC", 3L, 3L, 5L,   "FD", 2L, 2L, 4L,
  "MG", 2L, 2L, 3L,   "RI", 2L, 2L, 4L,   "SC", 2L, 3L, 5L
)
got <- dat %>%
  group_by(species, condition) %>%
  summarise(n = n_distinct(experiment), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = n, values_fill = 0L) %>%
  select(species, Pulsed, Continuous, Static) %>%
  mutate(species = as.character(species)) %>%
  arrange(species) %>% as.data.frame()

cat("===== experiments per species x condition (included only) =====\n")
print(got)

if (isTRUE(all.equal(got, as.data.frame(expected), check.attributes = FALSE))) {
  cat("\n*** Design MATCHES the study design table. Inputs are correct. ***\n")
} else {
  cat("\nEXPECTED:\n"); print(as.data.frame(expected))
  stop("Design does NOT match the study design table. Check the annotated master ",
       "files before trusting any downstream result.")
}

cat(sprintf("\nkinetic wells: %d included of %d   |   static wells: %d included of %d\n",
            nrow(kin), nrow(kin_all), nrow(sta), nrow(sta_all)))
cat(sprintf("excluded plates: %d (see `excluded_log`)\n\n", nrow(excluded_log)))
