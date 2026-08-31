# Per-species OD-to-cell calibration fits

suppressMessages({library(tidyverse)})
if (!exists("DATA_DIR")) DATA_DIR <- "."

cal <- read_csv(file.path(DATA_DIR, "calibration_metadata.csv"), show_col_types = FALSE) |>
  filter(use_for_fit == TRUE)

fit_one <- function(d) {
  m <- lm(cells_per_mL ~ target_OD, data = d)
  tibble(intercept   = coef(m)[1],
         slope       = coef(m)[2],
         R2          = summary(m)$r.squared,
         n_points    = nrow(d),
         n_OD        = n_distinct(d$target_OD),
         media       = paste(unique(d$media), collapse = "/"))
}
fits <- cal |> group_by(species) |> group_modify(~ fit_one(.x)) |> ungroup() |>
  arrange(desc(slope))

sc_slope <- fits$slope[fits$species == "SC"]
ratios <- fits |>
  transmute(species, media, slope,
            cells_per_OD_relative_to_SC = round(slope / sc_slope, 2),
            R2 = round(R2, 4))

write_csv(fits,   file.path(DATA_DIR, "08_calibration_fits.csv"))
write_csv(ratios, file.path(DATA_DIR, "08_crossspecies_ratios.csv"))

cat("=== per-species calibration fits (cells/mL ~ OD) ===\n")
fits |> mutate(slope = format(round(slope), big.mark=","),
               intercept = format(round(intercept), big.mark=","),
               R2 = round(R2,4)) |> as.data.frame() |> print(row.names = FALSE)

cat("\n=== cross-species: cells per OD, relative to SC ===\n")
as.data.frame(ratios) |> print(row.names = FALSE)

cat("\nWC-only spread (the headline, FD excluded as it is in YCFAG):\n")
wc <- ratios |> filter(media == "WC")
cat(sprintf("  highest %s (%.2fx)  ->  lowest %s (%.2fx)  =  %.1f-fold range\n",
            wc$species[which.max(wc$cells_per_OD_relative_to_SC)],
            max(wc$cells_per_OD_relative_to_SC),
            wc$species[which.min(wc$cells_per_OD_relative_to_SC)],
            min(wc$cells_per_OD_relative_to_SC),
            max(wc$cells_per_OD_relative_to_SC)/min(wc$cells_per_OD_relative_to_SC)))
