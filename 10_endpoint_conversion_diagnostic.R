# Endpoint calibration-extrapolation diagnostic (Fig S7)

suppressMessages({library(tidyverse); library(readxl)})
if (!exists("DATA_DIR")) DATA_DIR <- "."
NAVY <- "#1E2761"; ACCENT <- "#B85042"

fits <- read_csv(file.path(DATA_DIR,"08_calibration_fits.csv"), show_col_types=FALSE)

kin <- read_csv(file.path(DATA_DIR,"MASTER_kinetic_annotated.csv"), show_col_types=FALSE) |>
  filter(status=="Included")
endpoint_OD <- kin |> group_by(species) |> summarise(OD=mean(endpoint), .groups="drop")

pred <- endpoint_OD |> inner_join(fits |> select(species,slope,intercept), by="species") |>
  mutate(predicted_cells_mL = intercept + slope*OD,
         fold_beyond_calib   = round(OD/0.015, 0))

read_static <- function(sheet){
  s <- read_excel(file.path(DATA_DIR,"SC_RI_flow_cytometry_master.xlsx"),
                  sheet=sheet, skip=4)
  names(s)[1:9] <- c("Plate","Shaking_regime","Condition","n","mean_ev","sd_ev","cv","mean_cells","sd_cells")
  s |> filter(Condition=="Static") |>
    summarise(measured_cells_mL=mean(as.numeric(mean_cells), na.rm=TRUE)) |> pull()
}
measured <- tibble(
  species = c("SC","RI"),
  measured_cells_mL = c(read_static("SC_Plate_Summary"), read_static("RI_Plate_Summary")))

diag <- pred |> left_join(measured, by="species") |>
  mutate(ratio_pred_over_measured = round(predicted_cells_mL/measured_cells_mL, 2),
         ground_truth = ifelse(is.na(measured_cells_mL), "none (prediction only)", "measured"))

write_csv(diag |> select(species, OD, fold_beyond_calib, predicted_cells_mL,
                         measured_cells_mL, ratio_pred_over_measured, ground_truth),
          file.path(DATA_DIR,"10_conversion_diagnostic.csv"))

cat("=== ENDPOINT CONVERSION DIAGNOSTIC (not a result) ===\n")
diag |> mutate(across(c(predicted_cells_mL,measured_cells_mL), ~ifelse(is.na(.),NA,format(round(.),big.mark=",")))) |>
  select(species, endpoint_OD=OD, fold_beyond_calib, predicted_cells_mL,
         measured_cells_mL, ratio_pred_over_measured, ground_truth) |>
  as.data.frame() |> print(row.names=FALSE)

gt <- diag |> filter(ground_truth=="measured")

rawcal <- read_csv(file.path(DATA_DIR, "calibration_metadata.csv"), show_col_types = FALSE) |>
  filter(sample_type == "calibration",
         tolower(as.character(use_for_fit)) == "true",
         species %in% c("SC", "RI")) |>
  transmute(species, x = as.numeric(target_OD), y = as.numeric(cells_per_mL)) |>
  filter(is.finite(x), is.finite(y))

pi_stats <- rawcal |> group_by(species) |> group_modify(~{
  m <- lm(y ~ x, data = .x); n <- nrow(.x)
  tibble(sigma = summary(m)$sigma, n = n, dof = n - 2,
         xbar = mean(.x$x), Sxx = sum((.x$x - mean(.x$x))^2))
}) |> ungroup()

gt <- gt |> left_join(pi_stats, by = "species") |>
  mutate(pi_hw = qt(0.975, dof) * sigma * sqrt(1 + 1/n + (OD - xbar)^2 / Sxx),
         pi_lo = predicted_cells_mL - pi_hw,
         pi_hi = predicted_cells_mL + pi_hw,
         pi_pct = round(100 * pi_hw / predicted_cells_mL, 1))
cat("\n95% prediction interval at endpoint OD (as % of prediction):\n")
print(as.data.frame(gt[, c("species","predicted_cells_mL","pi_lo","pi_hi","pi_pct")]), row.names = FALSE)
fD <- ggplot(gt, aes(measured_cells_mL, predicted_cells_mL, colour=species)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey45") +
  geom_errorbar(aes(ymin=pi_lo, ymax=pi_hi), width=0.02*max(gt$measured_cells_mL), linewidth=0.6) +
  geom_point(size=4) +
  ggrepel::geom_text_repel(aes(label=sprintf("%s\n%.2fx", species, ratio_pred_over_measured)),
                           size=3.2, fontface="bold", box.padding=0.6) +
  scale_colour_manual(values=c(SC="#1E2761", RI="#E8833A"), guide="none") +
  scale_x_continuous(labels=scales::label_number(scale=1e-6, suffix=" M"),
                     limits=c(0, max(gt$measured_cells_mL,gt$predicted_cells_mL)*1.15)) +
  scale_y_continuous(labels=scales::label_number(scale=1e-6, suffix=" M"),
                     limits=c(0, max(gt$measured_cells_mL, gt$pi_hi, gt$predicted_cells_mL)*1.15)) +
  labs(title="Prediction vs measured endpoint",
       x="Measured endpoint cells/mL (flow cytometry)",
       y="Predicted cells/mL (calibration, extrapolated)",
       caption="Dashed = 1:1. Error bars = 95% prediction interval of the extrapolated calibration (~4-6%); both points miss the 1:1 line by far more, so the gap is not calibration scatter but OD-density non-linearity at high density. Extrapolated 15-36x beyond OD 0.015; diagnostic only, not used for conversion.") +
  theme_bw(11) + theme(panel.grid.minor=element_blank(),
                       plot.title=element_text(face="bold", colour=ACCENT),
                       plot.caption=element_text(colour="#6B7280", size=8, hjust=0))
ggsave(file.path(DATA_DIR,"FigDiag_predicted_vs_measured.png"), fD, width=7.0, height=6.2, dpi=300)
ggsave(file.path(DATA_DIR,"FigDiag_predicted_vs_measured.pdf"), fD, width=7.0, height=6.2)
cat("\nwrote 10_conversion_diagnostic.csv + FigDiag_predicted_vs_measured\n")
