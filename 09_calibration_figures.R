# Calibration figure (Fig 1)

suppressMessages({library(tidyverse); library(ggrepel)})
if (!exists("DATA_DIR")) DATA_DIR <- "."
NAVY <- "#1E2761"; ACCENT <- "#B85042"

SPP <- c(BH="#1D9E75", BT="#7F77DD", RI="#E8833A", SC="#1E2761", MG="#B85042", FD="#8A8D91")

theme_pub <- function(base=11) theme_bw(base) +
  theme(panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", colour=NAVY),
        plot.subtitle=element_text(colour="#6B7280", size=base-2),
        plot.caption=element_text(colour="#6B7280", size=base-3.5, hjust=0),
        legend.position="right")
save_fig <- function(p,n,w,h){ggsave(file.path(DATA_DIR,paste0(n,".png")),p,width=w,height=h,dpi=300)
  ggsave(file.path(DATA_DIR,paste0(n,".pdf")),p,width=w,height=h); cat("wrote",n,"\n")}

cal <- read_csv(file.path(DATA_DIR,"calibration_metadata.csv"), show_col_types=FALSE) |>
  filter(use_for_fit==TRUE)
fits <- read_csv(file.path(DATA_DIR,"08_calibration_fits.csv"), show_col_types=FALSE)
rat  <- read_csv(file.path(DATA_DIR,"08_crossspecies_ratios.csv"), show_col_types=FALSE)

pt <- cal |> group_by(species, media, target_OD) |>
  summarise(cells=mean(cells_per_mL), sd=sd(cells_per_mL), .groups="drop")

ord <- fits |> arrange(desc(slope)) |> pull(species)
pt   <- pt   |> mutate(species=factor(species, levels=ord))
fits <- fits |> mutate(species=factor(species, levels=ord),
                       lt=ifelse(media=="YCFAG","dashed","solid"))

lab <- rat |> mutate(species=factor(species, levels=ord),
                     txt=sprintf("%s  %.2fx%s", species, cells_per_OD_relative_to_SC,
                                 ifelse(media=="YCFAG","*","")))
ymax <- max(pt$cells)*1.02
lab  <- lab |> select(-any_of("slope")) |>
  left_join(fits |> select(species,slope,intercept), by="species") |>
  mutate(x=0.0155, y=intercept+slope*0.0155)

fMain <- ggplot(pt, aes(target_OD, cells, colour=species)) +
  geom_abline(data=fits, aes(slope=slope, intercept=intercept, colour=species, linetype=lt),
              linewidth=0.7, show.legend=FALSE) +
  geom_errorbar(aes(ymin=cells-sd, ymax=cells+sd), width=0.0003, alpha=0.7) +
  geom_point(size=2.3) +
  ggrepel::geom_text_repel(data=lab, aes(x=x, y=y, label=txt, colour=species),
                           hjust=0, direction="y", nudge_x=0.001, segment.size=0.3,
                           size=3.1, fontface="bold", xlim=c(0.016, NA)) +
  scale_colour_manual(values=SPP, guide="none") +
  scale_linetype_identity() +
  scale_x_continuous(limits=c(0, 0.024), breaks=seq(0,0.015,0.005),
                     labels=scales::label_number(accuracy=0.001)) +
  scale_y_continuous(labels=scales::label_number(scale=1e-6, suffix=" M")) +
  labs(x=expression("OD"[600]), y="Cell density (cells/mL)",
       caption="Labels: cells/OD relative to SC. *FD in YCFAG; others in WC. Spectrophotometer, OD 0.001-0.015. EC not calibrated.") +
  theme_pub()
save_fig(fMain, "FigCal_crossspecies", 8.8, 5.6)

r2lab <- fits |> mutate(txt=sprintf("R2 = %.3f\nslope = %.2e", R2, slope))
fSupp <- ggplot(pt, aes(target_OD, cells)) +
  geom_abline(data=fits, aes(slope=slope, intercept=intercept), colour=ACCENT,
              linewidth=0.6, linetype="dashed") +
  geom_errorbar(aes(ymin=cells-sd, ymax=cells+sd), width=0.0003, colour="grey55") +
  geom_point(size=2, colour=NAVY) +
  geom_text(data=r2lab, aes(x=0.001, y=Inf, label=txt), hjust=0, vjust=1.2,
            size=2.8, colour=NAVY, lineheight=0.95) +
  facet_wrap(~species, scales="free_y", nrow=2) +
  scale_x_continuous(breaks=seq(0,0.015,0.005)) +
  scale_y_continuous(labels=scales::label_number(scale=1e-6, suffix=" M")) +
  labs(x=expression("OD"[600]), y="Cell density (cells/mL)",
       caption="Points = replicate-well means +/- SD; dashed = linear fit. FD in YCFAG; others in WC.") +
  theme_pub()
save_fig(fSupp, "FigScal_curves", 9.5, 6.0)

cat("\ncalibration figures done\n")
