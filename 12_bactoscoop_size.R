# Single-cell size extraction (Fig S6)

suppressMessages({library(tidyverse)})
if (!exists("DATA_DIR")) DATA_DIR <- "."
NAVY <- "#1E2761"; ACCENT <- "#B85042"
PAL  <- c(Static="#888888", Pulsed="#7F77DD", Continuous="#1D9E75")
LEV  <- c("Static","Pulsed","Continuous")

sl <- read_csv(file.path(DATA_DIR,"microscopy_size_results_by_slide.csv"), show_col_types=FALSE) |>
  mutate(regime=factor(regime, levels=LEV),
         species=factor(species, levels=c("EC","BT")))

long <- sl |>
  select(species, slide, regime, n_cells,
         Width=median_width_um, Length=median_length_um, `Aspect ratio`=median_AR) |>
  pivot_longer(c(Width, Length, `Aspect ratio`), names_to="metric", values_to="value") |>
  mutate(metric=factor(metric,
           levels=c("Width","Length","Aspect ratio"),
           labels=c("Width (um)","Length (um)","Aspect ratio")))

theme_p <- function(b=11) theme_bw(b) + theme(
  panel.grid.minor=element_blank(),
  strip.background=element_rect(fill="#F4F6F8",colour=NA),
  strip.text=element_text(face="bold",colour=NAVY),
  plot.title=element_text(face="bold",colour=NAVY),
  plot.subtitle=element_text(colour="#6B7280",size=b-2),
  plot.caption=NULL,
  legend.position="bottom")

p <- ggplot(long, aes(regime, value, fill=regime)) +
  geom_point(aes(size=n_cells), shape=21, colour="grey20", stroke=0.5,
             position=position_jitter(width=0.12, height=0, seed=1), alpha=0.9) +
  facet_grid(metric ~ species, scales="free_y", switch="y") +
  scale_fill_manual(values=PAL, drop=FALSE, name=NULL) +
  scale_size_continuous(range=c(2.5,7), name="cells / slide") +
  labs(title="Single-cell size by shaking regime (EC, BT)",
       subtitle="Each point = one slide median | point size = cells measured",
       x=NULL, y=NULL,
       caption=NULL) +
  theme_p() + theme(strip.placement="outside")

ggsave(file.path(DATA_DIR,"FigS6_cell_size.png"), p, width=8.2, height=7.2, dpi=300)
ggsave(file.path(DATA_DIR,"FigS6_cell_size.pdf"), p, width=8.2, height=7.2)

cat("=== per-slide size summary ===\n")
print(as.data.frame(sl |> select(species,slide,regime,n_cells,median_width_um,median_length_um,median_AR)),
      row.names=FALSE)
cat("\nBT width by regime (slide medians):\n")
sl |> filter(species=="BT") |> group_by(regime) |>
  summarise(widths=paste(sort(median_width_um),collapse=", "), .groups="drop") |>
  as.data.frame() |> print(row.names=FALSE)
cat("\nwrote FigS6_cell_size\n")
