# Blinded morphology scoring and rater agreement (Fig 6)

suppressMessages({library(tidyverse); library(irrCAC); library(patchwork)})
if (!exists("DATA_DIR")) DATA_DIR <- "."
NAVY <- "#1E2761"; ACCENT <- "#B85042"
PAL  <- c(Static="#888888", Pulsed="#7F77DD", Continuous="#1D9E75")
LEV  <- c("Static","Pulsed","Continuous")

d <- read_csv(file.path(DATA_DIR,"microscopy_unblinded_scores.csv"), show_col_types=FALSE) |>
  mutate(condition=factor(condition, levels=LEV),
         species=factor(species, levels=c("BT","EC","FD","MG","RI","SC")))

qw <- function() {
  k <- 0:2; outer(k, k, function(a,b) 1 - ((a-b)/(max(k)-min(k)))^2)
}
ac2_one <- function(df, axis){
  w <- df |> select(field, rater, val=all_of(axis)) |>
    pivot_wider(names_from=rater, values_from=val) |>
    drop_na(Xingjian, Gaby)
  if (nrow(w) < 3) return(tibble(AC2=NA_real_, lo=NA_real_, hi=NA_real_))
  O <- table(factor(w$Xingjian, levels=0:2), factor(w$Gaby, levels=0:2))
  r <- gwet.ac1.table(as.matrix(O), weights = qw())
  ci <- r$coeff.ci
  nums <- as.numeric(unlist(regmatches(ci, gregexpr("[0-9.]+", ci))))
  tibble(AC2=r$coeff.val, lo=nums[1], hi=nums[2])
}
ac2 <- d |> group_by(species) |>
  summarise(agg=list(ac2_one(pick(everything()),"aggregation")),
            chn=list(ac2_one(pick(everything()),"chaining")),
            n_fields=n_distinct(field), .groups="drop") |>
  mutate(AC2_aggregation = map_dbl(agg,~.x$AC2), agg_lo=map_dbl(agg,~.x$lo), agg_hi=map_dbl(agg,~.x$hi),
         AC2_chaining    = map_dbl(chn,~.x$AC2), chn_lo=map_dbl(chn,~.x$lo), chn_hi=map_dbl(chn,~.x$hi)) |>
  select(-agg,-chn) |>
  mutate(across(starts_with("AC2")|ends_with("_lo")|ends_with("_hi"), ~round(.,2)))
write_csv(ac2, file.path(DATA_DIR,"11_AC2_agreement.csv"))

fld <- d |> group_by(species, field, slide, image, condition) |>
  summarise(aggregation=mean(aggregation, na.rm=TRUE),
            chaining   =mean(chaining,    na.rm=TRUE), .groups="drop")

cond_means <- fld |> group_by(species, condition) |>
  summarise(across(c(aggregation,chaining), list(mean=~round(mean(.),2), sd=~round(sd(.),2))),
            n_fields=n(), .groups="drop")
write_csv(cond_means, file.path(DATA_DIR,"11_condition_means.csv"))

rep_fields <- fld |> group_by(species, condition) |>
  mutate(cm_agg=mean(aggregation), cm_ch=mean(chaining),
         dist=sqrt((aggregation-cm_agg)^2 + (chaining-cm_ch)^2)) |>
  slice_min(dist, n=1, with_ties=FALSE) |>
  transmute(species, condition, field, slide, image,
            aggregation=round(aggregation,1), chaining=round(chaining,1)) |>
  ungroup()
write_csv(rep_fields, file.path(DATA_DIR,"11_representative_fields.csv"))

theme_p <- function(b=11) theme_bw(b) + theme(
  panel.grid.minor=element_blank(),
  strip.background=element_rect(fill="#F4F6F8",colour=NA),
  strip.text=element_text(face="bold",colour=NAVY),
  plot.title=element_text(face="bold",colour=NAVY),
  plot.subtitle=element_text(colour="#6B7280",size=b-2),
  legend.position="bottom")

pA <- ggplot(fld, aes(species, chaining)) +
  geom_boxplot(outlier.shape=NA, width=0.6, alpha=0.35, fill="#C9CDE0", colour="grey40") +
  geom_point(aes(colour=condition), position=position_jitter(width=0.18, height=0.05, seed=1), size=1.6, alpha=0.85) +
  scale_colour_manual(values=PAL, drop=FALSE) +
  scale_y_continuous(breaks=0:2, limits=c(-0.15,2.15)) +
  labs(title="A", subtitle="Chaining by species (0 single / 1 short / 2 long)",
       x=NULL, y="Chaining score", colour=NULL) +
  theme_p()

pB <- ggplot(fld, aes(condition, aggregation, fill=condition)) +
  geom_boxplot(outlier.shape=NA, width=0.65, alpha=0.55, colour="grey40") +
  geom_point(position=position_jitter(width=0.15, height=0.05, seed=1), size=1.3, alpha=0.7, colour="grey25") +
  facet_wrap(~species, nrow=1) +
  scale_fill_manual(values=PAL, drop=FALSE) +
  scale_y_continuous(breaks=0:2, limits=c(-0.15,2.15)) +
  labs(title="B", subtitle="Aggregation by regime (0 dispersed / 1 clustered / 2 aggregated)",
       x=NULL, y="Aggregation score", fill=NULL) +
  theme_p() + theme(axis.text.x=element_text(angle=45,hjust=1,size=7))

ac2_long <- bind_rows(
  ac2 |> transmute(species, axis="Aggregation", AC2=AC2_aggregation, lo=agg_lo, hi=agg_hi),
  ac2 |> transmute(species, axis="Chaining",    AC2=AC2_chaining,    lo=chn_lo, hi=chn_hi))
pC <- ggplot(ac2_long, aes(axis, species, fill=AC2)) +
  geom_tile(colour="white", linewidth=1) +
  geom_text(aes(label=sprintf("%.2f\n(%.2f-%.2f)",AC2,lo,hi), colour=AC2>0.80), lineheight=0.9,
            size=3.6, fontface="bold", show.legend=FALSE) +
  scale_fill_gradient2(low="#B85042", mid="#EAECEF", high="#1D9E75", midpoint=0.5,
                       limits=c(0,1), name="Gwet AC2") +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=NAVY)) +
  labs(title="C", subtitle="Inter-rater agreement, Gwet AC2 with 95% CI (2 blinded raters)", x=NULL, y=NULL) +
  theme_p() + theme(panel.grid=element_blank(), legend.position="right")

fig5 <- (pA | pC) / pB + plot_layout(heights=c(1, 1)) +
  plot_annotation(
    caption=NULL,
    theme=theme(plot.caption=NULL))

ggsave(file.path(DATA_DIR,"Fig6_morphology.png"), fig5, width=11, height=8.4, dpi=300)
ggsave(file.path(DATA_DIR,"Fig6_morphology.pdf"), fig5, width=11, height=8.4)

cat("=== AC2 agreement ===\n"); print(as.data.frame(ac2), row.names=FALSE)
cat("\n=== representative fields (for image pulls) ===\n"); print(as.data.frame(rep_fields), row.names=FALSE)
cat("\nwrote Fig6_morphology + 3 CSVs\n")
