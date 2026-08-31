# Shaking regime and inoculum OD in gut-bacterial growth experiments

Analysis code for the study of how the **hydrodynamic (shaking) regime** of a microplate
reader and the **optical-density-to-cell-number relationship** bias growth measurements of
human gut bacteria — two routine, usually-ignored settings that turn out to act
differently on different species.

Seven gut species — *Bacteroides thetaiotaomicron* (BT), *Escherichia coli* LF82 (EC),
*Faecalibacterium duncaniae* (FD), *Mediterraneibacter gnavus* (MG),
*Roseburia intestinalis* (RI), *Segatella copri* (SC) and *Blautia hydrogenotrophica* (BH) —
were grown under three regimes (static, pulsed, continuous shaking) in strictly anaerobic
96-well plates, with paired flow-cytometry, OD-to-cell calibration, and microscopy
sub-studies.

This repository contains the analysis code needed to reproduce the analyses and figures
from the raw data.

---

## Paper & data

- **Paper:** Adedamola Godwin Daodu, Pallabita Saha, Gabriela Bravo-Ruiseco, Xingjian Zhou,
  Karoline Faust. *(Journal and DOI to be added on acceptance.)*
- **Data:** the raw plate-reader, calibration, flow-cytometry, and microscopy data will be
  deposited on Zenodo and made openly available on publication *(DOI to be added)*.

> The raw data are **not** stored in this repository. On publication they will be available
> from the Zenodo record above; until then, please refer to the data availability statement
> in the paper.

---

## What's here

15 R scripts, one stage per file, run in numeric order. Scripts `01`–`05` do the analysis
and write intermediate result tables; `06`–`12` produce the figures.

| Script | Purpose |
|---|---|
| `00_load_and_filter.R` | Load raw data, apply quality-control filters |
| `00_run_all.R` | Run scripts 01–05 in order, with checks and logging |
| `01_design_and_identifiability.R` | What the design can and cannot estimate (regime vs batch) |
| `02_endpoint_OD.R` | Within-run endpoint-OD contrasts per regime |
| `03_PvsC_paired.R` | Pulsed-vs-continuous analysis and flow cytometry |
| `04_batch_effects.R` | Between-plate (run) variance components |
| `05_CoV_reproducibility.R` | Within-plate coefficient-of-variation reproducibility |
| `05_FD_growthcurve_startingOD.R` | *F. duncaniae* growth curves at two starting ODs |
| `06_figures.R` | Main figures |
| `07_supplementary_kinetics.R` | Supplementary kinetic figures |
| `08_calibration_fits.R` | Per-species OD-to-cell calibration fits |
| `09_calibration_figures.R` | Calibration figure |
| `10_endpoint_conversion_diagnostic.R` | Endpoint calibration-extrapolation diagnostic |
| `11_microscopy_scoring.R` | Blinded morphology scoring and rater agreement |
| `12_bactoscoop_size.R` | Single-cell size extraction |

A full data dictionary (every column of every raw file) will accompany the data on the
Zenodo record.

---

## Requirements

- **R** ≥ 4.5.0
- R packages: `tidyverse`, `lme4`, `lmerTest`, `emmeans`, `irrCAC`, `patchwork`,
  `ggrepel`, `broom`, `readxl`

Install the packages:

```r
install.packages(c("tidyverse","lme4","lmerTest","emmeans","irrCAC",
                   "patchwork","ggrepel","broom","readxl"))
```

Versions used for the published analysis: R 4.5.0; tidyverse 2.0.0, dplyr 1.1.4,
ggplot2 4.0.0, readr 2.1.5, readxl 1.5.0, lme4 2.0.1, lmerTest 3.2.1, emmeans 2.0.4,
irrCAC 1.4, patchwork 1.3.2, ggrepel 0.9.8, broom 1.0.10.

---

## Reproduce it

1. **Get the code** — clone or download this repository:
   ```bash
   git clone https://github.com/msysbio/shaking-regime-gut-bacteria.git
   ```
2. **Get the data** — the raw data files will be available from the Zenodo record on
   publication (DOI above). Download them and place them in the **same folder** as the
   scripts: the scripts and data must sit together in one flat directory.
3. **Set the folder path** at the top of a script (or your session):
   ```r
   DATA_DIR <- "/path/to/this/folder"
   ```
4. **Run the core analysis** (scripts 01–05, in order):
   ```r
   source(file.path(DATA_DIR, "00_run_all.R"))
   ```
5. **Generate the figures** (run as needed):
   ```r
   source(file.path(DATA_DIR, "06_figures.R"))                # main figures
   source(file.path(DATA_DIR, "07_supplementary_kinetics.R")) # supplementary kinetics
   source(file.path(DATA_DIR, "08_calibration_fits.R"))       # calibration fits
   source(file.path(DATA_DIR, "09_calibration_figures.R"))    # calibration figure
   source(file.path(DATA_DIR, "10_endpoint_conversion_diagnostic.R"))
   source(file.path(DATA_DIR, "11_microscopy_scoring.R"))
   source(file.path(DATA_DIR, "12_bactoscoop_size.R"))
   ```
   (Run `08` before `09` and `10`, which use its output.)

All scripts read and write to `DATA_DIR`; figures are written there as PNG and PDF.

---

## Notes

- **Flat folder.** Scripts assume code and data share one directory; moving files into
  subfolders requires editing the `DATA_DIR` paths.
- **Statistics.** Endpoint OD is modelled per species as `endpoint ~ regime + (1|run)`
  (REML; lme4/lmerTest), with pairwise contrasts via emmeans (Bonferroni-corrected within
  species) and rater agreement as Gwet's AC2 (irrCAC). Pulsed and continuous shaking are
  never compared directly — they never share a run, so each regime is compared only against
  its paired same-run static control.
- **Reproducibility.** `00_run_all.R` runs each script in a clean environment; scripts
  communicate only through files on disk.

---

## Citation

If you use this code or data, please cite the paper:

> Adedamola Godwin Daodu, Pallabita Saha, Gabriela Bravo-Ruiseco, Xingjian Zhou,
> Karoline Faust. *(Title, journal, and DOI to be added on acceptance.)*

and the Zenodo archive *(DOI to be added on publication)*.

## License

- Code: [MIT]
- Data (on Zenodo): [CC-BY-4.0 — or the lab's preferred licence]

## Contact

Adedamola Godwin Daodu — Laboratory of Molecular Bacteriology, KU Leuven.
Adedamolagodwin.daodu@kuleuven.be
