# PCBs in Northeast Atlantic common dolphins

Analysis code for:

> Fernández-Fernández, D. *et al.* Temporal trends and spatial patterns in PCB concentrations in common dolphins from the Northeast Atlantic: a Bayesian analysis. Submitted to *Scientific Reports*.

PCB concentrations in the blubber of 183 common dolphins (*Delphinus delphis*) stranded in France and Galicia (north-western Spain), 2001–2019, are modelled with Bayesian hierarchical models fitted by INLA, and projected forward to estimate when each group falls below two toxicity thresholds (17 and 41 mg/kg lipid weight).

## Structure

```
analysis.qmd        main analysis: models for males and females, projections
R/                  scripts for the immature model, tables and figures (run in order)
data/               input data
output/             everything the code produces (created on first run, not tracked)
```

## Requirements

R ≥ 4.5 with:

```r
install.packages(c("dplyr", "tidyr", "stringr", "ggplot2", "patchwork", "svglite",
                   "sf", "kableExtra", "knitr", "png", "psych", "vcd", "mgcv", "DHARMa"))
install.packages("INLA", repos = c(getOption("repos"),
                 INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

[Quarto](https://quarto.org) to render `analysis.qmd`, and [Inkscape](https://inkscape.org) for `R/10_figure_panels.R` (path set at the top of that script).

## How to run

Open `PCB DDE paper.Rproj` so the working directory is the project root, then:

1. Render `analysis.qmd`.
2. Run the scripts in `R/` in numerical order.

| Script | Produces |
|---|---|
| `analysis.qmd` | Table 2, base maps for Fig. 1, Table 3 (males, females), Fig. 4a, Supplementary Figs. S1–S2 |
| `R/01_model_selection_tables.R` | Supplementary Tables S2–S3 |
| `R/02_immature_model_selection.R` | Supplementary Table S4 |
| `R/03_year_effect_males.R` | Fig. 3a |
| `R/04_immature_year_effect.R` | Fig. 3b, Table 3 (immatures) |
| `R/05_immature_projection.R` | Immature projections and threshold years |
| `R/06_figure_projection_immatures.R` | Fig. 4b |
| `R/07_figure_probability_facets.R` | Supplementary Figs. S3–S4 |
| `R/08_figure_density.R` | Fig. 2 |
| `R/09_figure_threshold_timeline.R` | Fig. 5 |
| `R/10_figure_panels.R` | Assembles Figs. 3 and 4 from their panels |
| `R/11_table_temporal_distribution.R` | Table 1 |
| `R/12_check_immature_size.R`, `R/13_check_immature_length.R` | Body-size checks for the immature sex effect |

## Data

The raw data are not included in this repository. To run the code, place them in `data/PCBcorrectedDB.csv` (semicolon-separated, one row per individual) with these columns:

| Column | Description |
|---|---|
| `Country_Code` | Stranding country (4 = France, 5 = Spain) |
| `Lat`, `Lon` | Stranding location |
| `Year`, `Month` | Stranding date |
| `Sex`, `Length`, `Maturity`, `Pregnancy` | Biological variables; maturity inferred from body length |
| `CB149`, `CB187`, `ICES7`, `ICES9`, `sum18PCB` | PCB concentrations, mg/kg lipid weight |
| `ICES73` | Aroclor 1254-equivalent (ICES7 × 3): the response variable in all models |

`data/unique_models{m,f}.rds` hold the model-selection results used in the paper (WAIC, DIC and LCPO of every candidate model; no individual data). Refitting all 224 candidate models takes a while, so `analysis.qmd` reads these by default; set `run_selection <- TRUE` at the top to recompute them.

Threshold years come from Monte Carlo samples of the posterior (seed fixed at 1). Where the probability curve is flat around 0.8, as for France, another seed or INLA version can shift them by a year or more.

The map in `analysis.qmd` needs the GADM 3.6 country outlines (level 0) for Spain, France, the United Kingdom and Portugal. The GADM licence does not allow redistribution, so download them from [gadm.org](https://gadm.org) and place `gadm36_ESP_0.*`, `gadm36_FRA_0.*`, `gadm36_GBR_0.*` and `gadm36_PRT_0.*` in `data/shapefiles/`.

## Contact

Diego Fernández-Fernández — diego.fernandez@ieo.csic.es
