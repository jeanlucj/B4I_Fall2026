# Code

Save command-line scripts and shared R code here.

`figures/` holds scripts whose only job is to draw something for a talk or a
paper, as opposed to the pipeline scripts above them, which draw diagnostics as
a side effect of doing the analysis. They write their PNGs to
`output/figures/`, which is gitignored like the rest of `output/` — the script
is the thing worth versioning, the picture regenerates from it.

| script | what it draws |
|---|---|
| `figures/plot_pr_vs_as.R` | producer against associate effect, one panel per crop, in g/m2. The chain-averaged posterior mean, as opposed to `BGLR_{oat,pea}_Pr_vs_As.png`, which faces the four MCMC chains to show they agree |
| `figures/oat_pea_grid.R` | the oat x pea combination grid at 9% observed, for explaining how sparse the design is |
