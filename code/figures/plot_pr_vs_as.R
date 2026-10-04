# ============================================================
# PRODUCER VS ASSOCIATE EFFECT -- presentation figures
#
#   Rscript code/figures/plot_pr_vs_as.R
#
# One figure per crop, for talks and papers. The diagnostic version of this
# plot -- output/BGLR_{oat,pea}_Pr_vs_As.png, written by
# code/BGLR_multi_trait_model.R -- faces four panels, one per MCMC chain, to
# show that the chains agree. They do: oat r runs -0.470 to -0.479 across the
# four seeds. That is a convergence check, not a result, and four copies of the
# same scatter is the wrong object to put on a slide.
#
# SO THIS PLOTS THE POSTERIOR MEAN, averaged over the four chains, which is the
# BLUP every downstream script already uses (validation_functions.R averages
# the same way). Collapsing the chains loses nothing a reader needs: the
# disagreement between them is in the third decimal of r.
#
# Units are g/m2 on both axes. An accession's PRODUCER effect is what it does
# to its own yield; its ASSOCIATE effect is what it does to its partner's. Both
# are deviations in grain yield per square metre, so the negative slope is the
# trade-off the project exists to measure: the oats that yield best for
# themselves tend to be the ones that cost their pea most.
#
# Reads  output/BGLR_{oat,pea}_effects_all_seeds.csv
# Writes output/figures/pr_vs_as_{oat,pea}.png
# ============================================================

library(tidyverse)

here::i_am("code/figures/plot_pr_vs_as.R")

out_dir <- here::here("output", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Presentation settings, gathered so a slide deck's house style is one edit.
point_alpha <- 0.55
point_size  <- 2.2
base_size   <- 15
fig_width   <- 7
fig_height  <- 5.5
fig_dpi     <- 300

species <- list(
  oat = list(file = "BGLR_oat_effects_all_seeds.csv", id = "oatAcc",
             partner = "pea"),
  pea = list(file = "BGLR_pea_effects_all_seeds.csv", id = "peaAcc",
             partner = "oat")
)

pr_as_figure <- function(sp, cfg) {
  f <- here::here("output", cfg$file)
  if (!file.exists(f)) {
    stop("no effects file at ", f,
         "\n  Run code/BGLR_multi_trait_model.R first.", call. = FALSE)
  }

  # Average over chains: this is the BLUP, and it is what the pools, the power
  # calculation and the field book are all built from.
  eff <- readr::read_csv(f, show_col_types = FALSE) |>
    dplyr::group_by(acc = .data[[cfg$id]]) |>
    dplyr::summarise(Pr = mean(PrEff), As = mean(AsEff),
                     n_chains = dplyr::n(), .groups = "drop")

  r <- stats::cor(eff$Pr, eff$As)
  message(sp, ": ", nrow(eff), " accessions over ",
          unique(eff$n_chains), " chains, r = ", round(r, 3))

  p <- ggplot2::ggplot(eff, ggplot2::aes(Pr, As)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    ggplot2::geom_point(alpha = point_alpha, size = point_size,
                        colour = "grey20") +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "firebrick", linewidth = 1.1) +
    # Placed against the panel corner rather than at fixed data coordinates, so
    # the two crops' figures stay comparable although their ranges differ.
    ggplot2::annotate("text", x = Inf, y = Inf,
                      label = sprintf("r = %.2f", r),
                      hjust = 1.2, vjust = 1.6, size = 5) +
    ggplot2::labs(
      x = expression("Producer effect, own yield ("*g/m^2*")"),
      y = bquote("Associate effect, "*.(cfg$partner)*" yield ("*g/m^2*")")
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )

  path <- file.path(out_dir, paste0("pr_vs_as_", sp, ".png"))
  ggplot2::ggsave(path, p, width = fig_width, height = fig_height,
                  dpi = fig_dpi)
  message("  wrote ", path)
  invisible(p)
}

invisible(purrr::imap(species, \(cfg, sp) pr_as_figure(sp, cfg)))
