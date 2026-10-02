# ============================================================
# THE OAT x PEA DESIGN GRID, AT 9% OF COMBINATIONS EVALUATED
#
#   Rscript docs/figures/oat_pea_grid.R
#   Rscript docs/figures/oat_pea_grid.R --n-acc 30 --sparsity 0.09 --seed 7
#   Rscript docs/figures/oat_pea_grid.R --random        # unconstrained draw
#
# Oats on the rows, peas on the columns, one cell per possible pairing. A filled
# cell was evaluated; an empty one is a pairing nobody grew. The point of the
# picture is how much of the grid stays empty: at 9% observed, 819 of 900 cells
# carry no data, and every one of them is a combination the model has to predict
# without ever having seen it.
#
# WHICH CELLS GET FILLED. By default this uses the simulation's own sampler,
# sample_combinations() in code/sim_generate.R, so the pattern is the one the
# simulation actually draws rather than a generic scatter. That sampler builds
# the design from rounds of random oat-to-pea matching, which guarantees every
# accession appears a minimum number of times -- a real design would not leave an
# entry unevaluated.
#
# It needs `min_per_acc` lowered to 2 here. The default is
# SIM_MIN_PER_ACC = 3, which at 30 x 30 would need 3 x 30 = 90 observations, and
# 9% of 900 is 81 -- so the sampler refuses, correctly. The constraint binds
# because the figure uses a deliberately small panel; the simulation's own cells
# at this density run 200 or 400 accessions, where 9% is thousands of plots.
#
# --random draws 81 cells uniformly instead, which leaves some rows and columns
# completely empty. Worth looking at once: it shows what the minimum-per-
# accession constraint is buying.
#
# Output: docs/figures/oat_pea_grid.png
# ============================================================

library(tidyverse)

here::i_am("docs/figures/oat_pea_grid.R")

source(here::here("code", "sim_config.R"))
source(here::here("code", "sim_generate.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

n_acc       <- as.integer(arg_value("--n-acc", "30"))
sparsity    <- as.numeric(arg_value("--sparsity", "0.09"))
min_per_acc <- as.integer(arg_value("--min-per-acc", "2"))
seed        <- as.integer(arg_value("--seed", "7"))
random      <- "--random" %in% args

out_dir <- here::here("docs", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# The cells
#
# sample_combinations() returns a data frame of (oat, pea) index pairs. The
# `--random` branch is a plain uniform draw over the grid, for contrast.
# ------------------------------------------------------------

set.seed(seed)
n_cells <- n_acc * n_acc
n_obs   <- round(sparsity * n_cells)

obs <- if (random) {
  cells <- sample(n_cells, n_obs)
  tibble::tibble(oat = ((cells - 1L) %% n_acc) + 1L,
                 pea = ((cells - 1L) %/% n_acc) + 1L)
} else {
  tibble::as_tibble(sample_combinations(n_acc, n_acc, sparsity,
                                        min_per_acc = min_per_acc))
}

grid <- tidyr::expand_grid(oat = seq_len(n_acc), pea = seq_len(n_acc)) |>
  dplyr::left_join(dplyr::mutate(obs, evaluated = TRUE),
                   by = c("oat", "pea")) |>
  dplyr::mutate(evaluated = dplyr::coalesce(evaluated, FALSE))

n_eval <- sum(grid$evaluated)
message(n_eval, " of ", n_cells, " cells evaluated (",
        round(100 * n_eval / n_cells, 1), "%)")
message("oats with no plot: ", sum(!seq_len(n_acc) %in% obs$oat),
        "   peas with no plot: ", sum(!seq_len(n_acc) %in% obs$pea))

# ------------------------------------------------------------
# The figure
#
# Every cell is drawn, which is what makes the empty ones legible as absence
# rather than as background: the grid lines are the tile borders, so there is no
# separate grid layer to align. The evaluated cells are then redrawn with a
# darker edge -- a light fill on a light surface needs that relief to hold its
# shape at this cell size.
# ------------------------------------------------------------

SURFACE   <- "#fcfcfb"   # chart surface
GRID      <- "#e3e2de"   # recessive cell borders
FILL_EVAL <- "#9ec5f4"   # light blue, sequential-blue step 200
EDGE_EVAL <- "#2a78d6"   # same hue, step 450, for relief
INK       <- "#0b0b0b"
INK_SOFT  <- "#52514e"

brks <- unique(c(1, seq(5, n_acc, by = 5)))

p <- ggplot2::ggplot(grid, ggplot2::aes(x = pea, y = oat)) +
  # all cells: fill carries the state, border is the grid
  ggplot2::geom_tile(ggplot2::aes(fill = evaluated),
                     colour = GRID, linewidth = 0.15) +
  # evaluated cells again, for the darker edge
  ggplot2::geom_tile(data = dplyr::filter(grid, evaluated),
                     fill = NA, colour = EDGE_EVAL, linewidth = 0.3) +
  ggplot2::scale_fill_manual(
    values = c(`TRUE` = FILL_EVAL, `FALSE` = SURFACE),
    labels = c(`TRUE` = "evaluated", `FALSE` = "never grown"),
    breaks = c("TRUE", "FALSE"), name = NULL) +
  # Without this the "never grown" key is surface-on-surface and invisible.
  ggplot2::guides(fill = ggplot2::guide_legend(
    override.aes = list(colour = GRID, linewidth = 0.4))) +
  # oat 1 at the top, so the picture reads like the matrix it is
  ggplot2::scale_y_reverse(breaks = brks, expand = ggplot2::expansion(0)) +
  ggplot2::scale_x_continuous(breaks = brks, position = "top",
                              expand = ggplot2::expansion(0)) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "The oat × pea design grid",
    # Two lines deliberately: at this figure width a single line overflows the
    # right edge, and ggplot does not wrap subtitles.
    subtitle = sprintf(
      "%d of %d combinations evaluated (%.0f%%).\nEvery empty cell is a pairing the model must predict unseen.",
      n_eval, n_cells, 100 * n_eval / n_cells),
    x = "pea accession", y = "oat accession") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.background   = ggplot2::element_rect(fill = SURFACE, colour = NA),
    panel.background  = ggplot2::element_rect(fill = SURFACE, colour = NA),
    panel.grid        = ggplot2::element_blank(),
    plot.title        = ggplot2::element_text(colour = INK, face = "bold",
                                              size = 13),
    plot.subtitle     = ggplot2::element_text(colour = INK_SOFT, size = 9,
                                              margin = ggplot2::margin(b = 10)),
    axis.title        = ggplot2::element_text(colour = INK_SOFT, size = 9),
    axis.text         = ggplot2::element_text(colour = INK_SOFT, size = 8),
    axis.ticks        = ggplot2::element_blank(),
    legend.position   = "bottom",
    legend.justification = "left",
    legend.key.size   = ggplot2::unit(9, "pt"),
    legend.text       = ggplot2::element_text(colour = INK_SOFT, size = 9),
    plot.margin       = ggplot2::margin(12, 14, 10, 12))

out <- file.path(out_dir, "oat_pea_grid.png")
ggplot2::ggsave(out, p, width = 5.6, height = 6.1, dpi = 300, bg = SURFACE)
message("wrote ", out)
