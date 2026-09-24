# ============================================================
# LAY OUT THE VALIDATION TRIAL
#
# Turns the pools from code/validate_pool_selection.R into a field book: which
# oat meets which pea, at which location, in which block.
#
# The structure is a 2x2 factorial of oat pool x pea pool, equally represented
# at every location.  Both species are validated from the same plots, because
# an associate effect is read on the PARTNER's yield -- the oat contrast comes
# off pea yield and the pea contrast off oat yield -- so nothing has to be
# split between them.  Complete balance at every location also means losing a
# whole site, as AL was lost to crop failure in 2025, costs a fifth of the
# plots and confounds nothing.
#
# Two properties the analysis depends on, both checked here rather than
# assumed:
#
#   * Every focal accession meets both partner pools equally often, so partner
#     effects are orthogonal to the focal contrast and drop out of it.
#   * A handful of anchor combinations are sown twice at each location. Within-
#     location duplication is the only thing that separates plot error from
#     combination x location; repeating a combination across locations alone
#     confounds the two.
#
# Outputs: output/validation/<vintage>/field_book.csv
#          output/validation/<vintage>/design_check.txt
#          output/validation/<vintage>/design_balance.png
# ============================================================

library(tidyverse)

here::i_am("code/validate_design.R")

source(here::here("code", "validation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_root <- here::here("output", "validation")
vintage  <- format(Sys.Date(), "%Y-%m-%d")

# Shared with the other validate_* scripts; see VALIDATION_DEFAULTS in
# code/validation_functions.R. n_anchor_per_cell = 1 gives 4 cells x 2 plots =
# 8 plots per location, 40 across the trial, which is what a first estimate of
# specific-combination variance costs.
n_locations        <- validation_setting("n_locations")
plots_per_location <- validation_setting("plots_per_location")
n_anchor_per_cell  <- validation_setting("n_anchor_per_cell")
block_size         <- validation_setting("block_size")

design_seed <- 20260924L

# ============================================================
# Driver
# ============================================================

out_dir <- file.path(out_root, vintage)
pools_file <- file.path(out_dir, "pools.csv")

if (!file.exists(pools_file)) {
  stop("No pools for vintage ", vintage, ".\n  Run code/validate_pool_selection.R first.",
       call. = FALSE)
}

pools <- readr::read_csv(pools_file, show_col_types = FALSE)
oat_pools <- dplyr::filter(pools, species == "oat")
pea_pools <- dplyr::filter(pools, species == "pea")

message("pools: oat ", nrow(oat_pools), " accessions, pea ", nrow(pea_pools))

design <- make_validation_design(
  oat_pools, pea_pools,
  n_loc = n_locations, p_loc = plots_per_location,
  n_anchor = n_anchor_per_cell, block_size = block_size,
  seed = design_seed
)

readr::write_csv(design, file.path(out_dir, "field_book.csv"))

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------

chk <- check_validation_design(design, oat_pools, pea_pools)

report <- utils::capture.output({
  cat("VALIDATION TRIAL DESIGN CHECK\n")
  cat("vintage ", vintage, "  |  ", n_locations, " locations x ",
      plots_per_location, " plots\n\n", sep = "")

  cat("plots per location\n")
  print(as.data.frame(chk$plots_per_location), row.names = FALSE)

  cat("\ncells per location (all four should be equal)\n")
  print(as.data.frame(tidyr::pivot_wider(
    chk$cells_per_location, names_from = pea_pool, values_from = plots,
    names_prefix = "pea ")), row.names = FALSE)

  cat("\nplots per accession\n")
  cat("  oat: median ", stats::median(chk$oat_reps$plots),
      ", range ", paste(range(chk$oat_reps$plots), collapse = "-"), "\n", sep = "")
  cat("  pea: median ", stats::median(chk$pea_reps$plots),
      ", range ", paste(range(chk$pea_reps$plots), collapse = "-"), "\n", sep = "")

  cat("\ndistinct partners per oat accession: median ",
      stats::median(chk$partners_per_oat$partners),
      ", range ", paste(range(chk$partners_per_oat$partners), collapse = "-"),
      "\n", sep = "")

  cat("\nPARTNER BALANCE -- the property the contrast rests on.\n")
  cat("  Each focal accession should meet the two partner pools equally often;\n")
  cat("  imbalance is |plots against As+ minus plots against As-|.\n")
  cat("  worst oat imbalance: ", chk$max_oat_imbalance, "\n", sep = "")
  cat("  worst pea imbalance: ", chk$max_pea_imbalance, "\n", sep = "")

  cat("\nreplicated combinations (for plot error and SCA)\n")
  cat("  anchor plots: ", chk$anchor_plots, "\n", sep = "")
  cat("  combinations appearing more than once: ",
      nrow(chk$replicated_combinations), "\n", sep = "")
  if (nrow(chk$replicated_combinations) > 0) {
    cat("  replication of those: ",
        paste(sort(unique(chk$replicated_combinations$n)), collapse = ", "),
        "\n", sep = "")
  }

  cat("\nVERDICT\n")
  ok <- TRUE
  say <- function(pass, msg) {
    cat("  [", if (pass) "ok  " else "FAIL", "] ", msg, "\n", sep = "")
    pass
  }
  ok <- say(dplyr::n_distinct(chk$plots_per_location$plots) == 1,
            "every location has the same number of plots") && ok
  ok <- say(dplyr::n_distinct(chk$cells_per_location$plots) == 1,
            "all four cells equally represented at every location") && ok
  ok <- say(chk$max_oat_imbalance == 0,
            "every oat meets both pea pools equally often") && ok
  ok <- say(chk$max_pea_imbalance == 0,
            "every pea meets both oat pools equally often") && ok
  ok <- say(chk$anchor_plots == 2 * 4 * n_anchor_per_cell * n_locations,
            "anchors duplicated within every location") && ok
  ok <- say(nrow(chk$replicated_combinations) > 0,
            "some combinations are replicated, so plot error is estimable") && ok
  cat("\n", if (ok) "design passes every check" else
      "DESIGN HAS A PROBLEM -- see the FAIL lines above", "\n", sep = "")
})

writeLines(report, file.path(out_dir, "design_check.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ------------------------------------------------------------
# Figure: the balance, seen
# ------------------------------------------------------------

bal <- dplyr::bind_rows(
  chk$oat_partner_balance |>
    tidyr::pivot_longer(c(`As+`, `As-`), names_to = "partner_pool",
                        values_to = "plots") |>
    dplyr::mutate(species = "oat accessions"),
  chk$pea_partner_balance |>
    tidyr::pivot_longer(c(`As+`, `As-`), names_to = "partner_pool",
                        values_to = "plots") |>
    dplyr::mutate(species = "pea accessions")
)

p <- ggplot2::ggplot(bal, ggplot2::aes(plots, fill = partner_pool)) +
  ggplot2::geom_bar(position = ggplot2::position_dodge(preserve = "single")) +
  ggplot2::facet_wrap(~ species, scales = "free") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Each accession meets both partner pools equally often",
    subtitle = paste0("vintage ", vintage,
                      "; the two bars must coincide at every count, or ",
                      "partner effects leak into the contrast"),
    x = "plots against that partner pool", y = "accessions",
    fill = "partner pool"
  )

ggplot2::ggsave(file.path(out_dir, "design_balance.png"), p,
                width = 9, height = 4.5, dpi = 150)

message("\nwrote ", out_dir)
