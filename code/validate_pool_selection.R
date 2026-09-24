# ============================================================
# BUILD THE As+ / As- POOLS FOR THE VALIDATION TRIAL
#
# Two pools per species, both drawn from high-producer accessions and
# contrasted only on their associate effects.  Finding accessions with a LOW
# producer effect is not interesting -- poor lines are easy to find -- so the
# question the validation asks is whether, among good producers, we can tell
# apart the ones that help their partner from the ones that hurt it.
#
# The rule is a constrained optimisation: maximise mean(As+) - mean(As-)
# subject to the pools being disjoint, both pool means of Pr sitting above a
# threshold, and the two pool means of Pr being equal to within a tolerance.
# See build_pools() in code/validation_functions.R for why the obvious rule --
# top n on (Pr + As) against top n on (Pr - As) -- does not work.
#
# Pool membership is provisional until seed must be ordered.  Re-run this after
# every new trial is loaded; the diff against the previous vintage is written
# out, and churn in it is a result in its own right.
#
# Outputs: output/validation/<vintage>/pools.csv        both species, both pools
#          output/validation/<vintage>/pool_summary.csv what each pool achieved
#          output/validation/<vintage>/pool_diff.csv    vs the previous vintage
#          output/validation/<vintage>/pools_Pr_vs_As.png
# ============================================================

library(tidyverse)

here::i_am("code/validate_pool_selection.R")

source(here::here("code", "validation_functions.R"))

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

out_root <- here::here("output", "validation")

# Vintage folder. Dated, so a refresh never overwrites the pools a decision
# was taken on.
vintage <- format(Sys.Date(), "%Y-%m-%d")

# The settings that the three validate_* scripts must agree about live in
# VALIDATION_DEFAULTS, in code/validation_functions.R: power computed under a
# different partner filter would describe different pools than the ones
# selected here. Override a default only deliberately.
n_per_pool   <- validation_setting("n_per_pool")
pr_quantile  <- validation_setting("pr_quantile")
pr_tolerance <- validation_setting("pr_tolerance")
min_partners <- validation_setting("min_partners")

# ============================================================
# Driver
# ============================================================

out_dir <- file.path(out_root, vintage)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- validation_inputs(min_partners = min_partners)

cat("\n=== Vintage ===\n")
vint <- validation_vintage(inputs)
print(as.data.frame(vint), row.names = FALSE, digits = 4)

cat("\nReliability is what limits this experiment. PEV enters the contrast's\n",
    "standard error divided by the pool size, and nothing about the field\n",
    "design can reduce it -- only better estimates can.\n", sep = "")

# ------------------------------------------------------------
# Build
# ------------------------------------------------------------

built <- purrr::imap(inputs, \(inp, name) {
  n <- n_per_pool[[name]]
  message("\n--- ", name, ": ", n, " accessions per pool ---")
  bp <- build_pools(inp, n, pr_quantile = pr_quantile, tol = pr_tolerance)

  cat("\n=== ", toupper(name), " ===\n", sep = "")
  cat("eligible accessions (>= ", min_partners, " partners): ", inp$n_eligible,
      " of ", nrow(inp$accessions), "\n", sep = "")
  cat("candidates above the Pr threshold (", round(bp$pr_min, 2), "): ",
      bp$n_candidates, "\n", sep = "")
  cat("associate contrast  : ", round(bp$dAs, 2), " g/m2 on ", inp$response, "\n", sep = "")
  cat("producer difference : ", round(bp$dPr, 3), " g/m2 (tolerance ",
      pr_tolerance, ")\n", sep = "")
  cat("mean Pr, As+ / As-  : ", round(bp$mean_Pr_plus, 2), " / ",
      round(bp$mean_Pr_minus, 2),
      "   (population mean ", round(mean(inp$accessions$Pr), 2), ")\n", sep = "")
  cat("within-pool SD of true effects: ", round(sqrt(bp$sigma2_within), 2),
      " (PEV ", round(inp$PEV_As, 1), " + BLUP spread ",
      round(bp$var_blup_within, 1), ")\n", sep = "")

  stopifnot(
    "pools are not disjoint" =
      length(intersect(bp$pools$acc[bp$pools$pool == "As+"],
                       bp$pools$acc[bp$pools$pool == "As-"])) == 0,
    "a pool member falls below the partner filter" =
      all(bp$pools$n_partners >= min_partners)
  )
  bp
})

pools <- purrr::map(built, "pools") |> purrr::list_rbind()

pool_summary <- purrr::imap(built, \(bp, name) tibble::tibble(
  species = name, n_per_pool = bp$n, n_candidates = bp$n_candidates,
  pr_threshold = bp$pr_min, theta = bp$theta,
  dAs = bp$dAs, dPr = bp$dPr,
  mean_As_plus = bp$mean_As_plus, mean_As_minus = bp$mean_As_minus,
  mean_Pr_plus = bp$mean_Pr_plus, mean_Pr_minus = bp$mean_Pr_minus,
  sigma2_within = bp$sigma2_within, PEV = inputs[[name]]$PEV_As,
  reliability_As = inputs[[name]]$rel_As
)) |> purrr::list_rbind()

readr::write_csv(pools, file.path(out_dir, "pools.csv"))
readr::write_csv(pool_summary, file.path(out_dir, "pool_summary.csv"))
readr::write_csv(vint, file.path(out_dir, "vintage.csv"))

cat("\n=== Pool summary ===\n")
print(as.data.frame(pool_summary), row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# What changed since last time
#
# The pools are only worth validating if they are stable. If a couple of new
# trials reshuffle most of a pool, that says the effects are not settled, and
# it is much better to know that before seed is ordered than after.
# ------------------------------------------------------------

previous <- list.dirs(out_root, recursive = FALSE) |>
  setdiff(out_dir) |>
  (\(d) d[file.exists(file.path(d, "pools.csv"))])() |>
  sort()

if (length(previous) > 0) {
  prev_dir <- utils::tail(previous, 1)
  prev <- readr::read_csv(file.path(prev_dir, "pools.csv"), show_col_types = FALSE)
  dif <- pool_diff(prev, pools)

  readr::write_csv(dif$detail, file.path(out_dir, "pool_diff.csv"))

  cat("\n=== Change since ", basename(prev_dir), " ===\n", sep = "")
  print(as.data.frame(dif$summary), row.names = FALSE)

  churn <- dif$detail |>
    dplyr::filter(status != "retained") |>
    dplyr::count(species, name = "n_changed")
  retained <- dif$detail |>
    dplyr::filter(status == "retained") |>
    dplyr::count(species, name = "n_retained")

  dplyr::full_join(retained, churn, by = "species") |>
    dplyr::mutate(dplyr::across(c(n_retained, n_changed), \(x) tidyr::replace_na(x, 0L)),
                  pct_retained = round(100 * n_retained / (n_retained + n_changed))) |>
    as.data.frame() |>
    print(row.names = FALSE)

  cat("\nA pool that retains well under 70% of its members across a data\n",
      "vintage is a warning about the estimates, not about the design.\n", sep = "")
} else {
  cat("\nNo previous vintage to compare with; this is the baseline.\n")
}

# ------------------------------------------------------------
# Figure: where the pools sit in the Pr x As plane
# ------------------------------------------------------------

plot_data <- purrr::imap(inputs, \(inp, name) {
  inp$accessions |>
    dplyr::mutate(species = name) |>
    dplyr::left_join(dplyr::select(pools, species, acc, pool),
                     by = c("species", "acc")) |>
    dplyr::mutate(pool = tidyr::replace_na(pool, "not selected"))
}) |> purrr::list_rbind()

p <- plot_data |>
  dplyr::mutate(pool = factor(pool, levels = c("As+", "As-", "not selected"))) |>
  ggplot2::ggplot(ggplot2::aes(Pr, As, colour = pool, alpha = pool, size = pool)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
  ggplot2::geom_point() +
  ggplot2::scale_colour_manual(values = c("As+" = "#1b7837", "As-" = "#762a83",
                                          "not selected" = "grey60")) +
  ggplot2::scale_alpha_manual(values = c("As+" = 1, "As-" = 1, "not selected" = 0.35),
                              guide = "none") +
  ggplot2::scale_size_manual(values = c("As+" = 2.2, "As-" = 2.2, "not selected" = 1.2),
                             guide = "none") +
  ggplot2::facet_wrap(~ species, scales = "free") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    title    = "Validation pools: high producers, contrasted for associate effect",
    subtitle = paste0("vintage ", vintage, "; candidates restricted to >= ",
                      min_partners, " distinct partners and Pr above the ",
                      round(100 * pr_quantile), "th percentile"),
    x = "producer effect (own yield, g/m2)",
    y = "associate effect (partner yield, g/m2)", colour = NULL
  )

ggplot2::ggsave(file.path(out_dir, "pools_Pr_vs_As.png"), p,
                width = 10, height = 5, dpi = 150)

message("\nwrote ", out_dir)
