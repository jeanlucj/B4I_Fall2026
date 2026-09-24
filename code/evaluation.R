# ============================================================
# CONSOLE TOOLING FOR EVALUATING THIS PROJECT
#
# Driven from the console by EVALUATION.md (the analysis) and
# EVALUATION_SIMULATION.md (the simulation).  Nothing here is called by any
# pipeline script; sourcing it has no side effects beyond defining functions.
#
#   eval_load("analysis" | "simulation")  source the shared function files
#   eval_defs("code/<script>.R")          define ONE driver script's functions
#                                         and settings WITHOUT running it
#   eval_conflicts()                      which names two loaded files both define
#   eval_groups()                         the menu of debug() groups, fast -> slow
#   arm_evaluation(group) / disarm_evaluation()
#                                         debug() one module at a time
#   peek(x)                               one-line health summary of an object
#   check_grm(G) / check_matrix(Y) / check_alignment(...)
#                                         independent re-derivations, not the
#                                         pipeline's own opinion of itself
#
# The reason eval_defs() exists: most scripts in code/ define their functions
# and then immediately run a driver that talks to T3 or fits for twenty
# minutes.  Sourcing one to get at its functions runs the whole thing.
# eval_defs() evaluates only the top-level function definitions and the
# literal settings, and reports what it skipped.
#
# Sourced, not run.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------
# Module groups
#
# Ordered fast/offline -> slow/online within each pipeline.  One group is one
# `debug()` unit: arm exactly one, run the level's console lines, disarm.
# Group names are what EVALUATION.md's levels (A1...A12, S1...S11) refer to.
# ------------------------------------------------------------

EVAL_GROUPS <- list(

  # ---- analysis: offline, from data/ and an existing output/ ----------
  grm = c("read_grm", "collapse_grm", "grm_factor", "incidence"),

  curation_groups = c("identity_matrix", "identity_groups", "classify_groups",
                      "choose_representatives"),

  curation_family = c("pedigree_families", "family_pair_correlations",
                      "summarise_families", "pool_clonal_families",
                      "resolve_analysis_names"),

  matrix = c("scale_within_trial", "build_matrix", "trim_matrix",
             "grm_eigenvectors", "pea_blue", "pea_trait_blues"),

  megalmm = c("setup_megalmm_state", "run_megalmm", "megalmm_posterior",
              ".undecorate_rows"),

  cv = c("mask_cells", "margin_baselines", "cv_once"),

  bglr = c("covariance_components", "accession_effects", "rank_summary"),

  # ---- the validation trial (offline; reads the fit in output/) --------
  validation = c("validation_inputs", "validation_vintage", "build_pools",
                 ".pools_at", "pool_diff"),

  val_design = c("make_validation_design", "check_validation_design"),

  val_power  = c("contrast_se", "contrast_power", "power_grid",
                 "simulate_validation"),

  # ---- analysis: online (T3/Oat over BrAPI) ---------------------------
  t3 = c("connect_t3", "download_trial_observations", "download_observations",
         "fetch_observation_units", "pluck_chr"),

  trials = c("read_b4i_accessions", "candidate_trials", "count_b4i_accessions",
             "select_trials", "count_trait_coverage", "trial_metadata",
             "trait_availability", "clean_trait_name"),

  geno = c("read_accession_list", "resolve_accessions", "build_species_grm",
           "confirm_single_protocol", "marker_dosage", "fetch_pedigrees",
           "parent_names", "fetch_parent_records", "parent_genotype_status",
           "accession_vs_parents", "verify_parent_order"),

  # ---- simulation: all offline ----------------------------------------
  sim_design  = c("sim_grid", "sim_observed_parameters"),
  sim_panel   = c("sim_grms", "grm_factor", "draw_effect"),
  sim_cells   = c("sample_combinations"),
  sim_truth   = c("simulate_experiment"),
  sim_split   = c("standardize_within_env", "mask_observations",
                  "split_observations"),
  sim_basis   = c("grm_basis", "kron_basis"),
  sim_bglr    = c("fit_dge_ige", "run_scenario_bglr"),
  sim_megalmm = c("fit_megalmm", "run_scenario_megalmm"),
  sim_score   = c("interaction_part", "score_predictions",
                  "baseline_predictions")
)

# Fast/offline first, so a cheap bug surfaces before an expensive one.
EVAL_ORDER_ANALYSIS <- c("grm", "curation_groups", "curation_family", "matrix",
                         "megalmm", "cv", "bglr", "t3", "trials", "geno")

EVAL_ORDER_SIM <- c("sim_design", "sim_panel", "sim_cells", "sim_truth",
                    "sim_split", "sim_basis", "sim_bglr", "sim_megalmm",
                    "sim_score")

EVAL_ORDER_VALIDATION <- c("validation", "val_design", "val_power")

# The whole simulated-scenario chain: arm it and call run_scenario_bglr() to
# break once per stage.
EVAL_GROUPS$sim_pipeline <- unique(unlist(
  EVAL_GROUPS[c("sim_truth", "sim_split", "sim_basis", "sim_bglr", "sim_score")],
  use.names = FALSE))

# ------------------------------------------------------------
# Loading
#
# Which files are safe to source outright: the ones whose header says
# "Sourced, not run".  Everything else in code/ runs a driver on source.
# ------------------------------------------------------------

EVAL_SOURCES <- list(
  analysis = c("code/t3_functions.R",
               "code/curation_functions.R",
               "code/dge_ige_functions.R",
               "code/megalmm_setup.R"),
  simulation = c("code/dge_ige_functions.R",   # read_grm, grm_factor
                 "code/sim_config.R",
                 "code/sim_generate.R",
                 "code/sim_fit.R",
                 "code/megalmm_setup.R"),
  validation = c("code/validation_functions.R")
)

# Which driver scripts hold functions worth reaching with eval_defs().
EVAL_DEFS_FILES <- list(
  analysis = c("code/megalmm_build_inputs.R",   # the oat x pea matrix
               "code/megalmm_oat_pea.R",        # masking + baselines
               "code/find_trials_with_B4I_accessions.R",
               "code/create_GRMs_T3.R"),
  simulation = c("code/sim_run.R"),
  validation = c("code/validate_power.R",       # simulate_validation()
                 "code/validate_pool_selection.R",
                 "code/validate_design.R")
)

# name -> files that defined it, filled in by eval_load()/eval_defs()
.eval_defined <- new.env(parent = emptyenv())

.record_defs <- function(names, file) {
  for (nm in names) {
    prev <- get0(nm, envir = .eval_defined, ifnotfound = character(0))
    assign(nm, unique(c(prev, file)), envir = .eval_defined)
  }
}

#' Source the shared function files for one pipeline.
#'
#' @param which "analysis" or "simulation".
#' @param defs  Also pull in the driver scripts' own functions with
#'   eval_defs(), without running their drivers.
eval_load <- function(which = c("analysis", "simulation", "validation"),
                      defs = TRUE,
                      envir = globalenv()) {
  which <- match.arg(which)

  for (f in EVAL_SOURCES[[which]]) {
    before <- ls(envir, all.names = TRUE)
    sys.source(here::here(f), envir = envir, keep.source = TRUE)
    .record_defs(setdiff(ls(envir, all.names = TRUE), before), f)
    message("sourced ", f)
  }

  if (defs) for (f in EVAL_DEFS_FILES[[which]]) eval_defs(f, envir = envir)

  message("\n", which, " loaded. eval_groups(\"", which, "\") for the menu, ",
          "eval_conflicts() for names defined twice.")
  invisible(TRUE)
}

#' Evaluate a script's function definitions and literal settings, nothing else.
#'
#' A top-level assignment is evaluated when its right-hand side is a function
#' definition, a literal, or built only from the pure constructors in
#' `safe`.  Anything that would read a file, open a connection or fit a model
#' is skipped, so this is safe to run on a driver script.
eval_defs <- function(path, envir = globalenv(), quiet = FALSE) {
  full <- if (file.exists(path)) path else here::here(path)
  exprs <- parse(full, keep.source = TRUE)

  safe <- c("c", "list", "here::here", "here", "file.path", "paste", "paste0",
            "sprintf", "seq_len", "seq", "rep", "rev", "sort", "unique",
            "character", "numeric", "integer", "logical", "as.integer",
            "as.character", "as.numeric", "setdiff", "union", "intersect",
            "names", "unname", "modifyList", "rlang::set_names", "stats::setNames",
            "setNames", "[", "[[", "$", ":", "-", "+", "*", "/", "(")

  calls_in <- function(e) {
    if (!is.call(e)) return(character(0))
    c(deparse(e[[1]]), unlist(lapply(as.list(e)[-1], calls_in)))
  }

  n_fun <- n_set <- n_skip <- 0L
  fun_names <- set_names_ <- character(0)

  for (e in exprs) {
    ok <- FALSE
    if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=", "<<-") &&
        is.name(e[[2]])) {
      rhs <- e[[3]]
      is_fun <- is.call(rhs) && identical(as.character(rhs[[1]])[1], "function")
      if (is_fun) {
        ok <- TRUE
      } else {
        used <- setdiff(calls_in(rhs), "")
        ok <- length(used) == 0 || all(used %in% safe)
      }
      if (ok) {
        # A whitelisted call can still refer to an object the driver would
        # have built by now (`inputs[[1]]$Y`). Those are driver expressions
        # too; let them fail and count them as skipped.
        ok <- tryCatch({ eval(e, envir = envir); TRUE },
                       error = function(err) FALSE)
      }
      if (ok) {
        nm <- as.character(e[[2]])
        if (is_fun) { n_fun <- n_fun + 1L; fun_names <- c(fun_names, nm) }
        else        { n_set <- n_set + 1L; set_names_ <- c(set_names_, nm) }
        .record_defs(nm, path)
      }
    }
    if (!ok) n_skip <- n_skip + 1L
  }

  if (!quiet) {
    message(basename(path), ": ", n_fun, " function(s), ", n_set,
            " setting(s); skipped ", n_skip, " driver expression(s)")
    if (n_fun > 0) message("  ", paste(fun_names, collapse = ", "))
  }
  invisible(list(functions = fun_names, settings = set_names_, skipped = n_skip))
}

#' Names that more than one loaded file defines.
#'
#' Run it before trusting a walkthrough. `debug()` attaches to a *function
#' object*, so if two files define one name, arming the group debugs whichever
#' copy loaded last -- and if their bodies differ, stepping through one while
#' the pipeline calls the other is worse than not stepping at all.
#'
#' A clean repo answers "no name is defined by more than one loaded file".
#' Anything listed here is a refactor waiting to happen, not a setting.
eval_conflicts <- function(envir = globalenv()) {
  nms <- ls(.eval_defined)
  dup <- nms[vapply(nms, \(n) length(get(n, envir = .eval_defined)) > 1,
                    logical(1))]

  # A duplicated FUNCTION is a bug waiting to happen -- debug() would attach to
  # one copy while the pipeline calls the other.  A duplicated SETTING is just
  # two scripts each having their own settings block, which is the convention
  # here, so the two are reported apart.
  is_fun <- vapply(dup, \(n) exists(n, envir = envir, mode = "function"),
                   logical(1))

  show <- function(which, label) {
    if (!length(which)) return(invisible(NULL))
    cat(label, "\n", sep = "")
    for (n in sort(which))
      cat(sprintf("  %-22s %s\n", n,
                  paste(get(n, envir = .eval_defined), collapse = "  |  ")))
  }

  if (!any(is_fun)) {
    message("no FUNCTION is defined by more than one loaded file")
  } else {
    show(dup[is_fun], "Functions defined more than once -- fix these:")
  }
  show(dup[!is_fun], "\nSettings defined more than once (expected; one block per script):")

  invisible(list(functions = dup[is_fun], settings = dup[!is_fun]))
}

# ------------------------------------------------------------
# Arm / disarm
#
# debug(), not debugonce(): the flag persists until undebug()'d, and
# debugonce() can be neither cancelled nor seen by isdebugged().  The
# exists()/isdebugged() guards make both idempotent.
# ------------------------------------------------------------

arm_evaluation <- function(groups) {
  unknown <- setdiff(groups, names(EVAL_GROUPS))
  if (length(unknown)) {
    stop("unknown eval group(s): ", paste(unknown, collapse = ", "),
         "\n  known groups: ", paste(names(EVAL_GROUPS), collapse = ", "),
         call. = FALSE)
  }
  fns <- unique(unlist(EVAL_GROUPS[groups], use.names = FALSE))

  armed <- missing_fns <- character(0)
  for (fn in fns) {
    if (!exists(fn, mode = "function")) { missing_fns <- c(missing_fns, fn); next }
    if (!isdebugged(get(fn))) { debug(get(fn)); armed <- c(armed, fn) }
  }

  message("Armed debug() on ", length(armed), " function(s) in group(s): ",
          paste(groups, collapse = ", "), ".",
          "\n  In the debugger:  n = next line,  s = step into a call,",
          "  c = finish this function,  Q = quit.",
          "\n  disarm_evaluation() when done.")
  if (length(missing_fns)) {
    message("  not loaded (eval_load() / eval_defs() first): ",
            paste(missing_fns, collapse = ", "))
  }
  invisible(armed)
}

disarm_evaluation <- function(groups = names(EVAL_GROUPS)) {
  fns <- unique(unlist(EVAL_GROUPS[groups], use.names = FALSE))
  n <- 0L
  for (fn in fns)
    if (exists(fn, mode = "function") && isdebugged(get(fn))) {
      undebug(get(fn)); n <- n + 1L
    }
  message("Disarmed debug() on ", n, " function(s).")
  invisible(n)
}

#' The group menu, in fast -> slow order, with each group's functions and
#' whether they are currently loaded.
eval_groups <- function(which = c("all", "analysis", "simulation",
                                  "validation")) {
  which <- match.arg(which)
  known <- c(EVAL_ORDER_ANALYSIS, EVAL_ORDER_SIM, EVAL_ORDER_VALIDATION)
  order <- switch(which,
                  analysis   = EVAL_ORDER_ANALYSIS,
                  simulation = EVAL_ORDER_SIM,
                  validation = EVAL_ORDER_VALIDATION,
                  known)
  extra <- setdiff(names(EVAL_GROUPS), known)
  if (which == "all") order <- c(order, extra)

  for (g in order) {
    fns <- EVAL_GROUPS[[g]]
    have <- vapply(fns, \(f) exists(f, mode = "function"), logical(1))
    cat(sprintf("%-16s %s\n", g,
                paste0(fns, ifelse(have, "", "*"), collapse = ", ")))
  }
  cat("\n* = not loaded in this session\n")
  invisible(EVAL_GROUPS[order])
}

# ------------------------------------------------------------
# peek(): one line per object, tuned to the shapes that flow through this
# project, printing the failure signatures that otherwise pass silently.
# Returns its input invisibly, so it can sit in a pipe.
# ------------------------------------------------------------

peek <- function(x, label = deparse(substitute(x)), accessions = NULL) {
  lab <- paste0("[peek] ", label, ": ")

  # a build_grm() object, or a simulation
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(x$G) && is.matrix(x$G)) return(peek(x$G, paste0(label, "$G"), accessions))
    if (!is.null(x$obs) && !is.null(x$truth)) {
      s <- x$settings
      cat(lab, sprintf("simulation: %d x %d panel, %d obs (%.1f%% of cells), %d env, %d factor(s)\n",
                       s$n_oat, s$n_pea, nrow(x$obs),
                       100 * nrow(x$obs) / (s$n_oat * s$n_pea), s$n_envs,
                       s$n_factors), sep = "")
      cat("        V(true):", paste(sprintf("%s=%.3f", names(x$truth$V), x$truth$V),
                                    collapse = "  "), "\n")
      cat(sprintf("        obs per oat: %.1f  per pea: %.1f  (min %d / %d)\n",
                  nrow(x$obs) / s$n_oat, nrow(x$obs) / s$n_pea,
                  min(table(x$obs$oat)), min(table(x$obs$pea))))
      return(invisible(x))
    }
    if (!is.null(x$Y) && is.matrix(x$Y)) {
      cat(lab, "megalmm inputs, fill = ", x$fill %||% "?", "\n", sep = "")
      peek(x$Y, paste0(label, "$Y"))
      if (!is.null(x$X_Env)) peek(x$X_Env, paste0(label, "$X_Env"))
      return(invisible(x))
    }
  }

  # named numeric: effects, BLUEs, predictions
  if (is.numeric(x) && is.null(dim(x))) {
    fin <- x[is.finite(x)]
    cat(lab, sprintf("numeric[%d]  distinct=%d  NA=%d  range=[%s, %s]\n",
                     length(x), dplyr::n_distinct(x[!is.na(x)]), sum(is.na(x)),
                     if (length(fin)) formatC(min(fin), format = "g", digits = 4) else "NA",
                     if (length(fin)) formatC(max(fin), format = "g", digits = 4) else "NA"),
        sep = "")
    if (length(fin) && dplyr::n_distinct(fin) == 1)
      cat("        !! every finite value identical -- degenerate\n")
    if (!is.null(names(x)))
      cat("        names:", paste(utils::head(names(x), 4), collapse = ", "),
          if (length(x) > 4) "..." else "", "\n")
    return(invisible(x))
  }

  # matrix: a GRM, the oat x pea matrix, X_Env, a similarity matrix
  if (is.matrix(x)) {
    sq  <- nrow(x) == ncol(x)
    sym <- isTRUE(sq && isSymmetric(unname(x)))
    cat(lab, sprintf("matrix[%d x %d]  NA=%d (%.2f%% filled)",
                     nrow(x), ncol(x), sum(is.na(x)), 100 * mean(!is.na(x))),
        sep = "")
    if (sq) cat(sprintf("  square symmetric=%s  diag=[%.3g, %.3g]", sym,
                        min(diag(x), na.rm = TRUE), max(diag(x), na.rm = TRUE)))
    cat("\n")
    if (any(!is.na(x))) {
      rs <- rowSums(!is.na(x)); cs <- colSums(!is.na(x))
      if (!sq)
        cat(sprintf("        per row: min %d median %g | per col: min %d median %g\n",
                    min(rs), stats::median(rs), min(cs), stats::median(cs)))
      if (min(rs) == 0 || min(cs) == 0)
        cat("        !! an empty row or column -- MegaLMM's sampler goes to NaN on these\n")
    }
    if (!is.null(accessions) && !is.null(rownames(x))) {
      ov <- length(intersect(rownames(x), as.character(accessions)))
      cat(sprintf("        rowname overlap with accessions = %d / %d rows\n",
                  ov, nrow(x)))
      if (ov == 0 && nrow(x) > 0)
        cat("        !! zero overlap -- a name-mapping or collapse mismatch\n")
    }
    return(invisible(x))
  }

  # tibble: the plot table, observations, a grid
  if (is.data.frame(x)) {
    cat(lab, sprintf("tibble[%d x %d]\n", nrow(x), ncol(x)), sep = "")
    na <- vapply(x, \(col) sum(is.na(col)), integer(1))
    for (nm in names(x))
      cat(sprintf("        %-24s %-9s (NA=%d%s)\n", nm, class(x[[nm]])[1], na[[nm]],
                  if (na[[nm]] == nrow(x) && nrow(x) > 0) " -- ALL NA !!" else ""))
    return(invisible(x))
  }

  cat(lab, sprintf("%s  length=%d\n", paste(class(x), collapse = "/"), length(x)),
      sep = "")
  invisible(x)
}

# ------------------------------------------------------------
# Independent checks
#
# These re-derive a property rather than asking the pipeline whether it is
# happy.  Each prints and returns a tibble so it can be diffed across runs.
# ------------------------------------------------------------

#' Is this a usable relationship matrix, and how many rows carry no marker
#' information at all?
#'
#' build_grm() injects an ungenotyped accession as a prior-only line:
#' diagonal 1, every off-diagonal exactly 0.  Those lines are in the GRM and
#' in the models, related to nothing, and that is invisible unless counted.
check_grm <- function(G, label = deparse(substitute(G))) {
  if (is.list(G) && !is.null(G$G)) G <- G$G
  off <- rowSums(abs(G - diag(diag(G))))
  prior_only <- rownames(G)[off == 0]
  ev <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)

  out <- tibble::tibble(
    object = label, n = nrow(G),
    symmetric = isSymmetric(unname(G)),
    diag_min = min(diag(G)), diag_max = max(diag(G)),
    offdiag_min = min(G[upper.tri(G)]), offdiag_max = max(G[upper.tri(G)]),
    min_eigenvalue = ev, n_prior_only = length(prior_only)
  )
  print(as.data.frame(out), row.names = FALSE, digits = 4)
  if (length(prior_only))
    cat("  prior-only (no marker data):", length(prior_only), "->",
        paste(utils::head(prior_only, 6), collapse = ", "),
        if (length(prior_only) > 6) "..." else "", "\n")
  if (ev < -1e-6)
    cat("  !! negative eigenvalue: grm_factor() will drop directions here\n")
  invisible(dplyr::mutate(out, prior_only = list(prior_only)))
}

#' Shape and sparsity of an oat x pea matrix, and the trimming invariants.
check_matrix <- function(Y, min_per_col = 3, min_per_row = 2,
                         label = deparse(substitute(Y))) {
  rs <- rowSums(!is.na(Y)); cs <- colSums(!is.na(Y))
  out <- tibble::tibble(
    object = label, n_row = nrow(Y), n_col = ncol(Y),
    filled = sum(!is.na(Y)), pct_filled = 100 * mean(!is.na(Y)),
    min_per_row = min(rs), min_per_col = min(cs),
    rows_below = sum(rs < min_per_row), cols_below = sum(cs < min_per_col)
  )
  print(as.data.frame(out), row.names = FALSE, digits = 4)
  if (out$rows_below > 0 || out$cols_below > 0)
    cat("  !! trimming did not reach its own floor\n")
  invisible(out)
}

#' Do the phenotypes and the relationship matrices actually join?
#'
#' The failure this catches is silent: an accession renamed by curation on one
#' side and not the other leaves a level with no GRM row, and
#' `model.matrix()` drops the row rather than complaining.
check_alignment <- function(pheno, G_oat, G_pea,
                            oat_col = "germplasmName",
                            pea_col = "intercropGermplasmName") {
  if (is.list(G_oat) && !is.null(G_oat$G)) G_oat <- G_oat$G
  if (is.list(G_pea) && !is.null(G_pea$G)) G_pea <- G_pea$G

  oat <- unique(as.character(pheno[[oat_col]]))
  pea <- unique(as.character(pheno[[pea_col]]))
  miss_oat <- setdiff(oat, rownames(G_oat))
  miss_pea <- setdiff(pea, rownames(G_pea))

  out <- tibble::tibble(
    side = c("oat", "pea"),
    phenotyped = c(length(oat), length(pea)),
    in_grm = c(length(oat) - length(miss_oat), length(pea) - length(miss_pea)),
    missing = c(length(miss_oat), length(miss_pea)),
    grm_rows = c(nrow(G_oat), nrow(G_pea))
  )
  print(as.data.frame(out), row.names = FALSE)
  if (length(miss_oat))
    cat("  oat missing from G:", paste(utils::head(miss_oat, 8), collapse = ", "), "\n")
  if (length(miss_pea))
    cat("  pea missing from G:", paste(utils::head(miss_pea, 8), collapse = ", "), "\n")
  if (length(miss_oat) == 0 && length(miss_pea) == 0)
    cat("  every phenotyped accession has a GRM row\n")
  invisible(out)
}
