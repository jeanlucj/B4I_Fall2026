# ============================================================
# LOW-RANK DECOMPOSITION OF A FITTED OAT x PEA INTERACTION SURFACE
#
# Sourced, not run.
#
# WHAT THIS IS FOR. `dge_ige` predicts the interaction better than MegaLMM when
# combinations are sparse, and -- measured at full rank in
# output/simulation_kron_study.csv -- also when the interaction is higher-rank.
# But what it returns is a number per combination, not a score for an oat and a
# loading for a pea. MegaLMM returns the interpretable object and loses the
# accuracy. These functions recover the interpretable object from the better
# model.
#
# WHY NO DEREGRESSION STEP EXISTS. Fitting a second model to shrunken BLUPs
# would need deregression and inverse-variance weighting. Nothing is fitted
# here. The fitted surface is ALREADY an exact bilinear form -- see
# dge_ige_functions.R:392 -- so one SVD rewrites it in its own singular basis.
# That is an algebraic rotation of a term the model already estimated, not a
# secondary analysis: no information is created, none is double-counted, and
# there is no deregression to get wrong. Uncertainty comes from doing the same
# rotation inside each MCMC draw.
#
#   I_hat = A Beta B'                                  (A = U_a Lambda_a^(1/2))
#         = (Q_a U) D (Q_b V)'                         for svd(R_a Beta R_b')
#
# and because the two outer factors are orthonormal, sum(d_k^2) = ||I_hat||_F^2
# EXACTLY. So d_k^2 / sum(d_j^2) is an exact variance share, not an
# approximation. tests/test_decomp.R checks all of that as algebra.
#
# TWO THINGS THAT WILL MISLEAD IF IGNORED.
#
# 1. ROTATION. When two singular values are close, "component 2" is not a
#    well-defined object: the pair spans a plane and any rotation within it fits
#    equally well. No alignment scheme fixes that, so the rotation-INVARIANT
#    quantities -- cumulative share over the top k, and effective_rank() -- are
#    the ones to report. Per-component scores are secondary and must be read
#    with the `rot_diag` and `order_stable` diagnostics beside them.
#
# 2. CONVEXITY. Singular values are convex in the matrix, so
#    E[d_k(M)] >= d_k(E[M]): the average over draws of a spectrum is
#    systematically flatter and heavier-tailed than the spectrum of the average
#    surface. Both are reported and named separately. Reading the draw-wise
#    spread as though it were the point estimate would look exactly like the
#    prior-flattening bias the simulation is trying to measure.
#
# ON double_centre() VERSUS interaction_part(). They are the same operator.
# interaction_part() lives in sim_fit.R, which this file deliberately does not
# depend on -- the analysis pipeline sources this file without sim_fit.R. The
# duplication is guarded: tests/test_decomp.R loads both and asserts they agree.
# ============================================================

suppressPackageStartupMessages(library(tidyverse))

#' Sweep out both margins, leaving only the interaction.
#'
#' The same operator as `interaction_part()` in sim_fit.R; see the file header.
#' Scoring is done on the double-centred surface, so the decomposition is too --
#' otherwise the reported shares are shares of something other than what was
#' scored.
double_centre <- function(M) {
  M - rowMeans(M) - rep(colMeans(M), each = nrow(M)) + mean(M)
}

#' How many components the spectrum really has.
#'
#' The participation ratio `(sum d^2)^2 / sum d^4`: 1 for a rank-1 surface, r
#' for r equal singular values, and a continuous value in between. Preferred to
#' "count the ones above a threshold" because it needs no threshold and is
#' invariant to rotation within a degenerate block.
#'
#' `n90` is the blunter companion: how many components to reach 90% of the
#' variance.
effective_rank <- function(d) {
  d2 <- d^2
  tot <- sum(d2)
  if (!is.finite(tot) || tot <= 0) {
    return(list(participation = NA_real_, n90 = NA_integer_))
  }
  list(participation = tot^2 / sum(d2^2),
       n90 = which(cumsum(d2) / tot >= 0.9)[1])
}

#' Canonical correlations between two subspaces.
#'
#' Orthonormalises both and returns the singular values of the cross product,
#' which are the cosines of the principal angles: all 1 for identical spans, all
#' 0 for orthogonal ones. This is the right comparison when the truth is only
#' identified up to rotation -- which it is whenever `n_factors > 1`, so a
#' column-by-column correlation would understate recovery for reasons that have
#' nothing to do with the model.
subspace_cors <- function(X, Y) {
  if (is.null(X) || is.null(Y)) return(NA_real_)
  X <- as.matrix(X); Y <- as.matrix(Y)
  qx <- qr.Q(qr(X)); qy <- qr.Q(qr(Y))
  svd(crossprod(qx, qy), nu = 0, nv = 0)$d
}

#' The ceiling on recovery imposed by the truncated basis, with nothing fitted.
#'
#' At `kron_rank = q` the recovered scores are confined by construction to
#' `span(A)`, while the simulated truth has mass on every direction of G. So
#' there is a hard analytic limit on any subspace-recovery number, and it is
#' computable from the two spans alone. Recovery should be reported as a
#' FRACTION of this; without it a number at `kron_rank = 30` cannot be read, and
#' the apparent `n_acc` effect is partly 30-of-200 versus 30-of-400 directions
#' rather than anything about the model.
#'
#' Both spans are centred, because the interaction being decomposed is.
recovery_ceiling <- function(truth_U, A) {
  if (is.null(truth_U) || is.null(A)) return(NA_real_)
  ctr <- function(M) sweep(as.matrix(M), 2, colMeans(as.matrix(M)), "-")
  subspace_cors(ctr(truth_U), ctr(A))
}

#' Shared tail of the two decomposition paths: package an SVD into components.
.pack_components <- function(u, d, v, total_ss, rank = NULL, tol = 1e-10) {
  keep <- d > tol * max(c(d, .Machine$double.eps))
  if (!is.null(rank)) keep <- keep & (seq_along(d) <= rank)
  k <- max(sum(keep), 1L)
  d <- d[seq_len(k)]
  er <- effective_rank(d)
  share <- if (total_ss > 0) d^2 / total_ss else rep(NA_real_, k)
  list(d = d, share = share, cum_share = cumsum(share),
       scores = u[, seq_len(k), drop = FALSE],
       loadings = v[, seq_len(k), drop = FALSE],
       total_ss = total_ss, rank_kept = k,
       participation = er$participation, n90 = er$n90)
}

#' Decompose a fitted interaction surface. THE REFERENCE IMPLEMENTATION.
#'
#' Dense SVD of the double-centred surface. Slower than `decompose_bilinear()`
#' when the panel is large and the rank small, but it needs nothing but the
#' surface -- which is what makes the MegaLMM side free. MegaLMM's own factors
#' are NOT orthogonal (per-factor shares of its fitted interaction do not sum to
#' 1), so putting the two models on one scale means applying this same operator
#' to each model's surface rather than comparing a MegaLMM factor to a dge_ige
#' component.
#'
#' @return `d`, `share`, `cum_share`, `scores` (n_oat x k), `loadings`
#'   (n_pea x k), `total_ss`, `rank_kept`, `participation`, `n90`.
decompose_surface <- function(M, centre = TRUE, rank = NULL) {
  M <- as.matrix(M)
  S <- if (centre) double_centre(M) else M
  s <- svd(S)
  out <- .pack_components(s$u, s$d, s$v, sum(S^2), rank)
  rownames(out$scores)   <- rownames(M)
  rownames(out$loadings) <- colnames(M)
  out
}

#' Whiten the two bases once: the part of the bilinear route that is draw-free.
#'
#' `A` and `B` do not change between MCMC draws, so this is computed once and
#' reused -- which is not only cheaper but what makes the posterior summaries
#' comparable, since every draw is then expressed in the SAME outer coordinates
#' `Ua`, `Ub` and only the small `q x q` middle changes.
#'
#' Whitening is by SVD rather than QR on purpose: `qr()` may pivot, and a
#' pivoted `R` would silently break the `M = R_a Beta R_b'` identity. Centring
#' drops each basis's rank by one, which this handles -- the lost direction
#' arrives as a zero singular value.
.whiten <- function(A, B, centre = TRUE) {
  ctr <- function(M) if (centre) sweep(M, 2, colMeans(M), "-") else M
  sa <- svd(ctr(as.matrix(A))); sb <- svd(ctr(as.matrix(B)))
  list(Ua = sa$u, Ub = sb$u,
       Wa = sa$d * t(sa$v),                                  # D_a V_a'
       Wb = sb$v * rep(sb$d, each = nrow(sb$v)))             # V_b D_b
}

#' Decompose the same surface from its bilinear factors, without forming it.
#'
#' Centring is linear, so `double_centre(A Beta B') = (Ca A) Beta (Cb B)'` is
#' still bilinear and the decomposition reduces to an SVD of a `q x q` matrix.
#' Equivalent to `decompose_surface()` on the assembled surface --
#' tests/test_decomp.R asserts that -- and the route that keeps full rank
#' affordable on a large panel.
decompose_bilinear <- function(A, Beta, B, centre = TRUE, rank = NULL,
                               w = NULL) {
  A <- as.matrix(A); B <- as.matrix(B); Beta <- as.matrix(Beta)
  if (nrow(Beta) != ncol(A) || ncol(Beta) != ncol(B)) {
    stop("Beta is ", nrow(Beta), " x ", ncol(Beta), " but the bases are ",
         ncol(A), " and ", ncol(B), " columns wide: the reshape is wrong, ",
         "or a rank-deficient GRM returned fewer directions than asked for",
         call. = FALSE)
  }
  if (is.null(w)) w <- .whiten(A, B, centre = centre)
  M <- w$Wa %*% Beta %*% w$Wb
  s <- svd(M)
  out <- .pack_components(w$Ua %*% s$u, s$d, w$Ub %*% s$v, sum(s$d^2), rank)
  out$M <- M
  out$u_mid <- s$u
  out$v_mid <- s$v
  rownames(out$scores)   <- rownames(A)
  rownames(out$loadings) <- rownames(B)
  out
}

#' The one place the coefficient vector becomes a matrix.
#'
#' `byrow = TRUE` is not optional: column `(p-1)*kb + q` of the basis is
#' `A[,p] * B[,q]`, so the vector is row-major in the oat index. Getting it
#' wrong keeps the shape, raises no error, and still correlates well with the
#' truth -- see docs/specific-combination_kronecker.md, "Why it is fragile".
beta_from_vector <- function(b, ka, kb) {
  if (length(b) != ka * kb) {
    stop("coefficient vector is length ", length(b), ", expected ", ka * kb,
         call. = FALSE)
  }
  matrix(b, nrow = ka, ncol = kb, byrow = TRUE)
}

#' Read BGLR's streamed coefficient draws, dropping burn-in.
#'
#' THE BURN-IN TRAP. In `BGLR::Multitrait` the `saveEffects` write sits inside
#' `if (iter %% thin == 0)`, while the posterior mean accumulates under the
#' separate gate `(iter > burnIn) & (iter %% thin == 0)`. The header records
#' `nRow = nIter/thin`. So the file holds EVERY thinned iteration, burn-in
#' included, and the first `floor(burnIn/thin)` rows must go. Keeping them
#' inflates every posterior SD and flattens the spectrum -- which is
#' indistinguishable from the prior-flattening bias this machinery exists to
#' measure.
#'
#' The size assertion matters too: a killed job leaves a truncated file that
#' `readBinMatMultitrait()` returns as an NA-filled array with no error.
#'
#' @return `[n_kept, p, traits]`, burn-in removed.
read_beta_draws <- function(prefix, term = "ETA_G_mix", nIter, burnIn, thin,
                            p, traits, storage_mode = "double") {
  f <- paste0(prefix, term, "_beta.bin")
  if (!file.exists(f)) {
    stop("no streamed effects at ", f, ": was the fit run with ",
         "save_effects = TRUE?", call. = FALSE)
  }
  size     <- if (identical(storage_mode, "single")) 4L else 8L
  n_rows   <- nIter %/% thin
  expected <- 3 * size + n_rows * p * traits * size
  actual   <- file.size(f)
  if (actual != expected) {
    stop("effects file ", f, " is ", actual, " bytes, expected ", expected,
         " for ", n_rows, " x ", p, " x ", traits,
         ": a truncated file reads back as NA without erroring", call. = FALSE)
  }

  draws <- BGLR::readBinMatMultitrait(f, storageMode = storage_mode)
  if (anyNA(draws)) stop("streamed effects contain NA", call. = FALSE)

  drop_n <- burnIn %/% thin
  if (drop_n >= dim(draws)[1]) {
    stop("burn-in of ", burnIn, " at thin ", thin, " would drop all ",
         dim(draws)[1], " saved rows", call. = FALSE)
  }
  draws[seq.int(drop_n + 1L, dim(draws)[1]), , , drop = FALSE]
}

#' Match one draw's components to a reference, resolving sign and order.
#'
#' A draw's components carry no inherent labels: the joint sign of a
#' (score, loading) pair is unidentified, and the ORDER can swap between draws
#' whenever two singular values are close. Averaging unaligned scores across
#' draws therefore averages different things and shrinks them toward zero.
#'
#' Matching is greedy on `|t(scores) %*% ref_scores|`, which is exact enough at
#' these widths and needs no assignment-problem dependency. Within a block of
#' near-equal reference singular values no permutation is correct -- the truth
#' there is a continuous rotation -- so the block's rotation is measured and
#' reported rather than applied: `rot_diag` near 1 means the components are
#' individually stable, near 0 means only the block as a whole is meaningful.
#'
#' @param gap_tol Relative gap, as a fraction of `ref_d[1]`, below which two
#'   reference components count as degenerate.
#' @param k_stable How many leading components the `switched_top` diagnostic
#'   watches. Order stability over ALL components is uninformative: the trailing
#'   ones are noise and permute freely, so with enough of them at least one
#'   always swaps and the statistic reads 0 whatever the leading structure does.
align_to_reference <- function(scores, d, ref_scores, ref_d, gap_tol = 0.05,
                               k_stable = 3L) {
  k <- min(ncol(scores), ncol(ref_scores))
  C <- abs(crossprod(scores[, seq_len(k), drop = FALSE],
                     ref_scores[, seq_len(k), drop = FALSE]))

  perm <- integer(k); taken <- logical(k)
  for (j in seq_len(k)) {                     # reference component j
    cand <- which(!taken)
    i <- cand[which.max(C[cand, j])]
    perm[j] <- i; taken[i] <- TRUE
  }
  sgn <- vapply(seq_len(k), function(j)
    sign(crossprod(scores[, perm[j]], ref_scores[, j])[1, 1]), numeric(1))
  sgn[sgn == 0] <- 1

  # Degenerate blocks of the REFERENCE spectrum: consecutive components whose
  # relative gap is below tolerance belong to one rotationally-free block.
  gaps  <- if (k > 1) diff(ref_d[seq_len(k)]) / ref_d[1] else numeric(0)
  block <- cumsum(c(1L, as.integer(abs(gaps) >= gap_tol)))
  rot <- vapply(split(seq_len(k), block), function(ix) {
    if (length(ix) == 1L) return(1)
    s <- svd(crossprod(scores[, perm[ix], drop = FALSE],
                       ref_scores[, ix, drop = FALSE]))
    mean(abs(diag(s$u %*% t(s$v))))
  }, numeric(1))

  kt <- min(k_stable, k)
  list(perm = perm, sign = sgn,
       switched = !identical(perm, seq_len(k)),
       switched_top = !identical(perm[seq_len(kt)], seq_len(kt)),
       rot_diag = mean(rot[seq_len(min(length(rot), kt))]),
       gap_ratio = if (length(gaps)) min(abs(gaps)) else NA_real_,
       scores = sweep(scores[, perm, drop = FALSE], 2, sgn, "*"),
       d = d[perm])
}

#' Decompose every posterior draw and align them to a common reference.
#'
#' The reference is the decomposition of the POSTERIOR MEAN surface, not of the
#' first draw: draw 1 is one sample and would make the whole posterior summary
#' depend on it.
#'
#' @param draws `[n_draw, p, traits]` from `read_beta_draws()`.
#' @param trait Column index into the trait dimension.
#' @return A list with `reference` (the point estimate), `per_draw` (a tibble,
#'   one row per draw per component) and `diagnostics` (one row per draw).
#' @param k_proj How many leading reference subspaces to project each draw onto.
decompose_draws <- function(A, draws, B, trait, mean_beta = NULL,
                            centre = TRUE, rank = NULL, gap_tol = 0.05,
                            k_stable = 3L, k_proj = 5L) {
  ka <- ncol(A); kb <- ncol(B)
  n_draw <- dim(draws)[1]
  w <- .whiten(A, B, centre = centre)

  beta_bar <- if (is.null(mean_beta)) {
    beta_from_vector(colMeans(draws[, , trait, drop = FALSE][, , 1]), ka, kb)
  } else {
    beta_from_vector(mean_beta, ka, kb)
  }
  ref <- decompose_bilinear(A, beta_bar, B, centre = centre, rank = rank, w = w)
  kp <- min(k_proj, ncol(ref$u_mid), ncol(ref$v_mid))

  per <- vector("list", n_draw); diag_rows <- vector("list", n_draw)
  proj <- vector("list", n_draw)
  for (s in seq_len(n_draw)) {
    dc <- decompose_bilinear(A, beta_from_vector(draws[s, , trait], ka, kb), B,
                             centre = centre, rank = rank, w = w)
    al <- align_to_reference(dc$scores, dc$d, ref$scores, ref$d,
                             gap_tol = gap_tol, k_stable = k_stable)
    k <- length(al$d)
    tot <- sum(dc$d^2)
    per[[s]] <- tibble::tibble(
      draw = s, component = seq_len(k), d = al$d,
      share = if (tot > 0) al$d^2 / tot else NA_real_)

    # The draw's variance that lives in the REFERENCE's leading subspaces.
    #
    # This is the quantity whose point estimate is the reference's own
    # cumulative share, so an interval on it is an interval on the number being
    # reported. The draw-wise `share` above is NOT: each draw's spectrum is
    # flatter than the mean surface's (singular values are convex in the
    # matrix), so a draw's own component 1 carries systematically less than
    # component 1 of the mean surface. Quoting that interval beside the
    # reference value would put the point estimate outside its own interval.
    #
    # Both outer bases are shared across draws, so the projection reduces to
    # the middle: U_ref' M_draw V_ref.
    Pk <- vapply(seq_len(kp), function(j) {
      blk <- crossprod(ref$u_mid[, seq_len(j), drop = FALSE],
                       dc$M %*% ref$v_mid[, seq_len(j), drop = FALSE])
      if (tot > 0) sum(blk^2) / tot else NA_real_
    }, numeric(1))
    proj[[s]] <- tibble::tibble(draw = s, component = seq_len(kp),
                                cum_share_proj = Pk)

    diag_rows[[s]] <- tibble::tibble(
      draw = s, switched = al$switched, switched_top = al$switched_top,
      rot_diag = al$rot_diag,
      gap_ratio = al$gap_ratio, participation = dc$participation,
      total_ss = dc$total_ss)
  }

  list(reference = ref,
       per_draw = purrr::list_rbind(per),
       projected = purrr::list_rbind(proj),
       diagnostics = purrr::list_rbind(diag_rows))
}

#' Summarise the draws into one row per component, plus one row of diagnostics.
#'
#' THREE DIFFERENT NUMBERS, DELIBERATELY KEPT APART. Which one is the headline
#' matters, and the obvious choice is the wrong one.
#'
#'   share_draw_*     THE ESTIMAND. Each MCMC draw is a posterior sample of the
#'                    TRUE interaction surface, so the posterior of a draw's own
#'                    component k is the posterior of "component k's share of
#'                    the real interaction". This is the number to report, and
#'                    the only one of the three with an honest interval.
#'
#'   share_ref        component k's share of the POSTERIOR-MEAN surface. A
#'                    convenient summary but NOT an estimate of the above, and
#'                    biased towards concentration: singular values are convex
#'                    in the matrix, so averaging the draws cancels their
#'                    idiosyncratic directions and leaves a surface that looks
#'                    more low-rank than any draw. Measured on a dense rank-1
#'                    control: 0.691 against a posterior of 0.45 [0.37, 0.46].
#'                    Reporting it with the draw-wise interval beside it would
#'                    put the point estimate outside its own interval, which is
#'                    how this was caught.
#'
#'   cum_proj_*       how much of a draw's interaction lies in the reference's
#'                    leading-k subspace. Not a share of anything real -- an
#'                    IDENTIFIABILITY diagnostic: if a typical draw puts little
#'                    of its variance in the direction the mean surface picked
#'                    out, that direction is not a stable object to interpret.
component_summary <- function(dd, k_max = 10L) {
  ref <- dd$reference
  k <- min(k_max, ref$rank_kept)

  comp <- dd$per_draw |>
    dplyr::filter(component <= k) |>
    dplyr::group_by(component) |>
    dplyr::summarise(
      share_draw_mean = mean(share), share_draw_sd = stats::sd(share),
      share_draw_lo = stats::quantile(share, 0.025, names = FALSE),
      share_draw_hi = stats::quantile(share, 0.975, names = FALSE),
      d_mean = mean(d), .groups = "drop") |>
    dplyr::mutate(share_ref = ref$share[component],
                  cum_share_ref = ref$cum_share[component])

  pj <- dd$projected |>
    dplyr::filter(component <= k) |>
    dplyr::group_by(component) |>
    dplyr::summarise(
      cum_proj_mean = mean(cum_share_proj),
      cum_proj_lo = stats::quantile(cum_share_proj, 0.025, names = FALSE),
      cum_proj_hi = stats::quantile(cum_share_proj, 0.975, names = FALSE),
      .groups = "drop")
  comp <- dplyr::left_join(comp, pj, by = "component")

  dg <- dd$diagnostics
  list(components = comp,
       summary = tibble::tibble(
         n_draw = nrow(dg), rank_kept = ref$rank_kept,
         participation_ref = ref$participation, n90_ref = ref$n90,
         participation_draw_mean = mean(dg$participation, na.rm = TRUE),
         order_stable = mean(!dg$switched_top),
         order_stable_all = mean(!dg$switched),
         rot_diag = mean(dg$rot_diag, na.rm = TRUE),
         gap_ratio = mean(dg$gap_ratio, na.rm = TRUE),
         # the estimand, with its interval
         share1 = comp$share_draw_mean[1],
         share1_lo = comp$share_draw_lo[1],
         share1_hi = comp$share_draw_hi[1],
         # the mean-surface summary, which overstates concentration
         share1_meansurf = ref$share[1],
         cum_share3_meansurf = ref$cum_share[min(3L, ref$rank_kept)],
         # identifiability of the reported leading direction
         proj1_mean = comp$cum_proj_mean[1],
         proj1_lo = comp$cum_proj_lo[1], proj1_hi = comp$cum_proj_hi[1]))
}
