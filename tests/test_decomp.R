# test_decomp.R
#
# The low-rank decomposition of a fitted interaction surface
# (code/interaction_decomp.R).
#
# Every test here is an algebraic identity, because that is what the claim is:
# the fitted surface IS a scores-by-loadings object, exactly, and the SVD
# rewrites it without approximating. If any of these fail the claim is false,
# not merely imprecise.
#
# Three of the tests are NEGATIVE -- the wrong reshape, the wrong basis scaling
# and a rotation that a signed permutation cannot represent must all be
# detectably wrong. Without those, an oracle that only checks the right version
# is vacuous: code that returned the same thing for every input would pass.
#
# Panels are deliberately NON-square and the two bases of DIFFERENT rank, so no
# accidental transpose or index swap can conform by luck.
#
# Run: Rscript tests/test_decomp.R

library(tidyverse)
here::i_am("tests/test_decomp.R")
source(here::here("tests", "helper.R"))
load_code("dge_ige_functions.R", "interaction_decomp.R",
          "sim_config.R", "sim_generate.R", "sim_fit.R")

set.seed(20281001L)

N_OAT <- 9L; N_PEA <- 7L; KA <- 4L; KB <- 3L

A <- matrix(rnorm(N_OAT * KA), N_OAT, KA,
            dimnames = list(paste0("o", seq_len(N_OAT)), NULL))
B <- matrix(rnorm(N_PEA * KB), N_PEA, KB,
            dimnames = list(paste0("p", seq_len(N_PEA)), NULL))
Beta <- matrix(rnorm(KA * KB), KA, KB)
I <- A %*% Beta %*% t(B)

TOL <- 1e-10

# ------------------------------------------------------------
# 1. double_centre() is interaction_part(). The duplication is deliberate --
#    this file's functions must not depend on sim_fit.R -- so it is guarded.
# ------------------------------------------------------------

check_near(double_centre(I), interaction_part(I), TOL,
           "double_centre equals sim_fit.R's interaction_part")

# ------------------------------------------------------------
# 2. Reconstruction. The decomposition is exact, centred and uncentred.
# ------------------------------------------------------------

dc <- decompose_surface(I, centre = TRUE)
du <- decompose_surface(I, centre = FALSE)

rebuild <- function(x) x$scores %*% diag(x$d, nrow = length(x$d)) %*% t(x$loadings)

check_near(rebuild(dc), interaction_part(I), TOL,
           "centred decomposition rebuilds interaction_part(surface)")
check_near(rebuild(du), I, TOL,
           "uncentred decomposition rebuilds the surface itself")

# ------------------------------------------------------------
# 3. The factors are an orthonormal basis, and the centred scores are a PURE
#    interaction -- zero column means, so nothing additive leaked in.
# ------------------------------------------------------------

check_near(crossprod(dc$scores), diag(ncol(dc$scores)), TOL,
           "oat scores are orthonormal")
check_near(crossprod(dc$loadings), diag(ncol(dc$loadings)), TOL,
           "pea loadings are orthonormal")
check_near(colMeans(dc$scores), rep(0, ncol(dc$scores)), TOL,
           "centred scores have zero mean")
check_near(colMeans(dc$loadings), rep(0, ncol(dc$loadings)), TOL,
           "centred loadings have zero mean")

# ------------------------------------------------------------
# 4. The shares are EXACT variance shares. This is the whole point: it is what
#    licenses reading d_k^2 / sum(d^2) as "this component's fraction of the
#    interaction" rather than as an approximation to it.
# ------------------------------------------------------------

check_near(sum(dc$d^2), sum(interaction_part(I)^2), TOL,
           "sum of squared singular values equals the interaction sum of squares")
check_near(sum(dc$share), 1, TOL, "shares sum to one")
check_near(dc$cum_share[length(dc$cum_share)], 1, TOL,
           "cumulative share reaches one")

# ------------------------------------------------------------
# 5. The bilinear fast path equals the dense reference. Two independent routes
#    to the same object: SVD of a q x q matrix versus SVD of the n_oat x n_pea
#    surface.
# ------------------------------------------------------------

db <- decompose_bilinear(A, Beta, B, centre = TRUE)

check_near(db$d[seq_len(min(length(db$d), length(dc$d)))],
           dc$d[seq_len(min(length(db$d), length(dc$d)))], TOL,
           "bilinear path gives the same spectrum as the dense path")
check_near(rebuild(db), interaction_part(I), TOL,
           "bilinear path rebuilds the same surface")
check_near(abs(subspace_cors(db$scores[, 1, drop = FALSE],
                             dc$scores[, 1, drop = FALSE])), 1, 1e-8,
           "bilinear and dense leading score directions agree")

# ------------------------------------------------------------
# 6. Zero-padding invariance. fit_dge_ige() pads the surface back to the full
#    panel; the decomposition must see the same thing whether it is handed the
#    padded bases or the padded surface. The trap this guards is real: the
#    centring operator differs between the padded and unpadded panel, so
#    decomposing before padding would report shares of a different quantity
#    than the one that was scored.
# ------------------------------------------------------------

pad_rows <- function(M, extra) rbind(M, matrix(0, extra, ncol(M)))
A_pad <- pad_rows(A, 3L); B_pad <- pad_rows(B, 2L)
I_pad <- A_pad %*% Beta %*% t(B_pad)

check_near(decompose_bilinear(A_pad, Beta, B_pad)$d,
           decompose_surface(I_pad)$d, TOL,
           "padded bases and padded surface decompose identically")

# And padding really does change the answer, so the order is not a free choice:
# decomposing before padding would report shares of a different quantity than
# score_predictions() scores.
check(abs(decompose_surface(I_pad)$share[1] -
          decompose_surface(I)$share[1]) > 1e-4,
      "padding changes the centring, which is why it must happen first")

# ------------------------------------------------------------
# 7. Scale invariance. This is what licenses comparing a recovered score
#    against sim_generate.R's truth$U_oat, which is stored UNRESCALED while the
#    surface it came from was rescaled to a requested variance.
# ------------------------------------------------------------

sc <- decompose_surface(7.5 * I, centre = TRUE)
check_near(sc$d, 7.5 * dc$d, TOL, "singular values scale with the surface")
check_near(sc$share, dc$share, TOL, "shares are invariant to scale")
check_near(abs(subspace_cors(sc$scores[, 1, drop = FALSE],
                             dc$scores[, 1, drop = FALSE])), 1, 1e-8,
           "score directions are invariant to scale")

# ------------------------------------------------------------
# 8. Rank. A rank-r coefficient matrix gives exactly r components, and
#    effective_rank() reports r for r equal singular values.
# ------------------------------------------------------------

Beta1 <- outer(rnorm(KA), rnorm(KB))
d1 <- decompose_bilinear(A, Beta1, B)
check_near(d1$share[1], 1, 1e-8, "a rank-1 interaction is one component")
check_near(d1$participation, 1, 1e-6, "participation ratio is 1 at rank 1")
check(d1$n90 == 1L, "n90 is 1 at rank 1")

Beta2 <- outer(rnorm(KA), rnorm(KB)) + outer(rnorm(KA), rnorm(KB))
check(sum(decompose_bilinear(A, Beta2, B)$d > 1e-8) == 2L,
      "a rank-2 coefficient matrix gives exactly two non-null components")

check_near(effective_rank(rep(2, 5))$participation, 5, TOL,
           "participation ratio is r for r equal singular values")

# ------------------------------------------------------------
# 9. subspace_cors() and recovery_ceiling() behave at the extremes, and the
#    ceiling is BELOW one when the basis is truncated -- which is the whole
#    reason recovery has to be reported as a fraction of it.
# ------------------------------------------------------------

X <- qr.Q(qr(matrix(rnorm(N_OAT * 2), N_OAT, 2)))
check_near(subspace_cors(X, X), c(1, 1), 1e-8,
           "canonical correlations are 1 for identical spans")
Z <- qr.Q(qr(matrix(rnorm(N_OAT * N_OAT), N_OAT, N_OAT)))
check_near(max(subspace_cors(Z[, 1:2], Z[, 3:4])), 0, 1e-8,
           "canonical correlations are 0 for orthogonal spans")

truth_in <- A %*% matrix(rnorm(KA * 2), KA, 2)   # lies inside span(A)
check_near(min(recovery_ceiling(truth_in, A)), 1, 1e-8,
           "ceiling is 1 when the truth lies in the fitted span")
truth_out <- matrix(rnorm(N_OAT * 2), N_OAT, 2)  # generic, mostly outside
check(max(recovery_ceiling(truth_out, A)) < 0.999,
      "ceiling is below 1 when the truth leaves the truncated span")

# ------------------------------------------------------------
# 10. NEGATIVE: the reshape. byrow = FALSE keeps the shape, raises no error and
#     still produces a plausible surface -- it must be detectably wrong.
# ------------------------------------------------------------

b_vec  <- as.vector(t(Beta))                       # row-major, as BGLR returns
right  <- beta_from_vector(b_vec, KA, KB)
wrong  <- matrix(b_vec, nrow = KA, ncol = KB)      # byrow = FALSE
check_near(right, Beta, TOL, "beta_from_vector inverts the row-major flatten")
check(max(abs(decompose_bilinear(A, wrong, B)$share -
              decompose_bilinear(A, right, B)$share)) > 1e-3,
      "byrow = FALSE gives a detectably different spectrum")

# ------------------------------------------------------------
# 11. NEGATIVE: the guard on the reshape dimensions. grm_basis() caps its rank
#     at the number of non-null directions, so a rank-deficient GRM returns
#     fewer columns than asked for and the coefficient vector no longer
#     conforms. A guard that has stopped firing is worse than no guard.
# ------------------------------------------------------------

check_error(decompose_bilinear(A, matrix(0, KA + 1L, KB), B),
            "a Beta that does not match the bases is rejected")

# ------------------------------------------------------------
# 12. Alignment. A signed permutation must undo a sign flip and a swap exactly;
#     and it must NOT be able to undo a rotation inside a degenerate block --
#     rot_diag is what reports that, and a scheme that claimed success there
#     would be averaging different things across draws.
# ------------------------------------------------------------

Q <- qr.Q(qr(matrix(rnorm(N_OAT * 4), N_OAT, 4)))
ref_d_sep <- c(10, 6, 3, 1)                        # well separated

flipped <- Q; flipped[, 2] <- -flipped[, 2]
al <- align_to_reference(flipped[, c(2, 1, 3, 4)], ref_d_sep[c(2, 1, 3, 4)],
                         Q, ref_d_sep)
check(identical(al$perm, c(2L, 1L, 3L, 4L)), "alignment recovers the swap")
check_near(al$scores, Q, 1e-8, "alignment undoes the swap and the sign flip")
check_near(al$d, ref_d_sep, 1e-8, "alignment reorders the singular values")
check_near(al$rot_diag, 1, 1e-8, "rot_diag is 1 when components are separated")

# A rotation within the near-degenerate leading pair. No permutation represents
# it, so rot_diag must fall away from 1.
ref_d_deg <- c(10, 9.9, 3, 1)
th <- pi / 5
R  <- diag(4); R[1:2, 1:2] <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
al_deg <- align_to_reference(Q %*% R, ref_d_deg, Q, ref_d_deg)
check(al_deg$rot_diag < 0.95,
      "rot_diag detects a rotation inside a degenerate block")
check(al_deg$gap_ratio < 0.05,
      "gap_ratio flags the degenerate block")

# ------------------------------------------------------------
# 13. MegaLMM's factors are NOT orthogonal, so its raw per-factor shares do not
#     sum to one. This is the test that justifies applying the same operator to
#     both models rather than comparing a MegaLMM factor to a dge_ige
#     component.
# ------------------------------------------------------------

U_F <- matrix(rnorm(N_OAT * 3), N_OAT, 3)
Lam <- matrix(rnorm(3 * N_PEA), 3, N_PEA)
surf <- U_F %*% Lam
raw_shares <- vapply(1:3, function(k)
  sum(interaction_part(outer(U_F[, k], Lam[k, ]))^2) /
    sum(interaction_part(surf)^2), numeric(1))
check(abs(sum(raw_shares) - 1) > 1e-3,
      "MegaLMM's raw per-factor shares do not sum to one")
check_near(sum(decompose_surface(surf)$share), 1, TOL,
           "the same operator on MegaLMM's surface gives shares that do")

# ------------------------------------------------------------
# 14. The projected share is an interval on the RIGHT estimand. Fed a set of
#     "draws" that are all exactly the posterior mean, the projection onto the
#     reference's leading-k subspace must equal the reference's own cumulative
#     share. Without this the posterior interval reported beside a point
#     estimate is an interval on a different quantity -- measured on a dense
#     rank-1 control, the draw-wise interval 0.41-0.49 excludes its own point
#     estimate of 0.691.
# ------------------------------------------------------------

b_flat <- as.vector(t(Beta))
fake <- array(rep(b_flat, each = 4L), dim = c(4L, length(b_flat), 2L))
fake[, , 1] <- matrix(b_flat, nrow = 4L, ncol = length(b_flat), byrow = TRUE)
fake[, , 2] <- matrix(b_flat, nrow = 4L, ncol = length(b_flat), byrow = TRUE)

dd <- decompose_draws(A, fake, B, trait = 2L, mean_beta = b_flat)
cs <- component_summary(dd)
check_near(cs$components$cum_proj_mean, cs$components$cum_share_ref, 1e-8,
           "projected share equals the reference cumulative share when draws are the mean")
check_near(c(cs$summary$share1, cs$summary$share1_lo, cs$summary$share1_hi),
           rep(cs$summary$share1_meansurf, 3), 1e-8,
           "with no spread, the posterior share and the mean-surface share coincide")
check(all(diff(cs$components$cum_proj_mean) >= -1e-10),
      "projected cumulative shares are non-decreasing in k")
check_near(cs$summary$stable1, 1, 1e-12,
           "with no spread, the leading component is perfectly stable")
check_near(cs$components$comp_stable, rep(1, nrow(cs$components)), 1e-12,
           "and so is every component")

# ------------------------------------------------------------
# 15. The recovery FLOOR. The ceiling says how well the truncated basis could
#     do; this says how well it does by accident. It is not small, because the
#     truth is kinship-structured and so is every direction in the span -- so a
#     recovery number quoted without it, or a trait correlation quoted without
#     it, credits population structure as signal.
# ------------------------------------------------------------

truth_k <- A %*% matrix(rnorm(KA * 2), KA, 2)      # kinship-structured truth
nl <- recovery_null(truth_k, A, n_draw = 200L, seed = 42L)
check(nl[["median"]] > 0.05,
      sprintf("a random direction in the span already scores (%.2f)", nl[["median"]]))
check(nl[["p95"]] > nl[["median"]], "the null's 95th percentile exceeds its median")
check(nl[["p95"]] <= 1, "the null stays a correlation")
check(all(is.na(recovery_null(NULL, A))),
      "no truth means no floor, rather than a number")

finish("test_decomp")
