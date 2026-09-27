# test_kronecker.R
#
# The specific-combination basis and the reshape that turns fitted coefficients
# back into a full interaction surface. Derived in
# docs/specific-combination_kronecker.md.
#
# This is the most valuable file in the suite, because the failure mode is
# silent: get the index arithmetic wrong and the result keeps its shape, raises
# no error, and still correlates ~0.98 with the truth. On the simulation's
# square panels even a dimension check would pass. So the tests are algebraic
# identities, and two of them are NEGATIVE -- they assert that the wrong version
# is detectably wrong, which is what stops the oracle from being vacuous.
#
# Run: Rscript tests/test_kronecker.R

library(tidyverse)
here::i_am("tests/test_kronecker.R")
source(here::here("tests", "helper.R"))
load_code("dge_ige_functions.R")

set.seed(4)

# Deliberately NON-square (3 oats, 4 peas) and of different ranks (2 and 3), so
# that any accidental transpose or index swap cannot conform by luck. The
# simulation's panels are square, which is exactly why the test's are not.
A <- matrix(round(rnorm(3 * 2), 2), 3, 2, dimnames = list(paste0("o", 1:3), NULL))
B <- matrix(round(rnorm(4 * 3), 2), 4, 3, dimnames = list(paste0("p", 1:4), NULL))

# ------------------------------------------------------------
# 1. A row of the basis is the flattened outer product of the two coordinate
#    rows, in the order (a1b1, a1b2, ..., a2b1, ...)
# ------------------------------------------------------------

oi <- c(1L, 3L, 2L, 1L); pj <- c(2L, 4L, 1L, 3L)
X <- kron_basis(A, B, oi, pj)

check(identical(dim(X), c(4L, ncol(A) * ncol(B))),
      "kron_basis is n_obs x (rank_oat * rank_pea)")

for (k in seq_along(oi)) {
  expected <- as.vector(t(outer(A[oi[k], ], B[pj[k], ])))  # row-major flatten
  check_near(X[k, ], expected, tol = 1e-12,
             sprintf("basis row %d is the row-major outer product", k))
}

# column c = (p-1)*rank_pea + q must hold A[,p] * B[,q]
kb <- ncol(B)
for (p in seq_len(ncol(A))) for (q in seq_len(kb)) {
  c_idx <- (p - 1L) * kb + q
  check_near(X[, c_idx], A[oi, p] * B[pj, q], tol = 1e-12,
             sprintf("column %d is A[,%d] * B[,%d]", c_idx, p, q))
}

# ------------------------------------------------------------
# 2. X X' reproduces the Kronecker kernel
#
# This is the whole reason the basis exists: a BRR on X implies covariance X X',
# and that must equal G_oat[i,i'] * G_pea[j,j'] or the model being fitted is not
# the model claimed.
# ------------------------------------------------------------

Go <- tcrossprod(A); Gp <- tcrossprod(B)
check_near(X %*% t(X), Go[oi, oi] * Gp[pj, pj], tol = 1e-12,
           "X X' == G_oat[i,i'] * G_pea[j,j']")

# ------------------------------------------------------------
# 3. The reshape identity
#
# The prediction from the fitted basis must equal the prediction read off the
# reshaped surface, at every observed cell AND at cells never observed -- the
# latter being the point, since predicting unsown combinations cheaply is why
# the reshape exists.
# ------------------------------------------------------------

b <- rnorm(ncol(A) * ncol(B))
Beta_right <- matrix(b, nrow = ncol(A), ncol = ncol(B), byrow = TRUE)
surface <- A %*% Beta_right %*% t(B)

all_oi <- rep(seq_len(nrow(A)), each = nrow(B))
all_pj <- rep(seq_len(nrow(B)), times = nrow(A))
from_basis <- as.vector(kron_basis(A, B, all_oi, all_pj) %*% b)
from_surface <- surface[cbind(all_oi, all_pj)]

check_near(from_basis, from_surface, tol = 1e-12,
           "reshaped surface == basis prediction, at every cell")
check(identical(dim(surface), c(nrow(A), nrow(B))),
      "the surface is n_oat x n_pea")

# and the form the code actually uses: as.vector(t(Beta)) unrolls row-major to
# match the basis column order
check_near(as.vector(kron_basis(A, B, oi, pj) %*% as.vector(t(Beta_right))),
           surface[cbind(oi, pj)], tol = 1e-12,
           "as.vector(t(Beta)) matches the basis column order")

# ------------------------------------------------------------
# 4. NEGATIVE tests: the wrong versions must be detectably wrong
#
# Without these the identities above could pass for a trivial reason -- say if
# Beta happened to be symmetric -- and the oracle would be worthless.
# ------------------------------------------------------------

# byrow forgotten: Beta comes back transposed. With rank_oat != rank_pea the
# dimensions do not even conform, which is why the test matrices differ in rank;
# with equal ranks (the simulation's case) it conforms and is silently wrong.
Aq <- A[, 1:2]; Bq <- B[, 1:2]          # equal ranks, so it CAN conform
bq <- rnorm(4)
Beta_wrong <- matrix(bq, 2, 2)           # no byrow
Beta_ok    <- matrix(bq, 2, 2, byrow = TRUE)
basis_q <- as.vector(kron_basis(Aq, Bq, all_oi, all_pj) %*% bq)
surf_wrong <- (Aq %*% Beta_wrong %*% t(Bq))[cbind(all_oi, all_pj)]
surf_ok    <- (Aq %*% Beta_ok    %*% t(Bq))[cbind(all_oi, all_pj)]

check_near(surf_ok, basis_q, tol = 1e-12, "equal-rank case: byrow=TRUE is right")
check(max(abs(surf_wrong - basis_q)) > 1e-6,
      "dropping byrow gives a DIFFERENT surface (the silent bug is detectable)")
check(identical(dim(Aq %*% Beta_wrong %*% t(Bq)), c(nrow(Aq), nrow(Bq))),
      "and it still has the right shape, which is why a shape check would miss it")

# each/times swapped inside kron_basis: the same damage from the other end
kron_swapped <- function(A, B, oat_idx, pea_idx) {
  k <- ncol(A)
  A[oat_idx, rep(seq_len(k), times = k), drop = FALSE] *
    B[pea_idx, rep(seq_len(k), each = k), drop = FALSE]
}
X_sw <- kron_swapped(Aq, Bq, all_oi, all_pj)
check(max(abs(as.vector(X_sw %*% bq) - basis_q)) > 1e-6,
      "swapping each/times gives a DIFFERENT basis")
# it still reproduces the kernel, which is the trap: X X' is invariant to the
# column ORDER, so only the reshape identity catches this one
check_near(X_sw %*% t(X_sw),
           tcrossprod(Aq)[all_oi, all_oi] * tcrossprod(Bq)[all_pj, all_pj],
           tol = 1e-12,
           "the swapped basis STILL reproduces the kernel -- X X' cannot detect it")

finish("kronecker tests")
