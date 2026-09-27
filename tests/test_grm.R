# test_grm.R
#
# The relationship-matrix machinery in code/dge_ige_functions.R. Everything
# downstream joins on accession names and this is where names get rewritten, so
# these are the cheapest tests in the suite and the ones whose failure would be
# hardest to notice: a wrong collapse changes every variance component without
# changing any dimension.
#
# Run: Rscript tests/test_grm.R

library(tidyverse)
here::i_am("tests/test_grm.R")
source(here::here("tests", "helper.R"))
load_code("dge_ige_functions.R")

set.seed(1)

# ------------------------------------------------------------
# 1. grm_factor: L L' reproduces G
#
# The BGLR terms are a BRR on Z L, which is the RKHS model with kernel Z G Z'
# only if L L' == G. The production script asserts this at run time; here it is
# checked on matrices with known awkwardness.
# ------------------------------------------------------------

make_grm <- function(n, rank = n) {
  X <- matrix(rnorm(n * rank), n, rank)
  G <- tcrossprod(X) / rank
  dimnames(G) <- list(paste0("a", seq_len(n)), paste0("a", seq_len(n)))
  G
}

for (case in list(list(n = 20, r = 20), list(n = 20, r = 8), list(n = 5, r = 1))) {
  G <- make_grm(case$n, case$r)
  L <- grm_factor(G)
  check_near(tcrossprod(L), G, tol = 1e-10,
             sprintf("L L' == G (n=%d rank=%d)", case$n, case$r))
  check(identical(rownames(L), rownames(G)),
        sprintf("grm_factor keeps rownames (n=%d)", case$n))
  check(ncol(L) <= case$r + 1, "grm_factor drops null directions")
}

# A prior-only line -- diagonal 1, no covariance with anything -- is what
# build_grm() injects for an ungenotyped accession, and 8 oat / 20 pea of them
# are in the real GRMs. It must survive factoring.
G <- make_grm(10)
G[10, ] <- 0; G[, 10] <- 0; G[10, 10] <- 1
L <- grm_factor(G)
check_near(tcrossprod(L), G, tol = 1e-10, "L L' == G with a prior-only row")

# ------------------------------------------------------------
# 2. collapse_grm: averaging rows and columns of a group
#
# The oracle: collapsing a group of IDENTICAL lines must return the same matrix
# with one row per name, because averaging the relationships of identical lines
# is the same as taking the GRM of their averaged marker profiles. Built from a
# marker matrix so the answer is known independently of the function.
# ------------------------------------------------------------

n_mark <- 50
M <- matrix(rnorm(6 * n_mark), 6, n_mark)
M[2, ] <- M[1, ]           # a2 is a duplicate of a1
M[3, ] <- M[1, ]           # and so is a3
rownames(M) <- paste0("a", 1:6)

# NOT centred, deliberately. Centring uses the column means of whatever rows
# are present, so a 6-row matrix and its 4-row collapse would be centred
# differently and the oracle would compare two different quantities. Without
# centring, "average the GRM rows of a group" and "average the marker rows of a
# group" are the same operation and the identity is exact.
grm_of <- function(M) {
  G <- tcrossprod(M) / ncol(M)
  dimnames(G) <- list(rownames(M), rownames(M))
  G
}
G_full <- grm_of(M)

name_file <- tempfile(fileext = ".csv")
readr::write_csv(tibble::tibble(germplasmName = c("a2", "a3"),
                                analysis_name = c("a1", "a1")), name_file)

G_coll <- collapse_grm(G_full, name_file)
check(nrow(G_coll) == 4, "collapse_grm folds 6 lines into 4")
check(identical(rownames(G_coll), c("a1", "a4", "a5", "a6")),
      "collapse_grm keeps the surviving names in order")

# the independent answer: the GRM of the marker matrix with the duplicates
# averaged away (averaging identical rows changes nothing)
M_avg <- M[c(1, 4, 5, 6), ]
check_near(G_coll, grm_of(M_avg)[rownames(G_coll), rownames(G_coll)], tol = 1e-10,
           "collapsed GRM == GRM of the collapsed markers")

# a missing or malformed mapping must be a no-op, not an error
check(identical(collapse_grm(G_full, tempfile()), G_full),
      "collapse_grm no-ops on a missing file")
bad <- tempfile(fileext = ".csv")
readr::write_csv(tibble::tibble(wrong = 1, columns = 2), bad)
check(identical(collapse_grm(G_full, bad), G_full),
      "collapse_grm no-ops on a mapping without the expected columns")

# collapsing onto a name that ALREADY exists merges the two, which is the pea
# case (NDP170084G -> ND VICTORY): the result must lose exactly one row
nf2 <- tempfile(fileext = ".csv")
readr::write_csv(tibble::tibble(germplasmName = "a5", analysis_name = "a6"), nf2)
check(nrow(collapse_grm(G_full, nf2)) == 5,
      "collapse_grm merges onto an existing name")

# ------------------------------------------------------------
# 3. incidence and dummy_matrix
#
# model.matrix() DROPS rows with an NA factor level rather than complaining, so
# an accession absent from the GRM would silently shorten the design matrix and
# BGLR would fail somewhere unrelated. And it refuses a single-level factor
# outright, which is what fit_producer_associate() hit on a one-environment
# simulated experiment.
# ------------------------------------------------------------

lev <- paste0("a", 1:4)
Z <- incidence(c("a1", "a3", "a3", "a2"), lev)
check(identical(dim(Z), c(4L, 4L)), "incidence is n x levels")
check(identical(colnames(Z), lev), "incidence columns are the levels, in order")
check(all(rowSums(Z) == 1), "incidence has exactly one 1 per row")
check(Z[2, "a3"] == 1 && Z[3, "a3"] == 1, "incidence puts the 1 in the right place")

check(identical(dim(dummy_matrix(factor(rep("only", 5)))), c(5L, 1L)),
      "dummy_matrix handles a single-level factor")
check(all(dummy_matrix(factor(rep("only", 5))) == 1),
      "a single-level dummy is a column of ones")
check(identical(dim(dummy_matrix(factor(rep(c("a", "b", "c"), 2)))), c(6L, 3L)),
      "dummy_matrix handles several levels")
f <- factor(c("a", "b"), levels = c("a", "b", "unused"))
check(ncol(dummy_matrix(f)) == 2, "dummy_matrix drops unused levels")

# ------------------------------------------------------------
# 4. grm_basis: truncation by rank and by cumulative variance
# ------------------------------------------------------------

G <- make_grm(30)
check(ncol(grm_basis(G, rank = 5)) == 5, "grm_basis honours rank")
check(ncol(grm_basis(G, rank = 500)) <= 30,
      "grm_basis cannot return more directions than the matrix has")

# Truncating at a variance target does NOT reconstruct G exactly -- it drops the
# tail on purpose -- so the invariant is that the relative error is bounded by
# roughly the share dropped, and falls as the target rises.
rel_err <- function(v) {
  B <- grm_basis(G, variance = v)
  norm(tcrossprod(B) - G, "F") / norm(G, "F")
}
check(rel_err(0.999) < 0.01,
      "grm_basis at 99.9% variance reconstructs G to under 1% Frobenius error")
check(rel_err(0.999) < rel_err(0.5),
      "a higher variance target reconstructs G better")
check(ncol(grm_basis(G, variance = 0.5)) < ncol(grm_basis(G, variance = 0.999)),
      "a lower variance target keeps fewer directions")

finish("grm tests")
