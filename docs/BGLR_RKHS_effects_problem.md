# Why the old script's "producer" and "associate" effects were not accession effects

This explains the one bug in the original `code/BGLR_multi_trait_model.R` that
silently produced wrong answers. The variance components it estimated were
fine. The per-accession effects — the `GMA` rankings in
`BGLR_oat_rank_stability.csv`, and the points in the producer-vs-associate
scatterplots — were not.

---

## 1. What the script asked BGLR for

The genetic terms were set up as plot-level kernels:

```r
K_oat_plot <- incOatAcc %*% K_oat %*% t(incOatAcc)   # Z G Z'
ETA_pair <- list(..., G_oat = list(K = K_oat_plot, model = "RKHS"), ...)
```

`Z` is the plot × accession incidence matrix (508 columns, one per oat
accession) and `G` is the oat GRM. So `K_oat_plot` is `n × n`, where `n` is the
number of **plots**, not accessions.

This is a legitimate way to write the model. Saying

> the plots have covariance `Z G Z' ⊗ Σ`

is mathematically identical to saying

> each accession has an effect `g`, with `g ~ N(0, G ⊗ Σ)`, and plot *i* gets
> the effect of whichever accession it carries.

Both describe the same distribution over the data. The difference is what the
fitted object hands back to you.

## 2. What BGLR actually does with an RKHS term

`BGLR::Multitrait` does not sample a random effect per plot. It first
eigen-decomposes the kernel. From `BGLR:::setLT.RKHS_mt`:

```r
LT$EVD  <- eigen(LT$K, symmetric = TRUE)
keep    <- LT$EVD$values > 1e-10
LT$EVD$vectors <- LT$EVD$vectors[, keep]
LT$EVD$values  <- LT$EVD$values[keep]
LT$X    <- sweep(LT$EVD$vectors, MARGIN = 2, STATS = sqrt(LT$EVD$values), FUN = "*")
LT      <- setLT.BRR_mt(LT = LT, ...)     # <- then treated as a ridge regression
```

So an RKHS term is internally a **ridge regression on the eigenvectors of the
kernel**. Writing `K = V D V'`, it sets `X = V D^(1/2)` and samples
coefficients `b` with `b ~ N(0, I ⊗ Ω)`. The random effect on the plots is then
reconstructed at the end:

```r
ETA[[j]]$u <- ETA[[j]]$X %*% ETA[[j]]$post_beta
```

Two objects come out of this, and they are different things:

| object | dimension | meaning |
|---|---|---|
| `ETA[[j]]$beta` | rank(K) × 2 | coefficients **in the eigenvector basis of the plot kernel** |
| `ETA[[j]]$u` | n_plots × 2 | the genetic effect **per plot** |

Neither one is an effect per accession. The old script read `$beta`.

## 3. Where the accession identity went

Row *r* of `beta` is the loading on the *r*-th **eigenvector** of `Z G Z'`,
sorted by decreasing eigenvalue. An eigenvector is a weighted contrast across
all accessions at once — it is not "accession *r*". Eigenvector order is
determined by the spectrum of the kernel, which has nothing to do with the
alphabetical (or any other) order of `rownames(G)`.

The accession labels were in `colnames(Z)`, and the matrix product
`Z %*% G %*% t(Z)` throws them away: the result is indexed by plot, and the
eigendecomposition then re-indexes it by eigenvector. By the time `beta`
exists, there is no accession information left in it to recover.

## 4. Why nothing looked wrong

Three coincidences hid it:

1. **The row count is right.** `rank(Z G Z')` equals the number of distinct
   accessions (as long as `G` is full rank), because `Z` has exactly that many
   linearly independent columns. So `beta` had 508 rows for 508 oat accessions
   — exactly what you would expect from a correct extraction. The output CSV
   had the right number of lines.

2. **`beta` has no rownames, and `rownames_to_column()` invents them.** For a
   data frame with default row names, `tibble::rownames_to_column("oatAcc")`
   silently writes `"1"`, `"2"`, `"3"`, … So the `oatAcc` column was filled
   with eigenvector indices formatted as strings. Nothing errored, and unless
   you looked at the column you would not notice it held numbers rather than
   names like `IL24-2357`.

3. **The correlation in the scatterplot was roughly right.** The rows of `beta`
   are i.i.d. `N(0, Ω)` by construction, so the sample correlation between its
   two columns estimates the same producer–associate correlation that the
   accession-level effects have. The `r =` annotation on the plots was
   therefore approximately correct *even though not one plotted point was an
   accession*. This is the most misleading part: the summary statistic
   validated, so the plot looked trustworthy.

## 5. Numerical demonstration

Simulated data, 40 oat accessions across 1920 plots, with known true producer
and associate effects. The same data fitted both ways:

```
accessions: 40   plots: 1920   rank(ZGZ'): 40
dim(beta) from RKHS term: 40 2   rownames: NULL

OLD: cor(beta[,2], true Pr_oat)        = 0.012
OLD: cor(beta[,1], true As_oat->pea)   = 0.198
NEW: cor(L%*%beta[,2], true Pr_oat)    = 0.939
NEW: cor(L%*%beta[,1], true As->pea)   = 0.879

OLD: cor(beta cols) = 0.385  | NEW: cor(accession Pr, As) = 0.39  | simulated = 0.5
```

The old extraction correlates with the truth at 0.01 and 0.20 — it is noise.
The new one recovers the true effects at 0.94 and 0.88. And note the last line:
the *correlation* was estimated correctly by both (0.385 vs 0.39), which is
exactly why the diagnostic plots did not give the problem away.

Practically: the top of `BGLR_oat_rank_stability.csv` was a list of
eigenvectors in descending eigenvalue order, labelled 1, 2, 3, … Selecting oat
accession "1" from it would have meant selecting whichever accession happened
to sort first in the GRM.

## 6. What was *not* affected

- **Variance components.** `Ω` and the residual covariance were estimated from
  the correct model, because the plot-level kernel is a correct
  reparameterization. Re-fitting the same data both ways gives the same `Ω`
  (1.538 / 0.394 / 0.981 old vs 1.527 / 0.388 / 0.970 new — pure Monte Carlo
  difference).
- **The kernel structure.** Three RKHS terms plus an unstructured residual is
  the right translation of the proposal's model.
- **The trait ↔ role mapping.** With `Y = (peaYield, oatYield)`, the oat
  kernel's first column really is the associate effect on pea and the second
  the producer effect on oat. The old script's column naming was correct.

## 7. The fix, and why it is the same model

Factor the GRM as `G = L L'` (via its eigendecomposition, dropping numerically
null dimensions), and fit a `BRR` term on `X = Z L` instead of an `RKHS` term
on `Z G Z'`:

```r
grm_factor <- function(G, tol = 1e-8) {
  e <- eigen(G, symmetric = TRUE)
  keep <- e$values > tol * max(e$values)
  L <- sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "*")
  rownames(L) <- rownames(G)
  L
}

ETA$G_oat <- list(X = Z_oat %*% L_oat, model = "BRR")
```

With `b ~ N(0, I ⊗ Σ)`, the implied plot covariance is

```
Var(Z L b) = Z L L' Z' ⊗ Σ = Z G Z' ⊗ Σ
```

— identical to the RKHS term. But now the accession effects are one
multiplication away, and they keep their names:

```r
g <- L_oat %*% fit$ETA$G_oat$beta     # nAccessions × 2
rownames(g) <- rownames(L_oat)        # inherited from rownames(G_oat)
```

`Var(g) = L (I ⊗ Σ) L' = G ⊗ Σ`, which is exactly the distribution the proposal
specifies for `(Pr_oat, As_oat→pea)`.

It is also cheaper: BGLR no longer eigen-decomposes three
`n_plot × n_plot` matrices, only the much smaller GRMs, and it does so once
outside the sampler rather than inside `Multitrait`.

## 8. How to check it yourself

The one-line test on any future fit:

```r
rownames(L_oat %*% fit$ETA$G_oat$beta)   # must be accession names, not NULL
```

If the effects you are about to rank do not carry accession names that came
from `rownames(G)`, you are not looking at accession effects.
