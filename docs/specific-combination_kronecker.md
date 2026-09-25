# The specific-combination term: the Kronecker basis and the reshape

How the oat × pea interaction term in `code/sim_fit.R` is actually fitted, and
why the two lines that build it are the most fragile arithmetic in the
simulation. Written up because it is index arithmetic that nothing in the code
explains and nothing but level S6 of
[EVALUATION_SIMULATION.md](../EVALUATION_SIMULATION.md) checks.

For why the term is low-rank at all see
[SIMULATION.md](../SIMULATION.md#the-kronecker-term-is-low-rank-by-necessity);
for the model it belongs to see [BACKGROUND.md](../BACKGROUND.md). The code is
`grm_basis()` and `kron_basis()` in `code/sim_fit.R`, used by `fit_dge_ige()`.

---

## What we are trying to fit

The specific-combination effect `I_ij` of oat *i* with pea *j* is assumed to
have covariance

```
cov(I_ij , I_i'j')  =  G_oat[i,i'] × G_pea[j,j']
```

Two combinations are correlated to the extent that **both** their oats and
their peas are related. That is the Kronecker structure `G_oat ⊗ G_pea`.

There are two ways to fit it in BGLR.

- An `RKHS` term takes that kernel directly. But the kernel is
  n_obs × n_obs, and at 19,200 observations it is a 2.9 GB matrix to
  eigendecompose.
- A `BRR` term takes a design matrix `X` and implies covariance `X X'`. So a
  **skinny** `X` with `X X'` equal to that kernel gives the same model cheaply.

That `X` is the basis, and building it is what the rest of this document is
about.

---

## Step 1 — coordinates for each accession

`grm_basis()` eigendecomposes a relationship matrix, `G = U Λ U'`, and returns
`A = U Λ^½` truncated to the leading `rank` columns, so that

```
A A' = G_oat
```

Read `A[i,]` as **coordinates for oat *i*** in a small space, scaled so that the
dot product of two rows is exactly their relatedness. A worked example with
3 oats and rank 2:

```
        a1   a2
oat1  0.22 0.60
oat2 -0.54 1.64
oat3  0.89 0.69
```

The same for peas gives `B`, with `B B' = G_pea`:

```
        b1    b2
pea1 -1.28  0.57
pea2 -0.21  0.02
pea3  1.90  0.38
pea4  1.78 -0.05
```

---

## Step 2 — the identity, and the rows we actually need

```
G_oat ⊗ G_pea  =  (A A') ⊗ (B B')  =  (A ⊗ B)(A ⊗ B)'
```

So `A ⊗ B` is a design matrix with exactly the covariance we want. But
`A ⊗ B` carries a row for **every** oat × pea cell, and only some were sown. We
need only the rows we have.

**Row *k* of the design, for the combination (oat *i*, pea *j*), is the outer
product of `A[i,]` with `B[j,]`, flattened.** With rank 2 that is 2 × 2 = 4
numbers per row:

```
             a1*b1   a1*b2   a2*b1   a2*b2
oat1+pea2  -0.0462  0.0044 -0.1260  0.0120
oat3+pea4   1.5842 -0.0445  1.2282 -0.0345
oat2+pea1   0.6912 -0.3078 -2.0992  0.9348
```

That is all `kron_basis()` does:

```r
kron_basis <- function(A, B, oat_idx, pea_idx) {
  k <- ncol(A)
  A[oat_idx, rep(seq_len(k), each  = k), drop = FALSE] *
  B[pea_idx, rep(seq_len(k), times = k), drop = FALSE]
}
```

The two `rep()` calls choose which column of `A` and which of `B` multiply into
each output column:

| | value | picks |
|---|---|---|
| `rep(1:k, each = k)` | `1,1,2,2` | the oat direction *p* |
| `rep(1:k, times = k)` | `1,2,1,2` | the pea direction *q* |

so column `c = (p−1)·k + q` holds `A[,p] × B[,q]`. Formally this is the
**row-wise Kronecker** (or *Khatri–Rao*, or *face-splitting*) product — a
Kronecker product taken row by row rather than over the whole matrices.

**So the basis is every pairing of one oat eigen-direction with one pea
eigen-direction**: k² columns, spanning interaction surfaces of rank up to k².

It does what it was built to do — `X X'` reproduces the kernel exactly:

```
X %*% t(X)                              G_oat[i,i'] * G_pea[j,j']
          oat1+pea2 oat3+pea4 oat2+pea1          oat1    oat3    oat2
oat1+pea2    0.0182   -0.2286    0.2424        0.0182 -0.2286  0.2424
oat3+pea4   -0.2286    4.0213   -1.5018   =   -0.2286  4.0213 -1.5018
oat2+pea1    0.2424   -1.5018    5.8530        0.2424 -1.5018  5.8530
```

---

## Step 3 — the reshape

BGLR returns a coefficient vector `b` of length k². The fitted effect for an
observed combination is

```
Î_ij  =  Σ_c X[k,c]·b_c  =  Σ_p Σ_q  A[i,p] · b_{(p−1)k+q} · B[j,q]
```

Fold `b` into a k × k matrix with `Beta[p,q] = b_{(p−1)k+q}`, and that double
sum **is** a matrix product:

```
Î_ij  =  (A %*% Beta %*% t(B))[i,j]
```

which is what `fit_dge_ige()` returns as `interaction`.

That is the point of the reshape. It gives the **entire** n_oat × n_pea
surface — every cell, including combinations nobody grew — from three small
matrices, instead of building the full basis for all n_oat·n_pea cells. Cheap
prediction everywhere is the whole reason for going through the basis rather
than the kernel.

---

## Why it is fragile

`Beta[p,q] = b[(p−1)k+q]` is **row-major**. R's `matrix()` fills
**column-major**. Hence:

```r
Beta <- matrix(fit$ETA$int$b, nrow = ncol(A), ncol = ncol(B), byrow = TRUE)
```

Drop `byrow` and `Beta` comes back transposed:

```
correct (byrow = TRUE)      wrong (R's default)
    10   20                     10   30
    30   40                     20   40
```

and the surface becomes `A Beta' B'` instead of `A Beta B'`. In the toy example
the two differ by 33 units. But note what does **not** happen:

- no error is raised;
- the result is still 3 × 4, the right shape;
- it correlates **0.983** with the correct surface.

And in `code/sim_generate.R` the panels are square (`n_acc × n_acc`), so even a
dimension check would pass. The outcome is a plausible-looking interaction that
quietly mismatches oat coordinates with pea coordinates, making `dge_ige` look
slightly worse than it is — which is exactly the kind of error that would be
read as a result about the method rather than a bug.

Swapping `each` and `times` inside `kron_basis()` does the same damage from the
other end.

---

## The check

One line, in level S6 of [EVALUATION_SIMULATION.md](../EVALUATION_SIMULATION.md).
The prediction from the fitted basis must equal the prediction read off the
reshaped surface:

```r
lhs <- as.vector(kron_basis(A, B, oi, pj) %*% as.vector(t(Beta)))
rhs <- (A %*% Beta %*% t(B))[cbind(oi, pj)]
max(abs(lhs - rhs))        # ~7e-18
```

`as.vector(t(Beta))` unrolls row-major, matching the basis column order. It
agrees to machine precision, and it fails immediately if either `each`/`times`
or `byrow` is flipped. **Re-run it after touching either function.**

---

## The truncation, which is a real limitation

`SIM_KRON_RANK = 30` means `A A' ≈ G_oat`, not `=`. The interaction is fitted in
a 900-column subspace of the full Kronecker space, chosen as the leading
eigen-directions of each species.

MegaLMM's interaction is not truncated that way, so the head-to-head comparison
gives `dge_ige` a rank-limited interaction and MegaLMM a free one. When the two
are close, check whether raising `SIM_KRON_RANK` moves the result. Setting it to
`NA` builds the exact kernel instead, which is tractable only below
`SIM_KRON_EXACT_MAX = 3000` observations.

---

## Reproducing the examples

```r
library(tidyverse)
here::i_am("code/sim_fit.R")
source(here::here("code", "sim_fit.R"))

set.seed(4)
A <- matrix(round(rnorm(6), 2), 3, 2,
            dimnames = list(paste0("oat", 1:3), c("a1", "a2")))
B <- matrix(round(rnorm(8), 2), 4, 2,
            dimnames = list(paste0("pea", 1:4), c("b1", "b2")))

# the basis, for three observed combinations
oi <- c(1, 3, 2); pj <- c(2, 4, 1)
X <- kron_basis(A, B, oi, pj)
round(X %*% t(X), 4)                        # equals ...
round((A %*% t(A))[oi, oi] * (B %*% t(B))[pj, pj], 4)

# the reshape, and what getting it wrong costs
b <- c(10, 20, 30, 40)
Beta_right <- matrix(b, 2, 2, byrow = TRUE)
Beta_wrong <- matrix(b, 2, 2)
oi <- rep(1:3, each = 4); pj <- rep(1:4, times = 3)
basis <- as.vector(kron_basis(A, B, oi, pj) %*% b)
max(abs((A %*% Beta_right %*% t(B))[cbind(oi, pj)] - basis))   # ~7e-15
max(abs((A %*% Beta_wrong %*% t(B))[cbind(oi, pj)] - basis))   # 33.2
```
