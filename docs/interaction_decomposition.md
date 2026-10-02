# Reading the DGE-IGE interaction as scores and loadings

`dge_ige` predicts the oat × pea interaction better than MegaLMM when
combinations are sparse, and — measured at full rank in
`output/simulation_kron_study.csv` — also when the interaction is higher-rank.
But what it returns is a number per combination, `interaction$oat[i, j]`, not a
score for an oat and a loading for a pea. MegaLMM returns the interpretable
object and loses the accuracy.

This note is about recovering the interpretable object from the better model: the
algebra, what it is for biologically, and what the sweep measured. Code in
`code/interaction_decomp.R`, driven by `code/sim_decomp_run.R`, tested in
`tests/test_decomp.R` and `tests/test_fits.R`.

**Sections 1–3 build the algebra from the ground up.** Skip to
[§4](#4-what-the-sweep-measured) for results and
[§5](#5-what-any-of-this-means-biologically) for what they mean.

---

## 1. The algebra, from the ground up

### 1.1 A relationship matrix, and what factoring it buys

`G_oat` is 442 × 442 for the real panel: `G_oat[i, i']` is how related oat *i*
is to oat *i'*. We want a model where each oat has an effect and **related oats
have similar effects** — formally, effects `g` with `Var(g) = σ² G_oat`.

BGLR does not take a covariance matrix. It takes a design matrix and puts
independent, equally-shrunk coefficients on its columns (that is what `BRR`,
Bayesian ridge regression, means). So the covariance has to be smuggled in
through the design matrix. The trick is to **factor** `G`.

**Eigendecomposition.** Any symmetric matrix can be written

```
G = U Λ U'
```

- **`U`** is `n × n`, and its columns `u_1, u_2, …` are the **eigenvectors**.
  Each one is a *direction* in the space of accessions — a list of one weight
  per oat. `u_1` might be "northern material versus southern", `u_2` "hulled
  versus naked", and so on; they are whatever the marker data says the major
  axes of relatedness are. They are **orthonormal**: each has length 1, and any
  two are at right angles, so they describe genuinely separate things.
- **`Λ`** is diagonal, holding the **eigenvalues** `λ_1 ≥ λ_2 ≥ …`. Eigenvalue
  `λ_k` is how much genetic variance lies along direction `u_k`. Big λ means a
  major axis of population structure; λ near zero means a direction in which the
  panel barely varies.

This is principal components, exactly: `U`'s columns are the PCs of relatedness
and `Λ` says how important each is.

**Now the key step.** Define

```
A = U Λ^½          (each eigenvector scaled by the square root of its eigenvalue)
```

Then `A A' = U Λ^½ Λ^½ U' = U Λ U' = G`. So if we use `A` as the design matrix
and put independent coefficients `b ~ N(0, σ² I)` on its columns, the resulting
effects are

```
g = A b,    Var(g) = σ² A A' = σ² G
```

**exactly the kinship covariance we wanted.** `grm_basis()` and `grm_factor()`
in `code/dge_ige_functions.R` do this.

> **So your intuition is right, and here is the precise version.** Yes — every
> genetic term in the model is a Bayesian ridge regression on principal
> components of a relationship matrix. The `√λ` scaling is what makes it work:
> *equal shrinkage in the PC basis is the same thing as kinship covariance in
> the accession basis.* Without the scaling you would be shrinking every
> direction equally in the wrong metric, and major and minor axes of relatedness
> would be treated alike.

**Two words used throughout.** The **span** of a matrix `A`, written
`span(A)`, is every vector you can build as `A c` for some coefficient vector
`c` — everything reachable from those columns. A **basis** is a set of columns
whose span is the whole space of interest. If `A` has all 442 columns, `span(A)`
is the whole 442-dimensional accession space and `A` is a basis for it. **If `A`
has only 30 columns, `span(A)` is a 30-dimensional slice** — and any effect the
model estimates is confined to that slice. That single fact drives
[§3.3](#33-the-ceiling-and-the-floor) and most of the results.

### 1.2 The interaction, and the Kronecker product

Main effects are per-accession. The interaction is per-*combination*: `I[i, j]`
is how much oat *i* with pea *j* over- or under-performs relative to what oat
*i*'s and pea *j*'s own average contributions predict. That is specific
combining ability.

What covariance should two combinations have? Combination (*i*, *j*) should
resemble (*i'*, *j'*) to the extent that **the oats are related *and* the peas
are related** — a product, because if either pair is unrelated there is no
reason for their specific combining abilities to agree:

```
Cov[(i,j), (i',j')] = G_oat[i,i'] × G_pea[j,j']
```

Written over all pairs of cells at once, that product structure *is* the
**Kronecker product** `G_oat ⊗ G_pea`. It is a big matrix: one row and column
per cell of the grid, so 186,966 × 186,966 for the real panel.

### 1.3 Why the exact kernel cannot do the job

You could avoid the big matrix by building the kernel only over the
**observed** combinations — 2,063 × 2,063 for B4I, which is perfectly
tractable. `fit_producer_associate(kron_rank = NA)` does exactly that.

It fails, and not for want of compute. A design matrix built that way has **one
column per observed combination**. A combination nobody grew has no column, so
the model has no way to express its interaction — there is no coefficient to
look up. Measured: `r_int` comes back `NA`, because the returned surface has no
interaction component at the never-observed cells at all. Predicting a new pair
would need the cross-covariance between it and every training combination (a
BLUP extension), which is not implemented.

**And the never-observed cells are the entire question.** The programme wants to
know which *untried* oat × pea pairings will work.

### 1.4 The row-wise Kronecker basis, which predicts unobserved cells for free

Since `G_oat = A A'` and `G_pea = B B'`, a standard Kronecker identity gives

```
G_oat ⊗ G_pea = (A A') ⊗ (B B') = (A ⊗ B)(A ⊗ B)'
```

So `A ⊗ B` is a design matrix with exactly the covariance we want — and now
look at **one row of it**. The row for cell (*i*, *j*) is the outer product of
oat *i*'s coordinates `A[i, ]` and pea *j*'s coordinates `B[j, ]`, flattened:

```
row(i,j) = [ A[i,1]·B[j,1],  A[i,1]·B[j,2],  …,  A[i,q]·B[j,q] ]
```

It is built **entirely from the two accessions' own coordinates**. Nothing about
whether that pairing was ever grown enters. So a combination nobody planted has
a perfectly well-defined row, and the model can predict it. That is the whole
reason this construction is used rather than the exact kernel.

`kron_basis()` (`code/dge_ige_functions.R:168`) builds those rows for the
observed plots.

### 1.5 `kron_rank`: what is truncated, and what is not

`A ⊗ B` has `n_oat × n_pea` columns — 186,966 for the real panel. That is the
problem. So we keep only the **leading `q` eigenvectors of each species**, where
`q = kron_rank`, giving `q²` columns: at `q = 30`, **900 columns**.

> **Correcting the natural reading.** Three things are easy to get wrong here,
> and I had two of them wrong at first.
>
> 1. **It is `q` directions per species, not `q` components of the interaction.**
>    The interaction directions are *products* of an oat direction and a pea
>    direction, so `q` from each gives `q²` interaction directions. At
>    `q = n_acc` the products span the *entire* Kronecker space exactly —
>    nothing is lost — which is what `sim_kron_study.R` exploits.
> 2. **The main-effect terms are NOT truncated.** `grm_factor()` keeps every
>    direction with a non-negligible eigenvalue, so `G_oat` and `G_pea` enter at
>    full rank. Truncation is specific to the interaction term, because only
>    there does the column count multiply.
> 3. **Tractability is not the only reason, and not the main one.** The exact
>    kernel is tractable at B4I's 2,063 combinations — it is just useless for
>    prediction (§1.3). The low-rank basis is chosen because it *extends to
>    unobserved cells*; keeping the column count down is a bonus.
>
> And the sharp edge: the directions kept are chosen by **GRM eigenvalue** — by
> how much genetic variance each carries — **not by where the interaction
> actually lives.** Nothing guarantees the interaction aligns with major
> population structure. Rank 30 keeps 62.5% of the trace of the real oat GRM
> (508 genotyped accessions), and `30 × 30 = 900` of the 186,966 interaction
> directions the phenotyped grid spans (442 oats × 423 peas).

### 1.6 The `q × q` reshape, and what `I_hat` is

The 900 fitted coefficients are indexed by a **pair**: oat direction *a* paired
with pea direction *b*. `kron_basis()` lays them out row-major — all pea
directions for oat direction 1, then all for oat direction 2, and so on — so
column `(a−1)q + b` carries `A[,a] · B[,b]`.

Which means the coefficient *vector* is naturally a **`q × q` table**: rows
indexed by oat direction, columns by pea direction. That is `Beta`, and
`byrow = TRUE` in the reshape is what respects the row-major layout.
`beta_from_vector()` is the one place it happens; getting it wrong keeps the
shape, raises no error, and still correlates ~0.98 with the truth, which is why
`tests/test_kronecker.R` and `tests/test_decomp.R` both pin it.

Once it is a matrix, the whole fitted surface follows:

```
I_hat = kron_A %*% Beta %*% t(kron_B)
```

**`I_hat` is the fitted interaction surface**: `n_oat × n_pea`, one number per
cell — every oat × pea pairing, grown or not. A `q × q` table of coefficients
expands to the full grid. That is `code/dge_ige_functions.R:392`, and it is the
object everything below decomposes.

---

## 2. The decomposition

### 2.1 It is already a scores-and-loadings object

`I_hat = kron_A Beta kron_B'` is a product of three matrices, which makes it a
**bilinear form**: linear in the oat coordinates and linear in the pea
coordinates. Such a thing can always be rewritten as a sum of simple layers, by
the **singular value decomposition**.

Any matrix `M` can be written `M = U D V'` with `U` and `V` orthonormal and `D`
diagonal, non-negative and decreasing. Written out as a sum:

```
M = d_1 u_1 v_1'  +  d_2 u_2 v_2'  +  …
```

Each term is one **layer**: a single oat pattern `u_k` times a single pea
pattern `v_k`, scaled by `d_k`. Layer *k* says "oats with a high `u_k` do well
specifically with peas that have a high `v_k`". **`u_k` is the oat score, `v_k`
is the pea loading** — the object MegaLMM gives directly and `dge_ige` was
thought not to.

The layers are ordered: `d_1` is the largest, so layer 1 is the single best
rank-1 summary of the whole surface.

In code (`decompose_bilinear()`), after double-centring the two bases:

```
M = R_a Beta R_b'      (a small q × q matrix);   svd(M) = U D V'
oat scores   = Q_a %*% U        pea loadings = Q_b %*% V
```

No `n_oat × n_pea` matrix ever has to be formed. `decompose_surface()` does the
dense version and the two are tested to agree.

### 2.2 Why the shares are exact, not approximate

Because `U` and `V` are orthonormal, the layers do not overlap, so their
contributions add with no double counting:

```
||M||²_F = Σ_k d_k²
```

`||M||²_F`, the **Frobenius norm squared**, is just the sum of every squared
entry of the surface — its total size, and proportional to the variance of the
interaction across the grid. So

```
share of layer k  =  d_k² / Σ_j d_j²
```

is an **exact variance share**, not an approximation. Verified as algebra at
`tol = 1e-10`: reconstruction exact, both factor matrices orthonormal, centred
scores with zero column means, `Σd²` equal to the interaction sum of squares to
six figures, and a rank-1 `Beta` returning shares `1.000 0 0 …`.

### 2.3 Why there is no deregression step

The natural worry about analysing model output is that the inputs are shrunken
BLUPs: you would need to deregress them and weight each by its expected error
before fitting anything to them.

That worry applies to **fitting a new model**. Nothing is fitted here. The SVD
is a change of coordinates on a term the model already estimated — the same
surface, written in a different basis. No information is created, none is
double-counted, and there is no deregression step to get wrong. Uncertainty
comes from doing the rotation **inside each MCMC draw**, not from correcting a
point estimate.

Deregression becomes necessary one step later — regressing the recovered scores
on biological covariates ([§5.2](#52-turning-a-score-into-a-mechanism)). The
reliability and PEV machinery at `code/validation_functions.R:152-172` is what to
reuse then.

### 2.4 The participation ratio: how many layers are really there

Counting layers by eye needs a threshold. The **participation ratio** does not:

```
participation = (Σ d_k²)² / Σ d_k⁴  =  1 / Σ p_k²       where p_k = d_k²/Σd²
```

One dominant layer gives 1. Five equal layers give 5. Anything in between gives
a sensible fraction. It is a continuous *effective number of layers*.

> If that formula looks familiar: it is the **inverse Simpson index**. Ecologists
> use exactly this statistic to turn a list of species abundances into an
> "effective number of species". Here the abundances are the layers' variance
> shares, and the answer is the effective number of interaction layers. The
> companion `n90` is the blunter version — how many layers to reach 90% of the
> variance.

### 2.5 Three numbers, and why the obvious one is wrong

`component_summary()` returns three quantities per layer. Choosing the wrong one
as the headline is easy, and was caught only by noticing that a point estimate
sat outside its own interval.

| quantity | what it is |
|---|---|
| **`share1`** (+ `share1_lo/hi`) | **The estimand.** Each MCMC draw is a posterior sample of the *true* surface, so the posterior of a draw's own layer 1 is the posterior of layer 1's share of the real interaction. The only one of the three with an honest interval. |
| `share1_meansurf` | Layer 1's share of the posterior-**mean** surface. The only share comparable with MegaLMM, which has no streamed draws — but **not** an estimate of the above. |
| `proj1_mean` (+ interval) | How much of a typical draw's interaction lies in the direction the mean surface picked out. An identifiability diagnostic, not a share of anything. |

`share1_meansurf` runs systematically high because singular values are **convex**
in the matrix: averaging the draws cancels their idiosyncratic directions and
leaves a surface that looks more low-rank than anything real. In the sweep the
inflation is **+0.228 at 1.6% observed**, falling to +0.109 at 48% — worst
exactly where the data are thinnest.

### 2.6 Rotation: what can and cannot be named

When two singular values are close, "layer 2" is not a well-defined object — the
pair spans a plane and any rotation within it fits equally well. No alignment
scheme fixes that. So:

- **rotation-invariant quantities are primary**: cumulative share over the top
  *k*, and the participation ratio. Both are unchanged by rotation inside a
  degenerate block.
- per-layer scores are **secondary**, read with `stable1` (how often the leading
  layer really is the leading one) and `rot_diag` beside them.

`stable1` replaced an earlier `order_stable` that watched the top three layers
together. That was defective, and the sweep shows why
([§4.5](#45-one-of-the-instruments-was-defective)).

---

## 3. What a recovery number can and cannot mean

### 3.1 Compare subspaces, not columns

When the truth has more than one layer it is identified only **up to rotation**
within its own factor space — two different-looking pairs of score vectors can
describe the same surface. So the honest comparison is between *subspaces*, via
**canonical correlations** (the cosines of the principal angles): 1 for
identical spans, 0 for perpendicular ones. `subspace_cors()` does this. A
column-by-column correlation would understate recovery for reasons that have
nothing to do with the model.

A second reason: `sim_generate.R:278-290` rescales `I` but not the stored
`U`/`Lambda`, so `U %*% t(Lambda) ≠ I` and recovery is identified only up to
scale as well. `truth$U_oat`, `truth$Lambda_oat`, `truth$U_pea`,
`truth$Lambda_pea` were written and read nowhere in the repo until now.

### 3.2 The ceiling

At `kron_rank = q` the recovered scores live **by construction** in
`span(kron_A)`, a `q`-dimensional slice of the accession space, while the truth
has mass on every direction of `G`. So there is a hard analytic limit on any
recovery number, and `recovery_ceiling()` computes it **from the two spans with
nothing fitted**.

Recovery has to be quoted as a fraction of it. Otherwise the apparent `n_acc`
effect is partly 30-of-200 versus 30-of-400 directions rather than anything
about the model.

### 3.3 The floor, which is not small

The ceiling says how well the basis *could* do. `recovery_null()` says how well
it does **by accident** — and the gap between floor and ceiling is the only
range in which a recovery number carries information.

The floor is high here because the truth is itself kinship-structured, and so is
every direction in `span(kron_A)`. Measured on the **real 508-oat GRM at rank
30**: a heritable trait correlates **0.135** with a random direction in the span
on average, and **0.345** at the 95th percentile. Any vector assembled from
leading kinship eigenvectors correlates with any heritable quantity without
knowing a thing about it.

This is the number that governs whether a biological finding is real
([§5.3](#53-the-constraint-is-population-structure-not-sample-size)).

---

## 4. What the sweep measured

390 cells — 120 scenarios plus 10 null scenarios × 3 replicates — perfectly
balanced, 12 per `n_acc × sparsity × n_factors`, `kron_rank = 30`, 240 posterior
draws each. Results in `output/simulation_decomp/`. One validity check first:
`score_ceil1` for `dge_ige` reads 0.842/0.837/0.839/0.837/0.835 across the five
densities — flat, as a basis-only quantity must be.

### 4.1 The spectrum carries real rank information, and the nulls were essential

Participation ratio, draw-wise; truth is 1/3/5; the null cells have no
interaction at all:

| n_acc | observed | f1 | f3 | f5 | **null** |
|---|---|---|---|---|---|
| 200 | 1.6% | 8.73 | 8.69 | 8.52 | **8.68** |
| 200 | 9% | 8.21 | 8.19 | 8.47 | **9.59** |
| 200 | 48% | 3.86 | 5.40 | 6.49 | **10.5** |
| 400 | 1.6% | 7.34 | 8.23 | 8.36 | **9.18** |
| 400 | 9% | 4.14 | 6.88 | 7.32 | **10.4** |
| **400** | **48%** | **1.81** | **3.86** | **4.91** | **11.7** |

At n = 400 and 48% observed the participation ratio reads **1.81 / 3.86 / 4.91
against a truth of 1 / 3 / 5** — nearly unbiased. The slope on `n_factors` is
0.411 (SE 0.043, p < 2e-16).

**At 1.6% observed and n = 200 it is indistinguishable from the null** (8.73 vs
8.68). At n = 400 there is marginal separation for a rank-1 truth (7.34 vs
9.18).

Note the floor **rises** with density, 8.7 → 11.7: more data means more
well-estimated noise directions and a flatter spectrum. So there is no single
global floor, and comparing a cell against anything but its own
`n_acc × sparsity` null would be wrong. Without the null cells this table reads
as "the spectrum is uninformative everywhere", which is false.

### 4.2 Score recovery: the ceiling inverts the comparison

| observed | `dge_ige` raw | MegaLMM raw | `dge_ige` **of ceiling** | MegaLMM of ceiling |
|---|---|---|---|---|
| 1.6% | 0.658 | 0.145 | **0.778** | 0.145 |
| 4.8% | 0.724 | 0.369 | **0.862** | 0.369 |
| 9% | 0.752 | 0.748 | **0.896** | 0.748 |
| 16% | 0.789 | 0.902 | 0.942 | 0.902 |
| 48% | 0.814 | 0.971 | 0.975 | 0.971 |

`dge_ige` recovers 78–98% of what its truncated basis permits, **at every
density**. Its ceiling is ~0.84 and immovable at `kron_rank = 30`, which is why
MegaLMM overtakes it in absolute terms above about 9% observed. Raising the rank
is what lifts the cap; `sim_kron_study.R` measured that gain at 0.04–0.22 per
cell.

**The low-density number is signal, checked explicitly.** A random direction in
`span(kron_A)` scores 0.146 (n = 200) / 0.167 (n = 400) against the truth, 95th
percentile 0.34/0.41. Observed `dge_ige` recovery is 0.59–0.72. Note that
MegaLMM's 0.145 at 1.6% is *exactly* the random-direction level — it recovers
nothing there.

### 4.3 Whole-subspace recovery, and the orientation question

Mean canonical correlation over the top `n_factors` directions:

| n_factors | `dge_ige` scores | `dge_ige` loadings | MegaLMM scores | MegaLMM loadings |
|---|---|---|---|---|
| 1 | 0.678 | 0.644 | 0.648 | 0.656 |
| 3 | 0.609 | 0.570 | 0.515 | 0.521 |
| 5 | 0.576 | 0.543 | 0.417 | 0.427 |

`dge_ige` recovers the *whole* factor subspace better than MegaLMM in raw terms,
and the gap widens with rank — MegaLMM falls to 0.417 at rank 5 while `dge_ige`
holds 0.576.

And **both models recover the pea loadings about as well as the oat scores.**
That bears on [BOTH_ORIENTATIONS.md](../BOTH_ORIENTATIONS.md): the structural
asymmetry is real for main effects, but it does not appear in interaction
recovery for either model.

### 4.4 The cell closest to the real experiment

The real B4I data is **1.1% observed** — 2,063 combinations of 186,966 possible,
442 oats × 423 peas, 90.7% of combinations unreplicated. That is sparser than
the sparsest cell in the design. The nearest cell is n = 400, 1.6% observed:

| model | n_factors | raw `cor1` | ceiling | of ceiling | `share1` | participation |
|---|---|---|---|---|---|---|
| `dge_ige` | 1 | **0.673** | 0.794 | 0.846 | 0.301 | 7.34 |
| `dge_ige` | 3 | **0.701** | 0.846 | 0.826 | 0.245 | 8.23 |
| `dge_ige` | 5 | **0.773** | 0.872 | 0.887 | 0.236 | 8.36 |
| `megalmm` | 1 | 0.055 | 1.000 | 0.055 | — | — |
| `megalmm` | 3 | 0.122 | 1.000 | 0.122 | — | — |
| `megalmm` | 5 | 0.178 | 1.000 | 0.178 | — | — |
| *null* | — | — | — | — | *0.183* | *9.18* |

At the density B4I actually has, **`dge_ige` recovers the leading compatibility
axis at ρ ≈ 0.67–0.77 and MegaLMM recovers nothing** (0.055–0.178, at or below
the random-direction floor of 0.167).

### 4.5 One of the instruments was defective

`order_stable` at `n_factors = 1` sat flat at 0.16–0.23 across a 30-fold change
in density, while at `n_factors = 3` it rose 0.18 → 0.97 and at 5, 0.19 → 0.94.

The cause was `k_stable = 3L`: with a rank-1 truth, layers 2 and 3 are noise *by
construction* and permute freely between draws, so a top-three statistic reports
on them rather than on the layer anyone would interpret. A diagnostic that
cannot see density when the truth is rank-1, but can when it is rank-3, is
measuring the trailing layers.

`component_summary()` now returns `comp_stable` per layer and `stable1` for the
leading one. **On the rank-1 positive control the two read 1.00 and 0.09
respectively** — the leading layer was stable in every draw while the top-three
statistic called it unstable nine times out of ten. So the low `order_stable`
values in the sweep are very unlikely to mean the leading layer was unstable,
and the honest position is that **the sweep did not measure leading-layer
stability at all**; `stable1` has to come from a re-run. `rot_diag` was
unaffected and reads 0.72–1.00.

One cell was re-run under the fix, the B4I-relevant one — n = 200, 1.6%
observed, rank-1 truth, one replicate — and gives `stable1` = **0.62** against
`order_stable` = 0.06. So at the real density the leading layer is the leading
one in roughly two draws out of three: recoverable, but its *identity* as the
top layer is not certain draw to draw. That is a caveat on naming it, not on
using it — the direction itself recovers at 0.75 of ceiling in the same cell.

The same control now also prints the recovery floor beside the ceiling: 0.878
against a 0.901 ceiling and a 0.359 floor, so **0.96 of the usable
floor-to-ceiling range** rather than the less meaningful 0.98 of ceiling.

---

## 5. What any of this means biologically

The sections above are matrix algebra. Here is what each result is *for*.

### 5.1 The spectrum is an estimate of how many mechanisms there are

`I[i, j]` is specific combining ability — the part of a pairing's performance
that neither partner's own average explains. The decomposition says that surface
is a sum of layers, and **each layer is a candidate mechanism**:

> Layer 1 says: oats differ along one axis (`u_1`), peas differ along one axis
> (`v_1`), and a pairing does well when both are high (or both low). One axis is
> one story — a single trait mismatch or complementarity that you could name,
> measure and select on. Phenological synchrony, canopy height complementarity,
> N-fixation timing against N demand timing, rooting depth stratification.

So the **participation ratio is an estimate of the effective number of distinct
compatibility mechanisms**. That is a real biological question, not a numerical
one: is intercrop compatibility one thing with one explanation, or five things
at once?

Reading it requires the null, and the null is why §4.1 matters. A participation
ratio at the floor does **not** mean "many mechanisms" — it means *no
concentration is detectable*, and many weak mechanisms cannot be distinguished
from noise. What the sweep establishes is that the question is answerable **given
enough data**: at n = 400 and 48% observed, the ratio recovers the true number
of mechanisms to within about one.

**For B4I specifically, it is not answerable.** At 1.1% observed the spectrum
sits at or barely below its null. The programme will not be able to say how many
compatibility mechanisms there are. It can say where accessions sit on the
leading one — which is a different and more useful thing.

### 5.2 Turning a score into a mechanism

A recovered score is a column of numbers, one per oat. By itself it explains
nothing. The operational chain from score to biology is:

1. **Get the score and the loading.** `scores[, 1]` gives each oat a position on
   the leading compatibility axis; `loadings[, 1]` gives each pea its
   sensitivity. The prediction for an untried pairing is
   `d_1 × score[i] × loading[j]`.
2. **Correlate the score against traits you have measured.** Oat heading date,
   height, tillering, early vigour, straw strength. And the loading against pea
   traits: nodulation, growth habit, climbing ability, maturity. The hit — if
   there is one — *names the axis*.
3. **Read the pair together.** If `scores[,1]` tracks oat heading date and
   `loadings[,1]` tracks pea maturity, the mechanism is phenological synchrony,
   and it comes with a prescription a breeder can act on: pair early with early.
   If the score tracks oat height and the loading tracks pea climbing habit, it
   is structural support. The two sides constrain each other, which is what makes
   the pairing interpretable rather than just predictive.
4. **Allow for attenuation.** A score recovered at correlation ρ with the truth
   drags every downstream correlation toward zero by roughly that factor: an
   observed trait correlation ≈ ρ × the true one. At the B4I-like cell
   (ρ ≈ 0.67–0.77), a mechanism whose true correlation with the axis is 0.6 shows
   up at about 0.40–0.46.
5. **Test it against the kinship null, not against zero.** See §5.3 — this is the
   step that decides whether a finding is real.

**And the scores extend to accessions that were never grown.** Because a score is
a combination of kinship eigenvectors, it is defined for any genotyped
accession — so the model predicts which *untested* oat sits high on the
compatibility axis. That is precisely the property the validation experiment
needs, and it is why the truncated basis, for all its costs, is the right
construction.

### 5.3 The constraint is population structure, not sample size

This is the most important practical result in the document, and it is not
obvious from the algebra.

Suppose you have the recovered score and a measured trait, and you want to know
whether the correlation between them is evidence of a mechanism. The score is
defined for every genotyped accession, not just the phenotyped ones, so the
relevant *n* is the 508 oats in the GRM rather than the 442 with yield data.
Two different thresholds apply:

| | threshold on the TRUE trait correlation |
|---|---|
| Pure statistical power (n = 508, ρ = 0.72, 80% power) | **0.17** |
| Clearing the kinship null (95th percentile = 0.345) | **0.48** |

The second is nearly **three times** the first. The reason is §3.3: a heritable
trait correlates 0.135 on average, and up to 0.345, with a *random* direction in
`span(kron_A)`, because both are built out of the same population structure. So
a trait correlation of 0.3 with the recovered score is **not** evidence of a
mechanism, even though it is overwhelmingly significant against zero at n = 508.

**Practical consequences:**

- Test every trait association against `recovery_null()`'s distribution — random
  directions in the fitted span — rather than against zero. The full null
  distribution is far more powerful than the 95th-percentile threshold quoted
  above, so this is not as costly as the table suggests; but the naive *p*-value
  is simply wrong.
- Mechanisms with true correlation below ~0.5 to the compatibility axis are
  confounded with population structure at this panel size and cannot be
  separated. More accessions help only slowly, because the floor is set by the
  dimension of the span, not by *n*.
- The fix that would help most is not more accessions but a **less
  structure-aligned basis** — see §5.4.

### 5.4 What `kron_rank` means biologically

Rank 30 retains 62.5% of the real oat GRM's trace, and the directions kept are
the 30 largest axes of relatedness. So a recovered score is necessarily a
**smooth function of pedigree**: it varies between breeding groups and families
more readily than within them.

That is a biological limitation, not just an algebraic one:

> **The decomposition can only find compatibility mechanisms whose genetic
> architecture aligns with major population structure.** A mechanism driven by a
> single gene whose allele distribution cuts across pedigree groups — present in
> some northern and some southern material alike — lies mostly outside
> `span(kron_A)` and is invisible at rank 30.

Raising `kron_rank` makes finer, less structure-aligned mechanisms findable, and
simultaneously raises the ceiling on recovery and lowers the floor set by
structure. It costs `(rank/30)²` in basis columns. That is now the single most
consequential knob in the analysis, and it is why `sim_kron_study.R` matters
beyond settling the accuracy comparison.

### 5.5 Prediction arrives before explanation

The two headline results come apart, and the dissociation is itself the finding:

| | at B4I's density (~1.1% observed) |
|---|---|
| Can you locate the leading compatibility axis? | **Yes** — ρ ≈ 0.67–0.77, far above the 0.167 floor |
| Can you say how many axes there are? | **No** — the spectrum is at its null |

So the programme can rank oats and peas on a compatibility axis, predict untried
pairings, and select on it, **before** it can claim "there is one mechanism and
it is phenology". Breeding value first; biological explanation later, and only
with either much denser data or a higher `kron_rank`.

Stating it the other way round would be the error to avoid: a low participation
ratio at this density is not evidence of a single clean mechanism, because the
spectrum cannot tell that from noise yet.

---

## 6. Which model to use

**The short answer: for this programme, `dge_ige` for everything — and the
reasoning that got you to MegaLMM is right, it just points somewhere else.**

The evidence, by purpose:

| purpose | winner | where it holds |
|---|---|---|
| GMA / additive surface | `dge_ige` | 95–99% of cells (finding 2) |
| Interaction *prediction* (`r_int`) | `dge_ige` below ~9–16% observed; MegaLMM above. At **full rank** `dge_ige` also wins the rank-5 cells | `simulation_kron_study.csv` |
| Interaction *decomposition*, leading axis | `dge_ige` below ~9%; MegaLMM above | §4.2 |
| Interaction *decomposition*, whole subspace | **`dge_ige` everywhere**, widening with rank | §4.3 |
| At B4I's actual 1.1% density | **`dge_ige`, overwhelmingly** — 0.67–0.77 against MegaLMM's 0.055–0.178 | §4.4 |

So the claim "`dge_ige` is best for all purposes" is **correct for B4I**, and
would be too strong as a general statement: MegaLMM genuinely overtakes it on
absolute interaction recovery above roughly 9–16% observed. B4I is at 1.1%, well
below every density tested, and will stay sparse as long as the programme samples
broadly rather than deeply.

### Why the reduced-rank intuition was right but misdirected

Your reasoning was: the interaction is probably low-rank, MegaLMM is a
reduced-rank model, therefore use MegaLMM. The premise is right — the simulation
confirms the interaction *is* low-rank and that this is recoverable. What is
wrong is the inference, and the reason is instructive:

**MegaLMM imposes the low rank as a constraint during fitting. `dge_ige` fits the
interaction in a kinship kernel and the low rank is extracted afterwards.**

Those are not the same thing when data are sparse:

- MegaLMM's `K` factors have to be *estimated*, and estimating a factor needs
  combinations that share accessions — repeated structure. At 1.1% observed there
  is almost none, so the factors shrink to nothing and the model falls back on
  the column mean. That is why its recovery sits at the random-direction floor.
- `dge_ige`'s Kronecker kernel never estimates a factor. It borrows strength
  through kinship: an unobserved pairing is informed by every pairing involving
  *related* oats and *related* peas. That works at any density, which is why its
  recovery degrades so gently — 0.814 at 48% observed down to only 0.658 at 1.6%.

The general lesson, worth carrying to other problems: **estimate in the space
where the model is well-posed, then reduce; do not constrain the estimation to a
rank you guessed.** The reduction is free — it is a rotation (§2.3) — so there is
no cost to deferring it, and a large cost to imposing it early.

The one thing MegaLMM still has that `dge_ige` does not is an **untruncated**
interaction: its ceiling is 1.0 where `dge_ige`'s is 0.84. Closing that gap is
`kron_rank`, not a change of model.

---

## 7. The BGLR trap, in case this is ever re-plumbed

`BGLR::Multitrait` honours `ETA[[j]]$saveEffects = TRUE` and streams every draw
of that term's coefficients to `<saveAt><term>_beta.bin`. Three things about it:

1. **`saveAt` is a filename prefix, not a directory.** Pointing two array tasks
   at the same prefix makes them overwrite each other.
2. **The file contains the burn-in.** The write sits inside
   `if (iter %% thin == 0)` while the posterior mean accumulates under the
   separate gate `(iter > burnIn) & (iter %% thin == 0)`, and the header records
   `nRow = nIter/thin`. `read_beta_draws()` drops the first
   `floor(burnIn/thin)` rows. Keeping them moved the mean by 0.0994 against a
   3e-16 match when dropped — and an inflated posterior SD with a flattened
   spectrum is *exactly* what the prior-bias analysis is looking for, so this
   failure would have been read as a result.
3. **A truncated file reads back as NA without erroring**, which is what a killed
   job leaves. Hence the file-size assertion.

`tests/test_fits.R` pins all three with one assertion:
`colMeans(draws) == fit$ETA$G_mix$beta`. BGLR's returned `beta` *is* the running
posterior mean over exactly the post-burn-in thinned draws, so that equality can
only hold if the offset, the byte order, the trait-major layout and the storage
mode are all right.

Separately: MegaLMM's own factors are **not** orthogonal — per-factor shares of
`interaction_part(U_F %*% Lambda)` do not sum to one, asserted in
`tests/test_decomp.R` — so both models' surfaces go through the same
`decompose_surface()` rather than comparing a MegaLMM factor to a `dge_ige`
layer.

---

## 8. Running it

```
Rscript code/sim_decomp_run.R --check     # positive control; run first
Rscript code/sim_decomp_run.R --pilot     # two cheap cells, end to end
Rscript code/sim_decomp_run.R             # the design
sbatch -A <account> code/scinet/sim_decomp_array.sbatch
Rscript code/sim_decomp_run.R --combine   # after the array: one complete CSV
```

600 cells plus the nulls at 5 replicates, about 120 core-hours, roughly 6 h per
task over 20 tasks with a 24 h wall clock. Cells are cached per
scenario-replicate and the job is resumable. Outputs are
`output/simulation_decomp_results.csv` (one row per cell × trait × model) and
`output/simulation_decomp_spectrum.csv` (one row per component).

**`--combine` is the step it is easiest to forget.** Every array task writes both
CSVs from whatever was in the cache when *it* finished, so the files left behind
by a 20-task job are each a partial view. Running `--combine` on the login node
fits nothing, globs the whole cache and rewrites the two CSVs complete,
reprinting the summary tables. Do it before reading anything.

### Every flag

| flag | default | what it does |
|---|---|---|
| `--check` | — | The positive control: one dense rank-1 cell at `n_acc = 120`, 48% observed, short chain. Prints the diagnostics and `stop()`s if the draws do not average to BGLR's posterior mean, if the leading layer carries under 0.3, if the recovered score misses half its ceiling, or if MegaLMM's shares do not sum to 1. Run it first; a negative result from the sweep is only worth having if the pipeline can produce a positive one. |
| `--pilot` | — | Two cheap cells (`rep 1`, `n_acc = 200`, 48% observed), end to end. They are genuine design cells with design seeds, so the full run reuses their cache. |
| `--combine` | — | Rebuild both CSVs from the cache and reprint the tables, fitting nothing. See above. |
| `--refresh` | off | Ignore the cache and refit. Needed after any change to what is scored or stored — otherwise move `DECOMP_SCHEME`. |
| `--reps N` | `SIM_INT_REPS` (5) | Replicates per scenario. Replicate is the slow seed index, so raising it is additive: existing cells stay cached rather than being renumbered. |
| `--rank N` | `SIM_KRON_RANK` (30) | `kron_rank`, the truncation axis — and the flag most worth varying, for the reasons in §5.4. It sets the ceiling on recovery, the structural floor, and costs `(N/30)^2` in basis columns. |
| `--n-iter N` | 3000 | MCMC iterations. |
| `--burn-in N` | 600 | Burn-in. `read_beta_draws()` uses it to strip the leading rows BGLR streams but does not average. |
| `--thin N` | 10 | Thinning. Controls the fit and the reader together — they must agree, or the file-size assertion fires. Kept draws are `nIter/N − burnIn/N`. |
| `--filter EXPR` | — | An R expression over the design, e.g. `--filter "sparsity <= 0.09"`. |
| `--task N --ntasks M` | — | Array slicing: task `N` of `M`, taking every `M`th cell. Supplied by the sbatch wrapper. |
| `--no-null-cells` | nulls on | Drop the `interaction_pct = 0` scenarios. **Do not.** §4.1 is unreadable without them: they are the floor the participation ratio is measured against, and the floor moves with density. |
| `--no-megalmm` | on | Skip the MegaLMM pair. Saves about 32 of the 120 core-hours and is a reasonable first stage, but it drops the head-to-head. |
| `--drop-draws` | keep | Delete each cell's `.bin` once summarised. Saves roughly 1 GB over the sweep; the cost is that revisiting any decomposition choice — a different rank, another alignment scheme, `stable1` for cells fitted before it existed — then needs a refit rather than seconds. |

---

## 9. Scope: simulation only

`code/BGLR_multi_trait_model.R:94` sets `fit_mix_term <- FALSE`, and when TRUE it
builds the **exact** combination kernel inline and calls `BGLR::Multitrait`
directly rather than going through `fit_producer_associate`. So there is no
fitted interaction surface for the real experiment, the exact kernel is not
bilinear and this decomposition does not apply to it (§1.3), and `save_effects`
refuses outright when `kron_rank` is `NA` rather than silently producing
coefficients that do not reshape. The guard at `BGLR_multi_trait_model.R:233`
notes that over 90% of real combinations are unreplicated — confirmed above at
90.7% — so the exact term is near-inestimable there anyway.

Making this reach the real data means switching the production script to
`fit_producer_associate(kron_rank = ...)` with `save_effects = TRUE`. On the
evidence in §4.4 that is now worth doing: at the real density `dge_ige` recovers
the leading compatibility axis at ρ ≈ 0.67–0.77, which is enough to support the
trait-correlation exercise in §5.2 for any mechanism with a true correlation
above about 0.5.
