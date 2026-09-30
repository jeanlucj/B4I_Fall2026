# What moves MegaLMM's accuracy — an ANOVA of the simulation

Which of the three axes — what was simulated, how much data there was, and how
MegaLMM was configured — actually change its accuracy, and by how much.

Produced by `code/sim_anova.R`. For what the response variables mean see
[SIMULATION_GLOSSARY.md](../SIMULATION_GLOSSARY.md); for the design see
[SIMULATION.md](../SIMULATION.md).

Data: `megalmm` and `megalmm_U` rows only, **replicates 2–5** (800 runs, 1,600
rows). Replicate 1 is excluded because it holds two sets of cache files under
different seeds, both calling themselves rep 1 — see
`code/sim_filter_results.R`.

---

## Three things the analysis has to get right

### 1. Fisher z, not r

The responses run from about 0 to 0.95. `r` is bounded and heteroscedastic — its
sampling variance collapses as |r| → 1 — so a plain ANOVA gives the high cells
too much weight and the residuals fan. `atanh(r)` has variance ≈ 1/(n−3)
whatever `r` is. Everything is fitted on z and reported back on the r scale.

### 2. Eta and U are not two runs — it is a split plot

`megalmm` and `megalmm_U` are two **readouts of the same fit**: `Eta_mean` the
predicted phenotype, `U` the genetic value without the per-column intercept. They
are paired. Treating readout as an ordinary crossed factor would double the
apparent replication and shrink every standard error in the table.

The run is the whole plot, the readout is the subplot, and the two strata come
out exactly by analysing per-run means and differences:

```
mean_z = (z_Eta + z_U)/2    whole-plot stratum: everything except readout
diff_z =  z_Eta - z_U       subplot stratum: the INTERCEPT is the readout main
                            effect, and a main effect of A here IS readout × A
```

Each run contributes one mean and one difference, so this stays exact under the
unbalanced D-optimal fraction — which `aov(... + Error(run))` would not.

### 3. Effect size first, because everything is significant

With 800 runs almost every term clears p < 0.05. The table therefore leads with

- **partial ω²** — bias-corrected variance explained, and negative when a term
  explains less than its degrees of freedom would by chance;
- **Δr** — the spread the term produces **on the r scale**. For a main effect
  that is the range of its marginal means; for an interaction it is the range of
  the cell means **with the additive part removed**, the same
  `interaction_part()` split the rest of the project uses. Without that
  correction a weak interaction prints a larger Δr than the strong main effect
  inside it.

---

## Expected mean squares

For run *i*, readout *j*:

```
z_ij = mu + (whole-plot fixed effects)_i + (readout)_j
           + (readout x whole-plot)_ij + u_i + e_ij

u_i  ~ N(0, s2_run)   one MegaLMM fit versus another under the same settings
e_ij ~ N(0, s2_sub)   readout versus readout within one fit
```

With a fixed factor A of *a* levels, *n* runs per level and *N* runs:

| stratum | source | df | E(MS) |
|---|---|---|---|
| **whole plot** (on `mean_z`) | A | a−1 | σ²_run + σ²_sub/2 + n/(a−1) · Σα² |
| | whole-plot error | N−p | σ²_run + σ²_sub/2 |
| **subplot** (on `diff_z`) | readout | 1 | 2σ²_sub + N·δ² |
| | readout × A | a−1 | 2σ²_sub + 2n/(a−1) · Σ(δα)² |
| | subplot error | N−p | 2σ²_sub |

Each stratum is tested against its own residual, which is what the two `lm()`
fits do. The components follow from the two error mean squares:

```
s2_sub = MS_subplot_error / 2
s2_run = MS_wholeplot_error - MS_subplot_error / 4
```

**These coefficients are the balanced-design ones.** This design is a D-optimal
fraction and is not balanced, so the F tests use Type II sums of squares rather
than these coefficients. The EMS is here to show what is tested against what.

Fitted values:

| response | MS whole-plot err | MS subplot err | σ²_run | σ²_sub | run share of noise |
|---|---|---|---|---|---|
| `r_gma_oat` | 0.0111 | 0.0334 | 0.0027 | 0.0167 | 14% |
| `r_gma_pea` | 0.0099 | 0.0294 | 0.0025 | 0.0147 | 15% |
| `r_int_oat` | 0.0432 | 0.0033 | 0.0424 | 0.0016 | **96%** |
| `r_int_pea` | 0.0444 | 0.0040 | 0.0434 | 0.0020 | **96%** |

That contrast is itself a result. For **GMA** most of the noise is
readout-to-readout: Eta and U disagree sharply, because U structurally cannot
hold the associate effect that GMA contains. For the **interaction** they barely
disagree at all (96% of the noise is between runs), because the interaction lives
in the factor structure, which `U` does keep.

---

## What actually matters

Partial ω², whole-plot stratum, as planned (generative × data excluded):

| term | df | `r_gma_oat` | `r_gma_pea` | `r_int_oat` | `r_int_pea` |
|---|---|---|---|---|---|
| **sparsity** | 3 | **0.965** | **0.969** | **0.889** | **0.879** |
| **n_acc** | 1 | 0.561 | 0.611 | 0.209 | 0.222 |
| **environment** | 2 | 0.479 | 0.488 | 0.015 | 0.007 |
| **interaction** | 4 | 0.018 | 0.004 | **0.428** | **0.411** |
| sparsity × fixed_main_effect | 3 | **0.267** | 0.233 | 0.012 | 0.010 |
| fixed_main_effect | 1 | 0.135 | 0.151 | 0.022 | 0.014 |
| K | 1 | 0.061 | 0.101 | 0.002 | 0.003 |
| sparsity × K | 3 | 0.031 | 0.026 | −0.002 | −0.003 |
| eigen_variance | 1 | −0.001 | −0.001 | −0.001 | −0.001 |
| *the other 10 two-way terms* | | ≤0.001 | ≤0.015 | ≤0.006 | ≤0.008 |

Four readings:

1. **Sparsity is not one factor among several, it is the experiment.** ω² ≈ 0.97,
   and Δr = 0.76 on GMA: from r = 0.15 at 1.6% observed to 0.91 at 48%. Nothing
   else is within an order of magnitude.
2. **The generative factors separate cleanly by response, as they should.** The
   simulated interaction drives `r_int` (ω² ≈ 0.42) and is nearly irrelevant to
   `r_gma` (0.004–0.018). `environment` is the mirror: it matters for GMA (0.48)
   and not for the interaction (0.007–0.015). That is a sanity check passing —
   each response responds to the thing it is supposed to measure.
3. **Of the three analysis levers, one matters, one is marginal, one is
   inert.** `fixed_main_effect` is real (0.14–0.15 on GMA), `K` is small but
   consistent (0.06–0.10), and **`eigen_variance` does nothing at all** — partial
   ω² is *negative* for every one of the four responses, meaning it explains less
   than its single degree of freedom would by chance. 0.25 and 0.75 give 0.705
   and 0.704.
4. **Exactly one interaction is practically important**, and it is the one worth
   having: `sparsity × fixed_main_effect`.

### The one interaction that matters

Marginal means on the r scale, `r_gma_oat`:

| observed | `fixed_main_effect = FALSE` | `TRUE` | gain from pinning |
|---|---|---|---|
| 1.6% | 0.097 | **0.204** | **+0.107** |
| 4.8% | 0.423 | **0.624** | **+0.201** |
| 16% | 0.847 | 0.837 | −0.010 |
| 48% | 0.913 | 0.909 | −0.004 |

**Pin the main effect when the matrix is sparse; do not bother when it is
dense.** At 4.8% observed it doubles accuracy; by 16% it has stopped helping and
is very slightly harmful. B4I sits at about 3% observed, which is squarely in the
region where pinning is worth the most.

And it reverses for the interaction: on `r_int_oat`, pinning *costs* accuracy
(0.525 → 0.478). It buys the main effect at the price of some of the factor
structure. So the setting depends on what you want out of the model, not just on
how sparse the data is.

### The other levers, for the record

| lever | level | `r_gma_oat` | `r_int_oat` |
|---|---|---|---|
| K | 5 | 0.718 | 0.511 |
| | 10 | 0.691 | 0.493 |
| eigen_variance | 0.25 | 0.705 | 0.506 |
| | 0.75 | 0.704 | 0.498 |
| environment | one | 0.740 | 0.525 |
| | ten_stable | 0.735 | 0.503 |
| | ten_gxe | 0.628 | 0.476 |

`K = 5` beats `K = 10` slightly and consistently — more factors is not better
here. And note `environment`: ten *stable* environments cost essentially nothing
against one (0.735 vs 0.740), while ten with GxE cost 0.11. **It is the GxE that
hurts, not the number of environments** — which is the axis worth keeping.

### Eta versus U

Subplot stratum. The readout main effect is worth Δr = 0.30 on `r_gma_oat` and
0.21 on `r_gma_pea`, but only 0.03 on `r_int`. Eta is much better for general
mixing ability and no better for the interaction, exactly as the parameterisation
predicts: `U` omits the per-column intercept where the column species' associate
effect lives, and GMA contains that effect while the interaction part does not.

`readout × sparsity` is the largest subplot term (ω² = 0.886 on `r_gma_oat`,
Δr = 0.66): the Eta-over-U advantage is itself strongly sparsity-dependent.

### Size versus architecture of the simulated interaction

`interaction_pct` and `n_factors` are **nested, not crossed** — there is no "0%
with 5 factors" cell, which is why `sim_config.R` folds them into one composite
factor. Once the no-interaction level is dropped the remaining four levels *are*
a 2 × 2, and there they separate:

| term | df | `r_int_oat` ω² | `r_int_pea` ω² |
|---|---|---|---|
| architecture (1 vs 5 factors) | 1 | **0.343** | **0.325** |
| size (10% vs 20%) | 1 | 0.170 | 0.165 |
| size × architecture | 1 | 0.002 | 0.003 |

**Architecture matters about twice as much as size**, and they do not interact.
A rank-5 interaction is much harder to recover than a rank-1 one of the same
total variance — which is the more interesting finding, because rank is the thing
a real breeding programme has no control over.

---

## What was left out, and what it cost

The plan excluded **generative × data** interactions, on the grounds that they
are a question about the simulation rather than about MegaLMM. That is
defensible, but it is not free: those terms are real and large (F = 16.6 on
`r_gma_oat`), so their sum of squares goes into the whole-plot error and inflates
it by about 50% — from MS 0.0074 to 0.0111. **Every test in the table above is
therefore conservative.**

`Rscript code/sim_anova.R --all-two-way` fits them. Nothing in the conclusions
changes: the ranking of the analysis-model terms is identical, `eigen_variance`
stays inert, and `sparsity × fixed_main_effect` stays the dominant analysis
interaction (ω² rises 0.267 → 0.358). The one term that surfaces is
`environment × sparsity` (ω² = 0.319) — GxE hurts more when data is sparse, which
is worth knowing and is not a statement about MegaLMM's configuration.

`rep` was tested as a block and is negligible (F = 0.11, p = 0.74), so
replicates are exchangeable and are not modelled.

---

## Two things this cannot answer

**The fraction was built for main effects and all two-way interactions**, so the
three-way terms that would tell you whether the `sparsity × fixed_main_effect`
crossover itself moves with panel size or GxE are not estimable. That was the
right trade at 150 runs, but it is the boundary of what this design supports.

**Sparsity is fitted as a 4-level factor**, which spends 3 df on what is clearly
a monotone, saturating curve. That is right for an ANOVA whose job is to
apportion variance, but the practical statement — where the curve turns over — is
a question about shape, and `--extended` (8 sparsity levels) is the run that
answers it.
