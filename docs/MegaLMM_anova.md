# What moves MegaLMM's accuracy — an ANOVA of the simulation

Which of the three axes — what was simulated, how much data there was, and how
MegaLMM was configured — actually change its accuracy, and by how much.

Produced by `code/sim_anova.R`. For what the response variables mean see
[SIMULATION_GLOSSARY.md](../SIMULATION_GLOSSARY.md); for the design see
[SIMULATION.md](../SIMULATION.md).

> **These results are from the September 2026 Ceres run, which used the
> pre-October scoring scheme**: 20% of the observed cells held out, so the models
> were fitted at 0.8 × the labelled sparsity and the per-cell metrics were scored
> on the held-out cells rather than on the never-observed ones. Read every
> sparsity level as 0.8 × its label — the density gate quoted as "between 4.8%
> and 16%" is between **3.84% and 12.8%** actually fitted. The scheme has since
> changed (see [SIMULATION_GLOSSARY.md](../SIMULATION_GLOSSARY.md)); these
> analyses have not been rerun, because the conclusions they support do not turn
> on it.

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

## What partial ω² is, and why not η² or p

Three quantities get confused here, so they are worth separating.

**η² (eta squared)** is the obvious thing: the share of total variation a term
accounts for, `SS_effect / SS_total`. It has two problems. It depends on what
else is in the model — add a big factor and every other term's η² shrinks,
without anything about those terms having changed. And it is **biased upward**,
for a reason that matters here.

**The bias.** A term that does nothing at all still collects sum of squares,
because it soaks up noise along its degrees of freedom. In expectation a null
term with `df` degrees of freedom has

```
E(SS_effect) = df x MS_error
```

So η² is never zero even for a factor with no effect, and the more degrees of
freedom you give it the larger it looks. With 19 terms in this table, several of
them with 3 or 4 df, that is not a rounding problem.

**ω² (omega squared)** subtracts the noise a null term would have collected:

```
                    SS_effect - df_effect x MS_error
omega2_p  =  ----------------------------------------------------
             SS_effect - df_effect x MS_error + (df_error + 1) x MS_error
```

The **partial** in "partial ω²" is the denominator: it is this term plus the
error, *not* the total. So it answers "of the variation this term and the noise
between them could explain, how much is the term?" — which is comparable across
terms and across the four responses, and does not move when an unrelated factor
is added. The consequence is that **partial ω² values do not sum to 1**, and
should not be expected to.

**Why it can be negative.** If a term collects *less* sum of squares than its
degrees of freedom would collect by chance, the numerator goes negative. That is
not a defect, it is the informative case: it says the term is not merely small
but indistinguishable from — indeed below — what a null factor would produce.
η² cannot do this. It is bounded at zero, so a completely inert factor still
prints a small positive number and looks like *something*.

### Worked on the real numbers

`r_addsurf_oat`, whole-plot stratum. `MS_error = 0.01105` on 756 df, so a null term
with 1 df is expected to collect `SS ≈ 0.011` just from noise:

| term | df | SS | SS a null term would get | partial η² | partial ω² |
|---|---|---|---|---|---|
| `sparsity` | 3 | 233.4 | 0.033 | 0.965 | **0.965** |
| `fixed_main_effect` | 1 | 1.313 | 0.011 | 0.136 | **0.135** |
| `K` | 1 | 0.554 | 0.011 | 0.062 | **0.061** |
| `eigen_variance` | 1 | 0.0003 | 0.011 | 0.000 | **−0.001** |

For `sparsity` the correction is irrelevant — 233 against 0.03 of expected noise.
For `eigen_variance` it is the whole story: its sum of squares is **thirty times
smaller than a factor with no effect would be expected to produce**. Partial η²
rounds that to 0.000, which reads as "very small". Partial ω² makes it negative,
which reads correctly as "this lever does nothing, and we have enough data to say
so".

### Why Δr sits next to it

Partial ω² is a **share of variance, on the Fisher-z scale**. It says how much of
the explainable variation a term accounts for, and nothing about how big the
effect is in units anyone cares about. A term can own a large share of a small
amount of variation.

Δr is the other half: the spread the term produces **in correlation points**.
Read them together —

- **high ω², high Δr** — `sparsity` (0.965, 0.76). Dominant and consequential.
- **moderate ω², small Δr** — `K` (0.061, 0.027). Real and reproducible, worth
  2–3 correlation points. Statistically solid, practically marginal.
- **low ω², moderate Δr** would mean a large effect measured imprecisely. Nothing
  in this table is in that cell, which is itself worth knowing.
- **negative ω²** — `eigen_variance`. Nothing there.

Neither is a p-value, and deliberately so: with 800 runs, `K` sits at
p = 3 × 10⁻¹² and moves `r` by 0.027. The p-value tells you the effect is real.
It does not tell you it matters.

**A caution on rules of thumb.** The conventional 0.01/0.06/0.14 =
small/medium/large bands come from psychology, where effects are weak and noise
is large. They are useless here: a designed simulation with 800 runs and a factor
that spans a thirtyfold range in data volume puts `sparsity` at 0.97, which no
band accommodates. Compare terms against **each other** in this table, not
against an external convention.

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
| `r_addsurf_oat` | 0.0111 | 0.0334 | 0.0027 | 0.0167 | 14% |
| `r_addsurf_pea` | 0.0099 | 0.0294 | 0.0025 | 0.0147 | 15% |
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

| term | df | `r_addsurf_oat` | `r_addsurf_pea` | `r_int_oat` | `r_int_pea` |
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
   `r_addsurf` (0.004–0.018). `environment` is the mirror: it matters for GMA (0.48)
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

Marginal means on the r scale, `r_addsurf_oat`:

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

| lever | level | `r_addsurf_oat` | `r_int_oat` |
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

Subplot stratum. The readout main effect is worth Δr = 0.30 on `r_addsurf_oat` and
0.21 on `r_addsurf_pea`, but only 0.03 on `r_int`. Eta is much better for general
mixing ability and no better for the interaction, exactly as the parameterisation
predicts: `U` omits the per-column intercept where the column species' associate
effect lives, and GMA contains that effect while the interaction part does not.

`readout × sparsity` is the largest subplot term (ω² = 0.886 on `r_addsurf_oat`,
Δr = 0.66): the Eta-over-U advantage is itself strongly sparsity-dependent.

Mean `r_addsurf` by readout and sparsity (observed means, 192–208 runs per cell;
written by `code/sim_anova.R` to `output/simulation_anova_eta_vs_u.csv`):

| observed | Eta oat | U oat | **gap oat** | Eta pea | U pea | **gap pea** |
|---|---|---|---|---|---|---|
| 1.6% | 0.173 | 0.110 | **+0.063** | 0.161 | 0.184 | **−0.023** |
| 4.8% | 0.547 | 0.439 | **+0.108** | 0.505 | 0.515 | **−0.010** |
| 16% | 0.912 | 0.676 | **+0.236** | 0.900 | 0.726 | **+0.174** |
| 48% | 0.969 | 0.706 | **+0.263** | 0.964 | 0.764 | **+0.200** |

Two things to take from this.

**The gap widens with data, it does not close.** One might expect `U` to catch up
once there is enough information — it does the opposite. From 1.6% to 48% Eta
gains 0.80 while U gains only 0.60, because the per-column intercept that `U`
discards is precisely the part that more data estimates best. `U` is not a noisy
version of Eta; it is missing a component, and the missing component becomes more
recoverable, not less, as the matrix fills.

**At the sparse end for pea, U is very slightly ahead** (−0.023 and −0.010).
With 1.6–4.8% observed the per-column intercept is estimated from almost nothing,
so including it adds more noise than signal, and dropping it is marginally
better. The crossover is somewhere between 4.8% and 16% observed. This is small
and worth not over-reading, but it is the same story from the other side: the
intercept is an unshrunk fixed effect, so it helps exactly when it can be
estimated and hurts when it cannot.

Either way, **report Eta**. The one regime where `U` competes is the regime where
neither is usable.

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

The cell means make the size of that difference concrete. Mean `r_int` by
architecture and size (312–328 runs per cell; written to
`output/simulation_anova_arch_size.csv`):

| | size 10% | size 20% | **mean** |
|---|---|---|---|
| **1 factor**, oat | 0.401 | 0.481 | **0.441** |
| **5 factors**, oat | 0.258 | 0.320 | **0.289** |
| **1 factor**, pea | 0.387 | 0.472 | **0.430** |
| **5 factors**, pea | 0.246 | 0.307 | **0.276** |
| **mean**, oat | 0.330 | 0.400 | |
| **mean**, pea | 0.316 | 0.389 | |

Going from rank 1 to rank 5 costs about **0.15** in `r_int`; doubling the
interaction's share of variance from 10% to 20% buys about **0.07**. So rank
costs roughly twice what variance share buys, and the two are additive — the
absence of a size × architecture term means the penalty for rank 5 is the same
whether the interaction is large or small (0.143 at 10%, 0.161 at 20% for oat).

This is the more consequential finding of the two, because **the variance share
is something a breeding programme can partly influence and the rank is not**. A
rank-5 interaction of 20% (r = 0.320) is harder to recover than a rank-1
interaction of only 10% (r = 0.401): spreading the same specific-combination
variance across more independent directions hurts more than halving it. If real
oat × pea specific combining ability is high-rank, it will be hard to predict
however much of it there is.

---

## What was left out, and what it cost

The plan excluded **generative × data** interactions, on the grounds that they
are a question about the simulation rather than about MegaLMM. That is
defensible, but it is not free: those terms are real and large (F = 16.6 on
`r_addsurf_oat`), so their sum of squares goes into the whole-plot error and inflates
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
