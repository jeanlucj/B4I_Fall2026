# Why sparsity hurts MegaLMM more than a row average

For the models themselves see [BACKGROUND.md](../BACKGROUND.md); for the
simulation that produced these numbers see [SIMULATION.md](../SIMULATION.md).

## The puzzle

That a sparse matrix makes pea "environment" loadings hard to estimate is
obvious. The part that is not obvious:

> Why should that make MegaLMM worse at predicting **oats** than simply
> averaging each oat's own observations?

An oat's average uses only that oat's data. Nothing about pea loadings should
touch it. If MegaLMM cannot learn the pea structure, the worst one might expect
is that it *ignores* the pea dimension and falls back on something like a row
average.

It does not fall back on the row average. At 5% observed the row average
predicts held-out cells at r = 0.45 and MegaLMM at r = 0.10. Something
positively unhelpful is happening, and it is worth understanding because it
tells you what to fix.

## The one structural fact that explains most of it

**MegaLMM has no oat main effect.**

Its model of the matrix is

```
Y[oat i, pea j] = mu_j + (U_F Λ)[i,j] + U_R[i,j] + error
                  ^^^^
                  a per-PEA intercept. There is no per-OAT intercept.
```

The per-column intercept `mu_j` is the pea main effect, fitted directly and
cheaply. The oat main effect has no term of its own at all. To represent "this
oat is simply a better oat everywhere", MegaLMM must **discover a latent factor
whose loadings happen to be roughly constant across peas**, and estimate every
oat's score on it.

A row average *is* the oat main effect. It is handed to you for one parameter
per oat. MegaLMM has to reconstruct it out of `K` scores per oat and `K`
loadings per pea — and it can only do that if the data can tell it that the
peas rank oats alike.

That asymmetry is the whole story. Everything below is why the reconstruction
fails when the matrix is sparse, and where the model ends up instead.

## Mechanism 1: a factor model learns from overlap, and sparse designs have none

To learn that pea *j* and pea *j'* rank oats similarly — which is what a shared
factor asserts — you need oats grown with **both**. Two columns with no oat in
common carry no information about each other.

How many oats do two random pea columns share? With `n` oats and a fraction `p`
of the matrix observed, each column holds `n·p` oats, and the expected overlap
between two columns is **`n·p²`**:

| observed | obs per pea column | oats shared by two columns |
|---|---|---|
| **1.1%** (real B4I) | 2.2 | **0.02** |
| 5% | 10 | **0.5** |
| 15% | 30 | **4.5** |
| 45% | 90 | **40.5** |

That `p²` is the crux. Halving the density quarters the evidence about how
columns relate. At B4I's actual density, a randomly chosen pair of peas has a
one-in-fifty chance of sharing a single oat — so essentially no pair of pea
environments can be compared at all, and the factor structure the model exists
to estimate is simply not in the data.

At 45% every pair shares about 40 oats and the structure is easy. The crossover
we measured near 15% is where overlap first exceeds a handful.

This also explains something that looked odd: **a bigger panel does not help.**
400×400 at a given density has four times the observations, but the overlap
`n·p²` only doubles while the number of loadings to estimate also doubles. More
accessions is not more information per pea.

## Mechanism 2: there are more parameters than observations

Counting, for a 200×200 panel:

| observed | observations | parameters (K=10) | observations per parameter |
|---|---|---|---|
| 5% | 2,000 | 4,200 | **0.48** |
| 15% | 6,000 | 4,200 | 1.43 |
| 45% | 18,000 | 4,200 | 4.29 |

The parameters are 200 column intercepts + 200×K oat scores + 200×K pea
loadings. At 5% there are **twice as many parameters as data points**. A row
average, by contrast, has 200 parameters and 10 observations each.

This is also why K matters so much more than MegaLMM's shrinkage argument
suggests. At 5%, dropping K from 10 to 5 takes observations per parameter from
0.48 to 0.91 — from hopeless to merely bad — and we measured accuracy roughly
tripling as a result. The ARD prior is supposed to make surplus factors cheap,
and it does when there is data to shrink *with*; with one observation per
factor per oat there is almost nothing for it to work on.

## Mechanism 3: shrinkage sends it to the wrong margin

This is the part that answers the original question.

When the factor structure is unsupported, the ARD prior does the right thing
and shrinks `Λ` and `U_F` towards zero. But look at what the model reduces to
when it does:

```
Y[i,j] = mu_j + (something shrinking to 0) + error
       → mu_j
```

**MegaLMM's fallback is the pea column mean.** Not the oat row mean — the
*other* margin, the one with no oat information in it whatsoever. So as the
matrix gets sparser, MegaLMM's predictions do not degrade towards a row average
and stop there; they degrade towards a quantity that is uninformative about the
thing being predicted.

A row average degrades differently. With 10 observations it is noisy, but it is
still an estimate *of the right quantity*, so its accuracy sags gently — 0.45
at 5% against 0.53 at 45% in our runs. MegaLMM's accuracy collapses, because it
is not a noisier version of the row average; past a point it is a different
estimator altogether.

### The evidence

The prediction surface can be split into its row margin and its column margin
and each compared with what it should be:

| observed | MegaLMM row margin vs **training row means** | vs **true producer effect** | held-out: MegaLMM | held-out: row average |
|---|---|---|---|---|
| 5% | **0.36** | **0.29** | **0.10** | 0.45 |
| 15% | 0.95 | 0.94 | 0.87 | 0.50 |
| 45% | 0.99 | 0.98 | 0.97 | 0.53 |

At 5% MegaLMM's own row margin correlates only 0.36 with the oat averages
sitting in its input, and 0.29 with the true oat effect. It is not that the oat
effect was estimated and then diluted — it was never recovered. By 15% the same
quantity is at 0.94 and MegaLMM has overtaken the row average by a wide margin.

## Why the DGE-IGE model is immune

It makes a much stronger assumption and is rewarded for it when the data cannot
support anything weaker:

```
y[i,j] = Pr_i + As_j + error          (additive: every pea ranks oats identically)
```

- `Pr_i` is an **oat main effect with its own term**, estimated from every
  observation of oat *i* regardless of which peas those were. It never needs two
  columns to share an oat, because additivity asserts the answer that the factor
  model has to estimate.
- `Pr ~ N(0, G_oat σ²)` adds kinship borrowing on top, so an oat with ten
  observations also draws on its relatives. A row average has neither the term
  nor the borrowing, which is why DGE-IGE beats the row average too.

The cost of that rigidity is that it cannot represent interaction at all. That
is a real cost — it is why B4I's specific-combination term had to be dropped and
why the pilot, at 86% density, is the first fit in this project to include one.
But rigidity is not a defect at 1% density. There is nothing there to be
flexible about.

## This is not a criticism of MegaLMM

MegaLMM is built for multi-environment trials where a genotype is grown in a
useful fraction of the environments — 30% to 70% is typical — and where the
interesting question is precisely the genotype-by-environment structure that a
row average throws away. On the dense end of our own simulation it wins both
comparisons, r = 0.91 against DGE-IGE's 0.79 for total genetic value and 0.86
against 0.44 for the interaction.

Reading each pea accession as an environment is a legitimate and interesting
framing. It just puts the model in a regime two orders of magnitude sparser
than the one it was designed for, and asks it to reconstruct a main effect it
has no term for out of overlap that does not exist.

## What follows for the design

1. **The number to hit is overlap, not density.** Aim for `n·p²` in the tens,
   not for a percentage. At 400 peas that means roughly 15% of combinations
   observed; the arithmetic is above.
2. **Replicate specific pairs rather than spreading singletons.** This is the
   same conclusion the bivariate model reaches for `σ_PrAs`, from the other
   direction.
3. **Fewer peas, better covered, beats more peas thinly covered.** The pilot
   makes the point: 12 pea cultivars at 86% coverage supports a full
   interaction term, where B4I's 423 peas at 1.1% support none.
4. **Keep K small when the matrix is sparse.** K = 10 with 10 observations per
   column is asking for one factor per observation.
5. **Give the pea side a relationship matrix if you can.** MegaLMM's environment
   axis takes covariates but not kinship, which is a weaker channel than the
   `G_pea` the DGE-IGE model uses directly. See
   [BOTH_ORIENTATIONS.md](../BOTH_ORIENTATIONS.md).

## A fix: give the main effect a term of its own

If the diagnosis is that MegaLMM has no oat main effect, the remedy is to give
it one. MegaLMM already can: `setup_model_MegaLMM()` takes an undocumented
`Lambda_fixed` argument — present in the installed 0.9.5, used in neither
vignette — which pins the first *k* rows of `Λ` and samples only the rest.

Setting `Lambda_fixed = matrix(1, 1, p)` makes factor 1 an oat main effect:
every loading is 1, so its score `f_1` is estimated from **all** of an oat's
observations across every column, which is the row pooling a row average does.
It keeps the level-2 model `f_1 = U_F1 + E_F1` with its own `h²`, so it is a
main effect *with kinship borrowing* — close to DGE-IGE's `Pr_i ~ N(0, G_oat σ²)`
— and the free factors layer interaction on top. The fallback changes in the
right direction too: when the free factors shrink away, predictions go to
`mu_j + f_1i`, which has oat information in it.

It is implemented as `fixed_main_effect = TRUE` in `setup_megalmm_state()`, and
is a swept level in the simulation grid.

### Three settings it drags along

- **`scale_Y` must be FALSE.** Fixed loadings apply on the sampler's scale;
  with per-column standardisation "1" would mean "equal in each column's SD
  units", and those SDs are noisy when a column holds ten observations. `Y` is
  divided once by a global SD instead.
- **`tot_F_var` needs loosening for factor 1.** With loadings pinned, the
  main-effect variance *is* `var(f_1)`, and the default inverse-gamma prior
  concentrated near 1 would pin it near the unit scale whatever the data says.
  For free factors that is harmless because scale trades off against `Λ`; here
  there is no `Λ` to absorb it. A length-K vector is accepted, so factor 1 gets
  `V = 0.5, nu = 3` and the rest keep the default.
- **The saved `Λ` row 1 is constant but not 1.** `remove_nuisance_parameters`
  rescales `Λ` by `sqrt(var(F))` so factors have unit variance, so row 1 comes
  back at the main-effect SD. The test that the mechanism held is that its **sd
  across columns is zero** — which the code now checks and warns about, since
  this path is lightly exercised upstream and a silent failure would look
  exactly like the unconstrained model.

One trap worth recording. The first attempt died with `Wrong R type for mapped
matrix` from `record_sample_Posterior_array`. That is **not** a fixed-factor
bug: `Eta` is in `setup_model_MegaLMM`'s default save list but has no allocated
posterior array, and any run that does not override `posteriorSample_params`
hits it. Setting the save list explicitly, as the production code already did,
removes it.

### It works, and it helps exactly where predicted

One replicate, K = 10, covariates spanning 80% of pea genetic variance, rank-1
interaction at 20% of variance. `fixed_ok` was TRUE in every fixed run, so the
loadings really did stay constant.

| observed | | `r_total` | `r_interaction` | `f_1` vs true producer | row average vs true producer |
|---|---|---|---|---|---|
| **1.1%** (400 panel) | free | 0.053 | −0.053 | −0.013 | 0.610 |
| | **fixed** | **0.135** | **0.053** | **0.404** | 0.610 |
| **5%** | free | 0.096 | −0.027 | 0.448 | 0.801 |
| | **fixed** | **0.227** | −0.002 | **0.613** | 0.801 |
| 15% | free | **0.867** | **0.799** | 0.090 | 0.911 |
| | fixed | 0.815 | 0.698 | 0.773 | 0.911 |
| 45% | free | **0.966** | **0.946** | −0.008 | 0.976 |
| | fixed | 0.942 | 0.897 | 0.710 | 0.976 |

Sparse, it helps: `r_total` more than doubles at both 1.1% and 5%, and at 1.1%
the interaction correlation crosses from negative to positive. Across the full
6-setting K × covariate sweep at 5% it improved `r_total` in **all six**, from
+0.01 to +0.25.

Dense, it costs a little: −0.05 on `r_total` and −0.10 on the interaction at
15%, less at 45%. That is the expected price — one of K factors is spent on the
main effect, leaving K−1 for structure, and the model is constrained where it
did not need to be.

### But it does not close the gap

The validation criterion worth holding it to was that `f_1` should recover the
producer effect *at least as well as a row average*. **It does not, at any
density**: 0.40 against 0.61 at 1.1%, 0.61 against 0.80 at 5%, 0.77 against
0.91 at 15%, 0.71 against 0.98 at 45%.

So the fixed factor is a real improvement on the free model — compare 0.40
against −0.01 at 1.1% — and a real confirmation of the diagnosis, since forcing
the missing term back in is what recovered the signal. It is not yet a reason
to prefer MegaLMM over the additive model in a sparse design.

Two reasons to treat that comparison as a lower bound rather than a verdict.
`f_1` alone is not the whole fitted main effect: a free factor can also acquire
near-constant loadings, so part of the oat signal lives outside `f_1`, and the
row margin of the full prediction surface would be the fairer measure. And this
is one replicate at one K. The honest summary is that the mechanism is
confirmed and the remedy is directionally right, with the size of the remaining
gap not yet established.

## One honest caveat

These numbers come from a single replicate per configuration, at K = 10, with
the pea covariates set to span 80% of the pea genetic variance. The mechanism —
no oat main effect, overlap going as `p²`, shrinkage towards the wrong margin —
does not depend on those choices, and the row-margin diagnostic above is direct
evidence for it. But the exact crossover density does depend on K, and the full
grid has not been run since K became something we sweep. Expect the 15% figure
to move down with smaller K.
