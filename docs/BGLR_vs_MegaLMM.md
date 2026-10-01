# BGLR versus MegaLMM, with MegaLMM given its best shot

Two questions, and they have different answers:

- **GMA** — does analysing both yields jointly (BGLR `Multitrait`) beat two
  single-orientation factor models?
- **Interaction** — can the factor model decompose oat × pea specific combining
  ability better than a Kronecker kernel can, well enough that the decomposition
  might be biologically interpretable?

Produced by `code/sim_compare.R`. Column definitions:
[SIMULATION_GLOSSARY.md](../SIMULATION_GLOSSARY.md). The settings rule tested
here comes from [MegaLMM_anova.md](MegaLMM_anova.md).

**Both hypotheses hold. The second one holds conditionally, and the condition is
the problem.**

---

## How the comparison is made

**Paired, within scenario.** Every model is fitted to the *same* simulated data —
same scenario, same replicate, same seed, same train/held split — so the unit of
analysis is the within-scenario difference. Comparing marginal means would throw
the pairing away and drown a 0.07 contrast in the scenario-to-scenario variation
that sparsity alone drives from r = 0.15 to r = 0.91.

**Differences are taken on the Fisher-z scale** and reported back as a gap in
correlation points at the comparator's own mean, which is where the comparison
actually sits.

> **These results are from the September 2026 Ceres run, which used the
> pre-October scoring scheme**: 20% of the observed cells held out, so the models
> were fitted at 0.8 × the labelled sparsity and the per-cell metrics were scored
> on the held-out cells rather than on the never-observed ones. Read every
> sparsity level as 0.8 × its label — the density gate quoted as "between 4.8%
> and 16%" is between **3.84% and 12.8%** actually fitted. The scheme has since
> changed (see [SIMULATION_GLOSSARY.md](../SIMULATION_GLOSSARY.md)); these
> analyses have not been rerun, because the conclusions they support do not turn
> on it.

**Replicate 1 is excluded** — it holds two sets of cache files under different
seeds. 480 scenario-replicates remain, every one of which has both BGLR fits and
at least one MegaLMM run.

### Why "MegaLMM at its best" is estimated rather than subsetted

The natural approach — keep only rows with the preferred settings — does not
survive contact with the design. The run list is a **D-optimal fraction**, so
most scenarios were never given the preferred combination at all:

| filter | scenario-replicates kept |
|---|---|
| any MegaLMM run | 480 |
| `K = 5` | 352 (73%) |
| pin rule only | 384 (80%) |
| **both** | **192 (40%)** |

Worse, the survivors are unbalanced across exactly the axis the comparison is
about — 46% of the rank-1 cells survive against 33% of the rank-5 ones.

**And the averaging step is a no-op.** Of the 192 surviving scenario-replicates,
**4** have both `eigen_variance` levels. The fraction almost never assigns both
to the same scenario at a fixed `K` and pinning, so there is no replication there
to average over. (No loss: the ANOVA found `eigen_variance` inert, with negative
partial ω² on all four responses.)

So the settings go into the model as factors and the preferred configuration is
read off with `emmeans`, using all 480 and balanced by construction. The 192-row
subset is still computed as a check. The two agree in sign and ordering
everywhere, and differ by about 0.05 on the z scale — **the subset consistently
flatters MegaLMM**, which is the selection effect showing up exactly where it was
expected.

---

## 1. GMA: BGLR wins, and it is not close

MegaLMM minus BGLR, so negative favours BGLR.

| response | comparator | gap (r), subset | gap (r), all data | MegaLMM wins |
|---|---|---|---|---|
| `r_addsurf_oat` | `additive` | −0.072 | −0.088 | 5% |
| `r_addsurf_oat` | `dge_ige` | −0.073 | −0.089 | 3% |
| `r_addsurf_pea` | `additive` | −0.079 | −0.098 | 3% |
| `r_addsurf_pea` | `dge_ige` | −0.080 | −0.099 | 1% |

MegaLMM wins 1–5% of scenarios. All p < 10⁻¹⁶.

**The two BGLR models are indistinguishable from each other here** (−0.348 vs
−0.353 on z), which is the expected result: GMA is the additive part of the
surface, and adding a specific-combination term should not change it. It is a
sanity check passing.

The gap is strongly density-dependent, and closes:

| observed | gap (r), `dge_ige` vs MegaLMM |
|---|---|
| 1.6% | **−0.241** |
| 4.8% | −0.120 |
| 16% | −0.026 |
| 48% | −0.006 |

At 48% observed the two are near enough equivalent. The joint analysis earns its
advantage **where data is scarce** — which is the regime that matters, and is
where B4I sits.

This is the hypothesis as stated, and the mechanism is the one already
identified: MegaLMM has a per-column intercept and no per-row one, so one of the
two effects is always an unshrunk fixed effect. `Multitrait` puts kinship on all
four effects through `Σ ⊗ G`. When there is little information per accession,
shrinkage is worth a great deal.

---

## 2. Interaction: MegaLMM wins, but only above a density threshold

MegaLMM minus `dge_ige` (the only BGLR comparator that has an interaction term —
`additive` has none, so its `r_int` is `NA` by construction, not poor).

Overall: **+0.189** (oat) and **+0.162** (pea) on the r scale. But MegaLMM wins
only 51% and 45% of scenarios, which is the tell: the mean is not describing a
uniform advantage. Broken out, the structure is stark.

### Gap by interaction rank

| interaction | gap (r) | MegaLMM wins |
|---|---|---|
| rank 1, 10% variance | **+0.248** | 61% |
| rank 1, 20% variance | **+0.371** | 58% |
| rank 5, 10% variance | −0.050 | 31% |
| rank 5, 20% variance | +0.051 | 47% |

**The advantage is entirely at rank 1.** This is the hypothesis, confirmed: a
low-rank interaction surface is what a factor model is built to find, and at rank
5 MegaLMM has no advantage over a Kronecker kernel at all.

### But density gates it

| observed | gap (r) | MegaLMM wins |
|---|---|---|
| 1.6% | −0.195 | 12% |
| 4.8% | −0.244 | 6% |
| 16% | +0.247 | 75% |
| 48% | **+0.516** | **100%** |

### The two together — the table that decides it

`r_int_oat`, gap on the r scale; positive favours MegaLMM.

| | 1.6% | 4.8% | 16% | 48% |
|---|---|---|---|---|
| **rank 1** | −0.186 (12%) | −0.189 (10%) | **+0.373 (100%)** | **+0.572 (100%)** |
| **rank 5** | −0.205 (12%) | −0.315 (0%) | +0.035 (44%) | **+0.356 (100%)** |

`r_int_pea` is the same picture (−0.184, −0.267, +0.354, +0.558 at rank 1).

**Density is a gate, not a modifier.** Below about 16% observed, MegaLMM is worse
at recovering the interaction *at either rank* — it does not merely lose its
advantage, it goes clearly negative. Above the gate, rank decides the size of the
win: at 16% only the rank-1 case wins, and by 48% both do.

The crossover sits between 4.8% and 16% observed.

---

## What this means for B4I

**B4I is at about 3% of the oat × pea matrix observed.** That is below the gate,
at both ranks.

So the hoped-for outcome — use MegaLMM to decompose the interaction into
something biologically interpretable — is **not supported by the current data**,
and the simulation says so in the specific way that is most useful: it is not
that the method is wrong, it is that the design is too sparse for it.

And the two models fail differently, which matters:

| observed | `dge_ige` `r_int_oat` | `megalmm` `r_int_oat` |
|---|---|---|
| 1.6% | **0.205** | −0.002 |
| 4.8% | **0.285** | 0.029 |
| 16% | 0.416 | **0.549** |
| 48% | 0.479 | **0.852** |

The Kronecker kernel **degrades gracefully** — it still recovers the interaction
at r ≈ 0.2–0.3 where only a few percent of cells are observed, which is what a
kernel does: it borrows across relatives rather than needing to see structure
repeated. MegaLMM is at **zero** there, and then overtakes sharply once there is
enough to find factors in.

So below the gate the right thing to do is not "give up on the interaction", it is
"use the Kronecker kernel and do not expect axes you can interpret". The kernel
gives a covariance with real signal at B4I's density; what is unavailable is the
*decomposition*.

This lines up with what the real data already said from two directions: 91% of
observed oat × pea combinations occur in a single plot, `fit_mix_term` is off in
production because the term is weakly identified, and the fitted residual
correlation between the two yields cannot be split from specific-combination
covariance without replication.

**The actionable version:** getting to ~16% observed is what would make the
interaction question answerable, and it would be answerable by MegaLMM
specifically — with a 0.37 advantage over the Kronecker kernel if the underlying
interaction is low-rank. That is a statement about how to design the next round
of trials, not about which software to run on the current one. It also raises the
value of the validation trial's anchor plots, which are the only replicated
combinations in the programme.

For **GMA**, which is what selection decisions actually use, the answer is
unambiguous and does not depend on any of this: **use the BGLR `Multitrait`
model.** Its advantage is largest exactly where B4I sits.

---

## How much of this is the tuning bias?

The settings rule was chosen from the same data it is evaluated on, so it
flatters MegaLMM. Quantified three ways — on the r scale, MegaLMM minus BGLR:

| response | rule | average of all settings | per-cell oracle | rule − average | oracle − rule |
|---|---|---|---|---|---|
| `r_addsurf_oat` vs `dge_ige` | −0.073 | −0.099 | −0.089 | +0.026 | −0.016 |
| `r_addsurf_pea` vs `dge_ige` | −0.080 | −0.110 | −0.099 | +0.030 | −0.019 |
| `r_int_oat` vs `dge_ige` | +0.189 | +0.129 | +0.153 | +0.060 | −0.036 |
| `r_int_pea` vs `dge_ige` | +0.162 | +0.109 | +0.134 | +0.053 | −0.028 |

### What each column is

Every number is **MegaLMM minus BGLR on the r scale**, so negative favours BGLR
and positive favours MegaLMM. The three middle columns are three different
answers to "which MegaLMM settings do you get?", computed on the same
scenario-replicates so they are directly comparable.

| column | how MegaLMM's settings are chosen | is it achievable? |
|---|---|---|
| **rule** | the fixed rule under test: `K = 5`, pin the main effect when sparsity ≤ 4.8%. One choice, applied to every scenario. | **Yes** — it depends only on sparsity, which you know before fitting. |
| **average of all settings** | no choice at all: average over every setting the design happened to run for that scenario. | **Yes**, trivially — it is what you get by picking settings arbitrarily. |
| **per-cell oracle** | the best-performing setting *for each scenario separately*, chosen by looking at which one actually scored highest. | **No.** It uses the answer to pick the settings, so it cannot be done on real data where the truth is unknown. |

And the two derived columns:

- **rule − average** = what the rule buys you over not thinking about settings at
  all. This is the value of the tuning advice, and it is the part that could be
  inflated by having been chosen on this same data.
- **oracle − rule** = what is still on the table. It is the gap between the fixed
  rule and a per-scenario choice made with perfect foresight, so it is an **upper
  bound** on what any smarter rule could add. A small value means the fixed rule
  has captured nearly all of the available tuning benefit; a negative value means
  the rule is doing *better on this subset* than a per-cell choice would on
  average, which is what a rule fitted to these data looks like.

The reason for reporting all three is that "is the tuning biased?" is not
answerable in the abstract. What matters is whether the bias is large enough to
change a conclusion, and that needs the size of the tuning effect (rule −
average) set against the size of the effect being claimed.

### Reading them

The rule is worth **0.026–0.030** on GMA and **0.053–0.060** on the interaction,
against an arbitrary choice of settings. The per-cell oracle — the best setting
for each scenario, chosen *knowing the answer*, so not achievable — adds very
little beyond it, and on this subset actually sits below the rule, which is what
a rule tuned on these data looks like.

**The bias is real, it is in MegaLMM's favour, and it is far too small to change
either conclusion.** BGLR beats MegaLMM on GMA by 0.07–0.10; the entire tuning
advantage is 0.03. Giving MegaLMM its best shot does not close a third of that
gap. And on the interaction the advantage at rank 1 and 48% is 0.57 against a
tuning effect of 0.06.

So the instinct — that the bias would not be great — was right, and it is worth
having measured rather than assumed.

---

## Caveats

**The pinning rule is optimised for GMA and is mildly wrong for the
interaction.** `fixed_main_effect = TRUE` helps GMA when sparse but costs
interaction recovery at every density (0.525 → 0.478 overall). The rule pins at
1.6% and 4.8%, where it costs `r_int` about 0.005–0.018 — small, because almost
nothing is recoverable there anyway, but it is a cost, and it is paid in the
regime where MegaLMM already loses. Using a per-response rule would make
MegaLMM's interaction numbers slightly *better* at the sparse end and would not
move the gate.

**Only `dge_ige` can be compared on the interaction.** "Both BGLR models" is not
available for `r_int`; `additive` has no interaction surface.

**Sparsity is a 4-level factor here**, so "the gate is between 4.8% and 16%" is as
precise as this design can be — and the same applies to rank, where 1 wins and 5
does not, with nothing in between. Both boundaries are what
`code/sim_int_run.R` exists to locate: a **full factorial** over five sparsity
levels (adding 9%), three ranks (adding 3), and the pinning switch crossed rather
than fixed. See [code/sim_int_config.R](../code/sim_int_config.R) for what was
changed and why.

**This is one replicate structure.** Reps 2–5 of one grid; the gate's position is
estimated from 16–24 scenario-replicates per cell of the rank × density table.
The sign and ordering are unambiguous; the exact crossover is not.
