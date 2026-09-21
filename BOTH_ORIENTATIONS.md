# BOTH ORIENTATIONS

Ideas for using the oat-side and pea-side MegaLMM analyses together. Nothing
here is implemented; this is a design note to argue with. For the simulation
that motivates it see [SIMULATION.md](SIMULATION.md), and for the models
themselves [BACKGROUND.md](BACKGROUND.md).

## The asymmetry

Read from the oat's side, each pea is an environment. MegaLMM then gives the
**oat** axis a relationship matrix — `relmat = list(germplasmName = G_oat)` —
and the **pea** axis only a per-column intercept and K factor loadings, with
no kinship at all.

That asymmetry is structural, not a choice in how we called it. MegaLMM's
environment axis does not take a covariance matrix; its traits are traits, and
traits are not related to each other by pedigree. The pea covariates are the
designed substitute, entering as `X` in the Lambda prior so a pea's loadings
can be predicted from what kind of pea it is. But a prior on loadings is a
weaker channel than a covariance acting directly, which is what the DGE-IGE
model gives the pea side through `As ~ N(0, G_pea σ²)`.

So the oat-side analysis is well specified for oats and under-specified for
peas. Flipping it — peas as rows, oats as environments — is well specified for
peas and under-specified for oats. Neither orientation is wrong; each is half
an analysis.

## Why not just pick one

Because the quantity a breeder wants is symmetric. The producer–associate
framing treats oat and pea alike: each genotype has an effect on itself and an
effect on its partner. An analysis that models oat relatedness properly and pea
relatedness only through a prior will estimate the oat associate effect better
than the pea associate effect, for reasons that have nothing to do with
biology.

It also matters for the comparison in SIMULATION.md. Part of why MegaLMM trails
DGE-IGE on `r_total` is that the pea main effect sits in an unshrunk per-column
intercept rather than in a kinship-shrunk random effect. Running one
orientation and reporting it as "MegaLMM's accuracy" understates a framework
that was never given the pea relationships.

## Idea 1: covariates from the other orientation

The most direct version of what you suggested. Run the pea-side analysis,
extract something that characterises each pea, and hand it to the oat-side
analysis as an environmental covariate — replacing or supplementing the GRM
eigenvectors.

What to carry across is the real question. Candidates:

- **The pea's factor scores `U_F`** from the pea-side fit. These are the pea
  axis's genetic values on K latent dimensions, shrunk through `G_pea`. Using
  them as oat-side covariates says "predict how an oat performs with this pea
  from where that pea sits on the axes that organise pea performance."
- **The pea's fitted main effect** — its row mean of `Eta_mean` in the pea-side
  fit, which is a kinship-shrunk estimate of exactly the quantity the oat-side
  analysis keeps in an unshrunk intercept. This is the smallest, most targeted
  version and the one I would try first.
- **The pea's loadings on the oat-side factors**, which is circular on the
  first pass but not on later ones.

The appeal of `U_F` over the GRM eigenvectors is that eigenvectors are a
generic summary of relatedness, whereas factor scores are relatedness already
filtered through what actually predicted performance. On a flat spectrum like
the B4I pea panel — PC1 is 12% — that filtering should be worth a lot.

## Idea 2: iterate to a fixed point

The EM analogy you drew. Alternate:

```
initialise  pea covariates <- pea GRM eigenvectors
repeat
    oat-side fit   |  environment covariates = current pea summary
    pea summary    <- something extracted from the oat-side fit
    pea-side fit   |  environment covariates = current oat summary
    oat summary    <- something extracted from the pea-side fit
until the summaries stop moving
```

Attractive, and worth trying, but I would not expect it to behave like EM and
it is worth being clear about why.

EM converges because each step maximises the same objective. Here the two fits
maximise *different* objectives — one is a model of oat performance, the other
of pea performance — and neither is a conditional step of a single joint
likelihood. What this scheme is, formally, is closer to alternating conditional
estimation or to a Gibbs-flavoured backfitting, neither of which carries a
convergence guarantee. Two failure modes to watch:

- **Oscillation.** Each side chases the other's last estimate and they orbit
  rather than settle.
- **Variance collapse.** Each pass shrinks the summaries a little; after enough
  passes the covariates carry nothing and both analyses revert to their
  uninformed versions. This is the one I would bet on, because MegaLMM's ARD
  prior shrinks the covariate effects, and shrinking an already-shrunk quantity
  repeatedly is a contraction.

Mitigations if it does misbehave: damp the update (take a weighted average of
the old and new summary rather than replacing it), fix the summaries' scale
each pass, or simply stop after one round — which is Idea 1 and may capture
most of the available gain.

**What to measure.** Whether iterating helps is a cross-validation question,
not a convergence question. Held-out cell accuracy after 0, 1, 2, ... rounds
answers it directly, and if round 2 is no better than round 1 the fixed point
is irrelevant.

## Idea 3: average the two orientations

Cheaper and more robust than iterating. Fit both orientations once, predict
each held-out cell from each, and combine — a simple mean, or a weighted one
with weights from cross-validation.

Each orientation predicts the same cell, so they are directly combinable. And
each is strong exactly where the other is weak: the oat-side fit shrinks oat
effects properly, the pea-side fit shrinks pea effects properly, and a held-out
cell needs both. Model averaging with two complementary, individually-unbiased
predictors is usually hard to beat and cannot oscillate.

I would try this before iterating. It is a few lines, it has no convergence
question, and if it recovers most of the gap to DGE-IGE on `r_total` then the
asymmetry was the whole story and the more elaborate schemes are unnecessary.

## Idea 4: stop using MegaLMM for the main effects

The asymmetry only bites because we are asking a factor model to carry the
additive main effects as well as the interaction. It does not have to.

Fit the additive part once with the DGE-IGE model, which handles both species'
relatedness symmetrically and correctly. Take its residuals. Hand *those* to
MegaLMM, whose only job is then the interaction — the thing the simulation says
it is much better at (0.87 against 0.44 on a dense matrix, rank 1).

This gives up joint estimation, and the residuals carry the uncertainty of the
first stage as though it were zero, which is the standard two-stage criticism.
But it plays each framework to its strength rather than asking either to do the
other's job, and it sidesteps the orientation question entirely: the additive
stage is symmetric by construction, and the interaction stage has no main
effects left to mis-shrink.

Of the four, this is the one I think is most likely to produce a usable
analysis, and Idea 3 the one most likely to answer the narrower question of how
much the orientation asymmetry actually costs.

## What would settle it

All four are testable in the existing simulation framework, which already knows
the truth and already scores held-out cells. The comparison to run:

| variant | what it tests |
|---|---|
| oat-side only (current) | the baseline |
| pea-side only | is the asymmetry symmetric — does it cost the same both ways? |
| average of the two (Idea 3) | how much of the gap is orientation |
| one round of cross-covariates (Idea 1) | does a better pea summary beat GRM eigenvectors |
| iterate to convergence (Idea 2) | does more than one round buy anything |
| DGE-IGE residuals into MegaLMM (Idea 4) | does separating the jobs beat doing both badly |

Scored on `r_interaction` and `r_total` separately, because Idea 4 should win
the second without changing the first, and Idea 1 the reverse.
