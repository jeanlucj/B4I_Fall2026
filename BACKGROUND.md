# BACKGROUND

The theory behind the methods and the reasoning behind the decisions. For usage
see [README.md](README.md); for structure see [DESIGN.md](DESIGN.md); for how
accuracy is measured see [CROSS_VALIDATION.md](CROSS_VALIDATION.md); for when
either framework would work see [SIMULATION.md](SIMULATION.md); for what is
still open see [docs/B4I_followups.md](docs/B4I_followups.md).

## The breeding question

An intercrop plot contains two crops, and its value depends on both. A genotype
therefore has two distinct effects: what it does for **itself** (producer, or
direct genetic effect) and what it does to its **partner** (associate, or
indirect genetic effect). A breeder selecting oats for intercropping wants both,
and wants to know whether they trade off — whether the oat that yields best is
the one that suppresses the pea hardest.

The proposal in `docs/B4I_Proposal_Models.docx` sets out two related models, and
this repository fits the second of them plus a third framing.

## The bivariate DGE-IGE model

Conventionally the producer–associate model is fitted once per species, which
puts an oat's producer effect and its associate effect in two different
analyses — so the covariance between them, the parameter the breeder actually
wants, cannot be estimated. Stacking both species' yields into one bivariate
response fixes that: the oat kernel then carries a 2×2 covariance whose diagonal
is the producer and associate variances and whose off-diagonal is exactly the
covariance of interest. The same holds for pea.

### Why a BRR on `Z L` rather than an RKHS on `Z G Z'`

Both say the same thing about the data — `Var(Z L b) = Z G Z' ⊗ Σ` when
`G = L L'` — but not the same thing about what comes back. `BGLR::Multitrait`
fits an `RKHS` term by eigen-decomposing the kernel and running a ridge
regression on its eigenvectors, so the coefficients live in the eigenvector
basis of the *plot* kernel and carry no accession identity at all. Because
`rank(Z G Z')` equals the number of accessions, the result has exactly the shape
a per-accession effect would have, which is how the original script came to rank
eigenvectors and label them with accession names. Factoring `G` instead and
fitting a `BRR` on `Z L` keeps the coefficients on the accession scale:
`L %*% beta` is the accession effect, named, and `Var(L b) = G ⊗ Σ` is the
distribution the proposal specifies. `docs/BGLR_RKHS_effects_problem.md` has the
full account.

### Identifiability, and what the design will not support

Producer and associate effects separate only when accessions meet several
partners; a genotype grown with one partner has the two aliased. The
specific-combination term separates from the residual only when combinations are
replicated. In the B4I data 91% of oat–pea combinations occur in a single plot,
so that term is left out, and 99 oat and 84 pea accessions have a single
partner — almost all of them entries that appear only in the single harvested
2026 trial. Blocks are fitted only where a trial has more than one, since a
constant block within a trial is perfectly aliased with its fixed trial effect.

The fitted result: both producer–associate covariances are negative at the
posterior median but neither excludes zero (oat −12.3, 95% CrI [−69.8, +36.7];
pea −53.0, [−118.4, +1.8]). The residual covariance between the two yields on a
plot is clearly negative, −99.8 [−137.2, −64.5]. Note that the *raw* plot
correlation between the two yields is +0.50: the positive sign is entirely
trial and block, and the within-plot relationship reverses once environment is
accounted for. That is what competition looks like.

## Curation: which accessions are actually distinct

Genomic prediction assumes the rows of the GRM are different genotypes. Two
things break that here.

**Crossing failures.** A plant in the crossing nursery that self-pollinates
yields progeny recorded as several different crosses that are all the same
selfed line. The marker signature is unambiguous: a full-sib family segregates
and its members correlate around 0.74, while a family that never segregated sits
at 0.98–1.00. Grouping accessions by *both* parents and taking the mean pairwise
correlation within each family separates the two cleanly — the distribution is
bimodal with an almost empty valley between 0.88 and 0.97, which is also what
sets the 0.99 identity threshold, since within genuine families individual sib
pairs reach 0.97–0.99 legitimately.

Eleven of 74 oat families did not segregate. Where one female heads more than one
such family and the families cannot be told apart from each other — IL18-8562's
two families correlate 0.983 between and 0.985 within — the recorded pollen
parent is making no genetic difference and they are one selfed line, named
`<female>_self`.

**Duplicate entries.** Two names can be one genotype. This is the only thing
the pea side turns up: `ND VICTORY` and `NDP170084G` correlate at **r =
0.9996** across 6,164 markers, against 0.824 for the next-closest line and
0.701 for the closest of the other 49 `NDP17*` entries. So they are the same
line, not sibs.

**Why they are the same line is not known.** An earlier version of this
document called `NDP170084G` ND Victory's experimental designation. That was an
inference, not a reading of any record, and it is wrong: the cultivar
registration gives the designation as `NDP100144G`
([Plant Registrations 2023](https://acsess.onlinelibrary.wiley.com/doi/10.1002/plr2.20266)),
and `NDP100144G` is not in T3/Oat at all. T3 records no synonym linking
`NDP170084G` to `ND VICTORY` either — the only synonym on `ND VICTORY` is the
case variant `ND Victory`.

What the markers establish is that one of these two entries is not what its
name says. A seed or sample mix-up at genotyping, a maintenance reselection
issued under a new number, or a mislabelled T3 entry would all produce this,
and the marker data cannot tell them apart. The collapse is still right —
they are one genotype and belong in the model once — but which name is the
correct one is unresolved, and picking `ND VICTORY` as the representative was
a phenotype-count decision, not an identification.

Such accessions are **renamed, not dropped**. Dropping discards their phenotype
records; renaming keeps every plot and enters the genotype once, which is both
more data and a more honest representation. The GRM is collapsed the same way,
by averaging the rows and columns of a group — for identical lines that is the
GRM of their averaged marker profiles.

The pea side cannot be curated on pedigree at all: every B4I pea accession reads
`NA/NA` on T3/Oat. That is a gap rather than a defect, since the pea markers are
clean, but it means the family analysis simply does not run for pea.

What was actually collapsed, and under which thresholds, is reported in
[CURATION.md](CURATION.md), generated by `code/curation_report.R` from the settings
each curation run recorded for itself.

## The MegaLMM framing: pea accessions as environments

Read from the oat's side, each pea accession is an environment in which oats were
grown. 442 oats and 423 peas then define a genotype-by-environment matrix, and
MegaLMM's business is exactly such matrices: a low-rank factor model in which
factor loadings describe which environments rank genotypes alike.

The difficulty is density. Of 186,966 cells, 2,059 are filled — **1.1%**. A
multi-environment trial that MegaLMM is usually asked about might be 30–70%
filled. Whether the factor model can carry this is an empirical question, which
is why the analysis is built around cross-validation rather than around a fit.

### What goes in a cell

Three candidates, all implemented and compared rather than argued about:

| | what it is | why it might be right |
|---|---|---|
| `raw` | oat yield as measured | no assumptions imposed |
| `centered` | yield minus its trial mean | trials differ several-fold in mean yield; that difference is management and season, not the pea |
| `standardized` | centered, then divided by the trial SD | also equalises the spread, so no trial dominates the residual |

Centering and standardizing are done **within trial and before averaging**,
because the trial is the unit that carries the season and management effects.

The cross-validation settles it: `raw` is hopeless — even a simple oat mean
predicts held-out cells at only r = 0.02, because between-trial differences swamp
everything — while `centered` reaches 0.28 and `standardized` 0.21. Centering
removes the nuisance; standardizing removes a little real signal with it, since a
trial where genotypes genuinely spread out more is being flattened.

### Environmental covariates: why they are not optional here

After trimming (below) each pea environment holds about seven oats. Seven
observations cannot pin down a column's factor loadings. The extended MegaLMM
lets loadings be predicted from covariates describing the environment, and with
this little data per column that is not a refinement but the only thing giving
the model any purchase.

The covariates characterise the pea:

- **pea grain yield BLUE** across the whole experiment — how productive the
  partner is, adjusted for trial;
- **pea flowering date BLUE** — when it competes. Flowering and maturity BLUEs
  correlate 0.679, which is high enough that using both adds little; flowering
  is preferred because it was scored in all six trials while maturity was not
  scored in AL;
- **eigenvectors of the pea GRM**, enough to reach 80% of the pea genetic
  variance. This costs 58 eigenvectors, because the pea panel has weak
  population structure — PC1 is only 12% of the total. That flatness is worth
  noting: when a spectrum is this flat, principal components are a poor summary
  of a relationship matrix, and 58 covariates for 222 environments is a high
  ratio. The three covariate types are given separate ARD groups so the 58
  markers cannot swamp the two phenotypic covariates by sheer number.

### Trimming, and why it is not optional either

84 pea environments held a single oat and 169 held two or fewer. A column with
one observation has no residual variance to estimate, and MegaLMM's ARD sampler
does not object — it produces `NaN` and dies several frames deep in
`sample_Lambda_prec_ARD` with "missing value where TRUE/FALSE needed". Trimming
iteratively to at least 3 oats per pea and 2 peas per oat leaves a 237 × 222
matrix at 3.2% filled, which samples stably. The cross-validation masking carries
the same constraint: cells are held out only while their row and column can spare
one, because masking a row empty reproduces the same failure.

The missing-data map needs the opposite of the usual choice. For a dense matrix
the least-grouped map is safest; here it makes a group per column, one comes back
empty, and initialisation fails. The most-grouped map is both correct and fast,
since with seven observations per column there is no per-column missingness
pattern worth preserving.

## What the MegaLMM analysis found

It does not work on this data, and the failure is consistent across every
configuration tried. Mean correlation with held-out cells, over three
replicates at 20% masking:

| cell value | covariates | best MegaLMM | oat mean | oat + pea mean |
|---|---|---|---|---|
| centered | no | 0.031 | 0.259 | **0.279** |
| centered | **yes** | **0.072** | 0.259 | **0.279** |
| raw | no | 0.010 | 0.024 | 0.009 |
| raw | yes | 0.040 | 0.024 | 0.009 |
| standardized | no | 0.021 | 0.201 | 0.210 |
| standardized | yes | 0.002 | 0.201 | 0.210 |

All three of MegaLMM's prediction targets were scored — `Eta_mean` (the
predicted phenotype, the like-for-like comparison) and both genetic-value
matrices — so a weak result cannot be blamed on having picked the wrong one.

Two things stand out. The covariates **do** help: with `centered` values they
lift the model from 0.031 to 0.072, and they help in every filling. That is the
extended MegaLMM doing what it is supposed to. But the model does not come close
to a row mean. A factor model that cannot beat "this oat is a good oat" has not
found genotype-by-genotype structure, whatever its internal fit looks like.

That conclusion agrees with what the bivariate model already said from a
completely different direction: with 91% of oat–pea combinations appearing in
one plot, specific-combination effects are not estimable in this design. Two
methods, one answer. The sparsity is not a nuisance to be modelled around — it
is the finding.

None of this says the framing is wrong. It says this experiment does not yet
have the combination replication to support it. The design changes that would
follow — replicating specific oat×pea pairs rather than spreading singletons
thinly — are the same ones the bivariate model needs for `σ_PrAs`.
