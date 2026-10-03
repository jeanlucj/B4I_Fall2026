# tests/

Unit tests for the code under `code/`. Run them after any change:

```sh
Rscript tests/run_all.R          # fast tier, ~45 s
Rscript tests/run_all.R --all    # everything, adds ~10 s of real MCMC fits
Rscript tests/test_surface.R     # one file, when that is what you changed
```

Each file runs in its own process, prints `<name>: N passed, M failed` and exits
non-zero if anything failed. `run_all.R` aggregates and reports timings.

| file | what it pins | s |
|---|---|---|
| `test_grm.R` | GRM reading, collapsing onto analysis names, `grm_factor`, variance truncation | 1 |
| `test_kronecker.R` | the Khatri–Rao basis for the specific-combination term, including the `byrow` reshape | 1 |
| `test_curation.R` | identity groups, the selfed-female signature, representative choice, clonal families, analysis-name rules | 1 |
| `test_generator.R` | the bivariate simulator: exact variances and covariances, exact sparsity, interaction rank, GxE | 1 |
| `test_surface.R` | the additive/interaction decomposition, which margin is which effect, masking, scoring | 1 |
| `test_design.R` | grid collapses, the composite recoding, the D-optimal fraction, the cache key | 36 |
| `test_validation.R` | pool construction, the field design's balance, the three-term power formula and its Satterthwaite df, the named-only argument barrier | 3 |
| `test_fits.R` | **slow**: real BGLR Multitrait and MegaLMM chains — recovery, a null, and the streamed-draw identities behind the per-accession PEV | 12 |

## What makes these tests worth having

Not the framework — `tests/helper.R` is four assertions. It is the **oracles**.
Almost every test compares the code against a quantity known independently of
it:

- an **algebraic identity** (`M == additive_part(M) + interaction_part(M)`, to
  1e-12);
- a **planted answer** (a similarity matrix with a known duplicate chain in it,
  a family that was built not to segregate);
- a value that was **requested** (a scenario labelled 1.6% sparsity must observe
  exactly 1.6% of cells; `cor(oatPr, oatAs)` must come back at −0.065);
- a **negative control** (a two-way model on the raw grid coding must be
  singular; the naive Pr+As / Pr−As rule must produce overlapping pools).

A test that only pins today's output freezes bugs as well as behaviour, so there
are as few of those as possible.

The suite has already paid for itself three times:

- `test_kronecker.R` found `kron_basis()` using `ncol(A)` as the rank of *both*
  bases on its first run. Harmless while the two ranks agreed, wrong as soon as
  variance truncation made them differ.
- `test_validation.R` had **nine fully positional calls** to `contrast_power()`
  relying on the argument order `(dAs, sigma2_within, n, sigma2_e, P)`. Inserting
  any argument before `P` would have silently rebound all nine while the suite
  kept reporting PASS. The functions now put `...` ahead of everything but the
  first argument, so R refuses to match the rest positionally, and a
  `...length()` guard catches a misspelled name that would otherwise vanish into
  `...`. Two negative controls keep the barrier in place.
- **Nothing pinned `df == 2n - 2`**, so changing the degrees of freedom would
  have passed the whole suite in silence. It is pinned now, and the Satterthwaite
  correction has to reproduce it exactly as the zero-interaction boundary case.
- `test_fits.R` pins that per-accession PEV comes from the streamed draws and
  **not** from `L %*% SD.beta`, which ignores the posterior covariance among
  coefficients. That negative control exists so a future optimisation to the
  cheap-looking route fails loudly.
- `test_validation.R` found `make_validation_design()` taking its anchor plots
  from the *rotated* selection window, so the anchor combination differed at
  every location — the opposite of what anchors are for. Combination × location
  was not estimable in any design the script would actually have produced.
- `test_design.R` prompted the guard in `sim_design()` that names the parameter
  count instead of letting `optFederov` fail with "columns in expanded X".

Two other findings were about the *tests*, and are recorded in the files
themselves because they are easy to get wrong again: a near-zero variance
component is **not** a null (correlation is scale-free, so a 0.1%-variance
effect is still recovered at r ≈ 0.5), and one independent truth is not a null
either on a 40-accession panel (two independent draws reach |r| = 0.55 by
chance).

## Conventions

- `check(cond, msg)`, `check_near(actual, expected, tol, msg)`,
  `check_error(expr, msg)` — from `tests/helper.R`, sourced by every file.
- Tolerances say what kind of claim is being made: `1e-12` for an identity,
  `1e-9` for an empirically rescaled draw, loose bounds for anything Monte
  Carlo. A tight tolerance on a short MCMC chain is a test that fails for the
  wrong reason.
- Every file seeds its RNG at the top. A suite that fails one run in twenty
  teaches people to ignore it.
- Nothing reads `output/`, so the suite runs in a fresh clone with no
  credentials and no network.
