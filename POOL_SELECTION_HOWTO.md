# POOL SELECTION — how to run it, and how to iterate it

This is the operator's guide to `code/validate_pool_selection.R`: what to type,
what to read, what to type back, and how to tell when to stop. It is the one
step of the validation chain that cannot be automated, because choosing the
pools is a conversation with the seed store.

**Why the rule is what it is lives in
[VALIDATION_DESIGN.md](VALIDATION_DESIGN.md) §2** — the Pr-adjusted objective,
the two search modes, the feasibility arithmetic. This file does not repeat the
argument; it tells you how to run it. [README.md](README.md) §3 has the rest of
the pipeline.

---

## 1. Before you start

Pool selection is **offline**. It needs no T3 credentials; it reads the fit that
is already on disk:

| it reads | written by |
|---|---|
| `output/BGLR_{oat,pea}_effects_all_seeds.csv`, `BGLR_{oat,pea}_pev.csv` | `code/BGLR_multi_trait_model.R` |
| `output/BGLR_variance_components.csv` | same |
| `output/BGLR_fit_trials.csv`, `BGLR_fit_provenance.csv` | same — the fit's own record of what it saw |
| `output/trial_qc.csv` | `code/assemble_B4I_phenotypes.R` |
| `data/seed_search_ledger.csv` | **you**, round by round |
| the most recent `crossval_summary.csv` under `output/validation/` | `code/validate_crossval.R` |

`validation_inputs()` refuses to run if the fit's own record of what it saw
disagrees with the QC-filtered plot table, so a stale fit stops the run rather
than quietly selecting on the wrong accessions.

The cross-validation is the one soft dependency. λ is used for the power
**printed beside** the proposal and for nothing else: the search maximises a
λ-free score, and the chosen membership is invariant to λ across the whole
plausible range. If no `crossval_summary.csv` exists anywhere, the run warns,
falls back to a guess, and says so. If the newest one is from an older vintage,
the console prints `<- NOT this vintage`. Either way the accessions are the
same; only the number beside them moves.

Everything lands in `output/validation/<today>/`, dated, so a refresh never
overwrites the pools a decision was taken on.

---

## 2. Round 1 — the first proposal

```bash
Rscript code/validate_pool_selection.R
```

With an empty or all-`unsearched` ledger, `--search auto` resolves to
**unconstrained**: the optimum over all candidates, with nothing held fixed.
The console says which mode it chose and why, per species — a refresh must
never switch mode silently.

What to read in the console, in order:

1. **the vintage** — trials, plots, accessions, reliability. If this is not the
   fit you meant to select on, stop here.
2. **λ and its source vintage** — with the reminder that it moves the power and
   not the membership.
3. **the mode line** — `unconstrained` / `anchored`, with the usable anchor and
   exclusion counts that decided it.
4. **per species**: eligible accessions and the rule, candidates above the Pr
   threshold, ΔAs, the Pr-adjusted ΔAs with the slope `b`, ΔPr against its
   tolerance, and the **selection score** (`t` at λ = 1 with no interaction —
   *not* the `t` the trial will see, which is `t_reported` in
   `pool_summary.csv`, about half as large).
5. **`reached from theta = … in N move(s)`** — which θ-ranked start the answer
   was climbed from, and how far. `N = 0` means the θ scan was already at the
   optimum and the console says so outright; a small `N` (2 for oat and 3 for
   pea on this vintage) means the climb moved a couple of members for a gain of
   about 0.0015 in `t`. Either way an **unconstrained** proposal is essentially
   the θ solution — the search is here so that round 2 can be *constrained*.
6. **the search request** — one line per slot to fill.

### What round 1 writes

| file | what it is for |
|---|---|
| `pools.csv` | **the adopted design.** `validate_design.R` and `validate_power.R` read this |
| `pool_summary.csv` | what each pool achieved and at what constraint — 56 columns, the audit trail |
| `search_request.md` | **the deliverable: the printable list to take to the freezer** |
| `search_request.csv` | the same thing machine-readable, including every alternate |
| `search_input_buckets.csv` | every name the ledger supplied, accounted for in exactly one bucket |
| `search_curve.csv` | what each extra substitution would buy (anchored runs only) |
| `found_missing_balance.csv` | is the store biased in As or Pr? (once anything is searched) |
| `pool_diff.csv` | who entered and left each pool since the previous vintage |
| `pools_Pr_vs_As.png` | where the pools sit in the Pr × As plane |
| `vintage.csv` | the fit these pools came from |

and `data/seed_search_ledger.csv` gains one `unsearched` row per proposed
accession — 100 of them on round 1, 40 oat and 60 pea.

### Reading `search_request.md`

One block per species × pool, one **slot** per bag to find:

```
**Slot 3**
- [ ] **IL24-6924**  (Pr 16.7, As 11.4, first choice)
- [ ] IL23-621  (Pr 20.8, As -1.0, cost 0.1711)
- [ ] IL18-4663_IL18-2658_no_cross  (Pr 18.4, As -1.2, cost 0.1719)
      any one of: IL23-1206; IL23-1209; IL23-1225; IL23-1228; IL23-1236; IL23-1259
- [ ] IL19-5336  (Pr 20.6, As -1.1, cost 0.1782)
```

- **Find ONE accession per slot** — the first choice, or the first alternate you
  can lay hands on. A slot is a seat in the design, not a shopping list.
- `cost` is the drop in the contrast's `t` from taking that alternate instead of
  the first choice. Read it against λ's own standard error, which is worth about
  ±15 power points: a few thousandths is a genuine tie, and even round 1's
  alternates — the expensive ones, see §9, item 3 — cost about 3% of `t`.
- Alternates are **disjoint across slots**: each name is offered for one slot
  only, because one bag can fill one seat.
- An `_no_cross` name is a *collapsed analysis name* — one genotype the fit
  could not distinguish from its siblings. Any one of the listed entries
  satisfies it; write the one you actually have into `sowable_entry`.
  See [VALIDATION_DESIGN.md](VALIDATION_DESIGN.md) §7c.

---

## 3. The freezer trip — filling in the ledger

`data/seed_search_ledger.csv` is **hand-typed primary data**. It lives in
`data/` beside `trial_qc_manual.csv`, not in a regenerated vintage folder, it is
committed, and the driver only ever appends to it.

| column | written by | meaning |
|---|---|---|
| `species`, `acc` | the script | the accession, in the fit's naming |
| `round_proposed`, `pool_proposed` | the script | when and where it was first asked for |
| `status` | **you** | `unsearched` / `found` / `missing` / `missing_retry` |
| `searched_on` | **you** | the date you looked |
| `sowable_entry` | **you** | for a collapsed name, the entry you actually hold |
| `note` | **you** | anything the next round should know |

The three statuses do three different things, and the difference matters:

| status | what the search does with it |
|---|---|
| `found` | **anchors** it — the search must keep it (subject to the allowance) |
| `missing` | **excludes it, cumulatively across every later round.** A confirmed-absent accession is never proposed again |
| `missing_retry` | neither. It records a search that was not finished — "didn't look hard enough", "may be in the other freezer" — so the accession stays a candidate but does not anchor |

Use `missing_retry` honestly. `missing` is a permanent statement about the
store; `missing_retry` is a statement about the searcher. They are also kept
apart in `found_missing_balance.csv`, which tests whether the store is
systematically short of high- or low-As lines — a caveat on the As− pool's
interpretation if it is real, and one that `missing_retry` rows would otherwise
pollute.

Two rules the reader enforces, loudly:

- a status it does not recognise stops the run — no silent coercion to
  `unsearched`;
- two rows for the same `species` × `acc` stop the run. One accession, one
  status.

Rows you have filled in come back untouched on the next run. Only
`pool_proposed` on a still-`unsearched` row is refreshed, and the console
reports that as `n re-pooled`.

---

## 4. Round 2 and after — the anchored run

```bash
Rscript code/validate_pool_selection.R
```

The same command. `--search auto` switches to **anchored** for a species as soon
as its ledger has one usable `found` row, and says so. The search now holds the
found accessions, allows at most `--max-new` substitutions per pool (3 by
default), re-optimises, and asks for replacements only for the rest.

The round number is **derived, not incremented**: a round stays open while
anything it asked for is still `unsearched`. So re-running mid-round — which is
idempotent, and safe — does not advance the count, and the request you print
twice is the same round twice. `--round N` overrides it if you need to.

The anchored console adds three lines per species:

```
anchors: 34 of 34 usable used; 6 new accession(s) to find
allowance: 3 per pool (minimum feasible 3)  <- AT THE MINIMUM: every usable anchor must be kept, which costs power
t against the unconstrained optimum: 5.012 vs 5.2989 = 94.59%  (regret 0.2869 in t)
```

- **`n_new` vs `n_new_greedy`** — the parsimony stage hands slots back while the
  score stays within `score_slack_frac` (0.5% of `t`) of the optimum. When the
  unpruned optimum asked for more, the console says what those extra bags would
  have bought. It is usually nothing.
- **`AT THE MINIMUM`** is a warning, not a congratulation. See §9, item 1.
- **`regret`** is the unconstrained optimum *recomputed from this fit*, minus
  what the anchors allow: what the searching done so far has cost. It is printed
  every round including when it is zero.

Anchors are **one global list the search may reassign between pools.** An
accession you found for As+ can come back in As−. That is deliberate — a found
bag is found, and the As+/As− label is a property of the design, not of the seed
lot — but it will surprise you the first time. `anchor_pool_switches` counts it,
`previous_pool` in `search_request.csv` names them, and the console explains
itself whenever the count is non-zero.

### When the allowance is infeasible

Each pool must hold at least `n − max_new` anchors and the pools are disjoint,
so the allowance needs `2(n − max_new)` usable anchors: 34 of 40 for oat (85%
found) but **54 of 60 for pea (90%)**, because the allowance is an absolute
count while the pool is larger. Expect pea to hit this at 3 per pool.

The check runs **over every species before any pool is built**, by arithmetic
rather than by letting the first search throw — otherwise a run that fails on
oat never computes pea, and you learn "go and find more seed" one re-run at a
time. You get a table and then the run stops:

```
=== INFEASIBLE at --max-new 3 ===
 species n_per_pool n_usable_anchors max_new_per_pool anchors_required shortfall min_feasible_max_new new_accessions_at_min
     oat         20               30                3               34         4                    5                    10
```

Two ways out, both yours to choose — it is **never relaxed for you**, because
widening the allowance spends somebody's afternoon:

1. find `shortfall` more of the accessions already proposed, or
2. re-run with `--max-new min_feasible_max_new`, which asks for
   `new_accessions_at_min` new bags instead.

`search_curve.csv` prices option 2, and is written **before** the pools are
built, so it survives the run that stops.

---

## 5. Deciding whether to do another round

Four numbers, all in `pool_summary.csv`, all printed every round:

| number | the question it answers |
|---|---|
| `regret` = `score_unconstrained − score` | what the search done so far has cost |
| `pct_of_unconstrained` | the same thing as a fraction — how close to the free optimum you are |
| `score_if_remaining_found` | **what stopping now costs.** The anchored search re-run with the not-yet-searched accessions assumed found |
| `n_unsearched` | how many bags that assumption covers |

And `search_curve.csv`, one row per allowance from the minimum feasible up:

```
species, max_new_per_pool, feasible, n_anchors_usable, n_anchors_supplied, max_new_min,
score_unconstrained, power_unconstrained, delta_power_vs_unconstrained,
score_frac_of_unconstrained, n_new, n_new_greedy, score, power,
dAs, dPr, dAs_adj, at_min_budget, dpr_tol_relaxed
```

Read `power` against `n_new` down that table and "is another round worth it?"
becomes one monotone column.

**Calibrate every one of these against the noise floor.** λ is 0.759 ± 0.198
(oat) and 0.784 ± 0.140 (pea), i.e. about ±15 power points on the realised
contrast, and associate reliability is 0.29 / 0.36. A difference of a few tenths
of a power point is far below the sampling error of the quantity being
predicted. The index is **a screen against bad choices, not a fine-grained
preference order**, and an afternoon in a cold room to recover 0.3 points is not
a good trade.

---

## 6. Freezing the pools and re-issuing the design

Pool membership should be frozen **at the last moment compatible with seed
logistics, on the largest dataset available**. Everything before that is
provisional. When you freeze:

```bash
Rscript code/validate_power.R     # ~3.5 min; estimands read the adopted pools.csv
Rscript code/validate_design.R    # field book + six design checks
```

`validate_power.R`'s estimand section reads membership from `pools.csv` through
`pools_from_csv()` and recomputes every number from the fit, so
`estimand_power.csv` and `field_book.csv` describe **the same accessions**. An
adopted accession the fit no longer carries stops the run pointedly — that means
the fit moved under the adopted pools, which is a decision for a person.
`power_grid.csv` deliberately keeps re-deriving its own pools, because it sweeps
*n* and the adopted pools exist at one *n* only.

The whole chain, in the right order, is one command —
`Rscript code/validate_refresh.R` — and `--only validate` runs just the four
`validate_*` scripts. `validate_crossval.R` runs **first**, so the power printed
beside the pools comes from this vintage's λ rather than the last one's.

---

## 7. Flags

| flag | default | what it does |
|---|---|---|
| `--search MODE` | `auto` | `auto` is anchored as soon as a species has a usable `found` row, unconstrained otherwise, and prints which and why. `unconstrained` and `anchored` force it; `anchored` with nothing to anchor on is an error, not a silent fallback |
| `--max-new N` | 3 | substitutions allowed per pool — how many new bags you are willing to find. Needs `2(n − N)` usable anchors |
| `--ledger PATH` | `data/seed_search_ledger.csv` | the search record. Point it elsewhere to rehearse a round without touching the real one |
| `--round N` | derived | the round number. Only needed to override the "a round is open while anything it asked for is unsearched" rule |
| `--no-curve` | curve on | skip `search_curve.csv`. Don't: it is the output that answers "is another round worth it?" |

To rehearse a round, copy the ledger somewhere scratch, edit the copy, and run
with `--ledger`. That protects the real ledger — nothing else in the chain
writes to `data/` — but **not** the vintage folder: the rehearsal overwrites
`pools.csv` and everything beside it with its own answer. Re-run the real round
afterwards to put the vintage back, and check the console says `0 row(s) added`.

---

## 8. Worked trace — round 1 on the 2026-11-01 vintage

| | oat | pea |
|---|---|---|
| mode | unconstrained | unconstrained |
| reached from θ | 0.100, in 2 moves | 0.535, in 3 moves |
| ΔAs (g/m²) | 21.08 | 30.89 |
| ΔPr | 0.429 | −0.645 |
| Pr-adjusted ΔAs | 21.15 | 30.61 |
| selection score (`t` at λ = 1) | 5.299 | 5.100 |
| power at the measured λ | 0.727 | 0.889 |
| slots to find | 40 | 60 |
| regret | 0 | 0 |

The 2026-10-31 vintage was built by the θ scan alone, on the same fit, so the
diff against it isolates what the climb adds: `pool_diff.csv` shows 37 of 40 oat
and 58 of 60 pea members retained, 3 and 2 exchanged, for about 0.0015 in `t`.
That is the honest size of what the search adds to an **unconstrained**
proposal. It is not here for that; it is here so round 2 can be constrained.

A round 2 rehearsal — 34 of 40 oat and 54 of 60 pea marked found, the rest
missing, allowance 3 — asked for 6 new accessions per species and gave:

| | oat | pea |
|---|---|---|
| anchors used | 34 of 34 | 54 of 54 |
| ΔAs (g/m²) | 19.83 | 29.32 |
| selection score | 5.012 | 4.872 |
| % of the unconstrained optimum | 94.6 | 95.5 |
| regret (in `t`) | 0.287 | 0.228 |
| power | 0.716 | 0.876 |

So six missing bags cost oat 1.1 power points (0.727 → 0.716) and pea 1.3
(0.889 → 0.876) — and the allowance sat at its minimum in both cases, which is
the next section.

---

## 9. Five things that will surprise you

1. **The minimum feasible allowance costs power, unpredictably.** At the
   minimum every usable anchor *must* be kept, so the search cannot drop a
   badly-estimated one. Measured on this vintage at allowance 3 against
   allowance 4, over three random found-sets per species: oat at 34 of 40 found
   lost 0.1, 0.1 and 1.4 power points, and paid one or two extra bags to get
   them back; pea at 54 of 60 lost 0, 0 and 0.9, and got it back **for the same
   six bags**. So the penalty is real, often zero, occasionally a point, and its
   size and its price are **not** predictable from the counts — which accessions
   went missing decides both. If `at_min_budget` is `TRUE`, read
   `search_curve.csv` before you go anywhere: that is the whole reason it is
   written.
2. **A plain hill-climb on power asks for far more seed than it is worth.** The
   greedy optimum spends its entire allowance — at one point 11 new accessions
   instead of 8, to gain 0.03 power points. The parsimony stage is what stops
   that, and `n_new_greedy` is what it saved you.
3. **Round 1's alternates are the expensive ones.** Alternates are assigned
   disjointly, so the reservoir is thinnest when every slot is open: pea's 194
   candidates leave 134 non-members for 60 slots, and three alternates each
   would need 180. Every one of the 134 gets offered, the later slots get
   shorter lists, and round 1's alternates cost about 3% of `t` rather than the
   fraction of a percent a later round's do. Take the first choice in round 1
   if you can.
4. **An accession can change pool between rounds.** Global anchors; see §4.
5. **An anchored proposal need not be Pr-matched.** An adopted `|ΔPr|` above the
   1.0 g/m² tolerance triggers a warning, and it means the pre-registered `pool`
   estimand would have to condition on Pr — a `PREREG_MODEL` bump, which the
   driver will not make for you. Separately, if *no* feasible start satisfied the
   tolerance, the search relaxes it rather than returning nothing, says so
   loudly, and records `dpr_tol_relaxed = TRUE` — a column, not just a console
   line, because a warning in a long log gets missed.

---

## 10. Troubleshooting

| what you see | what it means |
|---|---|
| `unrecognised status in the ledger` | a typo in `status`. The four legal values are `unsearched`, `found`, `missing`, `missing_retry` |
| `duplicated species x acc` | two rows for one accession. Delete one |
| a name in `search_input_buckets.csv` with bucket `not_in_accession_table` | a naming problem, not a seed problem. The resolver tries the original `germplasmName` → `analysis_name` direction first; if it still fails, the name is not in the fit at all |
| bucket `eligible_but_below_Pr_threshold` | **reported, never quietly admitted.** The Pr threshold *is* the hypothesis — the pools are contrasted among *good producers* — so admitting it would change what the trial tests. The shortfall in g/m² is in the table; lowering `pr_quantile` has to be a deliberate, recorded act |
| bucket `ineligible` | below the reliability threshold, or below the two-partner floor, or both — the `reason` column says which, with the numbers. The floor is not negotiable: a single-partner accession has its producer and associate effects perfectly aliased, which is the one thing the pools exist to separate |
| `--search anchored was asked for, but a species has no usable anchor` | nothing is `found` yet for that species. Run unconstrained first |
| `no crossval_summary.csv anywhere` | λ is a guess, and the printed power with it. The accessions are unaffected. Run `validate_crossval.R` |
| `lambda ... <- NOT this vintage` | λ came from an older vintage. Same consequence as above; `lambda_source_vintage` records it |
| the round number did not advance | correct. A round is open while anything it asked for is `unsearched` |
| `\|dPr\| exceeds the tolerance` | see §9, item 5 |
