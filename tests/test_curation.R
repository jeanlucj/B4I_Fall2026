# test_curation.R
#
# The curation logic in code/curation_functions.R, tested against PLANTED
# answers rather than against T3. A synthetic similarity matrix with a known
# duplicate in it, and a synthetic pedigree with a known clonal family, are
# oracles: the right answer is decided before the code runs. That is also what
# makes these tests runnable with no network and no credentials.
#
# Run: Rscript tests/test_curation.R

library(tidyverse)
here::i_am("tests/test_curation.R")
source(here::here("tests", "helper.R"))
load_code("curation_functions.R")

set.seed(1)

# ------------------------------------------------------------
# 1. identity_groups: a planted duplicate, and a planted CHAIN
#
# Single linkage is deliberate -- a chain of near-identical accessions is one
# line however it is labelled -- so the chain is asserted rather than treated as
# a defect. acc3 and acc9 are never made similar to each other; they are one
# group only because acc7 bridges them.
# ------------------------------------------------------------

n <- 12
S <- diag(n)
rownames(S) <- colnames(S) <- paste0("acc", seq_len(n))
S[lower.tri(S)] <- runif(sum(lower.tri(S)), 0.5, 0.8)
S[upper.tri(S)] <- t(S)[upper.tri(S)]
S["acc3", "acc7"] <- S["acc7", "acc3"] <- 0.995
S["acc7", "acc9"] <- S["acc9", "acc7"] <- 0.993

grp <- identity_groups(S, 0.99)
check(nrow(grp) == 3, "identity_groups finds the three chained accessions")
check(setequal(grp$germplasmName, c("acc3", "acc7", "acc9")),
      "identity_groups finds the RIGHT three")
check(dplyr::n_distinct(grp$group) == 1, "they are one group, not three")
check(S["acc3", "acc9"] < 0.9,
      "and acc3/acc9 were never similar -- the chain is what joined them")

# nothing above the threshold means no groups at all, not an error
S_clean <- S
S_clean["acc3", "acc7"] <- S_clean["acc7", "acc3"] <- 0.7
S_clean["acc7", "acc9"] <- S_clean["acc9", "acc7"] <- 0.7
check(nrow(identity_groups(S_clean, 0.99)) == 0,
      "no near-identical pairs gives an empty result")

# raising the threshold above a pair's similarity must drop it
check(nrow(identity_groups(S, 0.996)) == 0,
      "a threshold above every pair leaves nothing grouped")

# ------------------------------------------------------------
# 2. classify_groups: the selfed-female signature
#
# One female, several recorded males, all genetically identical: that is a plant
# that self-pollinated in the crossing nursery. The oracle is that the flag
# fires on exactly that pattern and not on a group that is genuinely one cross.
# ------------------------------------------------------------

ped_self <- tibble::tibble(
  germplasmName = paste0("acc", c(3, 7, 9)),
  pedigree      = c("F/M1", "F/M2", "F/M3"),
  seed_parent   = "F",
  pollen_parent = c("M1", "M2", "M3"),
  both_parents  = TRUE
)
cls <- classify_groups(grp, ped_self, 0.5)
check(all(cls$selfing_suspect), "one female with three males flags as selfing")
check(all(cls$modal_seed == "F"), "the modal seed parent is identified")
check(all(cls$n_pollen_on_modal == 3), "all three pollen parents are counted")

ped_one_cross <- dplyr::mutate(ped_self, pedigree = "F/M1", pollen_parent = "M1")
cls2 <- classify_groups(grp, ped_one_cross, 0.5)
check(!any(cls2$selfing_suspect), "one cross does NOT flag as selfing")
check(all(cls2$same_cross), "it is recognised as a single cross instead")

ped_none <- dplyr::mutate(ped_self, pedigree = NA_character_,
                          seed_parent = NA_character_,
                          pollen_parent = NA_character_, both_parents = FALSE)
cls3 <- classify_groups(grp, ped_none, 0.5)
check(!any(cls3$selfing_suspect), "no pedigree cannot flag as selfing")
check(all(cls3$interpretation == "no pedigree recorded"),
      "and says so rather than guessing")

# ------------------------------------------------------------
# 3. choose_representatives: phenotype count, then alphabetical
#
# With no observations file every count is 0 and the representative is the
# alphabetically first name -- silently. That is worth pinning, because it is a
# tie-break masquerading as a decision.
# ------------------------------------------------------------

rep_no_obs <- choose_representatives(cls, NULL)
check(sum(rep_no_obs$action == "keep") == 1, "exactly one representative kept")
check(all(rep_no_obs$representative == "acc3"),
      "with no phenotype counts the alphabetically first name is kept")

obs_file <- tempfile(fileext = ".rds")
saveRDS(tibble::tibble(germplasmName = rep(c("acc3", "acc7", "acc9"),
                                          times = c(1, 50, 2)),
                       value = "1"), obs_file)
rep_obs <- choose_representatives(cls, obs_file)
check(all(rep_obs$representative == "acc7"),
      "the best-phenotyped member is kept when counts are available")
check(rep_obs$action[rep_obs$germplasmName == "acc7"] == "keep",
      "and it is the one marked keep")

# ------------------------------------------------------------
# 4. Full-sib families: a planted family that did not segregate
#
# A real cross segregates and its members correlate around 0.75; a family whose
# members all correlate above 0.985 never segregated. Both are planted here so
# the verdicts are known in advance.
# ------------------------------------------------------------

nm <- paste0("k", 1:6)
S2 <- matrix(0.75, 6, 6); dimnames(S2) <- list(nm, nm)
S2[1:3, 1:3] <- 0.995
diag(S2) <- 1

ped2 <- tibble::tibble(
  germplasmName = nm,
  seed_parent   = rep(c("A", "B"), each = 3),
  pollen_parent = rep(c("X", "Y"), each = 3),
  pedigree      = paste0(rep(c("A", "B"), each = 3), "/",
                         rep(c("X", "Y"), each = 3)),
  both_parents  = TRUE
)

fams <- pedigree_families(ped2, nm, 2)
check(dplyr::n_distinct(fams$family) == 2, "two families found")
fs <- summarise_families(family_pair_correlations(fams, S2), 0.985)
check(sum(fs$clonal_family) == 1, "exactly one family flagged as clonal")
check(fs$family[fs$clonal_family] == "A/X", "and it is the planted one")
check_near(fs$mean_r[fs$family == "B/Y"], 0.75, tol = 1e-9,
           "the segregating family's mean r is the planted 0.75")

# a family below the minimum size is not summarised at all
check(nrow(pedigree_families(ped2, nm, 4)) == 0,
      "min_family_size excludes families that are too small")

# ------------------------------------------------------------
# 5. resolve_analysis_names: the three rules, and their precedence
# ------------------------------------------------------------

an <- resolve_analysis_names(ped2, fs, S2, 0.99, 0.985, 0.01)
check(nrow(an) == 3, "only the clonal family is renamed")
check(all(an$analysis_name == "A_X_no_cross"),
      "a family that did not segregate becomes <seed>_<pollen>_no_cross")
check(all(stringr::str_detect(an$reason, "did not segregate")),
      "and the reason says why")

# Rule 3 overrides: a line identical to its GENOTYPED seed parent is named for
# the parent, discarding the _no_cross label. Tested by putting the parent in
# the similarity matrix at r = 1 with one of its progeny.
nm3 <- c(nm, "A")
S3 <- matrix(0.75, 7, 7); dimnames(S3) <- list(nm3, nm3)
S3[1:3, 1:3] <- 0.995
S3["k1", "A"] <- S3["A", "k1"] <- 0.999
diag(S3) <- 1
ped3 <- dplyr::bind_rows(ped2, tibble::tibble(
  germplasmName = "A", seed_parent = NA_character_,
  pollen_parent = NA_character_, pedigree = NA_character_, both_parents = FALSE))
fs3 <- summarise_families(
  family_pair_correlations(pedigree_families(ped3, nm3, 2), S3), 0.985)
an3 <- resolve_analysis_names(ped3, fs3, S3, 0.99, 0.985, 0.01)
check(an3$analysis_name[an3$germplasmName == "k1"] == "A",
      "a line confirmed identical to its genotyped seed parent takes its name")
check(an3$seed_parent_confirmed[an3$germplasmName == "k1"],
      "and is flagged as confirmed")

finish("curation tests")
