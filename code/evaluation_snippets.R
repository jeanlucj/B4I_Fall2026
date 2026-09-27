# evaluation_snippets.R -- the paste-along companion to
# EVALUATION_SIMULATION.md.
#
# NOT a script. Nothing here is meant to be run start to finish: it is the
# copy-pasteable form of the walkthrough, level by level, so that stepping
# through a module does not mean retyping its setup. Paste one block, read what
# lands in the Global Environment, disarm, move on.
#
# For the automated checks -- the identities, the planted answers, the requested
# variances -- see tests/ instead. These two do different jobs: the tests say a
# function computes what it claims, this file is for deciding whether it is the
# right thing to compute.

library(tidyverse)
here::i_am("code/evaluation.R")
source(here::here("code", "evaluation.R"))

eval_load("simulation")          # sim_config, sim_generate, sim_fit, megalmm_setup
eval_conflicts()                 # expect: no name defined twice
eval_groups("simulation")

panel <- sim_grms(100, seed = 100)
set.seed(7)
sim <- simulate_experiment(panel$G_oat, panel$G_pea, sparsity = 0.20,
                           n_factors = 1, interaction_pct = 0.20, n_envs = 1)
peek(sim)

# S1 Design
arm_evaluation("sim_design")
g <- sim_grid(SIM_LEVELS, n_reps = 1)
disarm_evaluation()

nrow(g); dplyr::n_distinct(g$scenario)
dplyr::count(g, n_factors, interaction_pct)
mm_grid <- tidyr::expand_grid(!!!SIM_MEGALMM_LEVELS)
nrow(mm_grid)

# S2 Panel
arm_evaluation("sim_panel")
panel <- sim_grms(200, seed = 200)
L <- grm_factor(panel$G_oat)
set.seed(1)
v <- draw_effect(L, 0.25)
disarm_evaluation()

max(abs(tcrossprod(L) - panel$G_oat)); var(v); mean(v)
# pea accessions that have no marker data and so are IID
sum(rowSums(abs(panel$G_pea - diag(diag(panel$G_pea)))) == 0)

# S3 Sparsity
arm_evaluation("sim_cells")
set.seed(2); cl <- sample_combinations(200, 200, 0.05)
disarm_evaluation()

nrow(cl)
round(0.05 * 200 * 200)
min(table(cl$oat)); min(table(cl$pea))
try(sample_combinations(200, 200, 0.002))      # must stop, not silently thin

# S4 Experiment
arm_evaluation("sim_truth")
set.seed(7)
sim <- simulate_experiment(panel$G_oat, panel$G_pea, sparsity = 0.20,
                           n_factors = 1, interaction_pct = 0.20, n_envs = 1)
disarm_evaluation()

sum(sim$truth$V)                                  # 1
var(sim$truth$producer);  sim$truth$V[["producer"]]
var(as.vector(sim$truth$interaction)); sim$truth$V[["interaction"]]
qr(sim$truth$interaction)$rank                    # == n_factors

# S5 Cross validation splits
arm_evaluation("sim_split")
a <- split_observations(sim, SIM_CV_FRACTION, seed = 7)
disarm_evaluation()
b <- split_observations(sim, SIM_CV_FRACTION, seed = 7)

identical(a$held$cell, b$held$cell)               # TRUE -- the load-bearing one
nrow(a$train); nrow(a$held)
min(table(a$train$oat)); min(table(a$train$pea))
length(intersect(a$train$cell, a$held$cell))      # 0

# S6 Kronecker basis
arm_evaluation("sim_basis")
A <- grm_basis(panel$G_oat, rank = 5)
B <- grm_basis(panel$G_pea, rank = 5)
disarm_evaluation()

set.seed(3); Beta <- matrix(rnorm(25), 5, 5)
oi <- c(1, 4, 9); pj <- c(2, 7, 11)
lhs <- as.vector(kron_basis(A, B, oi, pj) %*% as.vector(t(Beta)))
rhs <- (A %*% Beta %*% t(B))[cbind(oi, pj)]
max(abs(lhs - rhs))                               # 6.9e-18

# S7 BGLR
arm_evaluation("sim_bglr")
add <- fit_dge_ige(a$train, sim$G_oat, sim$G_pea, with_interaction = FALSE)
dge <- fit_dge_ige(a$train, sim$G_oat, sim$G_pea, with_interaction = TRUE)
disarm_evaluation()

dim(add$total); max(abs(add$interaction))         # 100 x 100; exactly 0
peek(as.vector(dge$interaction))

# S8 MegaLMM
arm_evaluation("sim_megalmm")
d <- file.path(tempdir(), "mmcheck"); dir.create(d, showWarnings = FALSE)
r <- run_scenario_megalmm(sim, K = SIM_MEGALMM_K,
                          eigen_variance = SIM_EIGEN_VARIANCE, seed = 7,
                          run_dir = d, fixed_main_effect = TRUE)
disarm_evaluation()
as.data.frame(r$scores)

