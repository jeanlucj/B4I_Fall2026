# test_trial_qc.R
#
# apply_trial_qc() decides which trials reach the production fit and the
# cross-validation. Its failure modes are both silent: drop a trial nobody meant
# to drop, or -- worse -- quietly keep every trial because the verdict file was
# missing, which looks exactly like "all trials passed".
#
# So the tests pin the three states apart: a verdict that drops something, a
# verdict that drops nothing, and no verdict at all.
#
# Run: Rscript tests/test_trial_qc.R

library(tidyverse)
here::i_am("tests/test_trial_qc.R")
source(here::here("tests", "helper.R"))
load_code("dge_ige_functions.R")

tmp <- tempfile(fileext = ".csv")
on.exit(unlink(tmp), add = TRUE)

plots <- tibble::tibble(
  studyName = rep(c("GOOD_A", "GOOD_B", "DROPME"), times = c(4, 3, 5)),
  oat_yield = seq_len(12))

# ------------------------------------------------------------
# 1. A verdict that drops a trial drops exactly that trial's plots.
# ------------------------------------------------------------

readr::write_csv(tibble::tibble(
  studyName = c("GOOD_A", "GOOD_B", "DROPME"),
  keep = c(TRUE, TRUE, FALSE)), tmp)

kept <- suppressMessages(apply_trial_qc(plots, qc_file = tmp))
check(nrow(kept) == 7, "the failing trial's plots are dropped")
check(!"DROPME" %in% kept$studyName, "and the trial itself is gone")
check(setequal(unique(kept$studyName), c("GOOD_A", "GOOD_B")),
      "while every passing trial survives intact")

# ------------------------------------------------------------
# 2. A verdict that drops nothing changes nothing. Not a tautology: an
#    off-by-one in the keep logic would show up here as an empty table.
# ------------------------------------------------------------

readr::write_csv(tibble::tibble(
  studyName = c("GOOD_A", "GOOD_B", "DROPME"), keep = TRUE), tmp)
check(nrow(suppressMessages(apply_trial_qc(plots, qc_file = tmp))) == 12,
      "an all-pass verdict leaves the table untouched")

# ------------------------------------------------------------
# 3. NO verdict file returns everything -- and must say so. An older analysis
#    has to keep working, but "the screen never ran" and "the screen passed
#    everything" are different facts and cannot look the same.
# ------------------------------------------------------------

msg <- utils::capture.output(
  invisible(out <- apply_trial_qc(plots,
                                  qc_file = file.path(tempdir(), "no_such.csv"))),
  type = "message")
check(nrow(out) == 12, "a missing verdict file keeps every trial")
check(any(grepl("no trial QC", msg)),
      "and says so, rather than passing silently")
check(length(utils::capture.output(
        invisible(apply_trial_qc(plots,
                                 qc_file = file.path(tempdir(), "no_such.csv"),
                                 quiet = TRUE)),
        type = "message")) == 0,
      "quiet = TRUE really is silent")

# ------------------------------------------------------------
# 4. It works off trialF too, which is the column name the model-side plot
#    table uses, and re-levels the factor so an empty level cannot reach lmer.
# ------------------------------------------------------------

readr::write_csv(tibble::tibble(
  studyName = c("GOOD_A", "GOOD_B", "DROPME"),
  keep = c(TRUE, TRUE, FALSE)), tmp)

model_side <- tibble::tibble(
  trialF = factor(plots$studyName), oatYield = seq_len(12))
m <- suppressMessages(apply_trial_qc(model_side, qc_file = tmp))
check(nrow(m) == 7, "trialF is recognised as the trial column")
check(!"DROPME" %in% levels(m$trialF),
      "and the dropped trial leaves no empty factor level behind")

finish("trial QC tests")
