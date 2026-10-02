# test_cache_freshness.R
#
# cache_is_fresh() decides whether a trial's downloaded observations are
# re-fetched. Getting it wrong is expensive in one direction and silent in the
# other, and it was silent: before this rule existed, a trial downloaded before
# its yields were uploaded kept that empty download forever, and the result was
# indistinguishable from "T3 has no yield data for this trial". That is exactly
# what happened to B4I_2026_ND and B4I_2026_NY on 2026-10-02.
#
# So the tests pin the decision at the boundary, and pin what happens when T3
# declines to say when a trial last changed.
#
# Run: Rscript tests/test_cache_freshness.R

library(tidyverse)
here::i_am("tests/test_cache_freshness.R")
source(here::here("tests", "helper.R"))
load_code("t3_functions.R")

cache <- tempfile(fileext = ".rds")
saveRDS(1, cache)
on.exit(unlink(cache), add = TRUE)

t_cache  <- file.mtime(cache)
before   <- t_cache - 3600   # trial changed an hour BEFORE the download
after    <- t_cache + 3600   # trial changed an hour AFTER it

# ------------------------------------------------------------
# 1. The rule itself: a cache is good only if it is newer than the trial.
# ------------------------------------------------------------

check(cache_is_fresh(cache, before),
      "a cache downloaded after the trial last changed is used")
check(!cache_is_fresh(cache, after),
      "a cache downloaded before the trial last changed is NOT used")

# ------------------------------------------------------------
# 2. The states that are not about time at all.
# ------------------------------------------------------------

check(!cache_is_fresh(tempfile(), before),
      "a cache file that does not exist is never fresh")
check(!cache_is_fresh(cache, before, refresh = TRUE),
      "refresh = TRUE overrides a cache that would otherwise be used")
check(!cache_is_fresh(NULL, before), "no cache path means no cache")

# ------------------------------------------------------------
# 3. An unknown modification time KEEPS the cache.
#
#    This is a deliberate choice and the riskiest one in the function, so it is
#    pinned rather than left to drift: re-downloading whenever T3 is silent
#    would re-fetch every cached trial on every run. The obligation that comes
#    with it -- saying so out loud -- lives in the caller.
# ------------------------------------------------------------

check(cache_is_fresh(cache, NA), "an NA modification time keeps the cache")
check(cache_is_fresh(cache, NULL), "a NULL modification time keeps the cache")
check(cache_is_fresh(cache, character(0)),
      "an empty modification time keeps the cache")
check(!cache_is_fresh(cache, NA, refresh = TRUE),
      "but refresh still wins over an unknown modification time")

# ------------------------------------------------------------
# 4. The real case, reconstructed. The September download of B4I_2026_NY sits
#    in the cache; the yields went up in October. The rule must re-fetch.
# ------------------------------------------------------------

sep_cache <- tempfile(fileext = ".rds")
saveRDS(1, sep_cache)
Sys.setFileTime(sep_cache, as.POSIXct("2026-09-19 13:02:00", tz = "UTC"))
on.exit(unlink(sep_cache), add = TRUE)

check(!cache_is_fresh(sep_cache, as.POSIXct("2026-10-01 09:00:00", tz = "UTC")),
      "the September cache is rejected once the trial changed in October")
check(cache_is_fresh(sep_cache, as.POSIXct("2026-08-01 09:00:00", tz = "UTC")),
      "and kept when the trial has not changed since")

# ------------------------------------------------------------
# 5. Character timestamps, which is how they arrive from BrAPI.
# ------------------------------------------------------------

check(!cache_is_fresh(sep_cache, "2026-10-01 09:00:00"),
      "a character timestamp is handled, not silently treated as unknown")

finish("cache freshness tests")
