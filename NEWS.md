# SpatialML 1.8.2 (2026-06-10)

## Packaging

* Added the `.Rbuildignore` file so that top-level development artefacts
  (the RStudio project file, the technical notes document, local tool
  settings and `cran-comments.md`) are excluded from the source tarball.
  The 1.8.1 release notes stated this exclusion was already in place,
  but the file itself was missing; `R CMD check --as-cran` consequently
  reported NOTEs about hidden and non-standard top-level files. Both
  NOTEs are now resolved.
* Added `cran-comments.md` (not shipped in the tarball) documenting the
  test environments and the resubmission context.

## DESCRIPTION metadata

* Added the maintainer's ORCID iD to `Authors@R`.
* Simplified the roles of Stefanos Georganos to `"aut"` (the previous
  `c("aut", "ctb")` was redundant).
* Added `Language: en-GB`, matching the British spelling used
  throughout the documentation.

## Documentation

* The `\donttest{}` examples of `grf()` and `predict.grf()` now set
  `nthreads = 1` explicitly. This was already CRAN-compliant (ranger
  defaults to two threads) and is purely a defensive clarification.

# SpatialML 1.8.1

## CRAN Resubmission

* Prepared the package for CRAN resubmission after archival. The known
  CRAN archival issue was an Rd cross-reference NOTE for unanchored
  external links to `ranger`; all such references now use
  `\link[ranger]{ranger}`.
* Version and date updated to 1.8.1 / 2026-06-09.
* Top-level development artefacts are excluded from `R CMD build`
  through `.Rbuildignore`.
* Removed "active development" wording from help pages to better match
  CRAN's publication-quality expectation.

## Dependency Reduction

* The package now depends only on `R` in `Depends:` and imports
  `ranger`, `stats`, `graphics` and `utils`.
* Removed the runtime dependencies on `caret` and `randomForest`.
  The two affected functions, `rf.mtry.optim()` and `grf.bw()`, are now
  implemented directly with `ranger` and base R helpers.
* Replaced broad namespace exports/imports with explicit
  `export()`, `S3method()` and `importFrom()` declarations.

## `grf()` Robustness

* `grf()` no longer repeatedly refits local forests when `ranger`
  returns a `NaN` out-of-bag prediction for the focal observation.
  The default fallback is now a leave-one-out local refit: the focal row
  is removed from the local training subset and predicted out-of-sample.
  This preserves finite local diagnostics without noisy console output.
* Added `oob.fallback = c("loo", "inbag")` to `grf()`. The default
  `"loo"` uses the new leave-one-out fallback; `"inbag"` reproduces the
  older in-sample fallback behaviour.
* The focal observation is now forced to be the first row of each local
  subset before reading `ranger`'s first training prediction. This avoids
  incorrect focal OOB extraction when two or more observations share the
  same coordinates.
* Degenerate zero-distance local kernels now receive equal weights
  instead of producing non-finite bi-square weights.
* The local-fit loop no longer rebuilds a full copy of `dframe` on every
  iteration. Local neighbourhoods are indexed by row position, reducing
  memory traffic for medium-sized data sets.
* Added validation for `coords`, `bw`, `ntree`, `kernel` and formula
  inputs, plus clearer fixed-bandwidth diagnostics when too few
  neighbours are available.
* `grf()` now returns an object of class `"grf"` so `predict()` dispatches
  to `predict.grf()`.
* `LocalModelSummary` now reports AIC and AICc for both OOB and
  predicted residuals.

## Bandwidth and mtry Selection

* `grf.bw()` no longer uses `caret::postResample()`. Its `Mixed` and
  `Low.Local` diagnostics are computed by an internal base R helper
  equivalent to squared correlation.
* `grf.bw()` returns `Best.BW = NA_real_` with an informative warning
  when every evaluated bandwidth produces a non-finite local R-squared,
  instead of returning `numeric(0)`.
* `grf.bw()` gains `verbose = TRUE`, allowing bandwidth-search progress
  messages to be suppressed in tests, vignettes and batch scripts.
* `grf.bw()` now caps an over-large `mtry` at the number of predictors
  and uses a one-predictor-safe default.
* `rf.mtry.optim()` gains the arguments `num.trees`, `cv.repeats`,
  `num.threads` and `verbose`, and supports three evaluation strategies
  via `cv.method`: `"oob"` (default), `"cv"`, `"repeatedcv"`.
* `rf.mtry.optim()` no longer returns a `caret` `train` object. The new
  return value is a list with `best.mtry`, `results`, `cv.method`,
  `num.trees` and `call`. Code that read `out$bestTune$mtry` should now
  read `out$best.mtry`.

## Documentation, Examples and Tests

* Added a `tests/testthat/` test suite (testthat 3rd edition) covering
  the public API of every exported function.
* Added a vignette, `vignettes/SpatialML.Rmd`, with a complete workflow:
  `mtry` tuning, bandwidth search, GRF fit, local-importance plotting
  and prediction at new locations.
* The vignette now uses the `knitr::knitr` engine with `markdown`
  output, avoiding the RStudio/Quarto/Pandoc probe path during local
  package checks.
* The `\donttest{}` examples in `grf.Rd`, `grf.bw.Rd`,
  `predict.grf.Rd` and `rf.mtry.optim.Rd` now use small, single-threaded
  examples suitable for CRAN checks. Heavier `Income` demonstrations are
  kept in `\dontrun{}` blocks.
* The maintainer's website (<https://stamatisgeoai.eu/>) is linked from
  the package help page and vignette in addition to `DESCRIPTION`.

# SpatialML 1.7.0

## CRAN compliance

* Fixed the missing package anchors in Rd cross-references that triggered
  the CRAN archival on 2025-11-28: every `\link{ranger}` is now
  `\link[ranger]{ranger}` (in `grf.Rd`, `grf.bw.Rd`, `predict.grf.Rd` and
  `SpatialML-package.Rd`).
* Moved the heavy run-time dependencies (`ranger`, `caret`,
  `randomForest`) from `Depends` to `Imports` in line with current CRAN
  best practice. `Depends` now contains only the minimum R version.
* Replaced `exportPattern("^[[:alpha:]]+")` in `NAMESPACE` with explicit
  `export()` statements; removed the unused `import(randomForest)` in
  favour of a narrow `importFrom()`.
* Added explicit `stats::` and `graphics::` namespacing to all internal
  calls (`predict`, `formula`, `terms`, `setNames`, `dist`, `sd`,
  `rnorm`, `rpois`, `runif`, `plot`).
* Added a package-level help page (`?SpatialML`).
* Switched URLs and DOIs in the documentation to the canonical
  `\doi{...}` macro.

## Bug fixes

* `grf()`: pre-allocation of the local-forest list previously read
  `as.list(rep(NA, length(ntrees)))`, producing a one-element list
  because `ntrees` is a scalar; fixed to `vector("list", Obs)`.
* `grf()`: replaced `subset(DataSetSorted, DNeighbour <= bw)` with
  bracket indexing to avoid the standard CRAN check NOTE about
  non-standard evaluation.
* `grf()`: the function now returns an object of class `"grf"` so that
  `predict(obj, ...)` dispatches correctly to `predict.grf()` without
  the user having to call the S3 method directly.
* `random.test.data()`: when `dep.var.dis = "normal"` the dependent
  variable is now drawn with `rnorm()` (was incorrectly drawn with
  `runif()`).
* `random.test.data()`: fixed the operator-precedence bug
  `for (i in 1:vars.no - 1)` (which evaluated to `0:(vars.no - 1)`); the
  predictor matrix is now filled vectorised and column-named
  `X1, X2, ...`.
* `predict.grf()`: replaced `object[[1]]` by `object$Global.Model`;
  validates that `object$Forests` exists and that `x.var.name` /
  `y.var.name` match columns of `new.data`; `...` is now propagated to
  both global and local `predict()` calls.
* `rf.mtry.optim()`: `caret::trainControl()` arguments are now passed by
  name (`cv.method` was previously matched positionally and could land
  in the wrong slot).

## New features and improvements

* `grf()`, `grf.bw()`, `random.test.data()`: added `match.arg()`-based
  validation for the `kernel` and `dep.var.dis` arguments and explicit
  validation of `coords`, `bw`, `ntree`, `dataset` and `step`.
* `grf()`: improved diagnostic stop messages when the fixed bandwidth
  yields too few neighbours.
* `grf()`: returned `LocalModelSummary` now also reports AIC and AICc
  for both OOB and predicted residuals.
* `rf.mtry.optim()`: new `plot.it = TRUE` argument lets the user disable
  plotting when running non-interactively or in vignettes.
* Removed every hard-coded `set.seed()` from inside the package
  functions (CRAN policy). Reproducibility is now under user control via
  a `set.seed()` call before each function.

## Documentation

* Fixed typos: "Radmom" -> "Random", "lcoal" -> "local",
  "Neightbours" -> "Neighbours", "tradition" -> "traditional".
* `grf.Rd`: the `\value{}` block was renamed `LocalModelSummary` to
  match the actual returned name (was incorrectly `lModelSummary`).
* `grf.bw.Rd`: the worked example used `weighted = TRUE` (a
  non-existent argument); replaced with `geo.weighted = TRUE`.
* `predict.grf.Rd`: examples now call `predict()` (which dispatches via
  the new `"grf"` class) instead of the explicit `predict.grf()`.
* Expanded the Description field of every Rd file with the bi-square
  weighting equation and a step-by-step description of the algorithm.
