# =============================================================================
# Optimal bandwidth selection for the Geographically Weighted Random Forest
#
# Exhaustive grid search: for every candidate bandwidth in [bw.min, bw.max]
# (inclusive, with stride 'step') a full grf() is fitted and its OOB local
# R-squared recorded. Two additional R-squared values are reported as
# diagnostics: one based on a 50/50 average of local and global predictions
# and one based on a 25/75 weighted average.
#
# The search is computationally expensive (one full grf() per bandwidth);
# for big data sets restrict the range with bw.min/bw.max and use a coarse
# 'step', then refine.
# =============================================================================

grf.bw <- function(formula, dataset, kernel = "adaptive", coords,
                   bw.min = NULL, bw.max = NULL, step = 1, trees = 500,
                   mtry = NULL, importance = "impurity", nthreads = 1,
                   forests = FALSE, geo.weighted = TRUE,
                   verbose = TRUE, ...) {

  # ---- Argument validation ----
  if (!is.data.frame(dataset))
    stop("'dataset' must be a data frame.", call. = FALSE)
  kernel <- match.arg(kernel, c("adaptive", "fixed"))
  if (missing(coords) || NCOL(coords) < 2L || NROW(coords) != nrow(dataset))
    stop("'coords' must be a 2-column object with the same number of rows ",
         "as 'dataset'.", call. = FALSE)
  if (!is.numeric(step) || length(step) != 1L || step <= 0)
    stop("'step' must be a positive number.", call. = FALSE)

  Obs <- nrow(dataset)
  f <- stats::formula(formula)
  RNames <- attr(stats::terms(f), "term.labels")
  ModelVarNo <- length(RNames)

  DepVarName <- row.names(attr(stats::terms(f), "factors"))[1]
  Y <- dataset[[DepVarName]]

  # Default search range: 5%-95% of Obs (or at least p+2 / >= 20 neighbours).
  if (is.null(bw.min)) bw.min <- max(round(Obs * 0.05, 0), ModelVarNo + 2L, 20L)
  if (is.null(bw.max)) bw.max <- max(round(Obs * 0.95, 0), ModelVarNo + 2L)
  if (bw.max < bw.min)
    stop("'bw.max' must be greater than or equal to 'bw.min'.", call. = FALSE)
  if (is.null(mtry)) mtry <- max(floor(ModelVarNo / 3), 1L)
  if (mtry > ModelVarNo) {
    if (isTRUE(verbose))
      message("'mtry' (", mtry, ") exceeds the number of predictors (",
              ModelVarNo, "); capping at ", ModelVarNo, ".")
    mtry <- ModelVarNo
  }

  # ---- Empty result table ----
  eval.bw.grf <- data.frame(Bandwidth = numeric(0),
                            Local     = numeric(0),
                            Mixed     = numeric(0),
                            Low.Local = numeric(0),
                            stringsAsFactors = FALSE)

  # ---- Main grid search ----
  count <- 1L
  for (abw in seq(from = bw.min, to = bw.max, by = step)) {

    eval.bw.grf[count, 1] <- abw

    grf.a <- grf(formula, dframe = dataset, bw = abw, kernel = kernel,
                 coords = coords, ntree = trees, mtry = mtry,
                 importance = importance, nthreads = nthreads,
                 forests = FALSE, geo.weighted = geo.weighted,
                 print.results = FALSE, progress = FALSE, ...)

    # Local R-squared (OOB) is the primary criterion.
    eval.bw.grf[count, 2] <- grf.a$LocalModelSummary$l.r.OOB

    if (isTRUE(verbose)) {
      message("Bandwidth: ", abw,
              " - R-squared local: ", round(eval.bw.grf[count, 2], 4))
    }

    # 50/50 blend of local OOB and global in-sample predictions.
    eval.bw.grf[count, 3] <- .cor_r2(
      pred = (grf.a$LGofFit$LM_yfitOOB + grf.a$Global.Model$predictions) / 2,
      obs  = Y)

    # 25/75 blend (mostly the global model, used as a stability check).
    eval.bw.grf[count, 4] <- .cor_r2(
      pred = (grf.a$LGofFit$LM_yfitOOB * 0.25) +
             (grf.a$Global.Model$predictions * 0.75),
      obs  = Y)

    count <- count + 1L
  }

  # ---- Pick the bandwidth that maximises the local R-squared ----
  # which.max() returns integer(0) if every value is NA / NaN, which would
  # leave Best.BW empty. Guard explicitly so callers always receive either
  # a single numeric value or an explicit NA they can test against.
  if (all(!is.finite(eval.bw.grf$Local))) {
    warning("All evaluated bandwidths produced non-finite local R-squared; ",
            "no optimal bandwidth could be selected. Try a wider bandwidth ",
            "range, increase 'trees', or set 'geo.weighted = FALSE'.",
            call. = FALSE)
    best.bw <- NA_real_
  } else {
    best.bw <- eval.bw.grf$Bandwidth[which.max(eval.bw.grf$Local)]
  }
  if (isTRUE(verbose))
    message("Best Bandwidth (Based on the Local Model): ", best.bw)

  list(tested.bandwidths = eval.bw.grf, Best.BW = best.bw)
}
