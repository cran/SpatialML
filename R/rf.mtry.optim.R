# =============================================================================
# Optimal mtry selection for a Random Forest
#
# Lightweight replacement for the previous caret::train(method = "rf")
# implementation. Three evaluation strategies are provided:
#
#   * "oob"        - one ranger fit per mtry; OOB error is the criterion.
#                    Fast, no resampling overhead, statistically valid for
#                    Random Forests, and the default.
#   * "cv"         - manual k-fold cross-validation (one repeat).
#   * "repeatedcv" - manual repeated k-fold CV (cv.repeats * cv.folds fits
#                    per mtry).
#
# For all three strategies the criterion to *minimise* is the average RMSE
# across folds (or the OOB RMSE in the "oob" case). Average R-squared is
# also recorded for diagnostics.
# =============================================================================

rf.mtry.optim <- function(formula, dataset,
                          min.mtry = NULL, max.mtry = NULL,
                          mtry.step = 1, num.trees = 500,
                          cv.method = c("oob", "repeatedcv", "cv"),
                          cv.folds = 10, cv.repeats = 5,
                          num.threads = NULL,
                          plot.it = TRUE, verbose = TRUE, ...) {

  # ---- Argument validation ----
  cv.method <- match.arg(cv.method)
  if (!is.data.frame(dataset))
    stop("'dataset' must be a data frame.", call. = FALSE)
  if (!is.numeric(mtry.step) || mtry.step <= 0)
    stop("'mtry.step' must be a positive number.", call. = FALSE)
  if (!is.numeric(num.trees) || num.trees < 1L)
    stop("'num.trees' must be a positive integer.", call. = FALSE)
  if (!is.numeric(cv.folds) || cv.folds < 2L)
    stop("'cv.folds' must be at least 2.", call. = FALSE)
  if (!is.numeric(cv.repeats) || cv.repeats < 1L)
    stop("'cv.repeats' must be at least 1.", call. = FALSE)

  # ---- Build the candidate mtry grid ----
  f <- stats::formula(formula)
  RNames <- attr(stats::terms(f), "term.labels")
  ModelVarNo <- length(RNames)

  if (is.null(min.mtry)) min.mtry <- 1L
  if (is.null(max.mtry)) max.mtry <- ModelVarNo
  if (max.mtry < min.mtry)
    stop("'max.mtry' must be greater than or equal to 'min.mtry'.",
         call. = FALSE)
  if (max.mtry > ModelVarNo) {
    if (verbose)
      message("'max.mtry' (", max.mtry, ") exceeds the number of predictors ",
              "(", ModelVarNo, "); capping at ", ModelVarNo, ".")
    max.mtry <- ModelVarNo
  }

  mtry.grid <- seq.int(from = as.integer(min.mtry),
                       to   = as.integer(max.mtry),
                       by   = as.integer(mtry.step))

  # Response, used only for CV-based RMSE / R-squared.
  DepVarName <- row.names(attr(stats::terms(f), "factors"))[1]
  Y <- dataset[[DepVarName]]
  N <- nrow(dataset)

  # ---- Evaluate each mtry ----
  results <- data.frame(mtry     = mtry.grid,
                        RMSE     = NA_real_,
                        Rsquared = NA_real_,
                        SDRMSE   = NA_real_,
                        SDRsq    = NA_real_)

  for (i in seq_along(mtry.grid)) {
    m <- mtry.grid[i]

    if (cv.method == "oob") {

      # Single ranger fit; ranger reports OOB MSE and OOB R-squared directly.
      fit <- ranger::ranger(formula = f, data = dataset,
                            num.trees = num.trees, mtry = m,
                            num.threads = num.threads, ...)
      results$RMSE[i]     <- sqrt(fit$prediction.error)
      results$Rsquared[i] <- fit$r.squared

    } else {

      # Manual (repeated) k-fold cross-validation. "cv" is just the special
      # case cv.repeats = 1 to keep the code path single.
      reps <- if (cv.method == "cv") 1L else as.integer(cv.repeats)

      rmse_runs <- numeric(reps * cv.folds)
      r2_runs   <- numeric(reps * cv.folds)
      k <- 1L

      for (rep_i in seq_len(reps)) {
        # Stratification is not needed for a continuous response; balanced
        # fold sizes are obtained by recycling 1:cv.folds across N rows.
        fold.id <- sample(rep(seq_len(cv.folds), length.out = N))

        for (fold in seq_len(cv.folds)) {
          test.idx  <- which(fold.id == fold)
          train.idx <- setdiff(seq_len(N), test.idx)

          fit <- ranger::ranger(formula = f,
                                data = dataset[train.idx, , drop = FALSE],
                                num.trees = num.trees, mtry = m,
                                num.threads = num.threads, ...)
          pred <- stats::predict(fit, dataset[test.idx, , drop = FALSE],
                                 num.threads = num.threads)$predictions

          rmse_runs[k] <- .rmse(pred, Y[test.idx])
          r2_runs[k]   <- .cor_r2(pred, Y[test.idx])
          k <- k + 1L
        }
      }

      results$RMSE[i]     <- mean(rmse_runs)
      results$Rsquared[i] <- mean(r2_runs, na.rm = TRUE)
      results$SDRMSE[i]   <- stats::sd(rmse_runs)
      results$SDRsq[i]    <- stats::sd(r2_runs, na.rm = TRUE)
    }

    if (verbose)
      message("mtry = ", m,
              " | RMSE = ", round(results$RMSE[i], 4),
              " | R-squared = ", round(results$Rsquared[i], 4))
  }

  # ---- Best mtry by lowest RMSE ----
  best.idx  <- which.min(results$RMSE)
  best.mtry <- results$mtry[best.idx]
  if (verbose)
    message("Best mtry: ", best.mtry,
            " (RMSE = ", round(results$RMSE[best.idx], 4), ")")

  # ---- Optional plot ----
  if (isTRUE(plot.it)) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    graphics::plot(results$mtry, results$RMSE, type = "b", pch = 19,
                   xlab = "mtry", ylab = "RMSE",
                   main = paste0("rf.mtry.optim (", cv.method, ")"))
    graphics::abline(v = best.mtry, lty = 2, col = "red")
  }

  list(best.mtry = best.mtry,
       results   = results,
       cv.method = cv.method,
       num.trees = num.trees,
       call      = match.call())
}
