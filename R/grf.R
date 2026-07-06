# =============================================================================
# Geographically Weighted Random Forest (GRF)
#
# For each observation i in space the function fits a Random Forest on the
# local neighbourhood of i (defined by the k nearest neighbours when the
# kernel is "adaptive", or by all observations within a fixed radius bw when
# the kernel is "fixed"). When geo.weighted = TRUE local observations are
# weighted by the bi-square kernel  w_ij = (1 - (d_ij / h)^2)^2,  where
# d_ij is the Euclidean distance between i and j and h is the local
# bandwidth (the largest distance retained in the neighbourhood).
#
# References:
#   Georganos et al. (2019) Geocarto International, doi:10.1080/10106049.2019.1595177
#   Georganos & Kalogirou (2022) ISPRS IJGI 11(9):471, doi:10.3390/ijgi11090471
# =============================================================================

grf <- function(formula, dframe, bw, kernel = c("adaptive", "fixed"),
                coords, ntree = 500, mtry = NULL, importance = "impurity",
                nthreads = NULL, forests = TRUE, geo.weighted = TRUE,
                print.results = TRUE, progress = print.results,
                oob.fallback = c("loo", "inbag"), ...) {

  # ---------------------------------------------------------------------------
  # 1. Argument validation
  # ---------------------------------------------------------------------------
  kernel <- match.arg(kernel)
  oob.fallback <- match.arg(oob.fallback)

  if (!is.data.frame(dframe))
    stop("'dframe' must be a data frame.", call. = FALSE)
  if (missing(coords) || is.null(coords))
    stop("'coords' must be supplied as a 2-column matrix or data frame.",
         call. = FALSE)
  if (NCOL(coords) < 2L)
    stop("'coords' must have at least two columns (X, Y).", call. = FALSE)
  if (NROW(coords) != nrow(dframe))
    stop("'coords' and 'dframe' must have the same number of rows.",
         call. = FALSE)
  if (missing(bw) || !is.numeric(bw) || length(bw) != 1L || bw <= 0)
    stop("'bw' must be a single positive numeric value.", call. = FALSE)
  if (!is.numeric(ntree) || ntree < 1L)
    stop("'ntree' must be a positive integer.", call. = FALSE)

  start.time <- Sys.time()

  # ---------------------------------------------------------------------------
  # 2. Parse the formula and pull the response from the data frame
  # ---------------------------------------------------------------------------
  f <- stats::formula(formula)
  RNames <- attr(stats::terms(f), "term.labels")
  ModelVarNo <- length(RNames)
  if (ModelVarNo < 1L)
    stop("The model formula must contain at least one predictor.",
         call. = FALSE)

  DepVarName <- row.names(attr(stats::terms(f), "factors"))[1]
  if (!(DepVarName %in% names(dframe)))
    stop("Dependent variable '", DepVarName, "' not found in 'dframe'.",
         call. = FALSE)
  Y <- dframe[[DepVarName]]

  K   <- ModelVarNo + 1L            # parameters used in AIC formula
  Obs <- nrow(dframe)

  # Default mtry follows the textbook regression heuristic floor(p / 3).
  if (is.null(mtry)) mtry <- max(floor(ModelVarNo / 3), 1L)

  # For an adaptive kernel, validate the neighbour count up front.
  Ne <- NA_integer_
  if (kernel == "adaptive") {
    if (bw < 2L || bw > Obs)
      stop("For an adaptive kernel, 'bw' (number of neighbours) must be ",
           "between 2 and the number of observations.", call. = FALSE)
    Ne <- as.integer(round(bw))
  }

  if (print.results) {
    message("\nNumber of Observations: ", Obs)
    message("Number of Independent Variables: ", ModelVarNo)
    if (kernel == "adaptive") {
      message("Kernel: Adaptive\nNeighbours: ", Ne)
    } else {
      message("Kernel: Fixed\nBandwidth: ", bw)
    }
  }

  # ---------------------------------------------------------------------------
  # 3. Global random forest (used as benchmark and for blended predictions)
  # ---------------------------------------------------------------------------
  Gl.Model <- ranger::ranger(formula = f, data = dframe, num.trees = ntree,
                             mtry = mtry, importance = importance,
                             num.threads = nthreads, ...)

  # In-sample predictions of the global model. These are *not* OOB; they are
  # used only for the optional blended prediction in predict.grf().
  Predict <- stats::predict(Gl.Model, dframe, num.threads = nthreads)
  yhat <- Predict$predictions

  if (print.results) {
    message("\n--------------- Global ML Model Summary ---------------\n")
    print(Gl.Model)

    message("\nImportance:\n")
    print(Gl.Model$variable.importance)

    g.RSS  <- sum((Y - yhat)^2)
    g.TSS  <- sum((Y - mean(Y))^2)
    g.r    <- 1 - (g.RSS / g.TSS)
    g.AIC  <- 2 * K + Obs * log(g.RSS / Obs)
    g.AICc <- g.AIC + ((2 * K * (K + 1)) / (Obs - K - 1))

    message("\nMean Square Error (Not OOB): ", round(g.RSS / Obs, 3))
    message("R-squared (Not OOB) %: ",        round(100 * g.r, 3))
    message("AIC (Not OOB): ",                round(g.AIC, 3))
    message("AICc (Not OOB): ",               round(g.AICc, 3))
  }

  # ---------------------------------------------------------------------------
  # 4. Pre-compute the full pairwise distance matrix
  #
  # Memory note: this is an N x N matrix of doubles (8 * N^2 bytes). For
  # N around 5000 this already costs 200 MB; very large datasets may need
  # an external nearest-neighbour search (e.g. RANN, FNN), but staying in
  # base R keeps the dependency footprint minimal.
  # ---------------------------------------------------------------------------
  Dij <- as.matrix(stats::dist(coords))

  # ---------------------------------------------------------------------------
  # 5. Allocate the result containers
  # ---------------------------------------------------------------------------
  LM_Forests <- if (isTRUE(forests)) vector("list", Obs) else NULL

  LM_LEst <- as.data.frame(stats::setNames(
    replicate(ModelVarNo, numeric(Obs), simplify = FALSE), RNames))

  LM_GofFit <- data.frame(y           = numeric(Obs),
                          LM_yfitOOB  = numeric(Obs),
                          LM_ResOOB   = numeric(Obs),
                          LM_yfitPred = numeric(Obs),
                          LM_ResPred  = numeric(Obs),
                          LM_MSE      = numeric(Obs),
                          LM_Rsq100   = numeric(Obs),
                          LPerm       = numeric(Obs))

  # Optional progress bar; only useful in interactive sessions.
  show.pb <- isTRUE(progress) && interactive() && Obs > 1L
  if (show.pb) {
    pb <- utils::txtProgressBar(min = 0, max = Obs, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  # ---------------------------------------------------------------------------
  # 6. Local random forests
  # ---------------------------------------------------------------------------
  for (m in seq_len(Obs)) {

    # 6a. Identify the local neighbourhood by index, *without* materialising
    # a per-iteration copy of dframe. The original implementation built
    # data.frame(dframe, DNeighbour = ...) on every iteration, which copied
    # the entire data frame N times.
    d   <- Dij[, m]
    # Keep the focal row first. This is important when two or more rows share
    # identical coordinates: ranger's training predictions are indexed in the
    # row order of SubSet, and the first prediction is treated as the focal one.
    ord <- c(m, setdiff(order(d), m))

    if (kernel == "adaptive") {
      keep <- ord[seq_len(Ne)]
      Kernel_H <- d[keep[Ne]]      # largest distance retained -> bandwidth
    } else {
      keep <- ord[d[ord] <= bw]
      Kernel_H <- bw
      if (length(keep) < (ModelVarNo + 1L))
        stop("Not enough neighbours within the fixed bandwidth at observation ",
             m, ". Increase 'bw'.", call. = FALSE)
    }

    SubSet <- dframe[keep, , drop = FALSE]

    # 6b. Bi-square geographical weights, zero outside the bandwidth. If all
    # selected observations have the same coordinates as the focal one, the
    # bandwidth is zero and the kernel is undefined; equal weights are the
    # least surprising fallback for that degenerate geometry.
    if (!is.finite(Kernel_H) || Kernel_H <= 0) {
      Wts <- rep(1, length(keep))
    } else {
      Wts <- (1 - (d[keep] / Kernel_H)^2)^2
      Wts[!is.finite(Wts)] <- 0
    }

    # 6c. Fit the local forest. Wrapped in a closure so the weighted and
    # unweighted paths share argument plumbing.
    fit_local <- function(weights = NULL) {
      ranger::ranger(formula = f, data = SubSet, num.trees = ntree,
                     mtry = mtry, importance = importance,
                     case.weights = weights,
                     num.threads = nthreads, ...)
    }

    if (isTRUE(geo.weighted)) {
      Lcl.Model <- fit_local(Wts)
    } else {
      Lcl.Model <- fit_local(NULL)
    }
    counter <- 1L

    if (isTRUE(forests)) LM_Forests[[m]] <- Lcl.Model

    # 6d. Record per-observation diagnostics. Indexing the importance vector
    # by name (rather than by position) keeps the result correct even if
    # ranger ever changes the order of variable.importance.
    local.importance <- Lcl.Model$variable.importance
    if (is.null(local.importance) || length(local.importance) == 0L) {
      local.importance <- stats::setNames(rep(NA_real_, ModelVarNo), RNames)
    }
    LM_LEst[m, ] <- local.importance[RNames]

    # In-sample prediction of the focal observation (re-used as the OOB
    # fallback below).
    l.predict <- stats::predict(Lcl.Model, dframe[m, , drop = FALSE],
                                num.threads = nthreads)
    oob.fit  <- Lcl.Model$predictions[[1]]
    pred.fit <- l.predict$predictions

    # Defensive OOB fallback. ranger returns NaN when no tree leaves the
    # focal row out-of-bag, which can happen with very small local windows or
    # highly concentrated case weights. Instead of retrying and printing a
    # warning, the default fits one extra local model without the focal row
    # and predicts it out-of-sample. If that is not possible, use the
    # in-sample prediction as a last resort so summary statistics remain
    # finite.
    if (length(oob.fit) != 1L || !is.finite(oob.fit)) {
      counter <- 2L
      loo.fit <- NA_real_
      if (oob.fallback == "loo" && nrow(SubSet) > 1L) {
        loo.weights <- if (isTRUE(geo.weighted)) Wts[-1] else NULL
        if (!is.null(loo.weights)) {
          loo.weights[!is.finite(loo.weights)] <- 0
          if (!any(loo.weights > 0)) loo.weights <- NULL
        }
        loo.fit <- tryCatch({
          loo.model <- ranger::ranger(
            formula = f, data = SubSet[-1, , drop = FALSE],
            num.trees = ntree, mtry = mtry, importance = importance,
            case.weights = loo.weights, num.threads = nthreads, ...)
          stats::predict(loo.model, SubSet[1, , drop = FALSE],
                         num.threads = nthreads)$predictions[1]
        }, error = function(e) NA_real_)
      }
      oob.fit <- if (length(loo.fit) == 1L && is.finite(loo.fit)) {
        loo.fit
      } else {
        pred.fit
      }
    }

    LM_GofFit[m, "y"]           <- Y[m]
    LM_GofFit[m, "LM_yfitOOB"]  <- oob.fit
    LM_GofFit[m, "LM_ResOOB"]   <- LM_GofFit[m, "y"] - oob.fit
    LM_GofFit[m, "LM_yfitPred"] <- pred.fit
    LM_GofFit[m, "LM_ResPred"]  <- LM_GofFit[m, "y"] - pred.fit
    LM_GofFit[m, "LM_MSE"]      <- Lcl.Model$prediction.error
    LM_GofFit[m, "LM_Rsq100"]   <- Lcl.Model$r.squared
    LM_GofFit[m, "LPerm"]       <- counter

    if (show.pb) utils::setTxtProgressBar(pb, m)
  }

  # ---------------------------------------------------------------------------
  # 7. Assemble the output object
  # ---------------------------------------------------------------------------
  grf.out <- list(Global.Model              = Gl.Model,
                  Locations                 = coords,
                  Local.Variable.Importance = LM_LEst,
                  LGofFit                   = LM_GofFit)
  if (isTRUE(forests)) grf.out$Forests <- LM_Forests

  # Local variable-importance summary statistics
  lvi <- data.frame(
    Min  = apply(grf.out$Local.Variable.Importance, 2, min),
    Max  = apply(grf.out$Local.Variable.Importance, 2, max),
    Mean = apply(grf.out$Local.Variable.Importance, 2, mean),
    StD  = apply(grf.out$Local.Variable.Importance, 2, stats::sd))

  # Aggregate goodness-of-fit
  l.RSS.OOB  <- sum(grf.out$LGofFit$LM_ResOOB^2)
  l.RSS.Pred <- sum(grf.out$LGofFit$LM_ResPred^2)
  TSS        <- sum((grf.out$LGofFit$y - mean(grf.out$LGofFit$y))^2)

  l.r.OOB     <- 1 - (l.RSS.OOB  / TSS)
  l.AIC.OOB   <- 2 * K + Obs * log(l.RSS.OOB  / Obs)
  l.AICc.OOB  <- l.AIC.OOB  + ((2 * K * (K + 1)) / (Obs - K - 1))

  l.r.Pred    <- 1 - (l.RSS.Pred / TSS)
  l.AIC.Pred  <- 2 * K + Obs * log(l.RSS.Pred / Obs)
  l.AICc.Pred <- l.AIC.Pred + ((2 * K * (K + 1)) / (Obs - K - 1))

  if (print.results) {
    message("\n--------------- Local Model Summary ---------------\n")
    message("\nResiduals OOB:\n")
    print(summary(grf.out$LGofFit$LM_ResOOB))

    message("\nResiduals Predicted (Not OOB):\n")
    print(summary(grf.out$LGofFit$LM_ResPred))

    message("\nLocal Variable Importance:\n")
    print(lvi)

    message("\nMean squared error (OOB): ",       round(l.RSS.OOB / Obs, 3))
    message("R-squared (OOB) %: ",                round(100 * l.r.OOB, 3))
    message("AIC (OOB): ",                        round(l.AIC.OOB, 3))
    message("AICc (OOB): ",                       round(l.AICc.OOB, 3))
    message("Mean squared error Predicted (Not OOB): ",
            round(l.RSS.Pred / Obs, 3))
    message("R-squared Predicted (Not OOB) %: ",  round(100 * l.r.Pred, 3))
    message("AIC Predicted (Not OOB): ",          round(l.AIC.Pred, 3))
    message("AICc Predicted (Not OOB): ",         round(l.AICc.Pred, 3))
  }

  grf.out$LocalModelSummary <- list(
    l.VariableImportance = lvi,
    l.MSE.OOB            = l.RSS.OOB  / Obs,
    l.r.OOB              = l.r.OOB,
    l.MSE.Pred           = l.RSS.Pred / Obs,
    l.r.Pred             = l.r.Pred,
    l.AIC.OOB            = l.AIC.OOB,
    l.AICc.OOB           = l.AICc.OOB,
    l.AIC.Pred           = l.AIC.Pred,
    l.AICc.Pred          = l.AICc.Pred
  )

  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  if (print.results)
    message("\nCalculation time (in seconds): ", round(time.taken, 4))

  # The S3 class is what makes predict(grf_obj, ...) dispatch to predict.grf().
  class(grf.out) <- "grf"
  grf.out
}
