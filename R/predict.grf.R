# =============================================================================
# Predict method for objects of class "grf"
#
# For each row in 'new.data':
#   1. compute the Euclidean distance to every training location stored in
#      object$Locations,
#   2. pick the local Random Forest fitted at the closest training location,
#   3. blend its prediction with the global Random Forest prediction
#      according to the weights local.w and global.w.
#
# 'object$Forests' must be present (i.e. grf() must have been called with
# forests = TRUE), otherwise the local forests are not available and the
# function stops with a clear error message.
# =============================================================================

predict.grf <- function(object, new.data, x.var.name, y.var.name,
                        local.w = 1, global.w = 0, ...) {

  # ---- Validate inputs ----
  if (!inherits(object, "grf"))
    stop("'object' must be of class 'grf'.", call. = FALSE)
  if (is.null(object$Forests))
    stop("'object' does not contain local forests; rerun grf() with ",
         "forests = TRUE.", call. = FALSE)
  if (!is.data.frame(new.data))
    stop("'new.data' must be a data frame.", call. = FALSE)
  if (!(x.var.name %in% names(new.data)))
    stop("Column '", x.var.name, "' not found in 'new.data'.", call. = FALSE)
  if (!(y.var.name %in% names(new.data)))
    stop("Column '", y.var.name, "' not found in 'new.data'.", call. = FALSE)
  if (!is.numeric(local.w) || !is.numeric(global.w))
    stop("'local.w' and 'global.w' must be numeric.", call. = FALSE)

  Obs <- nrow(new.data)
  predictions <- numeric(Obs)

  # Coerce coordinates once outside the loop. as.matrix() copes with
  # data.frame, matrix and tibble inputs.
  locations <- as.matrix(object$Locations)

  for (i in seq_len(Obs)) {

    x <- new.data[[x.var.name]][i]
    y <- new.data[[y.var.name]][i]

    # Euclidean distance from the new point to every training location.
    D <- sqrt((x - locations[, 1])^2 + (y - locations[, 2])^2)
    local.model.ID <- which.min(D)

    g.predict <- stats::predict(object$Global.Model,
                                new.data[i, , drop = FALSE], ...)
    g.prediction <- g.predict$predictions

    l.predict <- stats::predict(object$Forests[[local.model.ID]],
                                new.data[i, , drop = FALSE], ...)
    l.prediction <- l.predict$predictions

    # Linear blend; with the defaults (local.w = 1, global.w = 0) the result
    # is the local prediction only.
    predictions[i] <- global.w * g.prediction[1] + local.w * l.prediction[1]
  }

  predictions
}
