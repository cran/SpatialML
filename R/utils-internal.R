# Internal helpers (not exported).

# Squared correlation between predictions and observations, equivalent to
# the "Rsquared" element of caret::postResample(). Pairwise-complete to
# tolerate NA. Returns NA if fewer than two paired observations are usable.
.cor_r2 <- function(pred, obs) {
  ok <- !is.na(pred) & !is.na(obs)
  if (sum(ok) < 2L) return(NA_real_)
  r <- suppressWarnings(stats::cor(pred[ok], obs[ok]))
  if (is.na(r)) NA_real_ else r * r
}

# Root mean squared error.
.rmse <- function(pred, obs) {
  ok <- !is.na(pred) & !is.na(obs)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((pred[ok] - obs[ok])^2))
}
