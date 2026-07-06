# =============================================================================
# Random spatial test data generator
#
# Produces a data frame with one dependent variable (drawn from N(0, 1) or
# Poisson(7)), 'vars.no - 1' uniformly distributed predictors named
# X1, X2, ..., and optionally a regular nrows-by-ncols grid of (X, Y)
# coordinates.
#
# Mostly intended for the package examples and quick simulation studies.
# Call set.seed() before use if you need reproducible draws.
# =============================================================================

random.test.data <- function(nrows = 10, ncols = 10, vars.no = 3,
                             dep.var.dis = c("normal", "poisson"),
                             xycoords = TRUE) {

  dep.var.dis <- match.arg(dep.var.dis)
  if (!is.numeric(nrows)  || nrows  < 1L)
    stop("'nrows' must be a positive integer.",  call. = FALSE)
  if (!is.numeric(ncols)  || ncols  < 1L)
    stop("'ncols' must be a positive integer.",  call. = FALSE)
  if (!is.numeric(vars.no) || vars.no < 2L)
    stop("'vars.no' must be at least 2 (one dependent + one independent).",
         call. = FALSE)

  obs.no <- as.integer(nrows) * as.integer(ncols)

  # Dependent variable. Note: pre-1.7 the "normal" branch incorrectly used
  # runif(); fixed to rnorm().
  dep <- switch(dep.var.dis,
                normal  = stats::rnorm(obs.no),
                poisson = stats::rpois(obs.no, lambda = 7))

  # Predictor matrix, vectorised. Pre-1.7 this was a per-column loop with
  # an off-by-one bug (`for (i in 1:vars.no - 1)` evaluated to 0:(vars.no-1)).
  n.indep <- as.integer(vars.no) - 1L
  vars <- matrix(stats::runif(obs.no * n.indep),
                 nrow = obs.no, ncol = n.indep)
  colnames(vars) <- paste0("X", seq_len(n.indep))

  if (isTRUE(xycoords)) {
    X <- rep(seq_len(nrows), each = ncols)
    Y <- rep(seq_len(ncols), times = nrows)
    out <- data.frame(dep = dep, vars, X = X, Y = Y)
  } else {
    out <- data.frame(dep = dep, vars)
  }
  out
}
