# Tiny fixtures shared by the smoke tests. The data set is intentionally
# small (5 x 5 grid -> 25 observations) and the random forests are fitted
# with only 50 trees on a single thread to keep the whole test suite well
# under CRAN's typical check-time budget.

.ensure_spatialml_loaded <- function() {
  if (exists("random.test.data", mode = "function") &&
      exists("grf", mode = "function")) {
    return(invisible(TRUE))
  }

  candidates <- normalizePath(
    c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "..", "..")),
    mustWork = FALSE
  )
  is_pkg_root <- file.exists(file.path(candidates, "DESCRIPTION")) &
    dir.exists(file.path(candidates, "R"))
  pkg_root <- candidates[is_pkg_root][1]
  if (is.na(pkg_root)) {
    stop("Could not find the SpatialML package root. Run tests from the ",
         "package root or use testthat::test_dir('tests/testthat').",
         call. = FALSE)
  }

  r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$",
                        full.names = TRUE)
  for (r_file in r_files) sys.source(r_file, envir = globalenv())
  invisible(TRUE)
}

.ensure_spatialml_loaded()

.toy_data <- function(seed = 1, nrows = 5, ncols = 5, vars.no = 3) {
  set.seed(seed)
  random.test.data(nrows = nrows, ncols = ncols, vars.no = vars.no)
}

.toy_grf <- function(seed = 1, ..., bw = 8, kernel = "adaptive",
                     forests = TRUE, geo.weighted = TRUE) {
  df <- .toy_data(seed)
  set.seed(seed)
  grf(dep ~ X1 + X2, dframe = df, bw = bw, kernel = kernel,
      coords = df[, c("X", "Y")],
      ntree = 50, mtry = 1, nthreads = 1L,
      forests = forests, geo.weighted = geo.weighted,
      print.results = FALSE, progress = FALSE, ...)
}
