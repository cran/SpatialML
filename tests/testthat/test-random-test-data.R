test_that("random.test.data returns the expected shape with coordinates", {
  set.seed(1)
  df <- random.test.data(nrows = 5, ncols = 5, vars.no = 3)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 25L)
  expect_named(df, c("dep", "X1", "X2", "X", "Y"))
})

test_that("random.test.data without coordinates omits X and Y", {
  set.seed(1)
  df <- random.test.data(nrows = 3, ncols = 3, vars.no = 4,
                         xycoords = FALSE)
  expect_equal(nrow(df), 9L)
  expect_named(df, c("dep", "X1", "X2", "X3"))
})

test_that("random.test.data 'normal' draws come from rnorm", {
  set.seed(1)
  df <- random.test.data(nrows = 20, ncols = 20, vars.no = 2,
                         dep.var.dis = "normal")
  # rnorm has mean approx 0 and routinely produces negative values; the
  # earlier (buggy) version drew from runif and would fail both checks.
  expect_lt(abs(mean(df$dep)), 0.5)
  expect_true(any(df$dep < 0))
})

test_that("random.test.data 'poisson' draws are non-negative integers", {
  set.seed(1)
  df <- random.test.data(nrows = 5, ncols = 5, vars.no = 2,
                         dep.var.dis = "poisson")
  expect_true(all(df$dep >= 0))
  expect_true(all(df$dep == as.integer(df$dep)))
})

test_that("random.test.data validates inputs", {
  expect_error(random.test.data(vars.no = 1L), "vars.no")
  expect_error(random.test.data(nrows = 0L),   "nrows")
  expect_error(random.test.data(ncols = 0L),   "ncols")
  expect_error(random.test.data(dep.var.dis = "uniform"), "should be one of")
})
