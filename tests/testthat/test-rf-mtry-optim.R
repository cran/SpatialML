test_that("rf.mtry.optim with the OOB strategy returns the expected structure", {
  set.seed(1)
  df <- random.test.data(nrows = 5, ncols = 5, vars.no = 4)

  set.seed(1)
  res <- rf.mtry.optim(dep ~ X1 + X2 + X3, dataset = df,
                       num.trees = 50, num.threads = 1L,
                       cv.method = "oob",
                       plot.it = FALSE, verbose = FALSE)

  expect_named(res, c("best.mtry", "results", "cv.method",
                      "num.trees", "call"))
  expect_true(res$best.mtry >= 1L && res$best.mtry <= 3L)
  expect_s3_class(res$results, "data.frame")
  expect_named(res$results,
               c("mtry", "RMSE", "Rsquared", "SDRMSE", "SDRsq"))
  expect_equal(res$cv.method, "oob")
  expect_equal(res$num.trees, 50)
  # OOB does not perform resampling so the SDs are NA.
  expect_true(all(is.na(res$results$SDRMSE)))
  expect_true(all(is.na(res$results$SDRsq)))
})

test_that("rf.mtry.optim caps max.mtry at the number of predictors", {
  set.seed(1)
  df <- random.test.data(nrows = 5, ncols = 5, vars.no = 3)
  set.seed(1)
  res <- rf.mtry.optim(dep ~ X1 + X2, dataset = df,
                       max.mtry = 99, num.trees = 50,
                       num.threads = 1L, cv.method = "oob",
                       plot.it = FALSE, verbose = FALSE)
  expect_lte(max(res$results$mtry), 2L)
})

test_that("rf.mtry.optim validates inputs", {
  df <- .toy_data()
  expect_error(
    rf.mtry.optim(dep ~ X1 + X2, dataset = "not a df"),
    "data frame"
  )
  expect_error(
    rf.mtry.optim(dep ~ X1 + X2, dataset = df, cv.folds = 1),
    "cv.folds"
  )
  expect_error(
    rf.mtry.optim(dep ~ X1 + X2, dataset = df,
                  min.mtry = 5, max.mtry = 1),
    "greater than or equal"
  )
})

test_that("rf.mtry.optim 'cv' strategy runs (skipped on CRAN for time)", {
  skip_on_cran()

  set.seed(1)
  df <- random.test.data(nrows = 5, ncols = 5, vars.no = 3)

  set.seed(1)
  res <- rf.mtry.optim(dep ~ X1 + X2, dataset = df,
                       num.trees = 50, num.threads = 1L,
                       cv.method = "cv", cv.folds = 3,
                       plot.it = FALSE, verbose = FALSE)
  expect_equal(res$cv.method, "cv")
  expect_true(res$best.mtry %in% c(1L, 2L))
  # CV does report fold-to-fold standard deviations.
  expect_false(any(is.na(res$results$SDRMSE)))
})
