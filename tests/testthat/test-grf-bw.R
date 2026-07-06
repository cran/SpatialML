test_that("grf.bw returns the expected list shape", {
  df <- .toy_data()
  set.seed(1)
  # 100 trees and bandwidths >= 10 keep the stochastic local fits stable
  # while still exercising a real bandwidth grid.
  res <- grf.bw(dep ~ X1 + X2, dataset = df, kernel = "adaptive",
                coords = df[, c("X", "Y")],
                bw.min = 10, bw.max = 12, step = 1,
                trees = 100, mtry = 1, nthreads = 1L,
                forests = FALSE, verbose = FALSE)

  expect_named(res, c("tested.bandwidths", "Best.BW"))
  expect_s3_class(res$tested.bandwidths, "data.frame")
  expect_named(res$tested.bandwidths,
               c("Bandwidth", "Local", "Mixed", "Low.Local"))
  expect_equal(nrow(res$tested.bandwidths), 3L)
  expect_true(res$Best.BW %in% res$tested.bandwidths$Bandwidth)
  expect_true(is.numeric(res$Best.BW))
})

test_that("grf.bw validates inputs", {
  df <- .toy_data()
  expect_error(
    grf.bw(dep ~ X1 + X2, dataset = "not a df", kernel = "adaptive",
           coords = df[, c("X", "Y")]),
    "data frame"
  )
  expect_error(
    grf.bw(dep ~ X1 + X2, dataset = df, kernel = "adaptive",
           coords = df[, c("X", "Y")],
           bw.min = 20, bw.max = 10),
    "greater than or equal"
  )
})
