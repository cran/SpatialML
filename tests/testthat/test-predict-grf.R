test_that("predict() dispatches to predict.grf and returns a numeric vector", {
  m <- .toy_grf()

  set.seed(2)
  newd <- random.test.data(nrows = 2, ncols = 3, vars.no = 3)
  p <- predict(m, newd, x.var.name = "X", y.var.name = "Y")

  expect_type(p, "double")
  expect_length(p, 6L)
  expect_false(any(is.na(p)))
})

test_that("predict.grf errors when the model was fitted with forests = FALSE", {
  m <- .toy_grf(forests = FALSE)
  set.seed(2)
  newd <- random.test.data(nrows = 2, ncols = 2, vars.no = 3)
  expect_error(
    predict(m, newd, x.var.name = "X", y.var.name = "Y"),
    "forests = TRUE"
  )
})

test_that("predict.grf errors on a missing coordinate column", {
  m <- .toy_grf()
  set.seed(2)
  newd <- random.test.data(nrows = 2, ncols = 2, vars.no = 3)
  expect_error(
    predict(m, newd, x.var.name = "Lon", y.var.name = "Y"),
    "Lon"
  )
})

test_that("predict.grf blends local and global predictions correctly", {
  m <- .toy_grf()
  set.seed(2)
  newd <- random.test.data(nrows = 2, ncols = 2, vars.no = 3)

  p_local  <- predict(m, newd, x.var.name = "X", y.var.name = "Y",
                      local.w = 1, global.w = 0)
  p_global <- predict(m, newd, x.var.name = "X", y.var.name = "Y",
                      local.w = 0, global.w = 1)
  p_blend  <- predict(m, newd, x.var.name = "X", y.var.name = "Y",
                      local.w = 0.5, global.w = 0.5)

  # Linearity of the blend: 0.5 * local + 0.5 * global
  expect_equal(p_blend, 0.5 * p_local + 0.5 * p_global, tolerance = 1e-8)
})
