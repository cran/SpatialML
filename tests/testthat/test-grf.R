test_that("grf returns an object of class 'grf' with the expected slots", {
  m <- .toy_grf()

  expect_s3_class(m, "grf")
  expect_s3_class(m$Global.Model, "ranger")

  # All documented slots must be present.
  expected_slots <- c("Global.Model", "Locations",
                      "Local.Variable.Importance", "LGofFit",
                      "Forests", "LocalModelSummary")
  expect_true(all(expected_slots %in% names(m)))

  # Local diagnostics are populated for every observation.
  expect_equal(nrow(m$LGofFit), 25L)
  expect_equal(nrow(m$Local.Variable.Importance), 25L)
  expect_length(m$Forests, 25L)

  # LGofFit has the documented columns.
  expect_named(m$LGofFit,
               c("y", "LM_yfitOOB", "LM_ResOOB", "LM_yfitPred",
                 "LM_ResPred", "LM_MSE", "LM_Rsq100", "LPerm"))

  # LocalModelSummary contains the four headline statistics.
  expect_true(all(c("l.MSE.OOB", "l.r.OOB",
                    "l.MSE.Pred", "l.r.Pred") %in%
                  names(m$LocalModelSummary)))
})

test_that("grf with forests = FALSE drops the Forests slot", {
  m <- .toy_grf(forests = FALSE)
  expect_null(m$Forests)
  expect_s3_class(m, "grf")
})

test_that("grf with the fixed kernel runs end-to-end", {
  df <- .toy_data()
  set.seed(1)
  m <- grf(dep ~ X1 + X2, dframe = df, bw = 3,
           kernel = "fixed", coords = df[, c("X", "Y")],
           ntree = 50, mtry = 1, nthreads = 1L,
           print.results = FALSE, progress = FALSE)
  expect_s3_class(m, "grf")
})

test_that("grf with geo.weighted = FALSE runs end-to-end", {
  m <- .toy_grf(geo.weighted = FALSE)
  expect_s3_class(m, "grf")
})

test_that("grf handles small weighted windows without noisy OOB warnings", {
  df <- .toy_data()
  set.seed(1)
  expect_warning(
    m <- grf(dep ~ X1 + X2, dframe = df, bw = 4,
             kernel = "adaptive", coords = df[, c("X", "Y")],
             ntree = 20, mtry = 1, nthreads = 1L,
             print.results = FALSE, progress = FALSE),
    NA
  )
  expect_s3_class(m, "grf")
  expect_false(any(!is.finite(m$LGofFit$LM_yfitOOB)))
})

test_that("grf validates inputs", {
  df <- .toy_data()
  coords <- df[, c("X", "Y")]

  expect_error(
    grf(dep ~ X1 + X2, dframe = "not a data frame", bw = 8,
        kernel = "adaptive", coords = coords),
    "data frame"
  )
  expect_error(
    grf(dep ~ X1 + X2, dframe = df, bw = -1,
        kernel = "adaptive", coords = coords),
    "positive numeric"
  )
  expect_error(
    grf(dep ~ X1 + X2, dframe = df, bw = 8,
        kernel = "adaptive", coords = coords[1:10, , drop = FALSE]),
    "same number of rows"
  )
  expect_error(
    grf(dep ~ X1 + X2, dframe = df, bw = 1000,
        kernel = "adaptive", coords = coords),
    "between 2 and"
  )
})
