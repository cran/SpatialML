rf.mtry.optim <- function(formula, dataset, min.mtry=NULL, max.mtry=NULL, mtry.step=1, cv.method="repeatedcv", cv.folds=10, ...) {

  f <- formula(formula)
  RNames <- attr(terms(f), "term.labels")
  ModelVarNo <- length(RNames)

  if (is.null(min.mtry)) {min.mtry <- 1}
  if (is.null(max.mtry)) {max.mtry <-  ModelVarNo}

  if (cv.method == "repeatedcv") {
    control <- trainControl(cv.method, repeats=5, number=cv.folds, search="grid", ...)
  } else if (cv.method == "cv") {
    control <- trainControl(cv.method, number=cv.folds, search="grid", ...)
  } else {
    control <- trainControl(cv.method, number=cv.folds, search="grid", ...)
  }

  set.seed(123)

  tunegrid <- expand.grid(.mtry=seq(from=min.mtry, to=max.mtry, by=mtry.step))

  rf_gridsearch <- train(formula, data= dataset, method="rf", tuneGrid=tunegrid, trControl=control)

  print(rf_gridsearch)

  plot(rf_gridsearch)

  return(rf_gridsearch)
}


