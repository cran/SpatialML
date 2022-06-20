rf.mtry.optim <- function(formula, dataset) {

  control <- trainControl(method="repeatedcv",repeats=5, number=10, search="grid")

  set.seed(123)

  tunegrid <- expand.grid(.mtry=c(1:7))

  rf_gridsearch <- train(formula, data= dataset, method="rf", tuneGrid=tunegrid, trControl=control)

  print(rf_gridsearch)

  plot(rf_gridsearch)

  return(rf_gridsearch)
}


