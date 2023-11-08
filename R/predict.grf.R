#object = local Forests
predict.grf <- function(object, new.data, x.var.name, y.var.name, local.w=1, global.w=0, ...) {

  Obs <- nrow(new.data)

  predictions <- vector(mode="numeric", length=Obs)

  for(i in 1:Obs){

    x <- new.data[i, which(names(new.data)==x.var.name)]
    y <- new.data[i, which(names(new.data)==y.var.name)]

    locations <- object$Locations

    D <- sqrt((x-locations[,1])^2 + (y-locations[,2])^2)

    local.model.ID <- which.min(D)

    g.predict <- predict(object[[1]], new.data[i,], ...)
    g.prediction <- g.predict$predictions
    l.predict <- predict(object$Forests[[local.model.ID]], new.data[i,])
    l.prediction <- l.predict$predictions

    predictions[i] <- global.w * g.prediction[1] + local.w * l.prediction[1]
  }
return(predictions)
}
