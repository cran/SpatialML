#object = local Forests
predict.grf <- function(object, new.data, x.var.name, y.var.name, local.w=1, global.w=0,...) {

  Obs <- nrow(new.data)

  predictions <- vector(mode="numeric", length=Obs)

  for(i in 1:Obs){

    x <- new.data[i, which(names(new.data)==x.var.name)]
    y <- new.data[i, which(names(new.data)==y.var.name)]

    locations <- object$Locations

    D <- sqrt((x-locations$X)^2+(y-locations$Y)^2)

    local.model.ID <- which.min(D)

    g.prediction <- predict(object[[1]], new.data[i,])
    l.prediction <- predict(object$Forests[[local.model.ID]], new.data[i,])

    predictions[i] <- global.w * g.prediction[1] + local.w * l.prediction[1]
  }
return(predictions)
}
