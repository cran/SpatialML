grf <- function(formula, dframe, bw, kernel, coords, ntree=500, mtry=NULL, importance="impurity", nthreads = NULL, forests = TRUE, weighted = TRUE,  print.results=TRUE, ...)
  {
    start.time <- Sys.time()

    f <- formula(formula)

    RNames <- attr(terms(f), "term.labels")

    DepVarName <- row.names(attr(terms(f), "factors"))[1]

    Y.DF <- dframe[DepVarName]

    Y <- Y.DF[[1]]

    ModelVarNo <- length(RNames)
    K = ModelVarNo +1

    ntrees <- ntree

    Obs <- nrow(dframe)

    if (is.null(mtry)) {mtry= max(floor(ModelVarNo/3), 1)}

    if(print.results) {
      message("\nNumber of Observations: ", Obs)
      message("Number of Independent Variables: ", ModelVarNo)
    }


    if(kernel == 'adaptive')
    {
      Ne <- bw
      if(print.results) {message("Kernel: Adaptive\nNeightbours: ", Ne)}
    }
    else
    {
      if(kernel == 'fixed')
      {
        if(print.results) {message("Kernel: Fixed\nBandwidth: ", bw)}
      }
    }


    Gl.Model <- eval(substitute(ranger(formula, data = dframe, num.trees=ntree, mtry= mtry, importance=importance, num.threads = nthreads, ...)))
    Predict <- predict(Gl.Model, dframe, num.threads = nthreads)
    yhat <- Predict$predictions

    if(print.results) {
      message("\n--------------- Global Model Summary ---------------\n")
      print(Gl.Model)

      message("\nImportance:\n")
      print(Gl.Model$variable.importance)


      g.RSS <- sum((Y-yhat)^2)
      g.mean.y <- mean(Y)
      g.TSS<-sum((Y-g.mean.y)^2)

      g.r<-1-(g.RSS/g.TSS)

      g.AIC <- 2*K + Obs*log(g.RSS/Obs)

      g.AICc <- g.AIC + ((2*K*(K +1)) / (Obs - K - 1))

      message("\nMean Square Error (Not OBB): ", round(g.RSS/Obs,3))
      message("R-squared (Not OBB) %: ", round(100 * g.r,3))
      message("AIC (Not OBB): ", round(g.AIC,3))
      message("AICc (Not OBB): ", round(g.AICc,3))
    }

    DistanceT <- dist(coords)
    Dij <- as.matrix(DistanceT)

    if (forests == TRUE) {LM_Forests <- as.list(rep(NA, length(ntrees)))}

      LM_LEst <- as.data.frame(setNames(replicate(ModelVarNo, numeric(0), simplify = F), RNames[1:ModelVarNo]))

      LM_GofFit <- data.frame(y=numeric(0), LM_yfitOOB=numeric(0), LM_ResOOB=numeric(0), LM_yfitPred=numeric(0),
                              LM_ResPred=numeric(0), LM_MSE=numeric(0), LM_Rsq100=numeric(0), LPerm=numeric(0))

      for(m in 1:Obs){

        #Get the data
        DNeighbour <- Dij[,m]
        DataSet <- data.frame(dframe, DNeighbour = DNeighbour)

        #Sort by distance
        DataSetSorted <- DataSet[order(DataSet$DNeighbour),]

        if(kernel == 'adaptive')
        {
          #Keep Nearest Neighbours
          SubSet <- DataSetSorted[1:Ne,]
          Kernel_H <- max(SubSet$DNeighbour)
        }
        else
        {
          if(kernel == 'fixed')
          {
            SubSet <- subset(DataSetSorted, DNeighbour <= bw)
            Kernel_H <- bw
          }
        }

        #Bi-square weights
        Wts <- (1-(SubSet$DNeighbour/Kernel_H)^2)^2

        #Calculate WLM
        if (weighted == TRUE) {
          Lcl.Model<-eval(substitute(ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, case.weights=Wts, num.threads = nthreads, ...)))

          local.predicted.y <- Lcl.Model$predictions[[1]]
          counter <- 1
          while (is.nan(local.predicted.y)) {
            Lcl.Model<-eval(substitute(ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, case.weights=Wts, num.threads = nthreads, ...)))
            local.predicted.y <- Lcl.Model$predictions[[1]]
            counter <- counter + 1
           }
        } else
        {
          Lcl.Model<-eval(substitute(ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, num.threads = nthreads, ...)))
          counter <- 1
        }



        if (forests == TRUE) {LM_Forests[[m]] <- Lcl.Model}

        #Store in table
        #Importance
        for (j in 1:ModelVarNo) {
          LM_LEst[m,j] <- Lcl.Model$variable.importance[j]
        }

    #Observed y
    LM_GofFit[m,1] <- Y[m]
    LM_GofFit[m,2] <- Lcl.Model$predictions[[1]]
    LM_GofFit[m,3] <- LM_GofFit[m,1] - LM_GofFit[m,2]
    l.predict <- predict(Lcl.Model, dframe[m,], num.threads = nthreads)
    LM_GofFit[m,4] <- l.predict$predictions
    LM_GofFit[m,5] <- LM_GofFit[m,1] - LM_GofFit[m,4]
    LM_GofFit[m,6] <- Lcl.Model$prediction.error
    LM_GofFit[m,7] <- Lcl.Model$r.squared
    LM_GofFit[m,8] <- counter
  }
  if (forests == TRUE) {grf.out <- list(Global.Model=Gl.Model, Locations = coords, Local.Variable.Importance = LM_LEst, LGofFit=LM_GofFit, Forests=LM_Forests)}
  else {grf.out <- list(Global.Model=Gl.Model, Locations = coords, Local.Variable.Importance = LM_LEst, LGofFit=LM_GofFit)}

  if(print.results) {

    message("\n--------------- Local Model Summary ---------------\n")

    message("\nResiduals OOB:\n")
    print(summary(grf.out$LGofFit$LM_ResOOB))

    message("\nResiduals Predicted (Not OBB):\n")

    print(summary(grf.out$LGofFit$LM_ResPred))

  }
    lvi <- data.frame(Min = apply(grf.out$Local.Variable.Importance, 2, min), Max = apply(grf.out$Local.Variable.Importance, 2, max),
                     Mean = apply(grf.out$Local.Variable.Importance, 2, mean), StD = apply(grf.out$Local.Variable.Importance, 2, sd))


    l.RSS.OOB <- sum(grf.out$LGofFit$LM_ResOOB^2)
    l.RSS.Pred<-sum(grf.out$LGofFit$LM_ResPred^2)

    mean.y<-mean(grf.out$LGofFit$y)
    TSS<-sum((grf.out$LGofFit$y-mean.y)^2)

    l.r.OOB<-1-(l.RSS.OOB/TSS)
    g.AIC.OOB <- 2*K + Obs*log(l.RSS.OOB/Obs)
    g.AICc.OOB <- g.AIC.OOB + ((2*K*(K + 1)) / (Obs - K - 1))



   l.r.Pred<-1-(l.RSS.Pred/TSS)
   g.AIC.Pred <- 2*K + Obs*log(l.RSS.Pred/Obs)
   g.AICc.Pred <- g.AIC.Pred + ((2*K*(K +1)) / (Obs - K - 1))

  if(print.results) {

    message("\nLocal Variable Importance:\n")
    print(lvi)
    message("\nMean squared error (OOB): ", round(l.RSS.OOB/Obs,3))
    message("R-squared (OOB) %: ", round(100* l.r.OOB,3))
    message("AIC (OBB): ", round(g.AIC.OOB,3))
    message("AICc (OBB): ", round(g.AICc.OOB,3))
    message("Mean squared error Predicted (Not OBB): ", round(l.RSS.Pred/Obs,3))
    message("R-squared Predicted (Not OBB) %: ", round(100* l.r.Pred,3))
    message("AIC Predicted (Not OBB): ", round(g.AIC.Pred,3))
    message("AICc Predicted (Not OBB): ", round(g.AICc.Pred,3))
  }

   lModelSummary = list()
   lModelSummary$l.VariableImportance <- lvi
   lModelSummary$l.MSE.OOB <- l.RSS.OOB/Obs
   lModelSummary$l.r.OOB <- l.r.OOB
   lModelSummary$l.MSE.Pred <- l.RSS.Pred/Obs
   lModelSummary$l.r.Pred <- l.r.Pred

   grf.out$LocalModelSummary <- lModelSummary

   end.time <- Sys.time()
   time.taken <- end.time - start.time

   if(print.results) {message("\nCalculation time (in seconds): ", round(time.taken,4))}

  return(grf.out)
}
