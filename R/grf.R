# This function fits a geographical random forest model.
# Inputs:
# - formula: an object of class "formula" or one that can be coerced to that class
# - dframe: a data frame containing the variables in the model
# - bw: bandwidth, used for kernel density estimation
# - kernel: type of kernel to use ('adaptive' or 'fixed')
# - coords: coordinates for the geographical data
# - ntree: number of trees to grow in the forest
# - mtry: number of variables randomly sampled as candidates at each split
# - importance: type of importance measure ('impurity' or 'permutation')
# - nthreads: number of threads for parallel processing
# - forests: boolean indicating whether to save the local forests
# - geo.weighted: boolean indicating whether to use geographical weighting
# - print.results: boolean indicating whether to print the results
# - ...: additional arguments passed to the ranger function

grf <- function(formula, dframe, bw, kernel, coords, ntree=500, mtry=NULL, importance="impurity", nthreads = NULL, forests = TRUE, geo.weighted = TRUE,  print.results=TRUE, ...)
  {

    # Start timing the function execution
    start.time <- Sys.time()

    # Convert formula text to a formula object
    f <- formula(formula)

    # Extract variable names from the formula
    RNames <- attr(terms(f), "term.labels")

    # Get the name of the dependent variable
    DepVarName <- row.names(attr(terms(f), "factors"))[1]

    # Create a data frame for the dependent variable
    Y.DF <- dframe[DepVarName]

    # Convert the dependent variable data frame to a vector
    Y <- Y.DF[[1]]

    # Determine the number of independent variables and add 1 for degrees of freedom
    ModelVarNo <- length(RNames)
    K = ModelVarNo + 1

    # Set the number of trees in the model
    ntrees <- ntree

    # Count the number of observations in the data
    Obs <- nrow(dframe)

    # Define mtry if it is not provided [max(floor(Number of Variables/3), 1)]
    if (is.null(mtry)) {mtry= max(floor(ModelVarNo/3), 1)}

    # Print initial information if required
    if(print.results) {
      message("\nNumber of Observations: ", Obs)
      message("Number of Independent Variables: ", ModelVarNo)
    }

    # Configure the kernel type and its parameters
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

    # Fit the global random forest model using the ranger package
    Gl.Model <- eval(substitute(ranger(formula, data = dframe, num.trees=ntree, mtry= mtry, importance=importance, num.threads = nthreads, ...)))

    # Get predictions from the global model
    Predict <- predict(Gl.Model, dframe, num.threads = nthreads)

    yhat <- Predict$predictions

    # Print global model summary if required
    if(print.results) {
      message("\n--------------- Global ML Model Summary ---------------\n")
      print(Gl.Model)

      message("\nImportance:\n")
      print(Gl.Model$variable.importance)

      #calculate pseudoR2
      g.RSS <- sum((Y-yhat)^2)
      g.mean.y <- mean(Y)
      g.TSS<-sum((Y-g.mean.y)^2)

      g.r<-1-(g.RSS/g.TSS)

      g.AIC <- 2*K + Obs*log(g.RSS/Obs)

      g.AICc <- g.AIC + ((2*K*(K +1)) / (Obs - K - 1))

      message("\nMean Square Error (Not OOB): ", round(g.RSS/Obs,3))
      message("R-squared (Not OOB) %: ", round(100 * g.r,3))
      message("AIC (Not OOB): ", round(g.AIC,3))
      message("AICc (Not OOB): ", round(g.AICc,3))
    }

    # Calculate distances between observations based on coordinates
    DistanceT <- dist(coords)
    Dij <- as.matrix(DistanceT)

    # Initialize storage for local forests if required
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
        if (geo.weighted == TRUE) {
          Lcl.Model <- eval(substitute(ranger(formula, data = SubSet, num.trees=ntree, mtry= mtry, importance=importance, case.weights=Wts, num.threads = nthreads, ...)))

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

  # Compile outputs from the function
  if (forests == TRUE) {grf.out <- list(Global.Model=Gl.Model, Locations = coords, Local.Variable.Importance = LM_LEst, LGofFit=LM_GofFit, Forests=LM_Forests)}
  else {grf.out <- list(Global.Model=Gl.Model, Locations = coords, Local.Variable.Importance = LM_LEst, LGofFit=LM_GofFit)}

  if(print.results) {

    message("\n--------------- Local Model Summary ---------------\n")

    message("\nResiduals OOB:\n")
    print(summary(grf.out$LGofFit$LM_ResOOB))

    message("\nResiduals Predicted (Not OOB):\n")

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
    message("AIC (OOB): ", round(g.AIC.OOB,3))
    message("AICc (OOB): ", round(g.AICc.OOB,3))
    message("Mean squared error Predicted (Not OOB): ", round(l.RSS.Pred/Obs,3))
    message("R-squared Predicted (Not OOB) %: ", round(100* l.r.Pred,3))
    message("AIC Predicted (Not OOB): ", round(g.AIC.Pred,3))
    message("AICc Predicted (Not OOB): ", round(g.AICc.Pred,3))
  }

   lModelSummary = list()
   lModelSummary$l.VariableImportance <- lvi
   lModelSummary$l.MSE.OOB <- l.RSS.OOB/Obs
   lModelSummary$l.r.OOB <- l.r.OOB
   lModelSummary$l.MSE.Pred <- l.RSS.Pred/Obs
   lModelSummary$l.r.Pred <- l.r.Pred

   grf.out$LocalModelSummary <- lModelSummary

   # Calculate and print the time taken to run the function
   end.time <- Sys.time()
   time.taken <- end.time - start.time

   if(print.results) {message("\nCalculation time (in seconds): ", round(time.taken,4))}

  # Return the output list
  return(grf.out)
}
