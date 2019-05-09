grf <- function(formula, dframe, bw, kernel, coords, ntree=500, mtry=NULL, importance=TRUE, forests = TRUE)
  {
    f <- formula(formula)
    RNames <- attr(terms(f), "term.labels")
    ModelVarNo <- length(RNames)
    ntrees <- ntree
    Obs <- nrow(dframe)
    if (is.null(mtry)) {mtry= max(floor(ModelVarNo/3), 1)}

    message("\nNumber of Observations: ", Obs)
    message("Number of Independent Variables: ", ModelVarNo)

    if(kernel == 'adaptive')
    {
      Ne <- bw
      message("Kernel: Adaptive\nNeightbours: ", Ne)
    }
    else
    {
      if(kernel == 'fixed')
      {
        message("Kernel: Fixed\nBandwidth: ", bw)
      }
    }

    Gl.Model <- eval(substitute(randomForest(formula, data = dframe, ntree=ntree, mtry= mtry, importance=importance)))
    yhat<-predict(Gl.Model, dframe)

    message("Number of Variables: ", ModelVarNo)
    message("\n--------------- Global Model Summary ---------------\n")

    print(Gl.Model)

    message("\nImportance:\n")
    print(importance(Gl.Model))

    g.RSS<-sum((Gl.Model$y-yhat)^2)
    g.mean.y<-mean(Gl.Model$y)
    g.TSS<-sum((Gl.Model$y-g.mean.y)^2)

    g.r<-1-(g.RSS/g.TSS)

    message("\nResidual Sum of Squares (Predicted): ", round(g.RSS,3))
    message("R-squared (Predicted) %: ", round(100 * g.r,3))

    DistanceT <- dist(coords)
    Dij <- as.matrix(DistanceT)


    if (forests == TRUE) {LM_Forests <- as.list(rep(NA, length(ntrees)))}

      LM_LEst1 <- as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), RNames[1:ModelVarNo]))
      LM_LEst2 <- as.data.frame(setNames(replicate(ModelVarNo,numeric(0), simplify = F), RNames[1:ModelVarNo]))

      LM_GofFit <- data.frame(y=numeric(0), LM_yfitOOB=numeric(0), LM_ResOOB=numeric(0), LM_yfitPred=numeric(0), LM_ResPred=numeric(0), LM_MSR=numeric(0), LM_Rsq100=numeric(0))

      for(m in 1:Obs){

        #Get the data
        DNeighbour <- Dij[,m]
        DataSet <- data.frame(dframe, DNeighbour=DNeighbour)

        #Sort by distance
        DataSetSorted<- DataSet[order(DataSet$DNeighbour),]

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
        Wts<-(1-(SubSet$DNeighbour/Kernel_H)^2)^2

        #Calculate WLM
        Lcl.Model<-eval(substitute(randomForest(formula, data = SubSet, ntree=ntree, mtry= mtry, importance=importance)))

        if (forests == TRUE) {LM_Forests[[m]]<-Lcl.Model}

        #Store in table
        #Importance
        for (j in 1:ModelVarNo) {
          LM_LEst1[m,j] <- importance(Lcl.Model)[j,1]
          LM_LEst2[m,j] <- importance(Lcl.Model)[j,2]
        }

    #Observed y
    LM_GofFit[m,1]<-Gl.Model$y[m]
    LM_GofFit[m,2]<-Lcl.Model$predicted[[1]]
    LM_GofFit[m,3]<-LM_GofFit[m,1] - LM_GofFit[m,2]
    LM_GofFit[m,4]<-predict(Lcl.Model, dframe[m,])
    LM_GofFit[m,5]<-LM_GofFit[m,1] - LM_GofFit[m,4]
    LM_GofFit[m,6]<-Lcl.Model$mse[ntrees]
    LM_GofFit[m,7]<-100 * Lcl.Model$rsq[ntrees]

  }
  if (forests == TRUE) {grf.out<-list(Global.Model=Gl.Model, Locations = coords, Local.Pc.IncMSE= LM_LEst1, Local.IncNodePurity= LM_LEst2, LGofFit=LM_GofFit,Forests=LM_Forests)}
  else {grf.out<-list(Global.Model=Gl.Model, Locations = coords, Local.Pc.IncMSE= LM_LEst1, Local.IncNodePurity= LM_LEst2, LGofFit=LM_GofFit)}

   message("\n--------------- Local Model Summary ---------------\n")

   message("\nResiduals OOB:\n")
   print(summary(grf.out$LGofFit$LM_ResOOB))

   message("\nResiduals Predicted:\n")
   print(summary(grf.out$LGofFit$LM_ResPred))

   t1 <- data.frame(Min = apply(grf.out$Local.Pc.IncMSE, 2, min), Max = apply(grf.out$Local.Pc.IncMSE, 2, max),
                     Mean = apply(grf.out$Local.Pc.IncMSE, 2, mean), StD = apply(grf.out$Local.Pc.IncMSE, 2, sd))

   message("\n%IncMSE:\n")
   print(t1)

   t2 <- data.frame(Min = apply(grf.out$Local.IncNodePurity, 2, min), Max = apply(grf.out$Local.IncNodePurity, 2, max),
                     Mean = apply(grf.out$Local.IncNodePurity, 2, mean), StD = apply(grf.out$Local.IncNodePurity, 2, sd))

   message("\n%IncNodePurity: \n")
   print(t2)

   l.RSS.OOB<-sum(grf.out$LGofFit$LM_ResOOB^2)
   l.RSS.Pred<-sum(grf.out$LGofFit$LM_ResPred^2)

   mean.y<-mean(grf.out$LGofFit$y)
   TSS<-sum((grf.out$LGofFit$y-mean.y)^2)


   l.r.OOB<-1-(l.RSS.OOB/TSS)
   message("\nResidual Sum of Squares (OOB): ", round(l.RSS.OOB,3))
   message("R-squared (OOB) %: ", round(100* l.r.OOB,3))

   l.r.Pred<-1-(l.RSS.Pred/TSS)
   message("Residual Sum of Squares (Predicted): ", round(l.RSS.Pred,3))
   message("R-squared (Predicted) %: ", round(100* l.r.Pred,3))

   lModelSummary = list()
   lModelSummary$l.IncMSE<-t1
   lModelSummary$l.IncNodePurity<-t2
   lModelSummary$l.RSS.OOB<-l.RSS.OOB
   lModelSummary$l.r.OOB<-l.r.OOB
   lModelSummary$l.RSS.Pred<-l.RSS.Pred
   lModelSummary$l.r.Pred <-l.r.Pred


   grf.out$LocalModelSummary<-lModelSummary


  return(grf.out)
}
