grf.bw <- function(formula, dataset, kernel="adaptive", coords, bw.min = NULL, bw.max = NULL, step = 1, trees=500, mtry=NULL, importance="impurity", nthreads = 1, forests = FALSE, weighted = TRUE, ...) {

  Obs <- nrow(dataset)

  f <- formula(formula)
  RNames <- attr(terms(f), "term.labels")
  ModelVarNo <- length(RNames)

  DepVarName <- row.names(attr(terms(f), "factors"))[1]
  Y.DF <- dataset[DepVarName]
  Y <- Y.DF[[1]]



  if (is.null(bw.min)) {bw.min <- max(round(Obs*0.05,0), ModelVarNo+2, 20)}
  if (is.null(bw.max)) {bw.max <- max(round(Obs*0.95,0), ModelVarNo+2)}
  if (is.null(mtry)) {mtry= max(floor(ModelVarNo/3), 2)}

  #store goodness of fit statistics
  eval.bw.grf <- data.frame(Bandwidth=integer(),
                            Local=double(),
                            Mixed=double(),
                            Low.Local=double(),
                            stringsAsFactors=FALSE)
  set.seed(1234)
  count <- 1
  for(abw in seq(from= bw.min, to=bw.max, by=step)){

    eval.bw.grf[count,1] <- abw

    grf16.a <- eval(substitute(grf(formula, dframe=dataset, bw=abw, kernel, coords, ntree=trees, mtry = mtry, importance=importance, nthreads=nthreads, forests = FALSE, weighted = TRUE, print.results=FALSE, ...)))

    eval.bw.grf[count,2] <- grf16.a$LocalModelSummary$l.r.OOB

    message("Bandwidth: ", abw)
    message("R2 of Local Model: ", eval.bw.grf[count,2])

    metrics_bw_grf <- postResample(pred = (grf16.a$LGofFit$LM_yfitOOB + grf16.a$Global.Model$predictions)/2, obs = Y)
    eval.bw.grf[count,3] <- metrics_bw_grf[2]

    metrics_bw_grf <- postResample(pred = (grf16.a$LGofFit$LM_yfitOOB*0.25) +(grf16.a$Global.Model$predictions*0.75), obs = Y)
    eval.bw.grf[count,4] <- metrics_bw_grf[2]

    count <- count + 1
  }

  best.bw <- eval.bw.grf$Bandwidth[which(eval.bw.grf$Local == max(eval.bw.grf$Local))]

  message("Best Bandwidth (Based on the Local Model): ", best.bw)

  return(list(tested.bandwidths = eval.bw.grf, Best.BW = best.bw))
}




