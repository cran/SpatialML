## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4,
  warning  = FALSE,
  message  = FALSE
)

## -----------------------------------------------------------------------------
library(SpatialML)

## -----------------------------------------------------------------------------
set.seed(42)
df <- random.test.data(nrows = 8, ncols = 8, vars.no = 3)
head(df)

## -----------------------------------------------------------------------------
set.seed(1)
mtry.opt <- rf.mtry.optim(dep ~ X1 + X2, dataset = df,
                          cv.method = "oob", plot.it = FALSE,
                          verbose = FALSE)
mtry.opt$best.mtry
mtry.opt$results

## ----results = "hide"---------------------------------------------------------
set.seed(1)
bw.search <- grf.bw(dep ~ X1 + X2, dataset = df, kernel = "adaptive",
                    coords = df[, c("X", "Y")],
                    bw.min = 8, bw.max = 18, step = 2,
                    trees = 200, mtry = mtry.opt$best.mtry,
                    verbose = FALSE)

## -----------------------------------------------------------------------------
bw.search$tested.bandwidths
bw.search$Best.BW

## ----results = "hide"---------------------------------------------------------
set.seed(1)
m <- grf(dep ~ X1 + X2, dframe = df, bw = bw.search$Best.BW,
         kernel = "adaptive", coords = df[, c("X", "Y")],
         ntree = 200, mtry = mtry.opt$best.mtry,
         forests = TRUE, print.results = FALSE, progress = FALSE)

## -----------------------------------------------------------------------------
class(m)
m$LocalModelSummary$l.r.OOB              # local OOB R-squared
summary(m$LGofFit$LM_ResOOB)             # local OOB residuals
head(m$Local.Variable.Importance)        # local importance per observation

## ----fig.alt = "Local importance of predictor X1 across the synthetic 8 x 8 grid."----
imp <- m$Local.Variable.Importance$X1
plot(df$X, df$Y, pch = 19, cex = 2,
     col = grDevices::grey(1 - imp / max(imp)),
     xlab = "X", ylab = "Y",
     main = "Local importance of X1 (darker = more important)")

## -----------------------------------------------------------------------------
new.df <- random.test.data(2, 2, vars.no = 3)
predict(m, new.df, x.var.name = "X", y.var.name = "Y")

## -----------------------------------------------------------------------------
predict(m, new.df, x.var.name = "X", y.var.name = "Y",
        local.w = 0.5, global.w = 0.5)

## ----eval = FALSE-------------------------------------------------------------
# data(Income)
# Coords <- Income[, 1:2]
# 
# # 1. Search the optimal bandwidth (be patient)
# bw <- grf.bw(Income01 ~ UnemrT01 + PrSect01,
#              dataset = Income, kernel = "adaptive",
#              coords = Coords, bw.min = 30, bw.max = 80, step = 5)
# 
# # 2. Fit the final model
# m <- grf(Income01 ~ UnemrT01 + PrSect01, dframe = Income,
#          bw = bw$Best.BW, kernel = "adaptive", coords = Coords)
# 
# # 3. Local R-squared
# m$LocalModelSummary$l.r.OOB
# 
# # 4. Map the residuals
# plot(Coords[, 1], Coords[, 2], pch = 19,
#      col = ifelse(m$LGofFit$LM_ResOOB > 0, "red", "blue"),
#      xlab = "X", ylab = "Y", main = "GRF OOB residuals")

