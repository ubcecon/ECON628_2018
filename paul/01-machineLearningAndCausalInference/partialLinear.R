rm(list=ls())
p <- 4 # number of x's
mu.linear <- function(x) x%*%rep(1,p)
m.linear <- function(x) x%*%rep(2,p)
mu.step <- function(x) (x>0.5)%*%rep(1,p)
m.step <- function(x) (x>0.5)%*%rep(2,p)
theta <- 1

simulate <- function(n,p,mu,m) {
  theta <- 1
  x <- matrix(runif(n*p), ncol=p)
  d <- m(x) + rnorm(n)
  y <- theta*d + mu(x) + rnorm(n)
  data.frame(y=y,d=d,x=x)
}

library(grf)
library(hdm)
df <- simulate(100,p,mu.linear,m.linear)
mrfparams <- NULL
murfparams <- NULL
n.save <- NULL
partial.linear.rf <- function(df) {
  x.names <- names(df)[grep("x.",names(df))]
  if (is.null(mrfparams) || n.save!=nrow(df)) {
    # to save time, we only tune once per cluster worker and data set
    # size
    cat("tuning")
    m.rf  <- regression_forest(df[,x.names], df$d, num.trees=1000,
                               tune.parameters=TRUE)
    mrfparams <<- m.rf$tunable.params
    mu.rf  <- regression_forest(df[,x.names], df$y, num.trees=1000,
                                tune.parameters=TRUE)
    n.save <<- nrow(df)
    murfparams <<- mu.rf$tunable.params
  } else {
    cat("not tuning")
    m.rf  <- regression_forest(df[,x.names], df$d, num.trees=200,
                               tune.parameters=FALSE,
                               min.node.size =
                                 as.numeric(mrfparams["min.node.size"]),
                               alpha = as.numeric(mrfparams["alpha"]),
                               imbalance.penalty=as.numeric(mrfparams["imbalance.penalty"]),
                               sample.fraction = as.numeric(mrfparams["sample.fraction"]),
                               mtry=as.numeric(mrfparams["mtry"]))
    mu.rf  <- regression_forest(df[,x.names], df$y, num.trees=200,
                                tune.parameters=FALSE,
                               min.node.size =
                                 as.numeric(murfparams["min.node.size"]),
                               alpha = as.numeric(murfparams["alpha"]),
                               imbalance.penalty=as.numeric(murfparams["imbalance.penalty"]),
                               sample.fraction = as.numeric(murfparams["sample.fraction"]),
                               mtry=as.numeric(murfparams["mtry"]))

  }
  vhat <- df$d - predict(m.rf)$predictions
  ehat <- df$y - predict(mu.rf)$predictions
  lm(ehat ~ vhat)
}

## Manual sample splitting --- this turns out to be unneccessary. The
## default behavior of predict.regression_forest is to return
## predictions on the training data using only trees that were not fit
## on each observation. In other words, regression_forest already does
## the sample splitting for us.
##
## rf.tuneOnce <- function(x.names, y.name) {
##   parms <- NULL
##   function(df) {
##     if (is.null(parms)) {
##       rf  <- regression_forest(df[,x.names], df[,y.name], num.trees=500,
##                                  tune.parameters=TRUE)
##       parms <<- rf$tunable.params
##       rf
##     } else {
##       rf <- regression_forest(df[,x.names], df[,y.name], num.trees=200,
##                               tune.parameters=FALSE,
##                               honesty=FALSE,
##                               min.node.size =
##                                 as.numeric(parms["min.node.size"]),
##                               alpha = as.numeric(parms["alpha"]),
##                               imbalance.penalty=as.numeric(parms["imbalance.penalty"]),
##                               sample.fraction = as.numeric(parms["sample.fraction"]),
##                               mtry=as.numeric(parms["mtry"]))
##     }
##   }
## }
## n.save.split <- NULL
## m.hat.rf <- NULL
## mu.hat.rf  <- NULL
## partial.linear.split.rf <- function(df , splits=3) {
##   x.names <- names(df)[grep("x.",names(df))]
##   if (is.null(n.save.split) || n.save.split != nrow(df)) {
##     n.save.split <<- nrow(df)
##     m.hat.rf <<- rf.tuneOnce(x.names,"d")
##     mu.hat.rf <<- rf.tuneOnce(x.names,"y")
##   }
##   df$group <- sample(1:splits, nrow(df), replace=TRUE)
##   vhat <- df$d
##   ehat <- df$y
##   for(g in 1:splits) {
##     sdf <- subset(df, group!=g)
##     m <- m.hat.rf(sdf)
##     mu <- mu.hat.rf(sdf)
##     vhat[df$group==g] <- df$d[df$group==g] -
##       predict(m, newx=df[df$group==g,x.names])$predictions
##     ehat[df$group==g] <- df$y[df$group==g] -
##       predict(mu, newx=df[df$group==g,x.names])$predictions
##   }
##   lm(ehat ~ vhat)
## }


partial.linear.lasso <- function(df) {
  x.names <- names(df)[grep("x.",names(df))]
  fmla <- as.formula(paste(c("y ~ d",x.names), collapse=" + "))
  rlassoEffects(fmla, data=df, I = ~ d)
}

#summary(partial.linear.lasso(df))

# simulate a bunch of times in parallel
simulations <- 500 # number of simulations
library(parallel)
cl <- makeCluster(detectCores()/2)  # change as you see fit
clusterEvalQ(cl,library(hdm))
clusterEvalQ(cl,library(grf))

# R Socket cluster spawns new R sessions with empty environments, we
# need to make sure they load any needed libraries and have access to
# things from the main environment that they use
design <- c("linear") #,"step")
sim.df <- data.frame()
start.time <- Sys.time()
for (d in design) {
  if (d=="linear") {
    m <- m.linear
    mu <- mu.linear
  } else {
    m <- m.step
    mu <- mu.step
  }
  for (p in c(2,4,6,8)) {
    for (n in c(100, 200, 400, 800, 1600)) {
      clusterExport(cl,c("simulate","partial.linear.lasso",
                         "partial.linear.rf","p","mu","m",
                         "mrfparams","murfparams", "n.save"))
#                         "partial.linear.split.rf", "n.save.split",
#                         "m.hat.rf","mu.hat.rf","rf.tuneOnce"))

      thetas <- parSapply(cl, rep(n,simulations), function(n)
      {
        df <- simulate(n, p, mu, m)
        x.names <- names(df)[grep("x.",names(df))]
        fmla <- as.formula(paste(c("y ~ d",x.names), collapse=" + "))
        c(lm(fmla,data=df)$coef[2],
          partial.linear.rf(df)$coef[2],
          #partial.linear.split.rf(df)$coef[2],
          partial.linear.lasso(df)$coefficients)
      }
      )
      tmp <- (data.frame(t(thetas)) - 1)*sqrt(n)
      names(tmp) <- c("OLS","Random.Forest","Lasso")
      tmp$n <- n
      tmp$p <- p
      tmp$design <- d
      sim.df <- rbind(sim.df, tmp)
      cat("finished sample size ",n,"\n")
      cat("Elapsed time ", Sys.time()-start.time,"\n")
    }
    cat("finished p = ",p,"\n")
  }
  cat("finished design = ", d,"\n")
}
stopCluster(cl)
save(sim.df, file="partialLinearSim.Rdata")



