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

# There are many R packages for neural networks. I first tried neuralnet,
# but found it slow, then I tried RSNNS, and it seems fast
# enough. There might be a better option though
library(RSNNS)
df <- simulate(1000,p,mu.linear,m.linear)


partial.linear.net <- function(df, crossfits=5) {
  x.names <- names(df)[grep("x.",names(df))]
  fmla <- as.formula(paste(c("y ~ 1",x.names), collapse=" + "))
  n <- nrow(df)
  p <- ncol(df)
  n <- round(n*(1-1/pmax(1,crossfits)))
  rn <- floor(n^(1/(2*(1+1/(1+p)))))
  n <- nrow(df)
  ehat <- rep(NA, n)
  vhat <- rep(NA, n)
  if (crossfits>1) {
    ki <- sample(1:crossfits, n, replace=TRUE)
    for (k in 1:crossfits) {
      ynn <- mlp(df[ki!=k,x.names], df$y[ki!=k], size=c(rn),
                 linOut=TRUE, learnFunc="SCG")
      yhat <- predict(ynn, df[ki==k, x.names])
      ehat[ki==k] <- df$y[ki==k]-yhat
      dnn <-  mlp(df[ki!=k,x.names], df$d[ki!=k], size=c(rn),
                  linOut=TRUE, learnFunc="SCG")
      dhat <- predict(dnn, df[ki==k,x.names])
      vhat[ki==k] <- df$d[ki==k]-dhat
    }
  } else {
    ynn <- mlp(df[,x.names], df$y, size=c(rn), linOut=TRUE,
               learnFunc="SCG")
    yhat <- predict(ynn, )
    ehat <- df$y-yhat
    dnn <-  mlp(df[,x.names], df$d, size=c(rn), linOut=TRUE,
                learnFunc="SCG")
    dhat <- predict(dnn)
    vhat <- df$d-dhat
  }
  lm(ehat ~ vhat)
}

partial.linear.lasso <- function(df) {
  x.names <- names(df)[grep("x.",names(df))]
  fmla <- as.formula(paste(c("y ~ d",x.names), collapse=" + "))
  rlassoEffects(fmla, data=df, I = ~ d)
}

#summary(partial.linear.lasso(df))

# simulate a bunch of times in parallel
simulations <- 500 # number of simulations
library(parallel)
cl <- makeCluster(pmax(detectCores()-2,1))  # change as you see fit
clusterEvalQ(cl,library(hdm))
clusterEvalQ(cl,library(RSNNS))

# R Socket cluster spawns new R sessions with empty environments, we
# need to make sure they load any needed libraries and have access to
# things from the main environment that they use
design <- c("linear","step")
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
      clusterExport(cl,c("simulate",
                         "partial.linear.net","p","mu","m"))
#                         "partial.linear.split.rf", "n.save.split",
#                         "m.hat.rf","mu.hat.rf","rf.tuneOnce"))

      thetas <- parSapply(cl, rep(n,simulations), function(n)
      {
        df <- simulate(n, p, mu, m)
        x.names <- names(df)[grep("x.",names(df))]
        fmla <- as.formula(paste(c("y ~ d",x.names), collapse=" + "))
        c(lm(fmla,data=df)$coef[2],
          partial.linear.net(df,crossfits=5)$coef[2],
          partial.linear.net(df,crossfits=0)$coef[2])
          #partial.linear.split.rf(df)$coef[2],
          #partial.linear.lasso(df)$coefficients)
      }
      )
      tmp <- (data.frame(t(thetas)) - 1)*sqrt(n)
      names(tmp) <- c("OLS","Neural Network",
                      "Neural Network (no cross fitting)")
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
save(sim.df, file="partialLinearSimNet.Rdata")



