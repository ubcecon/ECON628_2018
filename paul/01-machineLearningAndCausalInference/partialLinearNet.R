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

library(RSNNS)
df <- simulate(1000,p,mu.linear,m.linear)


partial.linear.net <- function(df) {
  x.names <- names(df)[grep("x.",names(df))]
  fmla <- as.formula(paste(c("y ~ 1",x.names), collapse=" + "))
  n <- nrow(df)
  p <- ncol(df)
  rn <- floor(n^(1/(2*(1+1/(1+p)))))
  ynn <- mlp(df[,x.names], df$y, size=c(rn), linOut=TRUE)
  yhat <- predict(ynn)
  ehat <- df$y-yhat
  fmla <- as.formula(paste(c("d ~ 1",x.names), collapse=" + "))
  dnn <-  mlp(df[,x.names], df$d, size=c(rn), linOut=TRUE)
  dhat <- predict(dnn)
  vhat <- df$d-dhat
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
cl <- makeCluster(detectCores()/2)  # change as you see fit
clusterEvalQ(cl,library(hdm))
clusterEvalQ(cl,library(RSNNS))

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
          partial.linear.net(df)$coef[2])
          #partial.linear.split.rf(df)$coef[2],
          #partial.linear.lasso(df)$coefficients)
      }
      )
      tmp <- (data.frame(t(thetas)) - 1)*sqrt(n)
      names(tmp) <- c("OLS","Neural Network")
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



