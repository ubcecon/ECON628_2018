## Recreate Table 2 of QSW (2018)
if (!file.exists("qsw2018.Rdata")) {
  source("qsw-dataPrep.R",echo=TRUE)
}
load("qsw2018.Rdata")
library(lfe)
qsw$propre <- interaction(qsw$provincecheck, qsw$prefecturecheck)
summary(felm(pc1s ~ reform_cg2003 + reform_cg2003:commercial
             | newspaper_id + year | 0
             | propre, data=subset(qsw, !is.na(realfdi) & realfdi>0),
             na.action=na.exclude))



## try a lasso specification
# construct x matrix
m <- felm(pc1s ~ reform_cg2003 + reform_cg2003:commercial +
            poly(log(gdp), log(industrialshare), log(population),
                 log(realfdi), degree=5) +
            year:(log(gdp) + log(industrialshare) + log(population) +
                  log(realfdi)) +
            as.factor(year)
          | newspaper_id | 0 | propre,
          data=subset(qsw,realfdi>0), keepX=TRUE)
d <- m$X[,grep("reform_",colnames(m$X))]
x <- m$X[,!grepl("reform_",colnames(m$X))]
y <- subset(qsw,realfdi>0 & !is.na(pc1s))$pc1s
id <- as.factor(subset(qsw,realfdi>0 & !is.na(pc1s))$newspaper_id)
source("rlassoPanel.R")
# partial out x using lasso
eyx <- rlassoPanel(x, y=y, f=id)
ed1x <- rlassoPanel(x,y=d[,1],f=id)
ed2x <- rlassoPanel(x,y=d[,2],f=id)

# regress residuals of dependent variable on residuals of covariate(s)
# of interest
tmp <- data.frame(y=eyx$residuals, reform=ed1x$residuals,
                  reformC=ed2x$residuals, id=id)
summary(felm(y ~ reform + reformC | id | 0 | id, data=tmp))
