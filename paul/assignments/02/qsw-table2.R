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

