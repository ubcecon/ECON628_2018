# Data prep work for Qin, Stromberg, & Wu (2018)
#
# This should replicate the output of their 'rootdir'/Data/Code/*.do files

## Download data
if (!dir.exists("qsw2018")) {
  if (!file.exists("qsw.zip")) {
    download.file("https://www.aeaweb.org/doi/10.1257/aer.20170947.data",destfile="qsw.zip")
  }
  unzip("qsw.zip")
  ## why does MAC's zip create useless __MACOSX directory and  files?
  unlink("__MACOSX", recursive=TRUE) ## delete that garbage
  if (Sys.info()['sysname'] != "Linux") error("Use Linux or change lines 13-15 of qsw-dataPrep.R")
  system("mv data qsw2018") # you'll have to change this on windows,
                            # probably okay on mac, but remove the
                            # error check above
}

################################################################################
##do `datadir'/directory8111_en.do
library(readstata13)
direc <- read.dta13("qsw2018/Data/In/directory_81_11_raw.dta")
direc$general <-direc$category_content %in% c("daily", "evening", "metro")
direc$general_mixed <- direc$general |
  grepl("mixed",direc$category_content)

direc$type_g <- direc$general
direc$type_pd <- direc$supervisor_type=="party" &
  direc$category_content=="daily"
direc$type_pe <- direc$supervisor_type %in%  c("party","party_division") &
  direc$category_content %in% c("evening","metro")
direc$type_s <- direc$supervisor_type=="parent newspaper" & direc$general
direc$type_o <- with(direc, !(type_pd | type_pe | type_s))

vars <-  names(direc)[grep("type_",names(direc))]
direc[direc$supervisor_type=="" | direc$category_content=="",vars] <- NA

direc$ntype <- with(direc,ifelse(type_pd,"Party Daily",
                          ifelse(type_pe,"Party Evening",
                          ifelse(type_s, "Subsidiary",
                                 "Other"))))
direc$general[direc$type_o] <- FALSE
direc$type_g <- direc$general
direc <- direc[with(direc,order(provincecheck, prefecturecheck, year,
                                newspaper_ch)), ]
dir.create("qsw2018/Data/Out")
save(direc,file= "qsw2018/Data/Out/directory8111_en.Rdata")
################################################################################


################################################################################
##do `datadir'/coverage.do
rm(list=ls())
#' Create lags for panel data.
#'
#' This function creates lags (or leads) of panel data variables.
#' Input data should be sorted by i, t --- e.g.
#' df <- df[order(i,t),]
#' @param x Vector or matrix to get lags of.
#' @param i unit index
#' @param t time index
#' @param lag How many time periods to lag. Can be negative if leading
#' values are desired.
#' @return Lagged copy of x.
panel.lag <- function(x, i, t, lag=1, consecutive.time=TRUE) {
  if (!identical(order(i,t),1:length(i))) {
    stop("inputs not sorted.")
  }
  if (is.matrix(x)) {
    return(apply(x,MARGIN=2,FUN=function(c) { panel.lag(c,i,t,lag) }))
  }
  if (length(i) != length(x) || length(i) != length(t) ) {
    stop("Inputs not same length")
  }
  if (lag>0) {
    x.lag <- x[1:(length(x)-lag)]
    x.lag[i[1:(length(i)-lag)]!=i[(1+lag):length(i)] ] <- NA
    if (consecutive.time) x.lag[(t[1:(length(i)-lag)]+lag)!=t[(1+lag):length(i)] ] <- NA
    val <- (c(rep(NA,lag),x.lag))
  } else if (lag<0) {
    lag <- abs(lag)
    x.lag <- x[(1+lag):length(x)]
    x.lag[i[1:(length(i)-lag)]!=i[(1+lag):length(i)] ] <- NA
    if (consecutive.time) x.lag[(t[1:(length(i)-lag)]+lag) != t[(1+lag):length(i)] ] <- NA
    val <- (c(x.lag,rep(NA,lag)))
  } else { # lag=0
    return (x)
  }
  if (class(x)=="Date" & class(val)=="numeric") {
    stopifnot(0==as.numeric(as.Date("1970-01-01")))
    val <- as.Date(val, origin="1970-01-01")
  }
  return(val)
}
load("qsw2018/Data/Out/directory8111_en.Rdata")
df <- subset(direc, wisenews==1 & ntype !="Other")
df <- df[order(df$newspaper_id, df$year),]
#/* some newspapers enter the directory one year after they enter wisenews */
n <- nrow(df)
newid <- c(TRUE, df$newspaper_id[2:n] != df$newspaper_id[1:(n-1)])
newobs <- df[newid,]
newobs$year <- newobs$year-1
tmp <- rbind(df, newobs)
tmp <- tmp[order(tmp$newspaper_id, tmp$year),]
tmp <- tmp[,c("newspaper_en","newspaper_ch","newspaper_id","year")]

cvg <- read.dta13("qsw2018/Data/In/coverage.dta")

vars <- c("leader", "xinhua", "nreportepoch", "nreportxinhua",
          "ntcorruption", "disasters", "accidents", "sports",
          "crime3", "entertainment", sprintf("toppos_a%d",1:5),
          sprintf("county_village%d",2:4))
for(v in vars) {
  nv <- paste("s_",v,sep="")
  cvg[,nv] <- cvg[,v]/cvg$tarticles * 100
  cvg[cvg[,nv]>100,nv] <- NA
}
cvg$bias3 <- cvg$nreportepoch/
  (cvg$nreportxinhua + cvg$nreportepoch)*100
cvg$bias3[cvg$year<=2001] <- NA
cvg <- cvg[order(cvg$newspaper_en, cvg$year),]
for (L in 1:3) {
  bias3.lead <- panel.lag(cvg$bias3, cvg$newspaper_en, cvg$year, -L, consecutive.time=FALSE)
  year.lead <- panel.lag(cvg$year, cvg$newspaper_en,
                         cvg$year, -L,
                         consecutive.time=FALSE)
  rpl <- is.na(cvg$bias3) & (cvg$year==(2002-L) | (!is.na(year.lead) & year.lead==2002))
  cvg$bias3[rpl] <- bias3.lead[rpl]
}

coverage <- merge(cvg, tmp, by=c("newspaper_en","year"), all=FALSE)
coverage <-
  coverage[with(coverage,order(newspaper_en,newspaper_ch,newspaper_id,year)),]
save(coverage, file="qsw2018/Data/Out/coverage.Rdata")
################################################################################

################################################################################
##do `datadir'/pref_stat.do
library(openxlsx)
gdp <- read.xlsx("qsw2018/Data/In/GDP_IMF.xlsx", sheet="Sheet2",
                 colNames=TRUE, rows=1:33, cols=1:3)
names(gdp)[1] <- "year"

ads <- read.xlsx("qsw2018/Data/In/advertisement newspaper over years.xlsx",
                 sheet="Sheet1",rows=1:38,cols=1:5, colNames=TRUE)
ads <- subset(ads, year<=2011)
tmp <- merge(ads, gdp, by="year")
x <- tmp$newspaper_ad_sales1*1e4 /(tmp$GDP_RMB_current * 1e9)*100
library(Hmisc)
interp <- approxExtrap(tmp$year[!is.na(x)], x[!is.na(x)],
                       xout=tmp$year[is.na(x)], method="linear")
xi <- x
xi[is.na(x)] <- interp$y
tmp$npads_gdp <- xi
tmp$deflator <- tmp$GDP_RMB_current/tmp$GDP_RMB_cons
tmp$deflator <- tmp$deflator/tmp$deflator[nrow(tmp)]
tmp <- tmp[,c("npads_gdp","deflator","year")]

ctl <- read.dta13("qsw2018/Data/In/controls84to10.dta")

ps <- merge(ctl, tmp, by="year")
ps$el_ads <- with(ps, log10(gdp*1e4*npads_gdp/deflator/100))
ps <- ps[with(ps,order(provincecheck,prefecturecheck,year)),]
save(ps, file="qsw2018/Data/Out/pref_stat.Rdata")
################################################################################

################################################################################
##do `datadir'/ntypes.do
#/* Construct prefecture-year grid*/
ctl <- read.dta13("qsw2018/Data/In/controls84to10.dta")
tmp <- ctl[,c("provincecheck","prefecturecheck","year")]
foo <- subset(tmp, year==1984)
for (y in 1981:1983) {
  foo$year <- y
  tmp <- rbind(tmp, foo)
}
tmp <- tmp[with(tmp, order(provincecheck,prefecturecheck,year)),]

#/* Count number of Party Dailies and Commercial papers by prefecture and year */
load("qsw2018/Data/Out/directory8111_en.Rdata")
df <- subset(direc, ntype!="Other")
df$type_es <- df$type_pe | df$type_s
apd <- aggregate(type_pd ~ provincecheck + prefecturecheck + year +
                    admin_rank, data=df, FUN=sum)
aes <- aggregate(type_es ~ provincecheck + prefecturecheck + year +
                   admin_rank, data=df, FUN=sum)

foo <- merge(apd, aes, by=c("provincecheck","prefecturecheck","year","admin_rank"))
cnts <- merge(foo, tmp, by=c("provincecheck","prefecturecheck","year"),
             all=TRUE)
cnts$admin_rank[is.na(cnts$admin_rank)] <- "central"
cnts$type_pd[is.na(cnts$type_df)] <- 0
cnts$type_es[is.na(cnts$type_es)] <- 0

# /* Generate market structure strings */
## Skipping and hoping not needed

# /* controls */
ntypes1 <- merge(cnts, ps,
                 by=c("provincecheck","prefecturecheck","year"),
                 all.x=TRUE)
ntypes1 <- ntypes1[!is.na(ntypes1$year),]

# /* add constant characteristics */
pscons <- read.dta13("qsw2018/Data/In/pref_stat_cons.dta", nonint.factors = FALSE)
ntypes1 <- merge(ntypes1, pscons,
                 by=c("provincecheck","prefecturecheck"), all.x=TRUE)
for(v in c("long_march","pref_84","treaty_ports","npapers")) {
  ntypes1[is.na(ntypes1[,v]),v] <- 0
}

save(ntypes1,file="qsw2018/Data/Out/ntypes_1.Rdata")
################################################################################

################################################################################
##do `datadir'/WN_analysis.do
load("qsw2018/Data/Out/coverage.Rdata")
load("qsw2018/Data/Out/directory8111_en.Rdata")
vars <- c("newspaper_id", "provincecheck",
          "prefecturecheck", "year", "provincecheck_en",
          "prefecturecheck_en", "supervisor_type", "category_content",
          "general", "ntype", "admin_rank",
          names(direc)[grep("type_",names(direc))])
wn <- merge(coverage, direc[,vars], by=c("newspaper_id","year"))
wn <- wn[wn$ntype!="Other",]
wn <- merge(wn, ps, by=c("provincecheck","prefecturecheck","year"))
wn <- merge(wn, pscons)
for(v in c("long_march","pref_84","treaty_ports","npapers")) {
  wn[is.na(wn[,v]),v] <- 0
}

wn <- wn[order(wn$newspaper_id, wn$year),]
wn <- wn[wn$year!=2011 & !is.na(wn$newspaper_id) & wn$ntype!="Other",]
save(wn, file="qsw2018/Data/Out/WN_analysis.Rdata")
################################################################################

################################################################################
##do `datadir'/newspaper_ad_year.do
adm <- read.dta13("qsw2018/Data/In/newspaper_ad_month.dta")
adyr <- aggregate(ad_value_month ~ newspaper_id + year, data=adm,
                  FUN=mean)
names(adyr)[names(adyr)=="ad_value_month"] <- "ad_value"
save(adyr,file="qsw2018/Data/Out/newspaper_ad_year")
################################################################################

################################################################################
##do `datadir'/copy_in_out.do
# skip this because why bother?
################################################################################

################################################################################
# Tables_2_3_data.do
#/********* Load data and generate reform variable ***********/
load("qsw2018/Data/Out/directory8111_en.Rdata")
tmp <- subset(direc, admin_rank=="county" & general==TRUE &
                     year==2002)
tmp$num_cg_pref_2002 <- 1
tmp <- aggregate(num_cg_pref_2002 ~ provincecheck + prefecturecheck,
                 data=tmp, FUN=sum)

# Authors' Stata code says
# > merge 1:m provincecheck prefecturecheck using Working/pc1s_analysis,keep(2 3);
# But what is Working/pc1s_analysis? it's created by Table1_A2.do, so
# here we go...

load("qsw2018/Data/Out/WN_analysis.Rdata")
pca <- subset(wn, general & tarticles>250)
for(v in c("provincecheck","prefecturecheck")) {
  pca[,v] <- as.factor(pca[,v])
}
pca$propreyear <- interaction(pca$provincecheck, pca$prefecturecheck,
                              as.factor(pca$year))
varlist <- c("s_leader", "s_xinhua", "bias3", "s_ntcorruption",
           "s_disasters", "s_accidents", "s_sports", "s_entertainment",
           "s_crime3")
for(v in varlist) {
  pca[, sprintf("r%s",v)] <- pca[,v] - predict(lm(as.formula(paste(v," ~ propreyear")),
                                                  data=pca, na.action=na.exclude))
}
pc <-
  prcomp(pca[!is.na(pca$rbias3),sprintf("r%s",varlist)],scale=TRUE)
# Principal components on residuals, but then predict on non-residuals
# This mimics authors' Stata code, but given how opaque the Stata
# code makes this, I wonder whether this is intentional. Probably.
X <- as.matrix(pca[,varlist])
for(c in 1:ncol(X)) {
  X[,c]  <- X[,c]/pc$scale[c]
}
summary(X %*% pc$rotation)
pca$pc1 <- (X %*% pc$rotation)[,1]
pca$pc1s <- with(pca, (pc1-min(pc1,na.rm=TRUE))/(max(pc1,na.rm=TRUE) -
                                                 min(pc1,na.rm=TRUE)))
dir.create("qsw2018/Analysis/Working")
save(pca,file="qsw2018/Analysis/Working/pc1s_analysis")
#################################################################################

load("qsw2018/Analysis/Working/pc1s_analysis")
tmp <- merge(tmp, pca, by=c("provincecheck","prefecturecheck"), all.y=TRUE)
tmp$num_cg_pref_2002[is.na(tmp$num_cg_pref_2002)] <- 0

################################################################################
## WN_Markets.do

#/* Select the 79 newspapers in 30 prefectures in 2002. Add the directory data for their 30 prefectures. */
foo <- pca[,c("provincecheck","prefecturecheck","s")]
foo <- unique(foo)
foo$year <- 2002
vars <- c("provincecheck","prefecturecheck","year","provincecheck_en","prefecturecheck_en","admin_rank","newspaper_id",names(direc)[grep("type",names(direc))])
bar <- merge(foo, direc[,vars], all.x=TRUE)
#/* Count number of Party Dailies and Commercial papers by prefecture and year */
bar <- bar[bar$ntype!="Other",]
bar$type_es <- bar$type_pe | bar$type_s
foo <- pca[,c("newspaper_id","year")]
foo$wn <- TRUE
wnm <- merge(bar, foo, by=c("newspaper_id","year"),all.x=TRUE)
wnm$wn[is.na(wnm$wn)] <- FALSE

wnm$type_wn_pd <- wnm$type_pd & wnm$wn
wnm$type_wn_es <- wnm$type_es & wnm$wn;
wnml <- aggregate(. ~ provincecheck + prefecturecheck +
                    provincecheck_en + prefecturecheck_en +
                    s + year +
                   admin_rank,
                 FUN=sum,
                 data=wnm[,c("provincecheck","prefecturecheck","provincecheck_en","prefecturecheck_en"
                            ,"s","year","admin_rank",
                             "type_pd", "type_es", "type_wn_pd", "type_wn_es")])
wnm <- data.frame()
for (r in unique(wnml$admin_rank)) {
  df <- subset(wnml, admin_rank==r)
  df <- df[,names(df)!="admin_rank"]
  for (v in names(df)[grep("type",names(df))]) {
    names(df)[names(df)==v] <- paste(v,r,sep="")
  }
  wnm <- merge(wnm, df, all=TRUE)
}
for (v in names(wnm)[grep("type",names(wnm))]) {
  wnm[is.na(wnm[,v]),v] <- 0
  names(wnm)[names(wnm)==v] <- sub("type_","",v)
}
wnm <- wnm[,!(names(wnm) %in% c("wn_pdcounty","escounty","wn_escounty"))]

wnm$pdupper <- with(wnm,pdprefecture+pdprovince+pdcentral)
wnm$esupper <- with(wnm, escentral+esprovince+esprefecture)
wnm$wn_pdupper <- with(wnm, wn_pdprefecture+wn_pdprovince+wn_pdcentral)
wnm$wn_esupper  <- with(wnm, wn_escentral+wn_esprovince+wn_esprefecture)

wnm$tier <- NA
wnm$tier[with(wnm, (s==3) & pdupper>0 & esupper>0) ]  <- 1
## *In Chongqing, the upper-level government operates only one party daily.
## * There are party dailies owned by the prefectural level, but these do not cover the same markest as the county paper markets.;
wnm$tier[with(wnm, s!=3 & (pdupper==1 |  (prefecturecheck_en=="Chongqing" ) ))] <- 2
wnm$tier[with(wnm, is.na(tier)  & (esupper <4) & prefecturecheck_en!="Nanchang")]  <- 3
wnm$tier[is.na(wnm$tier)] <- 4

save(wnm, file="qsw2018/Analysis/Working/WN_Markets.Rdata")
################################################################################

df <- merge(tmp,
            wnm[,c("provincecheck","prefecturecheck","tier","esupper","pdcounty")],
            all=TRUE)
df$num_cg_pref_2002 <- df$pdcounty
df <- df[!is.na(df$pdcounty),]
df$commercial <- df$ntype %in% c("Party Evening","Subsidiary")
for(y in 2002:2004) {
  df[,sprintf("reform_cg%d",y)] <- with(df, num_cg_pref_2002*(year>=y))
}

qsw <- df
save(qsw,file="qsw2018.Rdata")
################################################################################
