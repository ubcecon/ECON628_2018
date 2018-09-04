## Download and run code used by Mullainathan & Spiess (2017)
if (!file.exists("4211.zip")) {
  download.file("https://www.aeaweb.org/articles/attachments?retrieve=hfL4yDm2epidrc9MWUOkEPFSp-mJmvum",
                destfile="4211.zip")
}
if (!file.exists("MullainathanSpiess/0_dataprep.R")) {
  unzip("4211.zip")
  ## why does MAC's zip create useless __MACOSX directory and  files?
  unlink("__MACOSX", recursive=TRUE) ## delete that garbage
}

if (!file.exists("ahs2011forjep.rdata")) {
  # download data
  if (!file.exists("tnewhouse.zip")) {
    download.file("http://www2.census.gov/programs-surveys/ahs/2011/AHS%202011%20National%20and%20Metropolitan%20PUF%20v1.4%20CSV.zip",
                  destfile="tnewhouse.zip")
  }
  unzip("tnewhouse.zip",files="tnewhouse.csv")

  # run M&S dataprep code
  source("MullainathanSpiess/0_dataprep.R")
}


## 1_loadpredictionstools uses relative paths, so we must change
## directories. Wrap in a tryCatch block so we don't end up in a
## different directory if there's an error
setwd("./MullainathanSpiess")
tryCatch(source("1_loadpredictiontools.R"),
         finally = setwd(".."))
# if there's errors, you likely need to install the missing
# packages,e.g. 'install.packages(c("xgboost","randomForest"))'

if (!files.exists("jepfittedmodels-ensemble")) {
  # Fit models: this takes awhile
  source("MullainathanSpiess/3_fitmodels.R")
}
source("MullainathanSpiess/4_table1.R")

fulldata$id <- 1:nrow(fulldata)
df <- melt(fulldata[,c(1:7,ncol(fulldata))],
           id=c("id","holdout","LOGVALUE"),
           variable_name="Method")
df$holdout[df$holdout==FALSE] <- "Training"
df$holdout[df$holdout==TRUE] <- "Holdout"
df$holdout <- as.factor(df$holdout)
names(df)[names(df)=="value"] <- "Prediction"
df$Error <- df$LOGVALUE - df$Prediction
fig <- ggplot(data=df, aes(x=Error, colour=Method)) +
  geom_density() + theme_minimal() +
  xlim(quantile(df$Error,c(0.02,0.98))) +
  facet_grid(holdout ~ .)

fig2 <- ggplot(data=df, aes(x=LOGVALUE, y=Prediction,
                            colour=Method, shape=holdout)) +
  geom_point(alpha=0.5) + theme_minimal() + geom_line(aes(y=LOGVALUE))

save(fig,fig2, file="jeperrorfig.RData")
