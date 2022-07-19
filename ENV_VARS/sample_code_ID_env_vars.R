library(corrplot) ## corrplot
library(Hmisc) ## rcorr
library(dplyr) ## sample_frac
library(pROC) ## AUC
library(MuMIn) ## stdize

##### MOCK DATA #####
## df is a data frame in which the first column is presence/absence (or abundance not including 0s) and the rest are the environmental variables.
pa <- c(rep(0,25), rep(1,25))
env <- matrix(runif(200,0,100), nrow=50, ncol=4)
df <- as.data.frame(cbind(pa, env))
vars <- c("env1","env2","env3","env4")
colnames(df) <- c("pa",vars)

## Standardise environmental variables 
df[,vars]  <- stdize(df[,vars])

####### CORRELATIONS #####
corrs <- cor(df[-1], method = c("spearman")) ## Save with separate name for each species. Probably not always necessary. Just instructive to look at to begin with.
jpeg(paste0(wd.out,"/corrs.jpg"), width=600, height=600)
corrplot(corrs)
dev.off()

#### Extract correlations > 0.7 and tabulate
corrs2 <- rcorr(as.matrix(df[,-1]))
corrs2$r[lower.tri(corrs2$r)] <- NA
a <- as.data.frame(which(abs(corrs2$r) >= 0.7 & abs(corrs2$r) != 1, arr.ind = TRUE))
##  correlation >= 0.7
r <- rownames(corrs2$r)[a$row]
c <- colnames(corrs2$r)[a$col]

corr.out <- data.frame(matrix(nrow=nrow(a), ncol=2, dimnames=list(c(), c("vars", "corrs"))), stringsAsFactors=F)
corr.out$vars <- paste0(r, " & ", c)
corr.out$corrs <- apply(a, MARGIN=1, FUN=function(x) { corrs2$r[x[1], x[2]] }) 
# write.csv(corr.out, file=paste0(wd.out, "/corrs.csv"), row.names=F) ## Save with seperate name for each species

##### DECIDE WHICH CORRELATED VARIABLES TO RETAIN #####
### Calibration and validation data
df$id <- 1:nrow(df)

dat.c <- sample_frac(df, size=0.7)
dat.v <- df[!(df$id %in% dat.c$id),]

### Variables to choose between - for variables that are measured in teh same seasons. 
vars <- paste0(c("env1", "env2", "env3"))

for (v in vars) {
  m <- glm(as.formula(paste0("pa ~ poly(",v,",2)")), family=binomial, dat.c) ## Make model with calibration data. May well want to use another method as well, e.g. random forest
  p <- predict(m, dat.v, re.form=NA, type="response") ## Predict model with validation data
  
  roc_obj <- roc(dat.v$pa, p)
  print(paste(v, round(auc(roc_obj),3))) ## May well want to use other measures of goodness of fit (e.g. deviance explained) or predictive poser (e.g. TSS)
}

## Tabulate predictive capacity / goodness of fit data for each native species and the invader composite metrics
## Summarise data across each species / summary metrics to identify the best performing environmental metrics for most species



