rm(list=ls())
# ----------------------
# Load data
# ---------------------

data = read.csv("../data/statesFourState.csv")
head(data)
dim(data)
str(data)


# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
dat_wo_U <- subset(data, state != "U")
dat_wo_U$state <- droplevels(dat_wo_U$state)

# ----------------------
### choice of variables
# ----------------------

selectedVars = c("annual_mean_temp", "pp_seasonality", "pp_warmest_quarter", "mean_diurnal_range","annual_pp", "mean_temperatre_wettest_quarter")

varCor2 = cor(dat_wo_U[, selectedVars])
varCor2

datSel = dat_wo_U[,c("state",selectedVars)]

rm(data, dat_wo_U)

# ----------------------
# models
# ----------------------


# ----------------------
# Calibration of the models
# ----------------------

# multimodal
# ----------------------

library(nnet)
SDM1 = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2), data = datSel, maxit =1500)
save(SDM1,file= "../data/Multinom_complete.rObj")

#load("../data/Multinom_complete.rObj")

# random Forest
# ----------------------

library(randomForest)
rs = runif(1,0,1)
set.seed(rs)
SDM2 = randomForest(state ~ . , data = datSel, ntree = 500)
save(SDM2,rs,file= "../data/RandomForest_complete.rObj")

#load("../data/RandomForest_complete.rObj")


# ----------------------
# Evaluation of the models
# ----------------------
# evaluation statistics
HK <- function (Pred, Obs)
{
	Misc = table(Pred, Obs)

    if (nrow(Misc)!=ncol(Misc)) stop("wrong misclassification table")
    Misc <- unclass(Misc)
    k  <- ncol(Misc)
    Nobs <- apply(Misc, 2, sum)
    Npred <- apply(Misc, 1, sum)
    N <- sum(Nobs)


   HK <- (sum(diag(Misc))/N - sum(as.numeric(Nobs)*as.numeric(Npred))/N/N ) / ( 1 - sum(as.numeric(Nobs)*as.numeric(Nobs))/N/N )

    return(HK)
}

# multinomial
# ----------------------
pred1 = predict(SDM1,new=datSel,"class")
(HK1 = HK(pred1, datSel$state)) # 0.202923

# random Forest
# ----------------------
set.seed(rs)
pred2 = predict(SDM2,new=datSel,"response", OOB=TRUE)
(HK2 = HK(pred2, datSel$state)) # 0.7468664


# ----------------------
# projection
# ---------------------

dataProj = read.csv("../data/transitionsFourState.csv")
head(dataProj)
dim(dataProj)
str(dataProj)

# multinomial
# ----------------------

proj1 = predict(SDM1,new=dataProj,"prob")
# sauvegarde
write.table(proj1, file = "../data/projection_multimod_complete.txt", quote=F, row.names=FALSE)


# random Forest
# ----------------------
set.seed(rs)
proj2 = predict(SDM2,new=dataProj,"prob", OOB=TRUE)
# sauvegarde
write.table(proj2, file = "../data/projection_rf_complete.txt", quote=F, row.names=FALSE)


