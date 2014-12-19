rm(list=ls())
# Load data

#data = read.csv("../data/statesFourState.csv")
data = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)

setwd("~/Documents/GitHub/STModel-Calibration/scripts")

# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
dat_wo_U <- subset(data, state != "U")
dat_wo_U$state <- droplevels(dat_wo_U$state)

selectedVars = c("state","annual_mean_temp", "pp_seasonality", "pp_warmest_quarter", "mean_diurnal_range","annual_pp", "mean_temperatre_wettest_quarter")


# ----------------------
# models
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



# cross validation - separation of the dataset
sampl = sample(1:nrow(datSel_wo_U), 2*nrow(datSel_wo_U)/3)
calib = datSel_wo_U[sampl,c("state",selectedVars_CCbio)]
valid = datSel_wo_U[-sampl,c("state",selectedVars_CCbio)]

# Run the models

# random Forest
# calib
library(randomForest)
rs = runif(1,0,1)
set.seed(rs)
SDM2 = randomForest(state ~ . , data = calib, ntree = 500)
save(SDM2,rs,sampl,file= "RandomForest_6vars.rObj")

# valid
set.seed(rs)
pred2 = predict(SDM2,new=datSel_wo_U,"response", OOB=TRUE)
(HK2 = HK(pred2[-sampl], valid$state))


# multimodal
#calib
library(nnet)


SDM1 = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2), data = calib, maxit =1500)

save(SDM1,sampl,file= "Multinom_6vars_version_new_class.rObj")
rm(SDM1)

summary(SDM1,Wald=TRUE)

pred1 = predict(SDM1,new=datSel_wo_U,"class")

HK1 = HK(pred1[-sampl], valid$state)