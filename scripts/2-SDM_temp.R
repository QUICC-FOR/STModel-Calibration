rm(list=ls())
# ----------------------
# load calibration data 
# ----------------------

data = read.csv("../data/statesFourState.csv")
#data = read.csv("../STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)



# ----------------------
### choice of variables
# ----------------------


selectedVars = c("annual_mean_temp",  "annual_pp")

datSel = data[,c("state",selectedVars)]

rm(data)


# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
str(datSel)
datSel_wo_U <- subset(datSel, state != "U")
datSel_wo_U$state <- droplevels(datSel_wo_U$state)

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
calib = datSel_wo_U[sampl,c("state",selectedVars)] 
valid = datSel_wo_U[-sampl,c("state",selectedVars)]

# Run the models

# random Forest
# calib
library(randomForest)
rs = runif(1,0,1)
set.seed(rs)
SDM2 = randomForest(state ~ . , data = calib, ntree = 500)
SDM2
save(SDM2,rs,sampl,file= "../data/RandomForest_temp.rObj")


# valid
set.seed(rs)
pred2 = predict(SDM2,new=datSel_wo_U,"response", OOB=TRUE)
(HK2 = HK(pred2[-sampl], valid$state)) 


# multimodal
#calib
library(nnet)
SDM1 = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) , data = calib, maxit =300)
summary(SDM1)
save(SDM1,file= "../data/Multinom_temp.rObj")

#valid
pred1 = predict(SDM1, new=valid,"class")
(HK1 = HK(pred1, valid$state)) 


# ----------------------
# projection 
# ----------------------
## ----recap data
load("../data/Multinom_temp.rObj")
load("../data/RandomForest_temp.rObj")
selectedVars = c("annual_mean_temp",  "annual_pp")
#---------------

dataProj = read.csv("../data/transitionsFourState.csv")
head(dataProj)

# means between climates of year of state 0 and year of state 1
dataProj$annual_mean_temp = apply(dataProj[,c("annual_mean_temp1", "annual_mean_temp2.")], 1, mean)
dataProj$annual_pp = apply(dataProj[,c("annual_pp1", "annual_pp2")], 1, mean)

datProjSel = dataProj[,selectedVars]

# projection
set.seed(rs)
projProba = predict(SDM2,new=datProjSel,"prob", OOB=TRUE)
head(projProba)


# sauvegarde
write.table(projProba, file = "../data/projection_neigbor_rf_temp.txt", quote=F, row.names=dataProj$X.plot)






