rm(list = ls())

#sdm = "rf"
sdm = "multinom"
propData = 0.33

#--
veget_pars = read.table(paste("../estimated_params/GenSA_initForFit_", sdm, "_",propData, ".txt", sep=""))
load(paste("initForFit_",sdm, "_", propData, ".RData", sep = ""))
load(paste("../estimated_params/GenSA_initForFit_", sdm, "_", propData, ".RData", sep = ""))
#---

load("../data/transitions_r1.RData")
dataProj = transitionData
select = unique(dataProj$plot[which(dataProj$annual_mean_temp<=10)])
dataProj_subset10 = dataProj[dataProj$plot %in% select,]
toremove = c(which(dataProj_subset10$state1 == "T" & dataProj_subset10$state2 == "B"), which(dataProj_subset10$state1 == "B" & dataProj_subset10$state2 == "T"))
dataProj_subset10= dataProj_subset10[-toremove,]
if(sdm == "rf") pred = read.table("../data/projection_rf_complete.txt", h=T)
if(sdm == "multinom") pred = read.table("../data/projection_multimod_complete.txt", h=T)
pred = pred[pred$plot %in% select,]
pred = pred[-toremove,]

#--
print(fit)
#--

#dim(datSel)
#dim(dataProj_subset10)
#length(select2)
datValid = cbind(dataProj_subset10[-select2,], pred[-select2,])
nrow(datValid)+nrow(datSel) == nrow(dataProj_subset10)
nrow(datValid)+nrow(datSel) == nrow(pred)

#--
load("scale_info.Robj")

#---
pars = as.numeric(veget_pars[,1])
names(pars) = rownames(veget_pars)
#--

ENV = data.frame(TP =datValid$annual_mean_temp , PP = datValid$tot_annual_pp)
ENV1 = (ENV$TP - vars.means["annual_mean_temp"])/ vars.sd["annual_mean_temp"]
ENV2 = (ENV$PP - vars.means["tot_annual_pp"])/ vars.sd["tot_annual_pp"]

e7 = 0
params=pars


    logit_alphab 	= params["ab0"] + params["ab1"]*ENV1 + params["ab2"]*ENV2 + params["ab3"]*ENV1^2 + params["ab4"]*ENV2^2 + params["ab5"]*ENV1^3 + params["ab6"]*ENV2^3
    logit_alphat 	= params["at0"] + params["at1"]*ENV1 + params["at2"]*ENV2 + params["at3"]*ENV1^2 + params["at4"]*ENV2^2 + params["at5"]*ENV1^3 + params["at6"]*ENV2^3
    logit_betab 	= params["bb0"] + params["bb1"]*ENV1 + params["bb2"]*ENV2 + params["bb3"]*ENV1^2 + params["bb4"]*ENV2^2 + params["bb5"]*ENV1^3 + params["bb6"]*ENV2^3
    logit_betat 	= params["bt0"] + params["bt1"]*ENV1 + params["bt2"]*ENV2 + params["bt3"]*ENV1^2 + params["bt4"]*ENV2^2 + params["bt5"]*ENV1^3 + params["bt6"]*ENV2^3
    logit_theta	= params["t0"] + params["t1"]*ENV1 + params["t2"]*ENV2 + params["t3"]*ENV1^2 + params["t4"]*ENV2^2 + params["t5"]*ENV1^3 + params["t6"]*ENV2^3
    logit_thetat	= params["tt0"] + params["tt1"]*ENV1 + params["tt2"]*ENV2 + params["tt3"]*ENV1^2 + params["tt4"]*ENV2^2 + params["tt5"]*ENV1^3 + params["tt6"]*ENV2^3
    logit_eps 	= params["e0"]  + params["e1"]*ENV1 + params["e2"]*ENV2  + params["e3"]*ENV1^2 + params["e4"]*ENV2^2 + params["e5"]*ENV1^3 + params["e6"]*ENV2^3 
    #e7*EB
 

## annual proba       
 logit_reverse <- function(x)
{
expx = ifelse(exp(x)==Inf, .Machine$double.xmax, exp(x))
expx/(1+expx)
}

itime = datValid$year2 - datValid$year1
    
macroPars = data.frame(alphab = 1-(1-logit_reverse(logit_alphab))^itime, 
alphat = 1-(1-logit_reverse(logit_alphat))^itime,
betab = 1-(1-logit_reverse(logit_betab))^itime,
betat = 1-(1-logit_reverse(logit_betat))^itime,
theta = 1-(1-logit_reverse(logit_theta))^itime,
thetat = 1-(1-logit_reverse(logit_thetat))^itime,
eps = 1-(1-logit_reverse(logit_eps))^itime)

summary(macroPars)


## proba transitions


pred = data.frame(matrix(NA, nrow(datValid), ncol = 4))
colnames(pred) = c("R", "T", "B", "M")
#
phib = macroPars$alphab*(datValid$M + datValid$B)*(1-macroPars$alphat*(datValid$T+datValid$M))
phit = macroPars$alphat*(datValid$M + datValid$T)*(1-macroPars$alphab*(datValid$B+datValid$M))
phim = macroPars$alphab*(datValid$M + datValid$B)*macroPars$alphat*(datValid$M + datValid$T)
pred[datValid$state1 == "R",] = data.frame(R=1 - phib - phit - phim, T=phit, B = phib, M = phim)[datValid$state1 == "R",]
#
pred[datValid$state1 == "T",] = data.frame(R=macroPars$eps, T=1- macroPars$eps - macroPars$betab*(datValid$B+datValid$M)*(1-macroPars$eps), B = rep(0, nrow(datValid)), M = macroPars$betab*(datValid$B+datValid$M)*(1-macroPars$eps))[datValid$state1 == "T",]
#
pred[datValid$state1 == "B",] = data.frame(R=macroPars$eps, T=rep(0, nrow(datValid)), B = 1- macroPars$eps - macroPars$betat*(datValid$T+datValid$M)*(1-macroPars$eps), M = macroPars$betat*(datValid$T+datValid$M)*(1-macroPars$eps))[datValid$state1 == "B",]
#
pred[datValid$state1 == "M",] = data.frame(R=macroPars$eps, T=macroPars$theta*macroPars$thetat*(1-macroPars$eps), B = macroPars$theta*(1-macroPars$thetat)*(1-macroPars$eps), M = (1 - macroPars$eps)*(1 - macroPars$theta))[datValid$state1 == "M",]
##

head(pred)
pred.state = factor(rep("R", nrow(pred)), levels= c("B","M","T","R"))
pred.state[1:nrow(pred)] = colnames(pred)[apply(pred, 1, which.max)]

##---------------------------

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

##---------------------------

HK(pred.state, datValid$state2)


