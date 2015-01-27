rm(list=ls())
# ----------------------
# Load data
# ---------------------

load("../data/stateData.RData")
head(stateData)
dim(stateData)
str(stateData)
#
#jpeg("../figures/plots_states.jpeg")
#colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
#plot(stateData$lat~stateData$lon, pch = 20, cex=.2, asp = 1, col = colo[stateData$state])
#dev.off()
#
#jpeg("../figures/plots_Rstates.jpeg")
#plot(stateData$lat~stateData$lon, pch = 20, cex=.2, asp = 1, col = ifelse(stateData$state=="R", 1, "grey"), main="R states")
#dev.off()

# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
dat_wo_U <- subset(stateData, state != "U")
#dat_wo_U$state <- droplevels(dat_wo_U$state)

# subset 10 degree
select = unique(dat_wo_U$plot[which(dat_wo_U$annual_mean_temp<=10)])
dat_subset10 = dat_wo_U[dat_wo_U$plot %in% select,]




# ----------------------
### choice of variables
# ----------------------
#
#varCor = cor(dat_subset10[,-c(1:9)]) # check correlation between variables
#
#nonCor = names(which(abs(varCor[,"annual_mean_temp"])<0.7))
#varCor = varCor[nonCor,nonCor]
#nonCor = names(which(abs(varCor[,"tot_annual_pp"])<0.7))
#varCor = varCor[nonCor,nonCor]
#nonCor = names(which(abs(varCor[,"mean_diurnal_range"])<0.7))
#varCor[nonCor,nonCor]


selectedVars = c("annual_mean_temp", "tot_annual_pp", "mean_diurnal_range", "pp_warmest_quarter", "pp_wettest_period", "mean_temp_wettest_quarter", "mean_temp_driest_quarter")


#library("ade4")
#var.pca = dudi.pca(dat_subset10[,-c(1:9)], scannf=FALSE, nf = 5)
#(var.pca$eig / sum(var.pca$eig))*100

#inertia.dudi(var.pca, row=F,col=T)

#pdf("../figures/PCA_selectedVariables.pdf")
#s.arrow(var.pca$co[selectedVars,], clab = .6, xlim = c(-2,2), sub = "axe 1: 47% ; axe 2: 30 %")
#dev.off()

varCor2 = cor(dat_subset10[, selectedVars])
varCor2


# ----------------------

datSel = dat_subset10[,c("state",selectedVars)]
datSel$state = as.factor(datSel$state)

rm(stateData, dat_wo_U)

# ----------------------
# scale transformation
# ----------------------
vars.means = apply(datSel[,2:8], 2, mean)
vars.sd = apply(datSel[,2:8], 2, sd)
save(vars.means, vars.sd, file = "scale_info.Robj")

datSel[,2:8] = scale(datSel[,2:8])

# ----------------------
# Calibration of the models
# ----------------------

# evaluation statistics (revised from Boulangeat et al. 2012 to handle large numbers)
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
# ----------------------
# subsample of the data
# ----------------------
source("subsample.r")
select = subsample.temp(dat_subset10$annual_mean_temp, .33)

jpeg("../figures/subsample_SDM_calib.jpeg", height=5000, width=5000, res=600)
plot(dat_subset10[,c("lon","lat")], pch = 20, cex=.2, col = "grey")
points(dat_subset10[select,c("lon","lat")], pch = 20, cex=.2, col = 1)
dev.off()

save(select, file = "subsample_SDM_calib.RData")

#-----------------------

datCal = datSel[select,]
datEval = datSel[-select,]

# ----------------------
# multimodal
# ----------------------

library(nnet)

##  complete
# ----------------------

SDM1 = multinom(state ~ (annual_mean_temp + tot_annual_pp  + pp_warmest_quarter + pp_wettest_period + mean_diurnal_range+ mean_temp_wettest_quarter +mean_temp_driest_quarter)^2 + I(tot_annual_pp^2) + I(annual_mean_temp^2) + I(pp_warmest_quarter^2) + I(pp_wettest_period^2) + I(mean_temp_driest_quarter^2) + I(mean_temp_wettest_quarter^2), data = datSel, maxit =1500)
save(SDM1,file= "../data/Multinom_complete.rObj")

# evaluation
pred1 = predict(SDM1,new=datSel,"class")
(HK1 = HK(pred1, datSel$state)) 
tt1 = table(pred1, datSel$state)
#error rate
1-sum(diag(tt1))/sum(tt1) 

# graph
# -----
dat_rc = data.frame(annual_mean_temp = seq(-5, 3, length.out = 200), tot_annual_pp = rep(0, 200), mean_diurnal_range=rep(0, 200),pp_warmest_quarter= rep(0, 200),mean_temp_wettest_quarter= rep(0, 200), mean_temp_driest_quarter=rep(0, 200),pp_wettest_period= rep(0,200) )

pred1_rc =  predict(SDM1,new=dat_rc,"probs")

##--- 
pdf("../figures/multinomial predictions.pdf")
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "multinomial")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
dev.off()
##--
jpeg("../figures/multinomial predictions_map.jpeg", height=5000, width=5000, res=600)
plot(dat_subset10[,"lon"], dat_subset10[,"lat"], cex = .1, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1, col = colo[as.character(pred1)])
title("multinomial")
dev.off()

#rm(SDM1)

## calib / eval
# ----------------------
SDM1_cal = multinom(state ~ (annual_mean_temp + tot_annual_pp  + pp_warmest_quarter + pp_wettest_period + mean_diurnal_range+ mean_temp_wettest_quarter +mean_temp_driest_quarter)^2 + I(tot_annual_pp^2) + I(annual_mean_temp^2) + I(pp_warmest_quarter^2) + I(pp_wettest_period^2) + I(mean_temp_driest_quarter^2) + I(mean_temp_wettest_quarter^2), data = datCal, maxit =1500)
save(SDM1_cal,file= "../data/Multinom_cal.rObj")

# evaluation
pred1 = predict(SDM1_cal,new=datEval,"class")
(HK1 = HK(pred1, datEval$state)) 
tt1 = table(pred1, datEval$state)
#error rate
1-sum(diag(tt1))/sum(tt1) 

# graph
# -----
pred1_rc =  predict(SDM1_cal,new=dat_rc,"probs")

##--- 
pdf("../figures/multinomial predictions_cal.pdf")
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "multinomial")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
dev.off()
##--



# ----------------------
# random Forest
# ----------------------

##  complete
# ----------------------

library(randomForest)
rs = runif(1,0,1)
set.seed(rs)
SDM2 = randomForest(state ~ . , data = datSel, ntree = 500)

save(SDM2,rs,file= "../data/RandomForest_complete.rObj")

(imp = importance(SDM2))

# evaluation random Forest
set.seed(rs)
pred2 = predict(SDM2,new=datSel,"response", OOB=TRUE)
(HK2 = HK(pred2, datSel$state)) # 0.7468664
#error rate
tt2 = table(pred2, datSel$state)
1-sum(diag(tt2))/sum(tt2)


#--

dat_rc = data.frame(annual_mean_temp = seq(-5, 3, length.out = 200), tot_annual_pp = rep(0, 200), mean_diurnal_range=rep(0, 200),pp_warmest_quarter= rep(0, 200),mean_temp_wettest_quarter= rep(0, 200), mean_temp_driest_quarter=rep(0, 200),pp_wettest_period= rep(0,200) )


set.seed(rs)
pred1_rc =  predict(SDM2,new=dat_rc,"prob")

##--- 
pdf("../figures/randomForest predictions.pdf", width = 8, height = 8)
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "random forest")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
dev.off()
##--
jpeg("../figures/randomForest predictions_map.jpeg", height=5000, width=5000, res=600)
plot(dat_subset10[,"lon"], dat_subset10[,"lat"], cex = .1, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1, col = colo[as.character(pred2)])
title("random forest")
dev.off()

#rm(SDM2)

## calib / eval
# ----------------------

set.seed(rs)
SDM2_cal = randomForest(state ~ . , data = datCal, ntree = 500)

save(SDM2_cal,rs,file= "../data/RandomForest_cal.rObj")

(imp = importance(SDM2_cal))

# evaluation random Forest
set.seed(rs)
pred2 = predict(SDM2_cal,new=datEval,"response", OOB=TRUE)
(HK2 = HK(pred2, datEval$state)) 
#error rate
tt2 = table(pred2, datEval$state)
1-sum(diag(tt2))/sum(tt2)

# graph
set.seed(rs)
pred1_rc =  predict(SDM2_cal,new=dat_rc,"prob")

##--- 
pdf("../figures/randomForest predictions_cal.pdf", width = 8, height = 8)
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "random forest")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
dev.off()
##--


# ----------------------
# projection
# ---------------------

load("../data/transitionData.RData")

dataProj = transitionData
head(dataProj)
dim(dataProj)
str(dataProj)

# subset 10 degree
select = unique(dataProj$plot[which(dataProj$annual_mean_temp<=10)])
dataProj_subset10 = dataProj[dataProj$plot %in% select,]

#pdf("../figures/transition plots.pdf")
#plot(dataProj_subset10[,"longitude"], dataProj_subset10[,"latitude"], cex = .5, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1)
#dev.off()

# rescale
dataRescaledProj = dataProj_subset10[,selectedVars]
dataRescaledProj = t(apply(dataRescaledProj, 1, function(x) {(x-vars.means)/vars.sd}))


# multinomial
# ----------------------
#load("../data/Multinom_complete.rObj")

proj1 = predict(SDM1,new=dataRescaledProj,"prob")
head(proj1)
summary(proj1)
# sauvegarde
write.table(proj1, file = "../data/projection_multimod_complete.txt", quote=F, row.names=FALSE)

# multinomial 2 states
# ----------------------
#load("../data/Multinom_complete_2steps.rObj")

projR =  predict(SDM1.R,new=data.frame(dataRescaledProj),"response")
proj1.2 = predict(SDM1.2,new=dataRescaledProj,"probs")
proj1.2 = data.frame(proj1.2)
proj1.2$R = projR
proj1.2 = t(apply(proj1.2, 1, function(x)x/sum(x)))
head(proj1.2)
summary(proj1.2)
# sauvegarde
write.table(proj1.2, file = "../data/projection_multimod_complete_2steps.txt", quote=F, row.names=FALSE)


# random Forest
# ----------------------
#load("../data/RandomForest_complete.rObj")

set.seed(rs)
proj2 = predict(SDM2,new=dataRescaledProj,"prob", OOB=TRUE)
head(proj2)
summary(proj2)
# sauvegarde
write.table(proj2, file = "../data/projection_rf_complete.txt", quote=F, row.names=FALSE)


