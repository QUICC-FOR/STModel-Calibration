rm(list=ls())
# ----------------------
# Load data
# ---------------------

load("../data/stateData.RData")
head(stateData)
dim(stateData)
str(stateData)

#jpeg("plots_states.jpeg")
#colo = c(T = "orange", B = "lightgreen", M = "blue", R =1  )
#plot(stateData$lat~stateData$lon, pch = 20, cex=.2, asp = 1, col = colo[stateData$state])
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
varCor = cor(dat_subset10[,-c(1:9)]) # check correlation between variables

nonCor = names(which(abs(varCor[,"annual_mean_temp"])<0.7))
varCor = varCor[nonCor,nonCor]
nonCor = names(which(abs(varCor[,"tot_annual_pp"])<0.7))
varCor = varCor[nonCor,nonCor]
nonCor = names(which(abs(varCor[,"mean_diurnal_range"])<0.7))
varCor[nonCor,nonCor]


selectedVars = c("annual_mean_temp", "tot_annual_pp", "mean_diurnal_range", "pp_warmest_quarter", "pp_wettest_period", "mean_temperatre_wettest_quarter", "mean_temp_driest_quarter")


library("ade4")
var.pca = dudi.pca(dat_subset10[,-c(1:5)], scannf=FALSE, nf = 5)
var.pca$eig / sum(var.pca$eig)

#inertia.dudi(var.pca, row=F,col=T)

pdf("../figures/PCA_selectedVariables.pdf")
s.arrow(var.pca$co[selectedVars,], clab = .6, xlim = c(-2,2))
dev.off()

varCor2 = cor(dat_subset10[, selectedVars])
varCor2


# ----------------------

datSel = dat_subset10[,c("state",selectedVars)]

rm(data, dat_wo_U)

# ----------------------
# scale transformation
# ----------------------
vars.means = apply(datSel[,2:7], 2, mean)
vars.sd = apply(datSel[,2:7], 2, sd)
save(vars.means, vars.sd, file = "scale_info.Robj")

datSel[,2:7] = scale(datSel[,2:7])

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


# multimodal
# ----------------------

library(nnet)

## stepwise

SDM1 = multinom(state ~ (annual_mean_temp + annual_pp  + pp_warmest_quarter + mean_diurnal_range+ mean_temperatre_wettest_quarter +mean_temp_driest_quarter)^2 + I(annual_pp^2) + I(annual_mean_temp^2) + I(pp_warmest_quarter^2)- mean_diurnal_range -pp_warmest_quarter:mean_diurnal_range + I(mean_temp_driest_quarter^2) + I(mean_temperatre_wettest_quarter^2), data = datSel, maxit =1500)

save(SDM1,file= "../data/Multinom_complete.rObj")


# evaluation

pred1 = predict(SDM1,new=datSel,"class")
(HK1 = HK(pred1, datSel$state)) # 0.33
tt1 = table(pred1, datSel$state)
#error rate
1-sum(diag(tt1))/sum(tt1) # 46%

dat_rc = data.frame(annual_mean_temp = seq(-5, 3, length.out = 200), annual_pp = rep(0, 200), mean_diurnal_range=rep(0, 200),pp_warmest_quarter= rep(0, 200),mean_temperatre_wettest_quarter= rep(0, 200), mean_temp_driest_quarter=rep(0, 200))

pred1_rc =  predict(SDM1,new=dat_rc,"probs")

##--- 
pdf("../figures/multinomial predictions.pdf")
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "multinomial")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
##--
plot(dat_subset10[,"longitude"], dat_subset10[,"latitude"], cex = .1, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1, col = colo[as.character(pred1)])
title("multinomial")
dev.off()

#rm(SDM1)

# multimodal two steps
# ----------------------

library(nnet)

datSel$rstate = ifelse(datSel$state =="R",  1, 0)

## manual stepwise (drop 2 square terms and some interaction terms)
#library(lmtest)
#mod0 = glm(rstate ~ 1 , data = datSel, family = "binomial")
#mod1 = glm(rstate ~ annual_pp , data = datSel, family = "binomial")
#mod2 = glm(rstate ~ annual_mean_temp + mean_temperatre_wettest_quarter + annual_pp + mean_temp_driest_quarter + pp_warmest_quarter + mean_diurnal_range, data = datSel, family = "binomial")
#mod3 = glm(rstate ~ (annual_mean_temp + mean_temperatre_wettest_quarter + annual_pp + mean_temp_driest_quarter + pp_warmest_quarter + mean_diurnal_range)^2, data = datSel, family = "binomial")
#mod4 = glm(rstate ~ (annual_mean_temp + mean_temperatre_wettest_quarter + annual_pp + mean_temp_driest_quarter + pp_warmest_quarter + mean_diurnal_range)^2 + I(annual_pp^2) + I(annual_mean_temp^2) + I(mean_diurnal_range^2), data = datSel, family = "binomial")
#lrtest(modAll, SDM1.R)

SDM1.R = glm(rstate ~ (annual_mean_temp + annual_pp  + pp_warmest_quarter + mean_diurnal_range+ mean_temperatre_wettest_quarter +mean_temp_driest_quarter)^2 + I(annual_pp^2) + I(annual_mean_temp^2) + I(pp_warmest_quarter^2)- mean_diurnal_range -pp_warmest_quarter:mean_diurnal_range + I(mean_temp_driest_quarter^2) + I(mean_temperatre_wettest_quarter^2), data = datSel, family = "binomial")

# step1 evaluation intermédiaire
predR = predict(SDM1.R,new=datSel,"response")
boxplot(predR~datSel$rstate)
title("model R")

source("BoulangeatEcoLet2012_fct.r")
cutoff = CutOff.optim(predR, datSel$rstate)$CutOff

SDM1.2 = multinom(state ~ (annual_mean_temp + annual_pp  + pp_warmest_quarter + mean_diurnal_range+ mean_temperatre_wettest_quarter +mean_temp_driest_quarter)^2 + I(annual_pp^2) + I(annual_mean_temp^2) + I(pp_warmest_quarter^2)- mean_diurnal_range -pp_warmest_quarter:mean_diurnal_range + I(mean_temp_driest_quarter^2) + I(mean_temperatre_wettest_quarter^2), data = datSel[which(datSel$state!="R"),], maxit =1500)

save(SDM1.R, SDM1.2, cutoff, file= "../data/Multinom_complete_2steps.rObj")

# evaluation

pred1.2 = predict(SDM1.2,new=datSel[predR<cutoff,],"class")
pred1.final = predR
pred1.final[predR>=cutoff] = "R"
pred1.final[predR<cutoff] = as.character(pred1.2)
pred1.final = as.factor(pred1.final)
(HK1.2 = HK(pred1.final, datSel$state)) # 0.31
tt1 = table(pred1.final, datSel$state)
#error rate
1-sum(diag(tt1))/sum(tt1) # 47%

dat_rc = data.frame(annual_mean_temp = seq(-5, 3, length.out = 200), annual_pp = rep(0, 200), mean_diurnal_range=rep(0, 200),pp_warmest_quarter= rep(0, 200),mean_temperatre_wettest_quarter= rep(0, 200), mean_temp_driest_quarter=rep(0, 200))

pred1_rcR =  predict(SDM1.R,new=dat_rc,"response")
pred1_rc = predict(SDM1.2,new=dat_rc,"probs")
pred1_rc = data.frame(pred1_rc)
pred1_rc$R = pred1_rcR
pred1_rc = t(apply(pred1_rc, 1, function(x)x/sum(x)))

##--- 
pdf("../figures/multinomial2steps predictions.pdf")
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "multinomial 2 steps")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
##--
plot(dat_subset10[,"longitude"], dat_subset10[,"latitude"], cex = .1, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1, col = colo[as.character(pred1.final)])
title("multinomial 2 steps model")
dev.off()

#rm(SDM1.R, SDM1.2)

# random Forest
# ----------------------

datSel = datSel[,1:7] ## careful not include rstate variable!

library(randomForest)
rs = runif(1,0,1)
set.seed(rs)
SDM2 = randomForest(state ~ . , data = datSel, ntree = 500)

save(SDM2,rs,file= "../data/RandomForest_complete.rObj")


# evaluation random Forest
set.seed(rs)
pred2 = predict(SDM2,new=datSel,"response", OOB=TRUE)
(HK2 = HK(pred2, datSel$state)) # 0.7468664
#error rate
tt2 = table(pred2, datSel$state)
1-sum(diag(tt2))/sum(tt2)


#--

dat_rc = data.frame(annual_mean_temp = seq(-5, 3, length.out = 200), annual_pp = rep(0, 200), mean_diurnal_range=rep(0, 200),pp_warmest_quarter= rep(0, 200),mean_temperatre_wettest_quarter= rep(0, 200), mean_temp_driest_quarter=rep(0, 200))

set.seed(rs)
pred1_rc =  predict(SDM2,new=dat_rc,"prob")

##--- 
pdf("../figures/randomForest predictions.pdf", width = 8, height = 8)
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "multinomial 2 steps")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
##--
plot(dat_subset10[,"longitude"], dat_subset10[,"latitude"], cex = .1, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1, col = colo[as.character(pred2)])
title("random forest")
dev.off()


#rm(SDM2)


# ----------------------
# projection
# ---------------------

dataProj = read.csv("../data/transitionsFourState_inf1.csv")
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


