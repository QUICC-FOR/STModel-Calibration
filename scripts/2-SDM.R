rm(list=ls())
# Load data
#dataProj = as.data.frame(read.table("../data/data_pairs_filter.txt"))
#dataProj$E = scale(dataProj$annual_mean_temp)
#dataProj$P = scale(dataProj$annual_pp)

data = read.csv("../data/statesFourState.csv")
head(data)
dim(data)

# ----------------------
### choice of variables
# ----------------------

varCor = cor(data[,-c(1:5)]) # check correlation between variables
varCor

# acp
library(ade4)
var.pca = dudi.pca(data[,-c(1:5)], scannf=FALSE, nf = 5)
var.pca$eig / sum(var.pca$eig)

s.corcircle(var.pca$co, clab = 0.5)
contrib = inertia.dudi(var.pca, row = FALSE, col = TRUE)$col.abs

nonCorVars = intersect(names(varCor[which(abs(varCor[,"annual_mean_temp"])<0.7),"annual_mean_temp"]), names(varCor[which(abs(varCor[,"pp_seasonality"])<0.7),"pp_seasonality"]))

contrib[nonCorVars,]

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"pp_warmest_quarter"])<0.7),"pp_warmest_quarter"]))

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"mean_diurnal_range"])<0.7),"mean_diurnal_range"]))

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"annual_pp"])<0.7),"annual_pp"]))

selectedVars = c("annual_mean_temp", "pp_seasonality", "pp_warmest_quarter", "mean_diurnal_range", "annual_pp", "mean_temperatre_wettest_quarter")

varCor2 = cor(data[, selectedVars])

datSel = data[,c("state",selectedVars)]

rm(data)


# ----------------------
# models 
# ----------------------
# evaluation statistics
source("BoulangeatEcoLet2012_fct.r") 

# cross validation - separation of the dataset
sampl = sample(1:nrow(datSel), 2*nrow(datSel)/3)
calib = datSel[sampl,c("state",selectedVars)] 
valid = datSel[-sampl,c("state",selectedVars)]

# Run the models

# random Forest
# calib
library(randomForest)
rs = runif(1,0,1)
set.seed(rs)
SDM2 = randomForest(state ~ . , data = calib, ntree = 500)

# valid
set.seed(rn)
pred2 = predict(SDM2,new=valid,"response", OOB=TRUE)
(HK2 = HK(pred2, valid$state)) 



# multimodal
#calib
library(nnet)
SDM1 = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2), data = calib, maxit =200)
summary(SDM1)

#valid
pred1 = predict(SDM1, new=valid,"class")
(HK1 = HK(pred1, valid$state)) 




#---------------------------------------------------------------------
# Get predictions
#newdat = expand.grid(E=seq(-1.7,7,0.1),P=seq(700, 1600, 10) )

# multinom
#library(nnet)
#SDM1 = multinom(state ~ poly(E,3,raw=TRUE) + poly(P,3,raw=TRUE) + E:P + lon* lat , data = data, maxit =200)
#
#pred_multi = predict(SDM1,new=dataProj,"prob")
#write.table(pred_multi,"../data/pred_states_multinom.txt")
#
#pred_class = predict(SDM1,new=dataProj,"class")
#(HK1 = HK(pred_class, dataProj$st0)) 
#(HK1 = HK(pred_class, dataProj$st1)) 
#
#
# random forest
#rn = runif(1,0,1)
#set.seed(rn)
#SDM2 = randomForest(state ~ . , data = data[, c("state", "E", "P","lat", "lon")], ntree = 500)
#
#set.seed(rn)
#pred_rf = predict(SDM2,new=dataProj,"prob")
#write.table(pred_rf,"../data/pred_states_randomForest.txt")
#
#pred_class = predict(SDM2,new=dataProj,"response")
#(HK2 = HK(pred_class, dataProj$st0)) 
#(HK2 = HK(pred_class, dataProj$st1)) 
#


#---------------------------------------------------------------------
#-- graphs

#------------------------

#for (i in 1:2)
#{
#if(i==1) {pred_proba = pred_rf; nam = "randomFor"}
#if(i==2) {pred_proba = pred_multi; nam = "multinom"}
#
#pdf(paste("../figures/SDM_temperature_",nam, ".pdf", sep = ""), width = 8, height= 6)
#
#par(mfrow = c(2,2))
#plot(lowess(cbind(dataProj$E,pred_proba[,"B"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
#plot(lowess(cbind(dataProj$E,pred_proba[,"T"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
#plot(lowess(cbind(dataProj$E,pred_proba[,"M"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
#plot(lowess(cbind(dataProj$E,pred_proba[,"R"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))
#
#par(mfrow = c(2,2))
#plot(dataProj$E,pred_proba[,"B"], pch=4, cex=.8, xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
#plot(dataProj$E,pred_proba[,"T"],pch=4, cex=.8,  xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
#plot(dataProj$E,pred_proba[,"M"], pch=4, cex=.8,  xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
#plot(dataProj$E,pred_proba[,"R"], pch=4, cex=.8, xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))
#
#
#par(mfrow = c(2,2))
##boxplot(pred_proba)
#boxplot(split(pred_proba[,"B"], ifelse(dataProj$st0=="B", "observed B", " observed other")), main="Boreal", ylab = "predicted B probability", notch=T)
#boxplot(split(pred_proba[,"T"], ifelse(dataProj$st0=="T", "observed T", " observed other")), main="Temperate", ylab = "predicted T probability", notch=T)
#boxplot(split(pred_proba[,"M"], ifelse(dataProj$st0=="M", "observed M", " observed other")), main="Mixed", ylab = "predicted M probability", notch=T)
#boxplot(split(pred_proba[,"R"], ifelse(dataProj$st0=="R", "observed R", " observed other")), main="Early succession (R)", ylab = "predicted R probability", notch=T)
#
#dev.off()
#}
##------------------------

# compare models
#pred_multi = predict(SDM1,new=dataProj,"prob")
#set.seed(rn)
#pred_rf = predict(SDM2,new=dataProj,"prob")
#
#par(mfrow = c(2,2))
#plot(lowess(cbind(pred_rf[,"B"], pred_multi[,"B"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Boreal")
#abline(0,1,lty=2)
#plot(lowess(cbind(pred_rf[,"T"], pred_multi[,"T"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Temperate")
#abline(0,1,lty=2)
#plot(lowess(cbind(pred_rf[,"M"], pred_multi[,"M"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Mixed")
#abline(0,1,lty=2)
#plot(lowess(cbind(pred_rf[,"R"], pred_multi[,"R"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Early succession")
#abline(0,1,lty=2)
#






