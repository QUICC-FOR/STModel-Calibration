rm(list=ls())
# Load data
dataProj = as.data.frame(read.table("../data/data_reshaped_RBTM.txt"))
dataProj$E = scale(dataProj$annual_mean_temp)
dataProj$P = scale(dataProj$annual_pp)

data = as.data.frame(read.table("../data/data_allyears_RBTM.txt"))
#data = data[-which(data$lat ==0), ] # remove wired localized plots
data$E = scale(data$annual_mean_temp)
data$P = scale(data$annual_pp)
data$state = data$class_final
str(data[,c("E", "P", "class_final")])
cor(data[, c("E", "P", "lat", "lon")]) # check correlation between variables

# cross validation
source("BoulangeatEcoLet2012_fct.r")
sampl = sample(1:nrow(data), 2*nrow(data)/3)
calib = data[sampl,] 
valid = data[-sampl,]

# Run the models
library(nnet)
SDM1 = multinom(state ~ poly(E,3,raw=TRUE) + poly(P,3,raw=TRUE) + E:P + lat*lon, data = calib, maxit =200)

pred1 = predict(SDM1, new=valid,"class")
(HK1 = HK(pred1, valid$state)) 

#--
library(randomForest)
set.seed(23)
SDM2 = randomForest(state ~ . , data = calib[, c("state", "E", "P","lat", "lon")], ntree = 500)

set.seed(23)
pred2 = predict(SDM2,new=valid,"response", OOB=TRUE)
(HK2 = HK(pred2, valid$state)) 


#---------------------------------------------------------------------
# Get predictions
#newdat = expand.grid(E=seq(-1.7,7,0.1),P=seq(700, 1600, 10) )

# multinom
library(nnet)
SDM1 = multinom(state ~ poly(E,3,raw=TRUE) + poly(P,3,raw=TRUE) + E:P + lon* lat , data = data, maxit =200)

pred_multi = predict(SDM1,new=dataProj,"prob")
write.table(pred_multi,"../data/data_pred_states_multinom.txt")

pred_class = predict(SDM1,new=dataProj,"class")
(HK1 = HK(pred_class, dataProj$st0)) 
(HK1 = HK(pred_class, dataProj$st1)) 


# random forest
set.seed(23)
#SDM2 = randomForest(state ~ . , data = data[, c("state", "E", "P","lat", "lon")], ntree = 500)

set.seed(23)
pred_rf = predict(SDM2,new=dataProj,"prob")
write.table(pred_rf,"../data/data_pred_states_randomForest.txt")

pred_class = predict(SDM2,new=dataProj,"response")
(HK2 = HK(pred_class, dataProj$st0)) 
(HK2 = HK(pred_class, dataProj$st1)) 



#---------------------------------------------------------------------
#-- graphs

#------------------------
# pred_proba =  pred_rf
#pred_proba =  pred_multi

pred_proba = (pred_rf + pred_multi)/2

par(mfrow = c(2,2))
plot(lowess(cbind(dataProj$E,pred_proba[,"B"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
plot(lowess(cbind(dataProj$E,pred_proba[,"T"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
plot(lowess(cbind(dataProj$E,pred_proba[,"M"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
plot(lowess(cbind(dataProj$E,pred_proba[,"R"]), f=1/10), type = "l", xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))

par(mfrow = c(2,2))
plot(dataProj$E,pred_proba[,"B"], pch=4, cex=.8, xlab = "Annual temperature", ylab="Boreal probability", ylim = c(0,1))
plot(dataProj$E,pred_proba[,"T"],pch=4, cex=.8,  xlab = "Annual temperature", ylab="Temperate probability", ylim = c(0,1))
plot(dataProj$E,pred_proba[,"M"], pch=4, cex=.8,  xlab = "Annual temperature", ylab="Mixed probability", ylim = c(0,1))
plot(dataProj$E,pred_proba[,"R"], pch=4, cex=.8, xlab = "Annual temperature", ylab="Early succession probability", ylim = c(0,1))

#dev.copy2pdf(file = "../figures/SDM_temperature_randomForest.pdf")
dev.copy2pdf(file = "../figures/SDM_temperature_multinom.pdf")

par(mfrow = c(2,2))
boxplot(split(pred_proba[,"B"], ifelse(dataProj$st0=="B", "State B", "other states")), main="Boreal", ylab = "predicted probability", notch=T)
boxplot(split(pred_proba[,"T"], ifelse(dataProj$st0=="T", "State T", "other states")), main="Temperate", ylab = "predicted probability", notch=T)
boxplot(split(pred_proba[,"M"], ifelse(dataProj$st0=="M", "State M", "other states")), main="Mixed", ylab = "predicted probability", notch=T)
boxplot(split(pred_proba[,"R"], ifelse(dataProj$st0=="R", "State R", "other states")), main="Early succession", ylab = "predicted probability", notch=T)

#dev.copy2pdf(file = "../figures/SDM_validation_randomForest.pdf")
dev.copy2pdf(file = "../figures/SDM_validation_multinom.pdf")

#------------------------

# compare models
pred_multi = predict(SDM1,new=dataProj,"prob")
set.seed(23)
pred_rf = predict(SDM2,new=dataProj,"prob")

par(mfrow = c(2,2))
plot(lowess(cbind(pred_rf[,"B"], pred_multi[,"B"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Boreal")
abline(0,1,lty=2)
plot(lowess(cbind(pred_rf[,"T"], pred_multi[,"T"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Temperate")
abline(0,1,lty=2)
plot(lowess(cbind(pred_rf[,"M"], pred_multi[,"M"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Mixed")
abline(0,1,lty=2)
plot(lowess(cbind(pred_rf[,"R"], pred_multi[,"R"])), type = "l", xlab = "RF", ylab = "multinom", xlim = c(0,1), ylim = c(0,1), main = "Early succession")
abline(0,1,lty=2)







