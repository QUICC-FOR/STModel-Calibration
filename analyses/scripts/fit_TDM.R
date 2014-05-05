setwd("/home/DominiqueGravel/Documents/Students/Supervision/Msc/Vissault_Msc13/multinom")

### FUNCTIONS ###
library(nnet)

# Run the models
fit_multinom = function(data) multiMod = multinom(t1 ~ ., data)
# Load data
data = read.table("states_class_dat.csv",header = T, sep= ",")
data = subset(data,data$t0!="Unclass" & data$t1!="Unclass")
data$transitions = factor(paste(data$t0, data$t1, sep=""))
data$E = data$av_annual_mean_tp
data$E2 = data$E^2

## subset data
dataC = data[which(data$t0=="C"),]
dataD = data[which(data$t0=="D"),]
dataM = data[which(data$t0=="M"),]
dataT = data[which(data$t0=="T"),]

subC = dataC[,c("t1","E","E2")]
subD = dataD[,c("t1","E","E2")]
subM = dataM[,c("t1","E","E2")]
subT = dataT[,c("t1","E","E2")]

modelC = fit_multinom(subC)
modelD = fit_multinom(subD)
modelM = fit_multinom(subM)
modelT = fit_multinom(subT)

Tgrad = seq(-2,6,0.01)

predC = predict(modelC, new = data.frame(E=Tgrad,E2=Tgrad^2),"probs")
predD = predict(modelD, new = data.frame(E=Tgrad,E2=Tgrad^2),"probs")
predM = predict(modelM, new = data.frame(E=Tgrad,E2=Tgrad^2),"probs")
predT = predict(modelT, new = data.frame(E=Tgrad,E2=Tgrad^2),"probs")

x11(height = 6, width = 8)
par(mar=c(5,5,2,1))
plot(Tgrad,predC[,3],type = "l",ylim=c(0,0.5),cex.axis = 1.25, cex.lab = 1.25, xlab = "Température moyenne annuelle", ylab = "Probabilité",lwd = 2,col = "darkgreen")
lines(Tgrad,predM[,2],col = "darkred",lwd = 2)
lines(Tgrad,predT[,2],col = "darkblue",lwd = 2)

plot(Tgrad,predC,3],type = "l")
lines(Tgrad,predM[,4])

