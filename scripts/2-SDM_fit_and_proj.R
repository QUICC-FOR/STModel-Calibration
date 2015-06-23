rm(list=ls())
# ----------------------
# Load data
# ---------------------

load("../data/transitions_r1.rdata")
head(stateData)
dim(stateData)
str(stateData)
#
#jpeg("../figures/plots_states.jpeg", height=5000, width=5000, res=600)
#colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
#plot(stateData$lat~stateData$lon, pch = 20, cex=.2, asp = 1, col = colo[stateData$state])
#dev.off()
#
#jpeg("../figures/plots_Rstates.jpeg", height=5000, width=5000, res=600)
#plot(stateData$lat~stateData$lon, pch = 20, cex=.2, asp = 1, col = ifelse(stateData$state=="R", 1, "grey"), main="R states")
#dev.off()


# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
dat_wo_U <- subset(stateData, state != "U")
#dat_wo_U$state <- droplevels(dat_wo_U$state)

#add data
#addData = read.csv("../data/stm_soil.csv")
#head(addData)
#summary(addData)
##
#dat = merge(dat_wo_U, addData, by = "plot_id", all.x=TRUE, all.y=FALSE)
#dim(dat)
#dim(dat_wo_U)

# check NA and ranges
summary(dat_wo_U)

# 
dat = dat_wo_U

# subset 10 degree and remove NA
select = unique(dat$plot_id[which(dat$annual_mean_temp<=10 & !is.na(dat$soil)& !is.na(dat$slp))])
dat_subset10 = dat[dat$plot_id %in% select,]

summary(dat_subset10)

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


selectedVars = c("annual_mean_temp", "tot_annual_pp", "mean_diurnal_range","ph_2cm", "slp")

#
#library("ade4")
#var.pca = dudi.pca(dat_subset10[,-c(1:9)], scannf=FALSE, nf = 5)
#round((var.pca$eig / sum(var.pca$eig))*100)
##
#inertia.dudi(var.pca, row=F,col=T)
##
###pdf("../figures/PCA_selectedVariables.pdf")
#s.arrow(var.pca$co[selectedVars,], clab = .6, xlim = c(-2,2), sub = "axe 1: 43% ; axe 2: 27 %")
#s.arrow(var.pca$co[selectedVars,], clab = .6, xlim = c(-2,2),yax=3, sub = "axe 1: 43% ; axe 3: 9 %")
#dev.off()

#varCor2 = cor(dat_subset10[, selectedVars])
#varCor2


# ----------------------

datSel = dat_subset10[,c("state",selectedVars, "lat","lon")]
datSel$state = as.factor(datSel$state)
#datSel$ph_2cm = scale(datSel$ph_2cm)
#datSel$lat = scale(datSel$lat)
#datSel$lon = scale(datSel$lon)

#rm(stateData, dat_wo_U)

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
select = subsample.temp(datSel$annual_mean_temp, .33)

jpeg("../figures/sdm/subsample_SDM_calib2.jpeg", height=5000, width=5000, res=600)
plot(dat_subset10[,c("lon","lat")], pch = 20, cex=.2, col = "grey")
points(dat_subset10[select,c("lon","lat")], pch = 20, cex=.2, col = 1)
dev.off()

save(select, file = "subsample_SDM_calib2.RData")
load(file = "subsample_SDM_calib2.RData")

#-----------------------

datCal = datSel[select,]
datEval = datSel[-select,]



# ----------------------
# random Forest
# ----------------------


library(randomForest)
rs = runif(1,0,1)

## calib / eval
# ----------------------

set.seed(rs)
#SDM2_cal = randomForest(state ~ annual_mean_temp + tot_annual_pp + mean_diurnal_range + ph_2cm + slp +lon +lat , data = datCal, ntree = 500)
SDM2_cal = randomForest(state ~ annual_mean_temp + tot_annual_pp + mean_diurnal_range + ph_2cm + slp , data = datCal, ntree = 500)
#save(SDM2_cal,rs,file= "../data/RandomForest_cal.RData")
save(SDM2_cal,rs,file= "../data/RandomForest_cal_woLatLon.RData")
(imp = importance(SDM2_cal))

# evaluation random Forest
set.seed(rs)
pred2 = predict(SDM2_cal,new=datEval,"response", OOB=TRUE)
(HK2 = HK(pred2, datEval$state)) 
#error rate
(tt2 = table(pred2, datEval$state))
1-sum(diag(tt2))/sum(tt2)

##  complete
# ----------------------

set.seed(rs)
#SDM2 = randomForest(state ~ . , data = datSel, ntree = 500)
SDM2 = randomForest(state ~ annual_mean_temp + tot_annual_pp + mean_diurnal_range + ph_2cm + slp , data = datSel, ntree = 500)

#save(SDM2,rs,file= "../data/RandomForest_complete.RData")
save(SDM2,rs,file= "../data/RandomForest_complete_woLatLon.RData")

(imp = importance(SDM2))
#
## evaluation random Forest
set.seed(rs)
pred2 = predict(SDM2,new=datSel,"response", OOB=TRUE)
(HK2 = HK(pred2, datSel$state)) #
#error rate
tt2 = table(pred2, datSel$state)
1-sum(diag(tt2))/sum(tt2)

#
##--
#
dat_rc = data.frame(annual_mean_temp = seq(-5, 3, length.out = 200), tot_annual_pp = rep(0, 200), mean_diurnal_range=rep(0, 200),ph_2cm= rep(0, 200),slp= rep(0, 200), lat=rep(0, 200),lon= rep(0,200) )
#
#
set.seed(rs)
pred1_rc =  predict(SDM2,new=dat_rc,"prob")
#
###--- 
pdf("../figures/sdm/randomForest predictions_woLatLon.pdf", width = 8, height = 8)
colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
plot(pred1_rc[,"R"]~seq(-5, 3, length.out = 200), type = "l", col = colo["R"], lwd = 1.5, ylim = c(0,1), xlab = "temperature gradient", ylab = "proportion of states", main = "random forest")
lines(pred1_rc[,"T"]~seq(-5, 3, length.out = 200), col = colo["T"], lwd = 1.5)
lines(pred1_rc[,"B"]~seq(-5, 3, length.out = 200), col = colo["B"], lwd = 1.5)
lines(pred1_rc[,"M"]~seq(-5, 3, length.out = 200), col = colo["M"], lwd = 1.5)
legend("topright", bty = "n", col = colo, legend = names(colo), lwd = 2)
dev.off()
###--
jpeg("../figures/sdm/randomForest predictions_map_woLatLon.jpeg", height=5000, width=5000, res=600)
plot(dat_subset10[,"lon"], dat_subset10[,"lat"], cex = .1, pch = 20, xlab = "longitude", ylab = "latitude", asp = 1, col = colo[as.character(pred2)])
title("random forest")
dev.off()
#

# ----------------------
# projection
# ---------------------


dataProj = transitionData
head(dataProj)
dim(dataProj)
str(dataProj)

dim(dataProj)


# rescale
dataRescaledProj = t(apply(dataProj[,c(selectedVars, "lat","lon")], 1, function(x) {(x-vars.means)/vars.sd}))
dataRescaledProj = cbind(dataProj[, c("plot","year1","year2","state1","state2")], dataRescaledProj)

head(dataRescaledProj)
summary(dataRescaledProj)
dim(dataRescaledProj)


# random Forest
# ----------------------
#load("../data/RandomForest_complete.RData")
#set.seed(rs)
#proj2 = predict(SDM2,new=dataRescaledProj,"prob", OOB=TRUE)
#proj2 = data.frame(cbind(proj2, dataRescaledProj))
#head(proj2)
#dim(proj2)
#summary(proj2)
### attention there are NAs !
## sauvegarde
#write.table(proj2, file = "../data/projection_rf_complete.txt", quote=F, row.names=FALSE)


load("../data/RandomForest_complete_woLatLon.RData")
set.seed(rs)
proj1 = predict(SDM2,new=dataRescaledProj,"prob", OOB=TRUE)
proj1 = data.frame(cbind(proj1, dataRescaledProj))
head(proj1)
dim(proj1)
summary(proj1)
## attention there are NAs !
# sauvegarde
write.table(proj1, file = "../data/projection_rf_complete_woLatLon.txt", quote=F, row.names=FALSE)




Temp.lim = c((-5-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], (10-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"])
Temp.ax = function(x)
{
temp = dataRescaledProj[,"annual_mean_temp"]*vars.sd["annual_mean_temp"]+vars.means["annual_mean_temp"]
axis(1, at = seq((round(min(temp))-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], (round(max(temp))-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], l=30), labels = seq(round(min(temp)), round(max(temp)), l = 30))
}

colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))

#jpeg("../figures/sdm/randomForest predictions_gradient.jpeg", height=5000, width=5000, res=600)
#
#par(mfrow = c(1,1))
#plot(proj2$B ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "B neighborhood", cex = .2, pch = 20, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col = colo["B"])
#Temp.ax()
#points(proj2$T ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "T neighborhood", cex = .2, pch = 20, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col=colo["T"])
#Temp.ax()
#points(proj2$M ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "M neighborhood", cex = .2, pch = 20, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col = colo["M"])
#Temp.ax()
#points(proj2$R ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "R neighborhood", cex = .2, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col = colo["R"])
#legend("topleft", bty = "n", col = colo, legend = names(colo), pch = 20)
##Temp.ax()
#dev.off()

jpeg("../figures/sdm/randomForest predictions_gradient_woLatLon.jpeg", height=5000, width=5000, res=600)

par(mfrow = c(1,1))
plot(proj1$B ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "B neighborhood", cex = .2, pch = 20, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col = colo["B"])
Temp.ax()
points(proj1$T ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "T neighborhood", cex = .2, pch = 20, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col=colo["T"])
Temp.ax()
points(proj1$M ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "M neighborhood", cex = .2, pch = 20, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col = colo["M"])
Temp.ax()
points(proj1$R ~ dataRescaledProj[,"annual_mean_temp"], xlab = "temperature", ylab = "R neighborhood", cex = .2, xaxt="n", xlim = Temp.lim, ylim = c(0,.6), col = colo["R"])
legend("topleft", bty = "n", col = colo, legend = names(colo), pch = 20)
#Temp.ax()
dev.off()

# random Forest - figure dans plan clim
# ----------------------
#load("../data/RandomForest_complete.RData")
load("../data/RandomForest_complete_woLatLon.RData")
load("scale_info.Robj")
## Temp -4 à 10
## PP 700 à 1300
Trange = c(-4, 10)
Tticks = max(Trange)-min(Trange)+1
PPrange = c(700, 1300)
Tbounds = (Trange - vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"]
PPbounds = (PPrange - vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"]
PPticks = (max(PPrange)-min(PPrange))/100 +1

tpseq=seq(Tbounds[1],Tbounds[2],l=200)
ppseq=seq(PPbounds[1],PPbounds[2],l=200)

scaled.axis<- function(){
temp = tpseq*vars.sd["annual_mean_temp"]+vars.means["annual_mean_temp"]
precip = ppseq*vars.sd["tot_annual_pp"]+vars.means["tot_annual_pp"]
axis(1, at = seq(Tbounds[1], Tbounds[2], l=Tticks), labels = seq(Trange[1], Trange[2], l = Tticks))
axis(2, at = seq(PPbounds[1],PPbounds[2], l=PPticks), labels = seq(PPrange[1], PPrange[2], l = PPticks))
}

ENV = expand.grid(TP =tpseq , PP = ppseq)

ENV.df = data.frame(annual_mean_temp = ENV$TP, tot_annual_pp = ENV$PP, mean_diurnal_range=rep(0, nrow(ENV)),ph_2cm= rep(0,  nrow(ENV)),slp= rep(0,  nrow(ENV)), lat=rep(0,  nrow(ENV)),lon= rep(0, nrow(ENV)) )

set.seed(rs)
proj2 = predict(SDM2,new=ENV.df,"prob", OOB=TRUE)

#---
tomap = as.factor(colnames(proj2)[apply(proj2, 1, which.max)])
colo = c(M = "lightgreen", B = rgb(44,133,113,maxColorValue=255), T = rgb(245,172,71,maxColorValue=255), R = rgb(218,78,48,maxColorValue=255))


jpeg("../figures/sdm/RF_neigborhood_plan_clim_woLatLon.jpeg", height=3000, width=5000, res=600)
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("",cex=2)
legend("center",legend = levels(tomap),fill = colo[levels(tomap)],bty = "n", cex = 0.8, ncol =2)
par(mar=c(5,5,0,2))

image(x=tpseq, y=ppseq, z = matrix(as.numeric(tomap), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Annual mean temperature (°C)", ylab = "Annual precipitations (mm)", col = colo[levels(tomap)], main = "", xaxt = "n", yaxt="n")
scaled.axis()

dev.off()



