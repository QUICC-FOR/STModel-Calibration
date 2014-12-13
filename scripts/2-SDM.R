rm(list=ls())
# Load data
#dataProj = as.data.frame(read.table("../data/data_pairs_filter.txt"))
#dataProj$E = scale(dataProj$annual_mean_temp)
#dataProj$P = scale(dataProj$annual_pp)

#data = read.csv("../data/statesFourState.csv")
data = read.csv("~/Documents/GitHub/STModel-Data/out_files/statesFourState.csv")
head(data)
dim(data)

setwd("~/Documents/GitHub/STModel-Calibration/scripts")


# ----------------------
### choice of variables
# ----------------------

varCor = cor(data[,-c(1:5)]) # check correlation between variables
varCor

# acp
library("ade4")
var.pca = dudi.pca(data[,-c(1:5)], scannf=FALSE, nf = 5)
var.pca$eig / sum(var.pca$eig)

s.corcircle(var.pca$co, clab = 0.5)
contrib = inertia.dudi(var.pca, row = FALSE, col = TRUE)$col.abs

nonCorVars = intersect(names(varCor[which(abs(varCor[,"annual_mean_temp"])<0.7),"annual_mean_temp"]), names(varCor[which(abs(varCor[,"pp_seasonality"])<0.7),"pp_seasonality"]))

contrib[nonCorVars,]

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"pp_warmest_quarter"])<0.7),"pp_warmest_quarter"]))

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"mean_diurnal_range"])<0.7),"mean_diurnal_range"]))

nonCorVars = intersect(nonCorVars, names(varCor[which(abs(varCor[,"annual_pp"])<0.7),"annual_pp"]))

selectedVars = c("annual_mean_temp", "pp_seasonality", "pp_warmest_quarter", "mean_diurnal_range","annual_pp", "mean_temperatre_wettest_quarter")

selectedVars_CCbio = c("annual_mean_temp", "pp_seasonality", "min_temp_coldest_period", "mean_temp_coldest_quarter", "mean_temperatre_wettest_quarter","annual_pp")

varCor2 = cor(data[, selectedVars])

datSel = data[,c("state",selectedVars_CCbio)]


# ----------------------
# Clean data
# ---------------------

# Clean Undefined state
str(datSel)
datSel_wo_U <- subset(datSel, state != "U")
datSel_wo_U$state <- droplevels(datSel_wo_U$state)

###########################################################################################################

# ----------------------
### Explo Data et PCA Steve
# ----------------------

require(reshape2)
require(ggplot2)

ggdata <- melt(datSel_wo_U,id=c("state"))
ggdata_all <- data[,5:34]

# Histograme
hists = ggplot(ggdata) + geom_histogram(aes(x=value,fill=state)) + facet_grid(state~variable,scales="free") + scale_fill_brewer(palette="Accent","State")
ggsave(hists,file="../figures/explo_hist_selVars_by_state.pdf",width=15,height=7)

#install.packages("GGally")
require(GGally)

ggplot <- function(...) ggplot2::ggplot(...) + scale_color_brewer(palette="Accent")

ggpairs(ggdata_all,
        columns=2:ncol(ggdata_all),
        lower = list(continuous = "density"), # data.frame with variables
        title="Climate exploration by state", # title of the plot
        colour = "state") # aesthetics, ggplot2 style

dev.copy2pdf(height=16,width=18,out.type="pdf")
dev.off()

### PCA
# See histogram figures

library(FactoMineR)

result.acp <- PCA(data[,-c(1:5)], scale.unit = TRUE)
summary(result.acp)

result.acp
result.acp$eig

plot(result.acp)
plot(result.acp, choix="var", select="contrib 10")
plot(result.acp, choix = "var", select="contrib 6")
plot(result.acp, choix = "var", select="contrib 3")

var.coord <- result.acp$var$xy.coords
eig <- result.acp$eig
var.contrib <- sweep(var.coord, 2, sqrt(eig[1:ncol(var.coord), 1]), FUN="/")


### Broken stick
(ev <- result.acp$eig[,1])
ev[ev > mean(ev)]

n <- length(ev)
bsm <- data.frame(j = seq(1:n), p = 0)
bsm$p[1] <- 1/n
for (i in 2:n) {
  bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
}
bsm$p <- 100 * bsm$p/n

layout(matrix(c(1:2), 1, 2))
barplot(result.acp$eig[,1],main="Eigenvalues",names.arg=1:nrow(result.acp$eig))
abline(h=mean(ev), col="red")
legend("topright", "Moyenne valeurs propres", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=TRUE,
    main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% eigenvalue", "Broken stick model"),
    pch=15, col=c("bisque",2), bty="n")

# ----------------------
### Discriminante analyze
# ----------------------

require(MASS)

dat.dfa <-datSel_wo_U[,c("state",selectedVars)]

state.dfa <- lda(state~., dat.dfa)
state.dfa

clim <- as.matrix(dat.dfa[,-1])
state <- dat.dfa[, 1]
test <- manova(clim ~ state)
summary(test, test = "Wilks")

plot(state.dfa, dimen=1, type="both", cex=1.2)
plot(state.dfa, abbrev=TRUE)

library(klaR)
partimat(state~., data = dat.dfa, method = "lda")
partimat(state ~ ., data = dat.dfa, method = "lda",     plot.matrix = TRUE, imageplot = FALSE)


###########################################################################################################

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

SDM1.a = multinom(state ~ . + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2), data = calib, maxit =1500)

save(SDM1.a,sampl,file= "Multinom_6vars_version_a.rObj")
rm(SDM1.a)


SDM1.b = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2), data = calib, maxit =1500)

save(SDM1.b,sampl,file= "Multinom_6vars_version_b.rObj")
rm(SDM1.b)

SDM1.ccbio = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(min_temp_coldest_period^2) + I(mean_temp_coldest_quarter^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2), data = calib, maxit =1500)

save(SDM1.ccbio,sampl,file= "Multinom_6vars_version_b.rObj")
rm(SDM1.ccbio)

SDM1.c = multinom(state ~ . + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2) + I(annual_mean_temp^3) + I(pp_seasonality^3) + I(pp_warmest_quarter^3) + I(mean_diurnal_range^3) +I(annual_pp^3) + I(mean_temperatre_wettest_quarter^3), data = calib, maxit =1500)

save(SDM1.c,sampl,file= "Multinom_6vars_version_c.rObj")
rm(SDM1.c)

SDM1.d = multinom(state ~ .^2 + I(annual_mean_temp^2) + I(pp_seasonality^2) + I(pp_warmest_quarter^2) + I(mean_diurnal_range^2) +I(annual_pp^2) + I(mean_temperatre_wettest_quarter^2) + I(annual_mean_temp^3) + I(pp_seasonality^3) + I(pp_warmest_quarter^3) + I(mean_diurnal_range^3) +I(annual_pp^3) + I(mean_temperatre_wettest_quarter^3), data = calib, maxit =1500)

save(SDM1.d,sampl,file= "Multinom_6vars_version_d.rObj")
rm(SDM1.d)

###### Evaluate models

load('../data/Multinom_6vars_version_a.rObj')
load('../data/Multinom_6vars_version_b.rObj')
load('../data/Multinom_6vars_version_c.rObj')
load('../data/Multinom_6vars_version_d.rObj')
load('../data/Multinom_6vars_version_d.rObj')
load('../data/RandomForest_7vars.rObj')

summary(SDM1.a,Wald=TRUE)
summary(SDM1.b,Wald=TRUE)
summary(SDM1.c,Wald=TRUE)
summary(SDM1.d,Wald=TRUE)


predA = predict(SDM1.a,new=datSel_wo_U,"class")
predB = predict(SDM1.b,new=datSel_wo_U,"class")
predC = predict(SDM1.c,new=datSel_wo_U,"class")
predD = predict(SDM1.d,new=datSel_wo_U,"class")
predccbio = predict(SDM1.ccbio,new=datSel_wo_U,"class")


(HKA = table(predccbio[-sampl], valid$state))
(HKB = HK(predB[-sampl], valid$state))
(HKC = HK(predC[-sampl], valid$state))
(HKD = HK(predD[-sampl], valid$state))
(HKD = HK(predccbio[-sampl], valid$state))

climate_grid <- read.csv("~/Documents/Maitrise/Analyse/dom_ouranos/dom_climate_grid.csv")
#climate_grid <- subset(climate_grid,lat>=rg_lat[1] & lat<=rg_lat[2])
#climate_grid <- subset(climate_grid,lon>=rg_lon[1] & lon<=rg_lon[2])

# Prob SDM
prob_RF = cbind(climate_grid,as.data.frame(predict(SDM2,new=climate_grid,"prob", OOB=TRUE)))
prob_Multi_a = cbind(climate_grid,as.data.frame(predict(SDM1.a,new=climate_grid,"prob")))
prob_Multi_b = cbind(climate_grid,as.data.frame(predict(SDM1.b,new=climate_grid,"prob")))
prob_Multi_b = cbind(climate_grid,as.data.frame(predict(SDM1.bstep,new=climate_grid,"prob")))
prob_Multi_c = cbind(climate_grid,as.data.frame(predict(SDM1.c,new=climate_grid,"prob")))
prob_Multi_d = cbind(climate_grid,as.data.frame(predict(SDM1.d,new=climate_grid,"prob")))

###########################################################################################################



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




