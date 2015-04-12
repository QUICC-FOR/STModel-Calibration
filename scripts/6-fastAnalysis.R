rm(list = ls())

sdm = "rf"
sdm = "cst"
propData = 0.337
option = "_less"
option = "5y"
option = "_less5y"
option = ""
option = "2"
option = "_less5y2"

#--
veget_pars = read.table(paste("../estimated_params/GenSA_initForFit_", sdm, "_",propData, option, ".txt", sep=""))
load(paste("initForFit_",sdm, "_", propData, ".RData", sep = ""))
load(paste("../estimated_params/GenSA_initForFit_", sdm, "_", propData, option, ".RData", sep = ""))
mat = read.table(paste("../estimated_params/traceMat_initForFit_", sdm, "_", propData, option, ".trMat", sep=""), h=T)
#--
print(fit)
#--
tpseq=seq(-4,2,l=100)
ppseq=seq(-2,5,l=100)
load("scale_info.Robj")

scaled.axis<- function(){
temp = tpseq*vars.sd["annual_mean_temp"]+vars.means["annual_mean_temp"]
precip = ppseq*vars.sd["tot_annual_pp"]+vars.means["tot_annual_pp"]
axis(1, at = seq((round(min(temp))-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], (round(max(temp))-vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"], l=18), labels = seq(round(min(temp)), round(max(temp)), l = 18))
axis(2, at = seq((700-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"], (1800-vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"], l=12), labels = seq(700, 1800, l = 12))
}

#---
pars = as.numeric(veget_pars[,2])
names(pars) = veget_pars[,1]
pars
estim.pars$value
#-----
# diagnostic convergence
#----
#jpeg(paste("../figures/diagnostic",sdm, "_",propData, option,".jpeg", sep=""), height=5000, width=5000, res=600)
par(mfrow = c(2,2), mar = c(4, 4, 2,1))
mat$number = 1:nrow(mat)
plot(nb.steps~number, data = mat, xlab = "anneal time", ylab = "nb.steps", cex=.2, type = "l")
plot(temperature~number, data = mat, cex=.2, type = "l")
plot(function.value~number, data = mat, cex=.2, type = "l")
plot(current.minimum~number, data = mat, cex = .2, type="l")
#dev.off()
#-----
# check estimated params and bounds
#----
parsBounds = data.frame(hi=unlist(par_hi),lo= unlist(par_lo))
parsBounds$estim = rep(0, nrow(parsBounds))
parsBounds[names(pars), "estim"] = pars
par(mfrow = c(1,1), mar = c(2, 2, 2,1))
plot(parsBounds$hi, ylim = c(min(parsBounds), max(parsBounds)), pch=15, cex=.5)
points(parsBounds$lo, col=1, pch=15, cex=.5)
points(parsBounds$estim, col=2, pch=19, cex=.5)
title("estimated params and given bounds")

pp = parsBounds
plot(parsBounds$estim, pp$estim)
#------
# look at transition probabilities
#------

ENV = expand.grid(TP =tpseq , PP = ppseq)
ENV1 = ENV$TP
ENV2 = ENV$PP

e7 = 0
params=parsBounds$estim
names(params) = rownames(parsBounds)

    logit_alphab 	= params["ab0"] + params["ab1"]*ENV1 + params["ab2"]*ENV2 + params["ab3"]*ENV1^2 + params["ab4"]*ENV2^2 + params["ab5"]*ENV1^3 + params["ab6"]*ENV2^3
    logit_alphat 	= params["at0"] + params["at1"]*ENV1 + params["at2"]*ENV2 + params["at3"]*ENV1^2 + params["at4"]*ENV2^2 + params["at5"]*ENV1^3 + params["at6"]*ENV2^3
    logit_betab 	= params["bb0"] + params["bb1"]*ENV1 + params["bb2"]*ENV2 + params["bb3"]*ENV1^2 + params["bb4"]*ENV2^2 + params["bb5"]*ENV1^3 + params["bb6"]*ENV2^3
    logit_betat 	= params["bt0"] + params["bt1"]*ENV1 + params["bt2"]*ENV2 + params["bt3"]*ENV1^2 + params["bt4"]*ENV2^2 + params["bt5"]*ENV1^3 + params["bt6"]*ENV2^3
    logit_theta	= params["th0"] + params["th1"]*ENV1 + params["th2"]*ENV2 + params["th3"]*ENV1^2 + params["th4"]*ENV2^2 + params["th5"]*ENV1^3 + params["th6"]*ENV2^3
    logit_thetat	= params["tt0"] + params["tt1"]*ENV1 + params["tt2"]*ENV2 + params["tt3"]*ENV1^2 + params["tt4"]*ENV2^2 + params["tt5"]*ENV1^3 + params["tt6"]*ENV2^3
    logit_eps 	= params["e0"]  + params["e1"]*ENV1 + params["e2"]*ENV2  + params["e3"]*ENV1^2 + params["e4"]*ENV2^2 + params["e5"]*ENV1^3 + params["e6"]*ENV2^3 
    #e7*EB
 
     logit_reverse <- function(x)
    {
    expx = ifelse(exp(x)==Inf, .Machine$double.xmax, exp(x))
    expx/(1+expx)
    }
    
macroPars = data.frame(alphab = logit_reverse(logit_alphab), 
alphat = logit_reverse(logit_alphat),
betab = logit_reverse(logit_betab),
betat = logit_reverse(logit_betat),
theta = logit_reverse(logit_theta),
thetat = logit_reverse(logit_thetat),
eps = logit_reverse(logit_eps))

summary(macroPars)
#---
#
#colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
pal = colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")


#jpeg(paste("../figures/estim_pars_",sdm, "_",propData, option,".jpeg", sep=""), height=3000, width=5000, res=600)
#
par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(macroPars))
{
image(x=tpseq, y=ppseq, z = matrix(macroPars[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(macroPars)[i], xaxt = "n", yaxt="n")
contour(x=tpseq, y=ppseq, z = matrix(macroPars[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
scaled.axis()
}

#dev.off()

### 
#------
# probabilites transition (no neighborhood restriction)
#------
head(ENV)
colnames(ENV) = c("annual_mean_temp", "tot_annual_pp")
ENV$mean_diurnal_range = rep(0, nrow(ENV))
ENV$pp_warmest_quarter = rep(0, nrow(ENV))
ENV$pp_wettest_period = rep(0, nrow(ENV))
ENV$mean_temp_wettest_quarter = rep(0, nrow(ENV))
ENV$mean_temp_driest_quarter = rep(0, nrow(ENV))

load("../data/RandomForest_complete.rObj")
library(randomForest)
set.seed(rs)
proj2 = predict(SDM2,new=ENV,"prob", OOB=TRUE)

nB = proj2["B"]
nT = proj2["T"]
nM = proj2["M"]

pTransitions = data.frame(pRT = macroPars$alphat*(nT+nM)*(1-macroPars$alphab*(nB+nM)),
pRB = macroPars$alphab*(nB+nM)*(1-macroPars$alphat*(nT+nM)),
pRM = macroPars$alphat*(nT+nM)*macroPars$alphab*(nB+nM),
pMT = macroPars$theta*macroPars$thetat,
pMB = macroPars$theta*(1-macroPars$thetat),
pTM = macroPars$betab,
pBM = macroPars$betat,
eps = macroPars$eps)

summary(pTransitions)

par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(pTransitions))
{
image(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(pTransitions)[i], xaxt = "n", yaxat="n")
contour(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
scaled.axis()
}

#------
# equilibrium analysis
##--------------------
pars= data.frame(
alphaB = macroPars$alphab,
alphaT = macroPars$alphat,
betaB = macroPars$betab,
betaT= macroPars$betat,
theta = macroPars$theta,
thetaT = macroPars$thetat,
eps = macroPars$eps)


model <- function(ti, states, parms)
{
with(as.list(c(states, parms)), {
R = 1 - T - B - M
# Differential equations describing the dynamics of the state variables
dB = theta*(1-thetaT)*(1-eps)*M + alphaB*(M+B)*(1-alphaT*(T+M))*R - betaT*(T+M)*(1-eps)*B - eps*B
dT = theta*thetaT*(1-eps)*M + alphaT*(M+T)*(1-alphaB*(B+M))*R - betaB*(B+M)*(1-eps)*T - eps*T
dM = betaT*(T+M)*(1-eps)*B + betaB*(B+M)*(1-eps)*T - theta*(1-eps)*M  - eps*M
return(list(c(dT, dB, dM)))
	})
}

#library(rootSolve)
#
#
#eq.winner = function(pars, model)
#{
#fparams = pars
#init = c(T = 0.33, B = 0.33, M=0.33)
#(eq = stode(y = init, func = model, parms = fparams, positive = TRUE)$y)
#return(names(which.max(eq)))
#}
#
##eq = t(apply(pars, 1, eq.winner, model =model))
##eq = as.factor(eq)
###-----
##jpeg(paste("../figures/equilibrium_map_", sdm, "_",propData, option, ".jpeg", sep=""), height=3000, width=3000, res = 300)
##
#colo = c(M = "orange", B = "darkgreen", T = "lightgreen")
#
#layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
#
#par(mar=c(0,0,0,0))
#plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#title("",cex=2)
#legend("center",legend = levels(eq),fill = colo[levels(eq)],bty = "n", cex = 0.8, ncol =2)
#par(mar=c(5,5,0,2))
#
#image(x=tpseq, y=ppseq, z = matrix(as.numeric(eq), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Annual mean temperature (°C)", ylab = "Annual precipitations (mm)", col = colo[levels(eq)], main = "", xaxt = "n", yaxt="n")
#scaled.axis()
##
##dev.off()


### 
#------
# invasibility
#------



eq.eigen = function(pars, model, eqState = "B")
{
fparams = pars
# Solve the model for the case where T/B and M are at 0
# and B/T = 1 - eps/alphaB/T (ie a l'équilibre si T/B=0 et M=0)
if(eqState=="B") init = c(T = 0, B = as.numeric(1-fparams["eps"]/fparams["alphaB"]), M=0)
if(eqState=="T") init = c(T = as.numeric(1-fparams["eps"]/fparams["alphaT"]), B=0, M=0)
if(is.infinite(sum(init))) init = c(T=0, B = 0, M=0)
(eq = stode(y = init, func = model, parms = fparams, positive = TRUE)$y)
(jacob = jacobian.full(y = eq, func= model, parms=fparams))
#print(jacob)
return(eigen(jacob)$values)
}

invT = t(apply(pars, 1, eq.eigen, model =model, eqState="B"))
invT.class = ifelse(invT[,3]<0, "stable", "invadeM")
invT.class = ifelse(invT[,1]>0, "invadeT", invT.class)
T.growth = apply(pars, 1, function(pars){1-pars["eps"]/pars["alphaT"]})
invT.class[T.growth<0] = "stable"
invT.class = as.factor(invT.class)

invB = t(apply(pars, 1, eq.eigen, model =model, eqState="T"))
invB.class = ifelse(invB[,3]<0, "stable", "invadeM")
invB.class = ifelse(invB[,2]>0, "invadeB", invB.class)
B.growth = apply(pars, 1, function(pars){1-pars["eps"]/pars["alphaB"]})
invB.class[B.growth<0] = "stable"
invB.class = as.factor(invB.class)

##--
coexist = numeric(length(invT))
coexist[invT.class=="stable" & invB.class=="stable"] = "AltSS"
coexist[invT.class=="stable" & invB.class=="invadeB"] = "B wins"
coexist[invT.class=="stable" & invB.class=="invadeM"] = "B dominates (with M)"

coexist[invT.class=="invadeM" & invB.class=="invadeM"] = "M wins (coexistence)"
coexist[invT.class=="invadeM" & invB.class=="invadeB"] = "B dominates (with M)"
coexist[invT.class=="invadeM" & invB.class=="stable"] = "T dominates (with M)"

coexist[invT.class=="invadeT" & invB.class=="stable"] = "T wins"
coexist[invT.class=="invadeT" & invB.class=="invadeM"] = "T dominates (with M)"
coexist[invT.class=="invadeT" & invB.class=="invadeB"] = "??"
table(coexist)
coexist = as.factor(coexist)
levels(coexist)
##---
colo = c(stable = "grey", invadeM = "blue", invadeT = "lightgreen", invadeB = "darkgreen")
par(mfrow = c(1,2), mar = c(4,4,1,1), cex=0.8)

image(x=tpseq, y=ppseq, z = matrix(as.numeric(invT.class), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = colo[levels(invT.class)], main = "invasibility of T (when B every where)", xaxt = "n", yaxt="n")
scaled.axis()
legend("topleft",legend = levels(invT.class),fill = colo[levels(invT.class)],bty = "n", cex = 0.8, ncol =1)

image(x=tpseq, y=ppseq, z = matrix(as.numeric(invB.class), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = colo[levels(invB.class)], main = "invasibility of B (when T everywhere)", xaxt = "n", yaxt="n")
scaled.axis()
legend("topleft",legend = levels(invB.class),fill = colo[levels(invB.class)],bty = "n", cex = 0.8, ncol =1)


##-----
#jpeg(paste("../figures/equilibrium_map_", sdm, "_",propData, option, ".jpeg", sep=""), height=3000, width=3000, res = 300)
#
colo = c(AltSS = "pink", 'M wins (coexistence)' = "lightblue", 'B wins' = "darkgreen", 'B dominates (with M)' ="grey30", 'T dominates (with M)' = "grey70", 'T wins' = "lightgreen", '??' = 1)

layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("",cex=2)
legend("center",legend = levels(coexist),fill = colo[levels(coexist)],bty = "n", cex = 0.8, ncol =2)
par(mar=c(5,5,0,2))

image(x=tpseq, y=ppseq, z = matrix(as.numeric(coexist), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Annual mean temperature (°C)", ylab = "Annual precipitations (mm)", col = colo[levels(coexist)], main = "", xaxt = "n", yaxt="n")
scaled.axis()
#
#dev.off()

estim.pars$value