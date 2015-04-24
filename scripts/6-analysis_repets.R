rm(list = ls())

sdm = "rf"
propData = 0.33
ordre = 3
step = 5

(name = paste(sdm,"_", propData, "_", ordre, "_",step, "y",sep=""))


#--
veget_pars = data.frame(matrix(NA, ncol = 9, nrow = 49))
logll = rep(0, 9)
for( i in 1:9)
{
estimatedPars = read.table(paste("../estimated_params/GenSA_", sdm, "_", propData, i, "_", ordre, "_",step, "y.txt", sep=""))
parnames = estimatedPars[,1]
veget_pars[,i] = as.vector(estimatedPars[,2])

load(paste("../estimated_params/GenSA_", sdm, "_", propData, i, "_", ordre, "_",step, "y.RData", sep=""))
logll[i]= estim.pars$value
}
rownames(veget_pars) = parnames
##-------

apply(veget_pars, 1, mean)
apply(veget_pars, 1, sd)
logll
which.min(logll)

#--
load("scale_info.Robj")
## Temp -4 à 10
## PP 700 à 1300
Trange = c(-4, 10)
Tticks = max(Trange)-min(Trange)+1
PPrange = c(700, 1300)
Tbounds = (Trange - vars.means["annual_mean_temp"])/vars.sd["annual_mean_temp"]
PPbounds = (PPrange - vars.means["tot_annual_pp"])/vars.sd["tot_annual_pp"]
PPticks = (max(PPrange)-min(PPrange))/100 +1

tpseq=seq(Tbounds[1],Tbounds[2],l=100)
ppseq=seq(PPbounds[1],PPbounds[2],l=100)

scaled.axis<- function(){
temp = tpseq*vars.sd["annual_mean_temp"]+vars.means["annual_mean_temp"]
precip = ppseq*vars.sd["tot_annual_pp"]+vars.means["tot_annual_pp"]
axis(1, at = seq(Tbounds[1], Tbounds[2], l=Tticks), labels = seq(Trange[1], Trange[2], l = Tticks))
axis(2, at = seq(PPbounds[1],PPbounds[2], l=PPticks), labels = seq(PPrange[1], PPrange[2], l = PPticks))
}

#---

#------
# look at transition probabilities
#------

ENV = expand.grid(TP =tpseq , PP = ppseq)
ENV1 = ENV$TP
ENV2 = ENV$PP


calcMacroPars = function(params, ENV)
{
names(params) = rownames(veget_pars)
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

return(macroPars)
}

macroPars_tab= apply(veget_pars, 2, calcMacroPars, ENV=ENV)
macroPars_tab = array(unlist(macroPars_tab), dim = c(nrow(macroPars_tab[[1]]), ncol(macroPars_tab[[1]]), length(macroPars_tab)))

# w mean
macroPars = data.frame(apply(macroPars_tab, c(1,2), function(x){sum(x/logll)/(sum(1/logll))}))

dim(macroPars)
head(macroPars)

colnames(macroPars) = c("alphab", "alphat", "betab", "betat", "theta", "thetat", "eps")
#---
#
pal = colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")

#jpeg(paste("../figures/ave_", ordre, "_", step, "y_params.jpeg", sep=""), height=3000, width=5000, res=600)

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

nB = proj2[,"B"]
nT = proj2[,"T"]
nM = proj2[,"M"]

pTransitions = data.frame(pRT = macroPars$alphat*(nT+nM)*(1-macroPars$alphab*(nB+nM)),
pRB = macroPars$alphab*(nB+nM)*(1-macroPars$alphat*(nT+nM)),
pRM = macroPars$alphat*(nT+nM)*macroPars$alphab*(nB+nM),
pMT = macroPars$theta*macroPars$thetat,
pMB = macroPars$theta*(1-macroPars$thetat),
pTM = macroPars$betab,
pBM = macroPars$betat,
eps = macroPars$eps)

summary(pTransitions)

#jpeg(paste("../figures/ave_", ordre, "_", step, "y_transitions.jpeg", sep=""), height=3000, width=5000, res=600)

par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(pTransitions))
{
image(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(pTransitions)[i], xaxt = "n", yaxt="n")
contour(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
scaled.axis()
}

#dev.off()

### 
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

library(rootSolve)


eq.winner = function(pars, model)
{
fparams = pars
init = c(T = 0.1, B = 0.1, M=0.8)
(eq = stode(y = init, func = model, parms = fparams, positive = TRUE)$y)
return(names(which.max(eq)))
}

eq = t(apply(pars, 1, eq.winner, model =model))
eq = as.factor(eq)
##-----
#jpeg(paste("../figures/ave_", ordre, "_", step, "y_equi.jpeg", sep=""), height=3000, width=5000, res=600)
#
colo = c(M = "lightgreen", B = rgb(44,133,113,maxColorValue=255), T = rgb(245,172,71,maxColorValue=255), R = rgb(218,78,48,maxColorValue=255))

layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("",cex=2)
legend("center",legend = levels(eq),fill = colo[levels(eq)],bty = "n", cex = 0.8, ncol =2)
par(mar=c(5,5,0,2))

image(x=tpseq, y=ppseq, z = matrix(as.numeric(eq), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Annual mean temperature (°C)", ylab = "Annual precipitations (mm)", col = colo[levels(eq)], main = "", xaxt = "n", yaxt="n")
scaled.axis()
#
#dev.off()

#------
## reactivite
#------
#
library(rootSolve)

reactivity = function(pars, model)
{
fparams = pars
init = c(T = 0.1, B = 0.1, M=0.8)
(eq = stode(y = init, func = model, parms = fparams, positive = TRUE)$y)
jacob = jacobian.full(y = eq, func= model, parms=fparams)
M = (jacob + t(jacob))/2
return(max(eigen(M)$values))
}

reac = apply(pars, 1, reactivity, model =model)
reac2 = reac
reac2[reac2<0] = 0

#jpeg(paste("../figures/ave_", ordre, "_", step, "y_reactivity.jpeg", sep=""), height=3000, width=5000, res=600)

par(mfrow = c(1,1), mar = c(4,4,4,4))
image(x=tpseq, y=ppseq, z = matrix(reac2, ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = c(rgb(0,0,0,.5),pal(12)), main = "reactivity (amplification)", xaxt = "n", yaxt="n")
contour(x=tpseq, y=ppseq, z = matrix(reac2, ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
scaled.axis()

#dev.off()

### 
#------
# invasibility
#------

library(rootSolve)

eq.eigen = function(pars, model, eqState = "B")
{
fparams = pars
# Solve the model for the case where T/B and M are at 0
# and B/T = 1 - eps/alphaB/T (ie a l'équilibre si T/B=0 et M=0)
if(eqState=="B") init = c(T = 0, B = as.numeric(1-fparams["eps"]/fparams["alphaB"]), M=0)
if(eqState=="T") init = c(T = as.numeric(1-fparams["eps"]/fparams["alphaT"]), B=0, M=0)
if(is.infinite(sum(init))) init = c(T=0, B = 0, M=0)
#(eq = stode(y = init, func = model, parms = fparams, positive = TRUE)$y)
(jacob = jacobian.full(y = init, func= model, parms=fparams))
#print(jacob)
return(eigen(jacob)$values)
}

 
 invT = t(apply(pars, 1, eq.eigen, model =model, eqState="B"))
 invT.class = apply(invT, 1, function(x){ifelse(max(x[-1])>0, "instable", "stable")})
 T.growth = apply(pars, 1, function(pars){1-pars["eps"]/pars["alphaT"]})
 invT.class[T.growth<0] = "stable"
 invT.class = as.factor(invT.class)
 table(invT.class)

 invB = t(apply(pars, 1, eq.eigen, model =model, eqState="T"))
 invB.class = apply(invB, 1, function(x){ifelse(max(x[-2])>0, "instable", "stable")})
 B.growth = apply(pars, 1, function(pars){1-pars["eps"]/pars["alphaB"]})
 invB.class[B.growth<0] = "stable"
 invB.class = as.factor(invB.class)
 table(invB.class)
 
##--
coexist = numeric(length(invT))
coexist[invT.class=="stable" & invB.class=="stable"] = "AltSS"
coexist[invT.class=="stable" & invB.class=="instable"] = "B dominates"
coexist[invT.class=="instable" & invB.class=="instable"] = "Mixture"
coexist[invT.class=="instable" & invB.class=="stable"] = "T dominates"
table(coexist)
coexist = as.factor(coexist)
levels(coexist)
##---
colo = c(M = "lightgreen", B = rgb(44,133,113,maxColorValue=255), T = rgb(245,172,71,maxColorValue=255), R = rgb(218,78,48,maxColorValue=255))
colo["AltSS"] = "pink"
colo['Mixture'] = colo["M"]
colo['B dominates'] = colo["B"]
colo['T dominates'] = colo["T"]
colo['??'] = 1

#jpeg(paste("../figures/ave_", ordre, "_", step, "y_invasibility.jpeg", sep=""), height=3000, width=5000, res=600)

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
