rm(list = ls())

sdm = "rf"
sdm = "multinom"
propData = 0.3

#--
veget_pars = data.frame(matrix(NA, ncol = 9, nrow = 49))
for( i in 1:9)
{
veget_pars[,i] = as.vector(read.table(paste("../estimated_params/GenSA_initForFit_", sdm, "_",propData, i,".txt", sep="")))
}
rownames(veget_pars) = rownames(as.vector(read.table(paste("../estimated_params/GenSA_initForFit_", sdm, "_",propData, i,".txt", sep=""))))
veget_pars

##-------

apply(veget_pars, 1, mean)
apply(veget_pars, 1, sd)
veget_pars[c("ab2","at0"),]

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



#--
load(paste("../estimated_params/GenSA_initForFit_", sdm, "_", propData, "1.RData", sep = ""))
estim.pars$value
#-----
# diagnostic convergence
#----
#jpeg(paste("../figures/diagnostic",fit,".jpeg", sep=""), height=5000, width=5000, res=600)
par(mfrow = c(2,2), mar = c(4, 4, 2,1))
mat = data.frame(estim.pars$trace.mat)
mat$number = 1:nrow(mat)
plot(nb.steps~number, data = mat, xlab = "anneal time", ylab = "nb.steps", cex=.2, type = "l")
plot(temperature~number, data = mat, cex=.2, type = "l")
plot(function.value~number, data = mat, cex=.2, type = "l")
plot(current.minimum~number, data = mat, cex = .2, type="l")
#dev.off()
#-----
# check estimated params and bounds
#----
#load(paste("initForFit_",sdm, "_", propData, "1.RData", sep = ""))

#------
# look at transition probabilities
#------

ENV = expand.grid(TP =tpseq , PP = ppseq)
ENV1 = ENV$TP
ENV2 = ENV$PP

e7 = 0
params=veget_pars[,1]
names(params) = rownames(veget_pars)

    logit_alphab 	= params["ab0"] + params["ab1"]*ENV1 + params["ab2"]*ENV2 + params["ab3"]*ENV1^2 + params["ab4"]*ENV2^2 + params["ab5"]*ENV1^3 + params["ab6"]*ENV2^3
    logit_alphat 	= params["at0"] + params["at1"]*ENV1 + params["at2"]*ENV2 + params["at3"]*ENV1^2 + params["at4"]*ENV2^2 + params["at5"]*ENV1^3 + params["at6"]*ENV2^3
    logit_betab 	= params["bb0"] + params["bb1"]*ENV1 + params["bb2"]*ENV2 + params["bb3"]*ENV1^2 + params["bb4"]*ENV2^2 + params["bb5"]*ENV1^3 + params["bb6"]*ENV2^3
    logit_betat 	= params["bt0"] + params["bt1"]*ENV1 + params["bt2"]*ENV2 + params["bt3"]*ENV1^2 + params["bt4"]*ENV2^2 + params["bt5"]*ENV1^3 + params["bt6"]*ENV2^3
    logit_theta	= params["t0"] + params["t1"]*ENV1 + params["t2"]*ENV2 + params["t3"]*ENV1^2 + params["t4"]*ENV2^2 + params["t5"]*ENV1^3 + params["t6"]*ENV2^3
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
pal = colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")

#jpeg(paste("../figures/estim_pars_", sdm, "_03x9.jpeg", sep=""), height=3000, width=5000, res=600)

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

pTransitions = data.frame(pRT = macroPars$alphat*0.5*(1-macroPars$alphab*0.5),
pRB = macroPars$alphab*0.5*(1-macroPars$alphat*0.5),
pRM = macroPars$alphat*0.5*macroPars$alphab*0.5,
pMT = macroPars$theta*macroPars$thetat*(1-macroPars$eps),
pMB = macroPars$theta*(1-macroPars$thetat)*(1-macroPars$eps),
pTM = macroPars$betab*0.5*(1-macroPars$eps),
pBM = macroPars$betat*0.5*(1-macroPars$eps),
eps = macroPars$eps)

#jpeg(paste("../figures/estim_transitions_", sdm, "_03x9.jpeg", sep=""), height=3000, width=5000, res=600)

par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(pTransitions))
{
image(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(pTransitions)[i], xaxt = "n", yaxat="n")
contour(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
scaled.axis()
}
#dev.off()

### 
#------
# invasibility
#------

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
return(eigen(jacob)$values)
}

invT = t(apply(pars, 1, eq.eigen, model =model, eqState="B"))
invT.class = ifelse(invT[,3]<0, "stable", "invadeM")
invT.class = ifelse(invT[,1]>0, "invadeT", invT.class)
invT.class = as.factor(invT.class)
T.growth = apply(pars, 1, function(pars){1-pars["eps"]/pars["alphaT"]})
invT.class[T.growth<0] = "stable"

invB = t(apply(pars, 1, eq.eigen, model =model, eqState="T"))
invB.class = ifelse(invB[,3]<0, "stable", "invadeM")
invB.class = ifelse(invB[,2]>0, "invadeB", invB.class)
invB.class = as.factor(invB.class)
B.growth = apply(pars, 1, function(pars){1-pars["eps"]/pars["alphaB"]})
invB.class[B.growth<0] = "stable"

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
par(mfrow = c(1,2), mar = c(4,4,1,1), cex=0.8)

image(x=tpseq, y=ppseq, z = matrix(as.numeric(invT.class), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = c("darkgreen","lightgreen", "grey"), main = "invasibility of T (when B every where)", xaxt = "n", yaxt="n")
scaled.axis()
legend("topleft",legend = levels(invT.class),fill = c("darkgreen","lightgreen", "grey"),bty = "n", cex = 0.8, ncol =1)

image(x=tpseq, y=ppseq, z = matrix(as.numeric(invB.class), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = c("darkgreen","lightgreen", "grey"), main = "invasibility of B (when T everywhere)", xaxt = "n", yaxt="n")
scaled.axis()
legend("topleft",legend = levels(invB.class),fill = c("darkgreen","lightgreen", "grey"),bty = "n", cex = 0.8, ncol =1)


##-----
jpeg(paste("../figures/equilibrium_map_", sdm, "_03x9.jpeg", sep=""), height=3000, width=3000, res = 300)
colo = c("pink","darkorange", "darkgreen", "orange3", "orange", "lightgreen")
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))

par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
title("",cex=2)
legend("center",legend = levels(coexist),fill = colo,bty = "n", cex = 0.8, ncol =2)
par(mar=c(5,5,0,2))

image(x=tpseq, y=ppseq, z = matrix(as.numeric(coexist), ncol = length(ppseq), nrow = length(tpseq)),xlab = "Annual mean temperature (°C)", ylab = "Annual precipitations (mm)", col = colo, main = "invasibility of B (when T everywhere)", xaxt = "n", yaxt="n")
scaled.axis()
dev.off()





#reactivite -> potentiel à amplifier les fluctuations env.
#M = (jacob + t(jacob))/2
#max(abs(eigen(M)$values))
#M.veg = (jacob.veg + t(jacob.veg))/2
#max(abs(eigen(M.veg)$values))
#



