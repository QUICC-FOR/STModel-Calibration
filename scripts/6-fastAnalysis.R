rm(list = ls())

veget_pars = read.table("../estimated_params/GenSA_initForFit_rf_0.05.txt")
load("initForFit_rf_0.05")
load("../estimated_params/GenSA_initForFit_rf_0.05.RData")
#--

load("scale_info.Robj")

#---
pars = as.numeric(veget_pars[,1])
names(pars) = rownames(veget_pars)
#pars = as.list(pars)

#-----
# diagnostic convergence
#----
jpeg(paste("../figures/diagnostic",fit,".jpeg", sep=""), height=5000, width=5000, res=600)
par(mfrow = c(2,2), mar = c(4, 4, 2,1))
mat = estim.pars$trace.mat
plot(1:nrow(mat),mat[,"nb.steps"], xlab = "anneal time", ylab = "nb.steps", cex=.2)
plot(temperature~nb.steps, data = mat, cex=.2)
plot(function.value~nb.steps, data = mat, cex=.2)
plot(current.minimum~nb.steps, data = mat, cex = .2)
dev.off()
#-----
# check estimated params and bounds
#----
parsBounds = cbind(unlist(par_lo), unlist(par_hi), unlist(pars))

par(mfrow = c(1,1), mar = c(2, 2, 2,1))

plot(parsBounds[,1], ylim = c(min(parsBounds), max(parsBounds)), pch=15, cex=.5)
points(parsBounds[,2], col=1, pch=15, cex=.5)
points(parsBounds[,3], col=2, pch=19, cex=.5)
title("estimated params and given bounds")

#------
# look at transition probabilities
#------

tpseq=seq(-4,2,l=100)
ppseq=seq(-2,5,l=100)
ENV = expand.grid(TP =tpseq , PP = ppseq)
ENV1 = ENV$TP
ENV2 = ENV$PP

e7 = 0
params=pars


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
#colo = c(R = rgb(.5,.5,.5,.5), T = rgb(1,0,0,.5), B = rgb(0.2,.8,.2,.5), M = rgb(0,0,1,.5))
pal = colorRampPalette(c("lightblue", "yellow", "orange"), space = "rgb")

#image(x=tpseq, y=ppseq, z = matrix(macroPars$alphab, ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = "alphab")
#contour(x=tpseq, y=ppseq, z = matrix(macroPars$alphab, ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
jpeg(paste("../figures/estim_pars_",fit,".jpeg", sep=""), height=3000, width=5000, res=600)

par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(macroPars))
{
image(x=tpseq, y=ppseq, z = matrix(macroPars[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(macroPars)[i])
contour(x=tpseq, y=ppseq, z = matrix(macroPars[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
}
dev.off()

### 
#------
# probabilites transition (no neighborhood restriction)
#------
nB = mean(datSel$EB)
nT = mean(datSel$ET)
nM = mean(datSel$EM)

pTransitions = data.frame(pRT = macroPars$alphat*(nT+nM)*(1-macroPars$alphab*(nB+nM)),
pRB = macroPars$alphab*(nB+nM)*(1-macroPars$alphat*(nT+nM)),
pRM = macroPars$alphat*(nT+nM)*macroPars$alphab*(nB+nM),
pMT = macroPars$theta*macroPars$thetat,
pMB = macroPars$theta*(1-macroPars$thetat),
pTM = macroPars$betab,
pBM = macroPars$betat,
eps = macroPars$eps)

par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(pTransitions))
{
image(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(pTransitions)[i])
contour(x=tpseq, y=ppseq, z = matrix(pTransitions[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
}


### 
#------
# invasibility: TODO
#------

alphaB = macroPars$alphab
alphaT = macroPars$alphat
betaB = macroPars$betab
betaT= macroPars$betat
theta = macroPars$theta
thetaT = macroPars$thetat
eps = macroPars$eps

# Compute the first? eigenvalues for C and D as invaders
source("invasibility_values.r")

#invT = invT2
invT = apply(cbind(invT1, invT2), 1, function(x){max(x, na.rm=TRUE)})
#invB = invB1
invB = apply(cbind(invB1, invB2), 1, function(x){max(x, na.rm=TRUE)})

# Interpret the invasability criterion
coexist = numeric(length(invT))
# unisp
#coexist[invT>0] = 3
#coexist[invT<0] = 4
#
#coexist[invB>0] = 2
#coexist[invB<0] = 4

###both
# Reciprocal resistance (alternative stable states)
coexist[invT<0 & (alphaB-eps)>0 & invB<0 & (alphaT-eps)>0] = 1

# Species B wins (instabilité au point B=0,T=kT + (ab-e)>0 et stabilité au point B=kB,T=0
# 
coexist[invB>0 & invT<0 & (alphaB-eps)>0] = 2
#coexist[invB>0 & invT<0] = 2

# Species T wins
coexist[invB<0 & invT>0 & (alphaT-eps)>0] = 3
#coexist[invB<0 & invT>0] = 3

# Reciprocal invasibility
#coexist[invB > 0 & invT > 0 & (alphaB-eps)>0 & (alphaT-eps)>0] = 4
coexist[invB > 0 & invT > 0] = 4

## deal with cases where alphaB or alphaT = 0
coexist[alphaB==0&alphaT>0] = 3
coexist[alphaT==0&alphaB>0] = 2



# instabilité vers crash
#coexist[] = 0
#

table(coexist)
# Plot the results
Z = matrix(coexist+1,nr = length(tpseq), nc = length(ppseq))
#quartz(width = 6, height = 6)
colo = c("white","pink", "darkgreen", "lightgreen", "orange")
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#title(title,cex=2)
legend("center",legend = c("other","AltSS","Boreal Wins","Temperate Wins","Coexistence"),fill = colo,bty = "n",horiz = TRUE,cex = 0.8)
par(mar=c(5,5,0,2))
image(tpseq,ppseq,Z,xlab = "Mean annual temperature", ylab = "Annual precipitation (mm)", cex.lab = 1.5, cex.axis = 1.25, col = colo, breaks = c(0:5))#grey(c(0:3)/3))

#dev.copy2pdf(file = "../figures/Coexistence_area_herbivores.pdf")
#dev.copy2pdf(file = "../figures/Coexistence_area_sansHerbivores.pdf")




