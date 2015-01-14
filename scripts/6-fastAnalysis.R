rm(list = ls())
#veget_pars = read.table("../estimated_params/GenSA_multinom_0.1.txt")
veget_pars = read.table("../estimated_params/GenSA_rf_0.1.txt")
pars = as.numeric(veget_pars[,1])
names(pars) = rownames(veget_pars)
pars = as.list(pars)
load("selectForFit_rf_0.1")

load("scale_info.Robj")

##----------------
neiborgh = "rf"
source("4-init_params.R")
##-----------------

#-----
# check estimated params and bounds
#----
parsBounds = cbind(unlist(par_lo), unlist(par_hi), unlist(pars))

#par(mfrow = c(1,1), mar = c(2, 2, 2,1))

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

    logit_alphab 	= pars$ab0 + pars$ab1*ENV1 + pars$ab2*ENV2 + pars$ab3*ENV1^2 + pars$ab4*ENV2^2 + pars$ab5*ENV1^3 + pars$ab6*ENV2^3
    logit_alphat 	= pars$at0 + pars$at1*ENV1 + pars$at2*ENV2 + pars$at3*ENV1^2 + pars$at4*ENV2^2 + pars$at5*ENV1^3 + pars$at6*ENV2^3
    logit_betab 	= pars$bb0 + pars$bb1*ENV1 + pars$bb2*ENV2 + pars$bb3*ENV1^2 + pars$bb4*ENV2^2 + pars$bb5*ENV1^3 + pars$bb6*ENV2^3
    logit_betat 	= pars$bt0 + pars$bt1*ENV1 + pars$bt2*ENV2 + pars$bt3*ENV1^2 + pars$bt4*ENV2^2 + pars$bt5*ENV1^3 + pars$bt6*ENV2^3
    logit_theta	= pars$t0 + pars$t1*ENV1 + pars$t2*ENV2 + pars$t3*ENV1^2 + pars$t4*ENV2^2 + pars$t5*ENV1^3 + pars$t6*ENV2^3
    logit_thetat	= pars$tt0 + pars$tt1*ENV1 + pars$tt2*ENV2 + pars$tt3*ENV1^2 + pars$tt4*ENV2^2 + pars$tt5*ENV1^3 + pars$tt6*ENV2^3
    logit_eps 	= pars$e0  + pars$e1*ENV1 + pars$e2*ENV2  + pars$e3*ENV1^2 + pars$e4*ENV2^2 + pars$e5*ENV1^3 + pars$e6*ENV2^3
 
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

par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)

for (i in 1:ncol(macroPars))
{
image(x=tpseq, y=ppseq, z = matrix(macroPars[,i], ncol = length(ppseq), nrow = length(tpseq)),xlab = "Temperature", ylab = "Precipitations", col = pal(12), main = colnames(macroPars)[i])
contour(x=tpseq, y=ppseq, z = matrix(macroPars[,i], ncol = length(ppseq), nrow = length(tpseq)), add=TRUE)
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
invT = apply(cbind(invT1, invT2), 1, function(x){x[which.max((x))]})
#invB = invB1
invB = apply(cbind(invB1, invB2), 1, function(x){x[which.max((x))]})

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
coexist[invT<0 & invB<0] = 1

# Species B wins (instabilité au point B=0,T=kT + (ab-e)>0 et stabilité au point B=kB,T=0
# 
coexist[invB>0 & invT<0 & (alphaB-eps)>0] = 2
#coexist[invB>0 & invT<0] = 2

# Species T wins
coexist[invB<0 & invT>0 & (alphaT-eps)>0] = 3
#coexist[invB<0 & invT>0] = 3

# Reciprocal invasibility
coexist[invB > 0 & invT > 0 & (alphaB-eps)>0 & (alphaT-eps)>0] = 4
#coexist[invB > 0 & invT > 0] = 4


# instabilité vers crash
#coexist[] = 0
#

table(coexist)
# Plot the results
Z = matrix(coexist+1,nr = length(tpseq), nc = length(ppseq))
quartz(width = 6, height = 6)
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




