rm(list = ls())
detach(params$best_par)

#load("../estimated_params/veget_m3_v2")
#source("4-init_params_v2.R")
#params = coarse_veget


#-- veget m3 -- esa ML=-3103.965
load("../estimated_params/coarse_veget_m3_step3")
source("4-init_params.R")
load("../estimated_params/coarse_veget_m3_step3_lim")
params = coarse_veget

##-- herbivores m3 - esa ML=-3096.31
load("../estimated_params/coarse_m3_step4")
source("4-init_params.R")
load("../estimated_params/coarse_m3_step4_lim")
#veget_pars = read.table("../estimated_params/par_herbivores_m3.txt")
params = coarse


##----------------

cbind(unlist(par_lo), unlist(par_hi), unlist(params$best_par))

#names(par_lo) = colnames(params$likeli)[-(1:3)]
#names(par_hi) = colnames(params$likeli)[-(1:3)]
head(params$likeli)
par(mfrow = c(7,7), mar = c(2, 2, 2,1))
for( i in colnames(params$likeli)[-(1:3)])
{
plot(params$likeli[,i], ylim = as.vector(unlist(c(par_lo[i], par_hi[i]))), main = i, pch=15)
points(y=params$best_par[i], x=20, col=2, pch=19)
}


attach(params$best_par)
tpseq=seq(min(data$ENV1),max(data$ENV1),length.out=100)
ppseq=seq(min(data$ENV2),max(data$ENV2),length.out=100)

#tpseq=seq(0,6,length.out=100)
#ppseq=seq(0.7,1.5,length.out=100)

ENV = expand.grid(TP =tpseq , PP = ppseq)
TP = ENV$TP
PP = ENV$PP
ENV1 = TP
ENV2 = PP



#v2
#  logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV2 + ab3*ENV1^2 + ab4*ENV2^2 +  ab5*ENV2*ENV1
#    logit_alphat 	= at0 + at1*ENV1 + at2*ENV2 + at3*ENV1^2 + at4*ENV2^2 + at5*ENV2*ENV1
#    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV2 + bb3*ENV1^2 + bb4*ENV2^2 + bb5*ENV2*ENV1
#    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV2 + bt3*ENV1^2 + bt4*ENV2^2 + bt5*ENV2*ENV1
#    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV2 + tb3*ENV1^2 + tb4*ENV2^2 + tb5*ENV2*ENV1
#    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV2 + tt3*ENV1^2 + tt4*ENV2^2 + tt5*ENV2*ENV1
#    logit_eps 	= e0  + e1*ENV1  + e2*ENV2 + e5*ENV2*ENV1
#

#v1
    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2 + tb5*ENV1^3 + tb6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3

alphab = exp(logit_alphab)/(1+exp(logit_alphab))
alphat = exp(logit_alphat)/(1+exp(logit_alphat))
betab = exp(logit_betab)/(1+exp(logit_betab))
betat = exp(logit_betat)/(1+exp(logit_betat))
thetab = exp(logit_thetab)/(1+exp(logit_thetab))
thetat = exp(logit_thetat)/(1+exp(logit_thetat))
eps = exp(logit_eps)/(1+exp(logit_eps))


par(mfrow = c(2,4), mar = c(4,4,1,1), cex=0.8)
plot(lowess(cbind(ENV1, alphab)), type ="l", xlab = "T", ylab = "alpha", ylim = c(0,1), col = "darkgreen")
legend("top", col = c("darkgreen","lightgreen"), legend=c("B", "T"), bty="n", lwd=1)
lines(lowess(cbind(ENV1, alphat)), col = "lightgreen")
plot(lowess(cbind(ENV1, betab)), type ="l", xlab = "T", ylab = "beta", ylim = c(0,1), col = "darkgreen")
lines(lowess(cbind(ENV1, betat)),col = "lightgreen")
plot(lowess(cbind(ENV1, thetab)), type ="l", xlab = "T", ylab = "theta", ylim = c(0,1), col = "darkgreen")
lines(lowess(cbind(ENV1, thetat)), col = "lightgreen")
plot(lowess(cbind(ENV1, eps)), type ="l", xlab = "T", ylab = "eps", ylim = c(0,1))
plot(lowess(cbind(ENV2, alphab)), type ="l", xlab = "PP", ylab = "alpha", ylim = c(0,1), col = "darkgreen")
lines(lowess(cbind(ENV2, alphat)), col = "lightgreen")
plot(lowess(cbind(ENV2, betab)), type ="l", xlab = "PP", ylab = "beta", ylim = c(0,1), col = "darkgreen")
lines(lowess(cbind(ENV2, betat)),col = "lightgreen")
plot(lowess(cbind(ENV2, thetab)), type ="l", xlab = "PP", ylab = "theta", ylim = c(0,1), col = "darkgreen")
lines(lowess(cbind(ENV2, thetat)), col = "lightgreen")
plot(lowess(cbind(ENV2, eps)), type ="l", xlab = "PP", ylab = "eps", ylim = c(0,1))



#dev.copy2pdf(file = "../figures/params_herbivores.pdf")
#dev.copy2pdf(file = "../figures/params_sansHerbivores.pdf")






### 

# correspondance
#aC = alphab
#aD = alphat
#bC = betab
#bD = betat
#sC = thetab
#sD = thetat
#e = eps
#invC = 1/2*((aD - e)*bC - (aD - e)*bD + aC*e - 2*aD*e - aD*sC - aD*sD + sqrt(aC^2*e^2 + aD^2*sC^2 + aD^2*sD^2 + (aD^2 - 2*aD*e + e^2)*bC^2 + (aD^2 - 2*aD*e + e^2)*bD^2 + 2*(aC*aD*e - aC*e^2)*bC + 2*(aC*aD*e - aC*e^2 + (aD^2 - 2*aD*e + e^2)*bC)*bD + 2*(aC*aD*e + (aD^2 - aD*e)*bC + (aD^2 - aD*e)*bD)*sC - 2*(2*aC*aD^2*e - (2*e^2 + e)*aC*aD - aD^2*sC + (aD^2 - aD*e)*bC + (aD^2 - aD*e)*bD)*sD))/aD
#invD = -1/2*((aC - e)*bC - (aC - e)*bD + 2*aC*e - aD*e + aC*sC + aC*sD - sqrt(aD^2*e^2 + aC^2*sC^2 + aC^2*sD^2 + 2*(aC*e - e^2)*aD*bC + (aC^2 - 2*aC*e + e^2)*bC^2 + (aC^2 - 2*aC*e + e^2)*bD^2 + 2*((aC*e - e^2)*aD + (aC^2 - 2*aC*e + e^2)*bC)*bD - 2*((2*aC^2*e - (2*e^2 + e)*aC)*aD + (aC^2 - aC*e)*bC + (aC^2 - aC*e)*bD)*sC + 2*(aC*aD*e + aC^2*sC + (aC^2 - aC*e)*bC + (aC^2 - aC*e)*bD)*sD))/aC
##----


ab = alphab
at = alphat
bb = betab
bt = betat
tb = thetab
tt = thetat
e = eps


# Compute the first? eigenvalues for C and D as invaders

invT1 = -(ab^2*bb - ab^2*bt + 2*ab^2*e + ab^2*tb + ab^2*tt - ab*at*e - ab*bb*e + ab*bt*e + sqrt(-4*ab^4*at*e*tb + ab^4*bb^2 + 2*ab^4*bb*bt - 2*ab^4*bb*tb + 2*ab^4*bb*tt + ab^4*bt^2 - 2*ab^4*bt*tb + 2*ab^4*bt*tt + ab^4*tb^2 + 2*ab^4*tb*tt + ab^4*tt^2 + 2*ab^3*at*bb*e + 2*ab^3*at*bt*e + 4*ab^3*at*e^2*tb + 2*ab^3*at*e*tb + 2*ab^3*at*e*tt - 2*ab^3*bb^2*e - 4*ab^3*bb*bt*e + 2*ab^3*bb*e*tb - 2*ab^3*bb*e*tt - 2*ab^3*bt^2*e + 2*ab^3*bt*e*tb - 2*ab^3*bt*e*tt + ab^2*at^2*e^2 - 2*ab^2*at*bb*e^2 - 2*ab^2*at*bt*e^2 + ab^2*bb^2*e^2 + 2*ab^2*bb*bt*e^2 + ab^2*bt^2*e^2))/(2*ab^2)

invT2 = -(ab^2*bb - ab^2*bt + 2*ab^2*e + ab^2*tb + ab^2*tt - ab*at*e - ab*bb*e + ab*bt*e - sqrt(-4*ab^4*at*e*tb + ab^4*bb^2 + 2*ab^4*bb*bt - 2*ab^4*bb*tb + 2*ab^4*bb*tt + ab^4*bt^2 - 2*ab^4*bt*tb + 2*ab^4*bt*tt + ab^4*tb^2 + 2*ab^4*tb*tt + ab^4*tt^2 + 2*ab^3*at*bb*e + 2*ab^3*at*bt*e + 4*ab^3*at*e^2*tb + 2*ab^3*at*e*tb + 2*ab^3*at*e*tt - 2*ab^3*bb^2*e - 4*ab^3*bb*bt*e + 2*ab^3*bb*e*tb - 2*ab^3*bb*e*tt - 2*ab^3*bt^2*e + 2*ab^3*bt*e*tb - 2*ab^3*bt*e*tt + ab^2*at^2*e^2 - 2*ab^2*at*bb*e^2 - 2*ab^2*at*bt*e^2 + ab^2*bb^2*e^2 + 2*ab^2*bb*bt*e^2 + ab^2*bt^2*e^2))/(2*ab^2)

invT3 =  -ab + e

invB1 = (ab*at*e + at^2*bb - at^2*bt - 2*at^2*e - at^2*tb - at^2*tt - at*bb*e + at*bt*e + sqrt(ab^2*at^2*e^2 - 4*ab*at^4*e*tt + 2*ab*at^3*bb*e + 2*ab*at^3*bt*e + 4*ab*at^3*e^2*tt + 2*ab*at^3*e*tb + 2*ab*at^3*e*tt - 2*ab*at^2*bb*e^2 - 2*ab*at^2*bt*e^2 + at^4*bb^2 + 2*at^4*bb*bt + 2*at^4*bb*tb - 2*at^4*bb*tt + at^4*bt^2 + 2*at^4*bt*tb - 2*at^4*bt*tt + at^4*tb^2 + 2*at^4*tb*tt + at^4*tt^2 - 2*at^3*bb^2*e - 4*at^3*bb*bt*e - 2*at^3*bb*e*tb + 2*at^3*bb*e*tt - 2*at^3*bt^2*e - 2*at^3*bt*e*tb + 2*at^3*bt*e*tt + at^2*bb^2*e^2 + 2*at^2*bb*bt*e^2 + at^2*bt^2*e^2))/(2*at^2)

invB2 = -at + e

invB3 = (ab*at*e + at^2*bb - at^2*bt - 2*at^2*e - at^2*tb - at^2*tt - at*bb*e + at*bt*e - sqrt(ab^2*at^2*e^2 - 4*ab*at^4*e*tt + 2*ab*at^3*bb*e + 2*ab*at^3*bt*e + 4*ab*at^3*e^2*tt + 2*ab*at^3*e*tb + 2*ab*at^3*e*tt - 2*ab*at^2*bb*e^2 - 2*ab*at^2*bt*e^2 + at^4*bb^2 + 2*at^4*bb*bt + 2*at^4*bb*tb - 2*at^4*bb*tt + at^4*bt^2 + 2*at^4*bt*tb - 2*at^4*bt*tt + at^4*tb^2 + 2*at^4*tb*tt + at^4*tt^2 - 2*at^3*bb^2*e - 4*at^3*bb*bt*e - 2*at^3*bb*e*tb + 2*at^3*bb*e*tt - 2*at^3*bt^2*e - 2*at^3*bt*e*tb + 2*at^3*bt*e*tt + at^2*bb^2*e^2 + 2*at^2*bb*bt*e^2 + at^2*bt^2*e^2))/(2*at^2)

#invT = invT2
invT = apply(cbind(invT1, invT2, invT3), 1, function(x){x[which.max((x))]})
#invB = invB1
invB = apply(cbind(invB1, invB2, invB3), 1, function(x){x[which.max((x))]})

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
coexist[invB>0 & invT<0 & (ab-e)>0] = 2
#coexist[invB>0 & invT<0] = 2

# Species T wins
coexist[invB<0 & invT>0 & (at-e)>0] = 3
#coexist[invB<0 & invT>0] = 3

# Reciprocal invasibility
coexist[invB > 0 & invT > 0 & (ab-e)>0 & (at-e)>0] = 4
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
image(tpseq,ppseq*1000,Z,xlab = "Mean annual temperature", ylab = "Annual precipitation (mm)", cex.lab = 1.5, cex.axis = 1.25, col = colo, breaks = c(0:5))#grey(c(0:3)/3))

#dev.copy2pdf(file = "../figures/Coexistence_area_herbivores.pdf")
#dev.copy2pdf(file = "../figures/Coexistence_area_sansHerbivores.pdf")




