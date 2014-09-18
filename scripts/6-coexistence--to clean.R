rm(list = ls())
#veget_pars = read.table("../estimated_params/par_herbivores_m3.txt")
veget_pars = read.table("../estimated_params/par_m3.txt")
#veget_pars = read.table("../estimated_params/par_v2.txt")

attach(veget_pars)
tpseq=seq(-2,2,0.01)
ppseq=seq(-2,2,0.01)
ENV = expand.grid(TP =tpseq , PP = ppseq)
TP = scale(ENV$TP)
PP = scale(ENV$PP)
ENV1 = TP
ENV2 = PP

    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV2 + ab3*ENV1^2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV2 + at3*ENV1^2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV2 + bb3*ENV1^2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV2 + bt3*ENV1^2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV2 + tb3*ENV1^2 + tb4*ENV2^2 + tb5*ENV1^3 + tb6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV2 + tt3*ENV1^2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1 + e2*ENV2  + e3*ENV1^2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3
#
#    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2 
#    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2 
#    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2 
#    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2 
#    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2 
#    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2 
#    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2
    
#logit_alphab 	= ac0 + act1*ENV1 + act2*ENV1^2 + acp1*ENV2 + acp2*ENV2^2
#logit_alphat 	= ad0 + adt1*ENV1 + adt2*ENV1^2 + adp1*ENV2 + adp2*ENV2^2
#logit_betab 	= bc0 + bct1*ENV1 + bct2*ENV1^2 + bcp1*ENV2 + bcp2*ENV2^2
#logit_betat 	= bd0 + bdt1*ENV1 + bdt2*ENV1^2 + bdp1*ENV2 + bdp2*ENV2^2
#logit_thetab	= tc0 + tct1*ENV1 + tct2*ENV1^2 + tcp1*ENV2 + tcp2*ENV2^2
#logit_thetat	= td0 + tdt1*ENV1 + tdt2*ENV1^2 + tdp1*ENV2 + tdp2*ENV2^2
#logit_eps 	= e0  + et1*ENV1  + et2*ENV1^2 + ep1*ENV2 + ep2*ENV2^2
 
alphab = exp(logit_alphab)/(1+exp(logit_alphab))
alphat = exp(logit_alphat)/(1+exp(logit_alphat))
betab = exp(logit_betab)/(1+exp(logit_betab))
betat = exp(logit_betat)/(1+exp(logit_betat))
thetab = exp(logit_thetab)/(1+exp(logit_thetab))
thetat = exp(logit_thetat)/(1+exp(logit_thetat))
eps = exp(logit_eps)/(1+exp(logit_eps))


# correspondance
aC = alphab
aD = alphat
bC = betab
bD = betat
sC = thetab
sD = thetat
e = eps


# Compute the first eigenvalues for C and D as invaders
invC = 1/2*((aD - e)*bC - (aD - e)*bD + aC*e - 2*aD*e - aD*sC - aD*sD + sqrt(aC^2*e^2 + aD^2*sC^2 + aD^2*sD^2 + (aD^2 - 2*aD*e + e^2)*bC^2 + (aD^2 - 2*aD*e + e^2)*bD^2 + 2*(aC*aD*e - aC*e^2)*bC + 2*(aC*aD*e - aC*e^2 + (aD^2 - 2*aD*e + e^2)*bC)*bD + 2*(aC*aD*e + (aD^2 - aD*e)*bC + (aD^2 - aD*e)*bD)*sC - 2*(2*aC*aD^2*e - (2*e^2 + e)*aC*aD - aD^2*sC + (aD^2 - aD*e)*bC + (aD^2 - aD*e)*bD)*sD))/aD

invD = -1/2*((aC - e)*bC - (aC - e)*bD + 2*aC*e - aD*e + aC*sC + aC*sD - sqrt(aD^2*e^2 + aC^2*sC^2 + aC^2*sD^2 + 2*(aC*e - e^2)*aD*bC + (aC^2 - 2*aC*e + e^2)*bC^2 + (aC^2 - 2*aC*e + e^2)*bD^2 + 2*((aC*e - e^2)*aD + (aC^2 - 2*aC*e + e^2)*bC)*bD - 2*((2*aC^2*e - (2*e^2 + e)*aC)*aD + (aC^2 - aC*e)*bC + (aC^2 - aC*e)*bD)*sC + 2*(aC*aD*e + aC^2*sC + (aC^2 - aC*e)*bC + (aC^2 - aC*e)*bD)*sD))/aC

# Interpret the invasability criterion
coexist = numeric(length(invC))

# Reciprocal resistance (alternative stable states)
coexist[invC<0 & invD<0] = 1

# Species C wins
coexist[invC>0 & invD<0] = 2

# Species D winds
coexist[invC<0 & invD>0] = 3

# Reciprocal invasibility
coexist[invC > 0 & invD > 0] = 4

table(coexist)
# Plot the results
Z = matrix(coexist,nr = length(tpseq), nc = length(ppseq))
quartz(width = 6, height = 6)
colo = c("pink", "darkgreen", "lightgreen", "orange")
layout(matrix(c(1,2),nr=2,nc=1,byrow=TRUE),heights = c(1,6))
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#title(title,cex=2)
legend("center",legend = c("AltSS","Boreal Wins","Temperate Wins","Coexistence"),fill = colo,bty = "n",horiz = TRUE,cex = 0.8)
par(mar=c(5,5,0,2))
image(tpseq,ppseq0,Z,xlab = "Mean annual temperature", ylab = "Annual precipitation", cex.lab = 1.5, cex.axis = 1.25, col = colo, breaks = c(0:4))#grey(c(0:3)/3))

#dev.copy2pdf(file = "Coexistence_area.pdf")




