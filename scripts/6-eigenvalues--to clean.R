rm(list = ls())
veget_pars = read.table("../estimated_params/par_herbivores_rf3_clim2.txt")
#veget_pars = read.table("../estimated_params/par_veget_m3.txt")
#veget_pars = read.table("../estimated_params/par_v2.txt")

attach(veget_pars)

ENV2 = 1.2
ENV1 = seq(-1,5,0.01) 

    logit_alphab 	= ab0 + ab1*ENV1 + ab2*ENV1^2 + ab3*ENV2 + ab4*ENV2^2 + ab5*ENV1^3 + ab6*ENV2^3
    logit_alphat 	= at0 + at1*ENV1 + at2*ENV1^2 + at3*ENV2 + at4*ENV2^2 + at5*ENV1^3 + at6*ENV2^3
    logit_betab 	= bb0 + bb1*ENV1 + bb2*ENV1^2 + bb3*ENV2 + bb4*ENV2^2 + bb5*ENV1^3 + bb6*ENV2^3
    logit_betat 	= bt0 + bt1*ENV1 + bt2*ENV1^2 + bt3*ENV2 + bt4*ENV2^2 + bt5*ENV1^3 + bt6*ENV2^3
    logit_thetab	= tb0 + tb1*ENV1 + tb2*ENV1^2 + tb3*ENV2 + tb4*ENV2^2 + tb5*ENV1^3 + tb6*ENV2^3
    logit_thetat	= tt0 + tt1*ENV1 + tt2*ENV1^2 + tt3*ENV2 + tt4*ENV2^2 + tt5*ENV1^3 + tt6*ENV2^3
    logit_eps 	= e0  + e1*ENV1  + e2*ENV1^2 + e3*ENV2 + e4*ENV2^2 + e5*ENV1^3 + e6*ENV2^3

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


## diagnostic plots

par(mfcol = c(2,3))

plot(alphab~ENV1, type ="l", ylim = c(0,1))
plot(alphat~ENV1, type ="l", ylim = c(0,1))
plot(betab~ENV1, type ="l", ylim = c(0,1))
plot(betat~ENV1, type ="l", ylim = c(0,1))
plot(thetab~ENV1, type ="l", ylim = c(0,1))
plot(thetat~ENV1, type ="l", ylim = c(0,1))



# correspondance
aC = alphab
aD = alphat
bC = betab
bD = betat
sC = thetab
sD = thetat
e = eps
TP = ENV1
PP = ENV2

# Compute the first eigenvalues for C and D as invaders
invC = 1/2*((aD - e)*bC - (aD - e)*bD + aC*e - 2*aD*e - aD*sC - aD*sD + sqrt(aC^2*e^2 + aD^2*sC^2 + aD^2*sD^2 + (aD^2 - 2*aD*e + e^2)*bC^2 + (aD^2 - 2*aD*e + e^2)*bD^2 + 2*(aC*aD*e - aC*e^2)*bC + 2*(aC*aD*e - aC*e^2 + (aD^2 - 2*aD*e + e^2)*bC)*bD + 2*(aC*aD*e + (aD^2 - aD*e)*bC + (aD^2 - aD*e)*bD)*sC - 2*(2*aC*aD^2*e - (2*e^2 + e)*aC*aD - aD^2*sC + (aD^2 - aD*e)*bC + (aD^2 - aD*e)*bD)*sD))/aD

invD = -1/2*((aC - e)*bC - (aC - e)*bD + 2*aC*e - aD*e + aC*sC + aC*sD - sqrt(aD^2*e^2 + aC^2*sC^2 + aC^2*sD^2 + 2*(aC*e - e^2)*aD*bC + (aC^2 - 2*aC*e + e^2)*bC^2 + (aC^2 - 2*aC*e + e^2)*bD^2 + 2*((aC*e - e^2)*aD + (aC^2 - 2*aC*e + e^2)*bC)*bD - 2*((2*aC^2*e - (2*e^2 + e)*aC)*aD + (aC^2 - aC*e)*bC + (aC^2 - aC*e)*bD)*sC + 2*(aC*aD*e + aC^2*sC + (aC^2 - aC*e)*bC + (aC^2 - aC*e)*bD)*sD))/aC

# Plot the eigenvalues
#quartz(height = 5, width = 6)
par(mar = c(5,6,2,1), mfrow = c(1,1))
plot(TP, invD, type = "l", col = "lightgreen", ylim = range(invD,invC), xlab = "Temperature", ylab = "Per capita growth rate \nas invader",lwd = 1.5,cex.lab = 1.5, cex.axis = 1.25)
lines(TP,invC, col = "darkgreen",lwd = 1.5)
abline(h = 0,lwd = 2)
legend("topright",lty = 1, col = c("lightgreen","darkgreen"), legend = c("Temperate", "Boreal"),bty = "n",lwd = 2)


#dev.copy2pdf(file = "Invasion_rate")


