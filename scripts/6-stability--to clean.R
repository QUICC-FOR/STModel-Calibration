library(rootSolve)

# Define the model
model = function(t,y,pars) {
	with(as.list(c(y,pars)), {
		# Fraction of empty patches converted into the different states following a disturbance
		pC = aC*(C+M)
		pD = aD*(D+M)
		pM = pC*pD
		pC_ = pC*(1-pD)
		pD_ = pD*(1-pC)

		# Regeneration state
		R = 1 - C - D - M

		# Differential equations describing the dynamics of the state variables
		dCdt = pC_*R + sC*M - bD*(D+M)*C - e*C
		dDdt = pD_*R + sD*M - bC*(C+M)*D - e*D
		dMdt = pM*R + bC*(C+M)*D + bD*(D+M)*C - sC*M - sD*M - e*M		
		list(c(dCdt,dDdt,dMdt))
		})
	}

# Wrapper to compute equilibrium and stability
get_eigen = function(TP,PP,par) {
	with(par, {	
	
	# Compute the logit
	logit_alphad 	= ad0 + adt1*TP + adt2*TP^2 + adp1*PP + adp2*PP^2
	logit_alphac 	= ac0 + act1*TP + act2*TP^2 + acp1*PP + acp2*PP^2
	logit_betad 	= bd0 + bdt1*TP + bdt2*TP^2 + bdp1*PP + bdp2*PP^2
	logit_betac 	= bc0 + bct1*TP + bct2*TP^2 + bcp1*PP + bcp2*PP^2
	logit_thetad	= td0 + tdt1*TP + tdt2*TP^2 + tdp1*PP + tdp2*PP^2
	logit_thetac	= tc0 + tct1*TP + tct2*TP^2 + tcp1*PP + tcp2*PP^2
	logit_eps 		= e0 + et1*TP + et2*TP^2 + ep1*PP + ep2*PP^2

	# Back transform into probabilities
	alphac = exp(logit_alphac)/(1+exp(logit_alphac))
	alphad = exp(logit_alphad)/(1+exp(logit_alphad))
	betac = exp(logit_betac)/(1+exp(logit_betac))
	betad = exp(logit_betad)/(1+exp(logit_betad))
	thetac = exp(logit_thetac)/(1+exp(logit_thetac))
	thetad = exp(logit_thetad)/(1+exp(logit_thetad))
	eps = exp(logit_eps)/(1+exp(logit_eps))

	# Vector of parameters to feed stode
	pars = c(aC = alphac, aD = alphad, bC = betac, bD = betad, sC = thetac, sD = thetad, e = eps)	

	# Initial conditions
	y = c(C = 0.2,  D = 0.2, M = 0.2)

 
	# Get the equilibrium
#	eq = stode(y=y, func=model, parms=pars, pos=TRUE)[[1]]
	eq = runsteady(y=y, func=model, parms=pars,times = c(0, 1000))[[1]]

	# Compute the Jacobian
	J = jacobian.full(y=eq,func=model,parms=pars)
	
	# Compute the largest eigen value
	MaxEig = max(as.numeric(eigen(J)$values))

	c(eq,MaxEig)
	})
}

# Paramètres
# Compute the logit
setwd("/Users/DGravel/Documents/Projects_On_Going/Maple_migration/transition_maple/analyses/v3")

par = read.table("par_v2.txt")

Tgrad = seq(1,6,0.01)

res = matrix(nr = length(Tgrad), nc = 4)

for(i in 1:length(Tgrad)) res[i,] = get_eigen(TP = Tgrad[i], PP = 1.1, par)

# Plot the results (EQUILIBRIUM)
quartz(height = 5, width = 6)
par(mar=c(5,6,2,1))
plot(Tgrad,res[,1],type ="l",ylim=c(0,1.1),cex.axis = 1.25, cex.lab = 1.5, xlab = "Température moyenne annuelle", ylab = "Proportion du paysage",lwd = 2,col = "darkcyan" )
lines(Tgrad,res[,2],col="orange",lwd = 2)
lines(Tgrad,res[,3],col ="palegreen3",lwd = 2)
lines(Tgrad,1-res[,1]-res[,2]-res[,3],col ="darkred",lwd = 2)

legend("topleft",bty = "n", col = c("darkcyan","orange","palegreen3","darkred"),legend = c("Boréal","Tempéré","Mixte","Régénération"),lty = 1,lwd=3)

# Plot the results (STABILITY)
quartz(height = 5, width = 6)
par(mar=c(5,6,2,1))
plot(Tgrad,-res[,4],type ="l",cex.axis = 1.25, cex.lab = 1.5, xlab = "Température moyenne annuelle", ylab = "Résilience",lwd = 2)

#dev.copy2pdf(file = "Stability.pdf")





