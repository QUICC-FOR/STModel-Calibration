
# SUPPOSITIONS
# Densité constante (zero sum) --> permet de simplifier considérablement l'équation de base parce que l'on a pas à se soucier de l'origine d'un espace vide (sinon il faudrait ajouter tout plein de termes de proba)
# Remplacement immédiat
# e: probabilité de mortalité
# Pour l'instant, seulement l'érable à sucre est sensible à la présence du sapin
# Dispersion globale à l'intérieur du peuplement
# Régénération de type lotterie, similaire à ce qui était étudié dans le papier de 2006

F = seq(0,1,0.01)
M = 1-F
a0 = 1
a1 = -1
eM = 0.01
eF = 0.01

lF = 0.5
lM = a0 + a1*F
pMM = (lM*M)/(lM*M+lF*F)

plot(F,pMM,type="l")
abline(a=1,b=-1,lty=3)



%%%%%%%%%%%%%%
syms M F
syms a0 a1 lF eM eF

pMM = a0*M/(a0*M + lF*F)
pFM = 1 - pMM

pMF = (a0+a1)*M/((a0+a1)*M+lF*F)
pFF = 1 - pMF

dM = pMM*(eM*M) + pMF*(eF*F) - eM*M
dF = pFM*(eM*M) + pFF*(eF*F) - eF*F


%%%%%%%
solve(dM,M)

solve(dF,M)


%%%%%%%%%%%
rm(list = ls())
eqM = function(F) -(F*a0*eF*lF - F*eM*lF^2 + F*a1*eF*lF)/(a0^2*eF + a0*a1*eF - a0*eM*lF - a1*eM*lF)
 
eqF = function(F) -(F*a0*eF*lF - F*eM*lF^2 + F*a1*eF*lF)/(a0^2*eF + a0*a1*eF - a0*eM*lF - a1*eM*lF)
 
F = seq(0,1,0.01)
a0 = 1
a1 = -0.85
lF = 0.5
eM = 0.01
eF = 0.01
  
plot(F,eqM(F),ylim = c(0,1), type="l",xlab = "Fir", ylab = "Maple")
lines(F,eqF(F))


%%%%%%%%%%%


dT = 1
nT = 100000
res = numeric(nT/dT)
a0 = 1
a1 = 0
lF = 1.4
eM = 100^-1
eF = 75^-1

M = 0.5
F = 1-M

for(i in 1:(nT/dT)) {
	pMM = a0*M/(a0*M + lF*F)
	pFM = 1 - pMM

	pMF = (a0+a1)*M/((a0+a1)*M+lF*F)
	pFF = 1 - pMF

	dM = pMM*(eM*M) + pMF*(eF*F) - eM*M
	dF = pFM*(eM*M) + pFF*(eF*F) - eF*F

	M = M + dM*dT
	F = F + dF*dT
	
	res[i] = M
	}
plot(1:(nT/dT),res,type="l")


### PROBLÈME: IL N'Y A PAS DE COEXISTENCE! La rétroaction est tellement forte
# Solution: il faut revenir à un cas où le recrutement à faible abondance sera plus fort que sous une lotterie simple avec différence de fécondité
# En d'autres mots, il faut de la densité dépendance négative
# La fonction de calcul du taux de recrutement n'est pas bonne



