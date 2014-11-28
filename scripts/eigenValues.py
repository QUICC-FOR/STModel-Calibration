from sympy import *
from scipy import linalg
import numpy as np
from sympy import init_printing

T, B, M, R = symbols('T B M R')
theta, thetaT, eps, alphaB, alphaT, betaT, betaB = symbols('theta thetaT eps alphaB alphaT betaT betaB')

# Empty state
R = 1 - T - B - M


# Differential equations describing the dynamics of the state variables
dBdt = theta*(1-thetaT)*(1-eps)*M + alphaB*(M+B)*(1-alphaT*(T+M))*R - betaT*(T+M)*(1-eps)*B - eps*B
dTdt = theta*thetaT*(1-eps)*M + alphaT*(M+T)*(1-alphaB*(B+M))*R - betaB*(B+M)*(1-eps)*T - eps*T
dMdt = betaT*(T+M)*(1-eps)*B + betaB*(B+M)*(1-eps)*T - theta*(1-eps)*M  - eps*M

# Compute the Jacobian
systemEq = Matrix([dTdt, dBdt, dMdt])
stateVar = Matrix([T, B, M])

# equilibrium = solve(systemEq, T,B,M)

Jacob = systemEq.jacobian(stateVar)


# Solve the model for the case where T and M are at 0
# and B = 1 - eps/alphaB (ie a l'équilibre si T=0 et M=0)
J_B = Jacob.subs([(T,0),(B,1-eps/alphaB) , (M,0)])
eigs_J_B = J_B.eigenvals().keys()
invT = eigs_J_B

# Solve the model for the case where B and M are at 0
# and T = 1 - e/at (ie a l'équilibre si B=0 et M=0)
J_T = Jacob.subs([(T,1-eps/alphaT), (B,0), (M,0)])
eigs_J_T = J_T.eigenvals().keys()
invB = eigs_J_T

##  solve eq
mu1, alphaB, B, alphaT, mu2, eps, M = symbols('mu1 alphaB B alphaT mu2 eps M')
equ = mu1*M + alphaB*(M+B)*(1-alphaT*M)*(1-M-B)-mu2*M*B-eps*B

equilibrium = solve(equ,M)


