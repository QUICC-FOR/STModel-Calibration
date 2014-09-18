# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------
from sympy import *
from scipy import linalg
import numpy as np
from sympy import init_printing
init_printing()
init_session()

alphat, alphab, M, T, B, trRT, trRB, trRM = symbols('alphat alphab M T B trRT trRB trRM')
eq1 = alphat*(M+T)*alphab*(M+B)
eq2 = alphab*(M+B)*(1-alphat*(M+T))
res = solve([Eq(eq1, trRM), Eq(eq2, trRB)], [alphat, alphab])
simplify(res)


# carrying capacity
I, m, p, z, c = symbols('I m p z c')
K = (c*I - (p+q*exp(-z*c*I)))*(1-m) - m
Ik = solve(Eq(K, 0), I)
simplify(Ik)

