"""
Cálculo da probabilidade de falha para a concentração de cloretos
Mauro Real
Versão com D0 como variável aleatória normal
09/02/2023
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import norm
from scipy.stats import lognorm
from scipy.stats import beta
from scipy.special import erf
from realpy import *




#
# Step 0 - Column: g(Y, Z) = Y*Z-1140 = 0
#


def gfunction(x, d):

     #
     # Função estado limite de despassivação
     #
     tsl = d[0]
     t0 = d[1]
     D0 = x[0]
     Ccr = x[1]
     Cs = x[2]
     alpha = x[3]
     xd = x[4]
     tlim = 30.
     eta = (t0/tlim)**alpha
     g = Ccr - Cs*(1.00 - erf((0.001*xd)/(2*np.sqrt(D0*eta*tsl))))

     return g


#
# Data input
#

# Tempo de vida útil requerido
tsl= 120 # vida útil (service life) em anos

# Dados de entrada determinísticos

t0 =float(28./365.) #t0 a idade do concreto quando exposto aos íons [anos]

# Coeficiente de difusão de cloretos em m2/anos - distribuição lognormal
mediaD0 = 3.e-12 * 3600 * 24 * 365 # m2/ano 
desvioD0 = 0.20*mediaD0
#
# Concentração crítica de cloretos - distribuição beta
lower = 0.45
upper = 1.25 
q = 0.22
r = 0.36

#
# Concentração superficial de cloretos - distribuição lognormal
mediaCs=5.40
desvioCs=0.82

#
# alpha = fator de envelhecimento do concreto - distribuição normal
mediaalpha=0.47
desvioalpha=0.029

#
# xd = cobrimento da armadura - distribuição normal
mediaxd = 72.20 # mm
desvioxd = 5.3 # mm


# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'D0', 'vardist': 'lognormal', 'varmean': mediaD0, 'varstd': desvioD0 },
    {'varname': 'Ccr', 'vardist': 'beta', 'parameter1': lower, 'parameter2': upper, 'parameter3': q, 'parameter4': r},
    {'varname': 'Cs', 'vardist': 'lognormal', 'varmean': mediaCs, 'varstd': desvioCs },
    {'varname': 'alpha', 'vardist': 'normal', 'varmean': mediaalpha, 'varstd': desvioalpha },
    {'varname': 'xd', 'vardist': 'normal', 'varmean': mediaxd, 'varstd': desvioxd },    
]

# Design variables

dvar = [
    {'varname': 'tsl', 'varvalue': tsl},
    {'varname': 't0', 'varvalue': t0}
]

#
# FORM method
#
abaco = Reliability(xvar, dvar, gfunction, None)
abaco.form2(iHLRF=True, toler=1.e-3)
#







