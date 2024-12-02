"""
Cálculo da probabilidade de falha de um grupo de estacas

Created on Mon Dec 02 10:30:00 2024
Reliabilty Analysis
Example 6-4 - Probability of failure of a group of piles
@author: MVREAL
"""
from realpy import *


def gfunction(x, d):

    g = d[0]*x[0] -d[1]*x[1]
    return g

#
# Data input
#

# Dados de entrada
mu_r = 5000.00
sigma_r= 1500.00
sm = 3000.00
delta_s = 0.33
ns = 100_000
mu_s = sm * np.sqrt(delta_s ** 2 + 1)
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'R', 'vardist': 'normal', 'varmean': mu_r, 'varstd': sigma_r},
    {'varname': 'S', 'vardist': 'lognormal', 'varmean': mu_s, 'varcov': delta_s},
]

# Design variables

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

#
# MCS method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.mc2(100, 100000, 0.05, 1.00)
#