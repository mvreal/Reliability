"""
Cálculo da probabilidade de falha de um grupo de estacas

Created on Mon Dec 02 10:30:00 2024
Reliabilty Analysis
Example 6-5 - Horizontal displacement of a platform
@author: MVREAL
"""
from realpy import *


def gfunction(x, d):

    F = x[0]
    A = x[1]
    B = x[2]
    epsilon = x[3]
    desloc_lim = d[0]

    g = desloc_lim - (A*F + B*F**2 + epsilon)
    return g

#
# Data input
#

# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'F', 'vardist': 'lognormal', 'varmean': 25.00, 'varcov': 0.23},
    {'varname': 'A', 'vardist': 'normal', 'varmean': 0.0113, 'varcov': 0.30},
    {'varname': 'B', 'vardist': 'normal', 'varmean': 0.0006, 'varcov': 0.30},
    {'varname': 'epsilon', 'vardist': 'normal', 'varmean': 0.0000, 'varstd': 0.10},
    
]

# Design variables

dvar = [
    {'varname': 'desloc_lim', 'varvalue': 1.00},
]

#
# FORM method
#
platform = Reliability(xvar, dvar, gfunction, None, None)
platform.form(iHLRF=True, toler=1.e-6)
#