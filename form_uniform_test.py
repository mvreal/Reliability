"""
Created on Mon Oct 18 18:50:00 2021
Reliabilty Analysis
Test of uniform distribution in FORM
Beck, A.T. e Ávila da Silva Júnior, C.R.
Strategies for finding the design point under bounded random variables
Structural Safety 58 (2016) 79-93
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Beam: g(D, P) = (Pb.L/4 - e) - P(D - D^2/L)
#


def gfunction(x, d):
    L = d[0]
    epsilon = d[1]
    Pb = d[2]
    D = x[0]
    P = x[1]

    g = (Pb * L / 4 - epsilon) - P * (D - D ** 2 / L)
    return g


# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'D', 'vardist': 'uniform', 'parameter1': 0.00, 'parameter2': 0.60},
    {'varname': 'P', 'vardist': 'uniform', 'parameter1': 0.60, 'parameter2': 1.00},
    ]
# Design variables

dvar = [
    {'varname': 'L', 'varvalue': 1.00},
    {'varname': 'epsilon', 'varvalue': 0.05},
    {'varname': 'Pb', 'varvalue': 1.00}
]
# Correlation matrix

corrmatrix = [[1.00, 0.00],
              [0.00, 1.00]]
#
# FORM method
#
x0 = [0.30, 0.80]
beam = Reliability(xvar, dvar, gfunction, x0, corrmatrix)
beam.form(iHLRF=True, toler=1.e-6)
#