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
    {'varname': 'D', 'vardist': 'uniform', 'parameter1': 0.00, 'parameter2': 1.00},
    {'varname': 'P', 'vardist': 'uniform', 'parameter1': 0.60, 'parameter2': 1.00},
    ]
# Design variables

dvar = [
    {'varname': 'L', 'varvalue': 1.00},
    {'varname': 'epsilon', 'varvalue': 0.050},
    {'varname': 'Pb', 'varvalue': 1.00}
]
#
# MCS method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.mc(100, 100_000, 0.03, 1.00)
#