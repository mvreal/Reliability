"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.8 - Nonlinear limit state function with non-normal independent variables
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Beam: g(Y, Z, M) = Y*Z-M = 0
#


def gfunction(x, d):

    g = d[0] - d[1] * x[0] / (x[1] * x[2] ** 4)
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'w', 'vardist': 'normal', 'varmean': 10.00, 'varstd': 0.40},
    {'varname': 'E', 'vardist': 'normal', 'varmean': 2.e7, 'varstd': 0.5e7},
    {'varname': 'h', 'vardist': 'normal', 'varmean': 0.40, 'varstd': 0.01}
]

# Design variables

dvar = [
    {'varname': 'gamma1', 'varvalue': 7. / 360.},
    {'varname': 'gamma2', 'varvalue': 360.}
]
#
# FORM method
#
beam = Reliability(xvar, dvar, gfunction)
beam.form(iHLRF=True, toler=1.e-6)
#