"""
Created on Mon Sep 17 17:00:00 2021
Reliabilty Analysis
Example 6.6 - Nonlinear limit state function with non-normal variables
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Beam: g(Y, Z, M) = Y*Z-M = 0
#


def gfunction(x, d):

    g = d[0]*x[0]*x[1]-d[1]*x[2]
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 38.00, 'varcov': 0.10},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 60.00, 'varcov': 0.05},
    {'varname': 'M', 'vardist': 'frechet', 'varmean': 1000.00, 'varcov': 0.30}
]
# Design variables

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

# Correlation matrix

#
# FORM method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.form(iHLRF=True, toler=1.e-6)
#