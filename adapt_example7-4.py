""""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.4 - Nonlinear limit state function with  normal independent variables
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
    {'varname': 'Y', 'vardist': 'normal', 'varmean': 40.00, 'varcov': 0.125},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 50.00, 'varcov': 0.05},
    {'varname': 'M', 'vardist': 'normal', 'varmean': 1000.00, 'varcov': 0.20}
]

# Design variables

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

# MC-IS adaptative method
#
beam = Reliability(xvar, dvar, gfunction)
beam.adaptive(100, 5000, 0.03, 1.50)
#