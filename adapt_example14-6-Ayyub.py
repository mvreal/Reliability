"""
Created on Mon Oct 18 15:15:00 2021
Reliabilty Analysis
Example 14-6 - Nonlinear limit state function with nonnormal uncorrelated variables
Ayyub, B.M.; McCuen, R.H. Probability, Statistics and Reliability for Engineers and Scientists, 2nd ed.
Boca Raton, Chapman & Hall, 2003.
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Column: g(Y, S, W, L) = YS-WL^2/4 = 0
#


def gfunction(x, d):
    d[0] = 1
    g = x[0] * x[1] - x[2] * x[3] ** 2 / 4
    return g


# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Y', 'vardist': 'normal', 'varmean': 38.00, 'varcov': 0.05},
    {'varname': 'S', 'vardist': 'lognormal', 'varmean': 100.00, 'varcov': 0.05},
    {'varname': 'W', 'vardist': 'Gumbel', 'varmean': 0.30, 'varcov': 0.25},
    {'varname': 'L', 'vardist': 'normal', 'varmean': 180.00, 'varcov': 0.05}
]
# Design variables

dvar = [
    {'varname': 'factor1', 'varvalue': 1.00},
]

#
# MC-IS adaptative method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.adaptive(100, 1000, 0.05, 1.50)
#