"""
Created on Mon Mar 15 15:26:00 2021
Reliabilty Analysis
Example 4.8 - Choi et al. Reliability based structural design p. 132.
@author: MVREAL
"""
import numpy as np
from realpy import *

#
# Step 0 - Function: g(x1, x2) = x1**4+2*x2**4-20 = 0
#


def gfunction(x, d):

    g = d[0] * x[0] ** 4 + 2 * d[1] * x[1] ** 4 - 20.00
    return g

#
# Data input
#

# Design variables

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'X1', 'vardist': 'normal', 'varmean': 10.00, 'varcov': 0.50},
    {'varname': 'X2', 'vardist': 'normal', 'varmean': 10.00, 'varcov': 0.50}
]

#
# SORM method
#
test = Reliability(xvar, dvar, gfunction, None, None)
test.sorm(iHLRF=True, toler=1.e-6)
#