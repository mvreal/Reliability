"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.9 - Nonlinear limit state function with non-normal correlated variables
@author: MVREAL
"""
import numpy as np
from realpy import *
import time

#
# Step 0 - Column: g(Y, Z) = Y*Z-1140 = 0
#


def gfunction(x):

    g = x[0]*x[1]-1140
    return g


ti = time.time()
#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'Y', 'vardist': 'lognormal', 'varmean': 38.00, 'varcov': 0.10},
    {'varname': 'Z', 'vardist': 'normal', 'varmean': 54.00, 'varcov': 0.05}
]


# SORM method
#
beam = Reliability(xvar, gfunction, None, None)
beam.sorm()
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#