"""
Created on Mon Mar 15 15:26:00 2021
Reliabilty Analysis
Example 4.8 - Choi et al. Reliability based structural design p. 132.
@author: MVREAL
"""
import numpy as np
from realpy import *
import time

#
# Step 0 - Function: g(x1, x2) = x1**4+2*x2**4-20 = 0
#


def gfunction(x):

    g = x[0] ** 4 + 2 * x[1] ** 4 - 20.00
    return g


ti = time.time()
#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'X1', 'vardist': 'normal', 'varmean': 10.00, 'varcov': 0.50},
    {'varname': 'X2', 'vardist': 'normal', 'varmean': 10.00, 'varcov': 0.50}
]

#
# SORM method
#
test = Reliability(xvar, gfunction, None, None)
test.sorm()
tf = time.time()
ttotal = tf - ti
print(f'Processing time = {ttotal}')
#