"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example Melchers - 1989
MELCHERS, R.E. Importance sampling in structural systems.
Structural Safety, 6(1989) 3-10
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Limit state funcion
#


def gfunction(x):

    g = 1.44e5 * np.pi * x[0] * x[4] - 24.1344 * x[1] * x[2] ** 2 * x[3]
    return g

#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'X1', 'vardist': 'lognormal', 'varmean': 0.0988, 'varhmean': 0.095759, 'varcov': 0.08666},
    {'varname': 'X2', 'vardist': 'lognormal', 'varmean': 10.6489, 'varhmean': 11.871692, 'varcov': 0.18652},
    {'varname': 'X3', 'vardist': 'frechet', 'varmean': 63.0, 'varhmean': 89.110336, 'varcov': 0.16},
    {'varname': 'X4', 'vardist': 'normal', 'varmean': 0.672,'varhmean': 0.695869, 'varcov': 0.10},
    {'varname': 'X5', 'vardist': 'normal', 'varmean': 38.0, 'varhmean': 36.546130, 'varcov': 0.10}
]
#
# FORM method
#
test = Reliability(xvar, gfunction)
test.mc(100, 5000, 0.02, 1.00)
#