"""
Created on Tue Apr 12 14:29:00 2021
Reliabilty Analysis
Example Melchers - 1990
MELCHERS, R.E. Search based importance sampling.
Structural Safety, 9(1990) 117-128
@author: MVREAL
"""
from class_reliability import *

#
# Step 0 - Limit state funcion
#


def gfunction(x, ng, ns):

    g = np.zeros((ns, ng))
    g[:, 0] = x[:, 0] + x[:, 1] - x[:, 2] - x[:, 3] + 6.00
    g[:, 1] = 3 * x[:, 0] * x[:, 3] - x[:, 1] ** 2 * x[:, 2] + 11.00

    return g

#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'X1', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
    {'varname': 'X2', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
    {'varname': 'X3', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
    {'varname': 'X4', 'vardist': 'normal', 'varmean': 0.00, 'varstd': 1.00},
]
#
# FORM method
#
test = Reliability(xvar, gfunction)
test.multig_mc(2, 100, 5000, 0.02, 1.00)
#