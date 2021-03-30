"""
Created on Mon Feb 08 17:12:00 2021
Reliabilty Analysis
Example 7.6 - Nonlinear limit state function with normal correlated variables
@author: MVREAL
"""
from class_reliability import *

#
# Step 0 - Beam: g(Y, Z, M) = Y*Z-M = 0
#


def gfunction(x):

    g = x[0]*x[1]-x[2]
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


# Correlation matrix

corrmatrix = [[1.00, 0.40, 0.00],
              [0.40, 1.00, 0.00],
              [0.00, 0.00, 1.00]]
#
# MCS method
#
beam = Reliability(xvar, gfunction, None, corrmatrix)
beam.mc(100, 5000, 0.03, 1.00)
#