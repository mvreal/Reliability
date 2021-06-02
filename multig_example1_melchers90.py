"""
Created on Tue Apr 12 14:29:00 2021
Reliabilty Analysis
Example Melchers - 1990
MELCHERS, R.E. Search based importance sampling.
Structural Safety, 9(1990) 117-128
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Limit state funcion
#


def gfunction1(x, d):

    g = d[0] * x[0] * x[3] + d[1] * x[1] ** 2 * x[2] + d[2]
    return g


def gfunction2(x, d):

    g = d[3] * x[0] + d[4] * x[1] + d[5] * x[2] + d[6] * x[3] + d[7]
    return g


#
# Data input
#


# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'X1', 'vardist': 'normal', 'varmean': 0.001, 'varstd': 1.00},
    {'varname': 'X2', 'vardist': 'normal', 'varmean': 0.001, 'varstd': 1.00},
    {'varname': 'X3', 'vardist': 'normal', 'varmean': 0.001, 'varstd': 1.00},
    {'varname': 'X4', 'vardist': 'normal', 'varmean': 0.001, 'varstd': 1.00},
]

dvar = [
    {'varname': 'factor1b', 'varvalue': 3.00},
    {'varname': 'factor2b', 'varvalue': -1.00},
    {'varname': 'constb', 'varvalue': 11.00},
    {'varname': 'factor1a', 'varvalue': 1.00},
    {'varname': 'factor2a', 'varvalue': 1.00},
    {'varname': 'factor3a', 'varvalue': -1.00},
    {'varname': 'factor4a', 'varvalue': -1.00},
    {'varname': 'consta', 'varvalue': 6.00}
]

#
# multig method
#
glist = [gfunction1, gfunction2]
test = Reliability(xvar, dvar, glist, None, None)
test.multig(xvar, dvar, glist)