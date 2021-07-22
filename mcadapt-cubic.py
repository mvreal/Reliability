"""
Created on Wed Jul 21 11:25:00 2021
Reliabilty Analysis
Example with a cubic function
Ghalehnovi, M., Rashki, M., & Ameryan, A. (2020).
First order control variates algorithm for reliability analysis of engineering structures. Applied Mathematical Modelling, 77, 829–847.
https://doi.org/10.1016/j.apm.2019.07.049
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Function: g(X) = X1**3+X1**2*X2+X2**3-18 = 0
#


def gfunction(x, d):
    d[0]=1.0
    d[1]=1.0
    g =x[0] ** 3 + x[0] ** 2 * x[1] + x[1] ** 3 - 18.00
    return g


#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation

xvar = [
    {'varname': 'X1', 'vardist': 'normal', 'varmean': 10.00, 'varstd': 5.00},
    {'varname': 'X2', 'vardist': 'normal', 'varmean': 9.90, 'varstd': 5.00}
]

# Design variables

dvar = [
    {'varname': 'gamma1', 'varvalue': 1.00},
    {'varname': 'gamma2', 'varvalue': 1.00}
]

# Correlation matrix


#
# MC-IS method
#
beam = Reliability(xvar, dvar, gfunction, None, None)
beam.adaptive(100, 5000, 0.02, 1.50)
#