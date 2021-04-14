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


def gfunction1(x):

    g = x[0] + x[1] - x[2] - x[3] + 6.00
    return g


def gfunction2(x):

    g = 3 * x[0] * x[3] - x[1] ** 2 * x[2] + 11.00
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
#
# FORM method for gfunction1
#
x0 = [0.00, 0.00, 0.00, 0.00]
test = Reliability(xvar, gfunction1, x0, None)
beta, x0, alpha, normgradyk, sigmaxneqk = test.form(iHLRF=False, toler=1.e-3)
#
# MC method for gfunction1
#


#
# MC method
#
test = Reliability(xvar, gfunction1, x0)
beta1, pf1, cov_pf1, nsimul1, ttotal1 = test.adaptive(40, 200, 0.05, 1.50, igraph=False)
print('Beta1 =', beta1,'pf1 =', pf1,'cov_pf1 =', cov_pf1, 'nsimul1 =', nsimul1, 'ttotal1 =', ttotal1)

#

#
# FORM method for gfunction2
#
x0 = [1.00, 0.00, 0.00, -1.00]
test = Reliability(xvar, gfunction2, x0, None)
beta, x0, alpha, normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)
#
# MC method for gfunction1
#

#
# MC method
#
test = Reliability(xvar, gfunction2, x0)
beta2, pf2, cov_pf2, nsimul2, ttotal2 = test.adaptive(40, 200, 0.05, 1.50, igraph=False)
print('Beta2 =', beta2,'pf2 =', pf1,'cov_pf2 =', cov_pf2, 'nsimul2 =', nsimul2, 'ttotal2 =', ttotal2)

#
pff = 0.5*pf1 + 0.5*pf2
print('pff =', pff)