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

    g = d[0] * x[0] + d[2] * x[2] + d[3] * x[3] + d[4] * x[4] + d[5] * x[5] + d[6] * x[6]

    return g

#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'X1', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10},
    {'varname': 'X2', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10},
    {'varname': 'X3', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10},
    {'varname': 'X4', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10},
    {'varname': 'X5', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10},
    {'varname': 'X6', 'vardist': 'gumbel', 'varmean': 20.00, 'varcov': 0.30},
    {'varname': 'X7', 'vardist': 'gumbel', 'varmean': 25.00, 'varcov': 0.30},
]

# Design variables

dvar = [
    {'varname': 'factor0', 'varvalue': 1.00},
    {'varname': 'factor1', 'varvalue': 0.00},
    {'varname': 'factor2', 'varvalue': 2.00},
    {'varname': 'factor3', 'varvalue': 2.00},
    {'varname': 'factor4', 'varvalue': 1.00},
    {'varname': 'factor5', 'varvalue': -5.00},
    {'varname': 'factor6', 'varvalue': -5.00}
]
#
# FORM method for gfunction1
#
x0 = [60.00, 60.00, 60.00, 60.00, 60.00, 20.00, 25.00]
test = Reliability(xvar, dvar, gfunction1, x0, None)
beta, x0, alpha, normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)
#
# MC method
#
test = Reliability(xvar, dvar, gfunction1, x0, None)
beta1, pf1, cov_pf1, nsimul1, ttotal1 = test.bucher(100, 5000, 0.05, 2.00, igraph=True)
print('Beta1 =', beta1, 'pf1 =', pf1, 'cov_pf1 =', cov_pf1, 'nsimul1 =', nsimul1, 'ttotal1 =', ttotal1)

#
print('Final results:')
print('pf1 =', pf1)
print('cov_pf1 =', cov_pf1)
print('nsimul1 =', nsimul1)
