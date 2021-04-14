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
    g[:, 0] = x[:, 0] + 2 * x[:, 2] + 2 * x[:, 3] + x[:, 4] - 5 * x[:, 5] - 5 * x[:, 6]
    g[:, 1] = x[:, 0] + 2 * x[:, 1] + x[:, 3] + x[:, 4] - 5 * x[:, 5]
    g[:, 2] = x[:, 1] + 2 * x[:, 2] + x[:, 3] - 5 * x[:, 6]

    return g

#
# Data input
#
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'X1', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10, 'varhmean': 60.00},
    {'varname': 'X2', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10, 'varhmean': 60.00},
    {'varname': 'X3', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10, 'varhmean': 60.00},
    {'varname': 'X4', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10, 'varhmean': 60.00},
    {'varname': 'X5', 'vardist': 'lognormal', 'varmean': 60.00, 'varcov': 0.10, 'varhmean': 60.00},
    {'varname': 'X6', 'vardist': 'gumbel', 'varmean': 20.00, 'varcov': 0.30, 'varhmean': 20.00},
    {'varname': 'X7', 'vardist': 'gumbel', 'varmean': 25.00, 'varcov': 0.30, 'varhmean': 25.00},
]
#
# MC method
#
test = Reliability(xvar, gfunction)
test.multig_mc(3, 100, 5000, 0.02, 1.00)
#