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

    g = x[1] + 2 * x[2] + x[3] - 5 * x[6]

    return g


def gfunction2(x):

    g = x[0] + 2 * x[2] + 2 * x[3] + x[4] - 5 * x[5] - 5 * x[6]

    return g


def gfunction3(x):

    g = x[0] + 2 * x[1] + x[3] + x[4] - 5 * x[5]

    return g


#
# Data input
#
alpha = np.zeros((3, 7))
beta = np.zeros(3)
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
#
# FORM method for gfunction1
#
test = Reliability(xvar, gfunction1, None, None)
beta[0], x0, alpha[0, :], normgradyk, sigmaxneqk = test.form(iHLRF=False, toler=1.e-3)

#

#
# FORM method for gfunction2
#
test = Reliability(xvar, gfunction2, None, None)
beta[1], x0, alpha[1, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)

#

#
# FORM method for gfunction3
#
test = Reliability(xvar, gfunction3, None, None)
beta[2], x0, alpha[2, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)

pf = norm.cdf(-beta)
#
print('Final results:')
print('pf1 =', pf[0])
print('pf2 =', pf[1])
print('pf3 =', pf[2])


pfinf = pf.max()
pfsup = pf.sum()
print('pfinf =', pfinf)
print('pfsup =', pfsup)

print('beta =', beta)

print('alpha =', alpha)
alpha_sign = np.sign(alpha)
alpha2 = alpha_sign * alpha ** 2
print(('alpha2 =', alpha2))

ro = np.dot(alpha, alpha.T)
print('ro =', ro)

pa = np.zeros((3, 3))
pb = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        if i != j:
            pa[i, j] = norm.cdf(-beta[i]) * norm.cdf(-((beta[j]-ro[i, j]*beta[i])/np.sqrt(1.-ro[i, j]**2)))
            pb[i, j] = norm.cdf(-beta[j]) * norm.cdf(-((beta[i]-ro[i, j]*beta[j])/np.sqrt(1.-ro[i, j]**2)))

print('pa =', pa)
print('pb =', pb)

pfij_inf = np.zeros((3, 3))
pfij_sup = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        if i != j:
            if ro[i, j] >= 0.00:
                pfij_inf[i, j] = pa[i, j] + pb[i, j]
                pfij_sup[i, j] = np.max([pa[i, j], pb[i, j]])
            else:
                pfij_inf[i, j] = np.min([pa[i, j], pb[i, j]])
                pfij_sup[i, j] = 0.00


pf_inf = pf[0] + np.max([0.00, pf[1] - pfij_inf[1, 0]]) + np.max([0.00, pf[2] - pfij_inf[2, 0] - pfij_inf[2, 1]])
pf_sup = sum(pf) - pfij_sup[1, 0] - np.max([pfij_sup[2, 0], pfij_sup[2, 1]])
print('pf_inf =', pf_inf)
print('pf_sup =', pf_sup)