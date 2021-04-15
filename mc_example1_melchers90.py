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

    g = 3 * x[0] * x[3] - x[1] ** 2 * x[2] + 11.00
    return g

def gfunction2(x):

    g = x[0] + x[1] - x[2] - x[3] + 6.00
    return g


#
# Data input
#
alpha = np.zeros((2, 4))
beta = np.zeros(2)

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
x0 = [-1, 0, 0, 1]
test = Reliability(xvar, gfunction1, x0, None)
beta[0], x0, alpha[0, :], normgradyk, sigmaxneqk = test.form(iHLRF=False, toler=1.e-3)

#

#
# FORM method for gfunction2
#
test = Reliability(xvar, gfunction2, None, None)
beta[1], x0, alpha[1, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)

#

pf = norm.cdf(-beta)
#
print('Final results:')
print('pf1 =', pf[0])
print('pf2 =', pf[1])


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

pa = np.zeros((2, 2))
pb = np.zeros((2, 2))
for i in range(2):
    for j in range(2):
        if i != j:
            pa[i, j] = norm.cdf(-beta[i]) * norm.cdf(-((beta[j]-ro[i, j]*beta[i])/np.sqrt(1.-ro[i, j]**2)))
            pb[i, j] = norm.cdf(-beta[j]) * norm.cdf(-((beta[i]-ro[i, j]*beta[j])/np.sqrt(1.-ro[i, j]**2)))

print('pa =', pa)
print('pb =', pb)

pfij_inf = np.zeros((2, 2))
pfij_sup = np.zeros((2, 2))
for i in range(2):
    for j in range(2):
        if i != j:
            if ro[i, j] >= 0.00:
                pfij_inf[i, j] = pa[i, j] + pb[i, j]
                pfij_sup[i, j] = np.max([pa[i, j], pb[i, j]])
            else:
                pfij_inf[i, j] = np.min([pa[i, j], pb[i, j]])
                pfij_sup[i, j] = 0.00


pf_inf = pf[0] + np.max([0.00, pf[1] - pfij_inf[1, 0]])
pf_sup = sum(pf) - pfij_sup[1, 0]
print('pf_inf =', pf_inf)
print('pf_sup =', pf_sup)