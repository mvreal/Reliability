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

    g = x[0] + 2 * x[2] + 2 * x[3] + x[4] - 5 * x[5] - 5 * x[6]

    return g

def gfunction2(x):

    g = x[1] + 2 * x[2] + x[3] - 5 * x[6]

    return g

def gfunction3(x):

    g = x[0] + 2 * x[1] + x[3] + x[4] - 5 * x[5]

    return g


#
# Data input
#
alpha = np.zeros((3, 7))
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
x0 = [60.00, 60.00, 60.00, 60.00, 60.00, 20.00, 25.00]
test = Reliability(xvar, gfunction1, x0, None)
beta1, x0, alpha[0, :], normgradyk, sigmaxneqk = test.form(iHLRF=False, toler=1.e-3)
#
# MC method
#
test = Reliability(xvar, gfunction1, x0, None)
beta, pf1, cov_pf1, nsimul1, ttotal1 = test.bucher(100, 5000, 0.05, 1.50, igraph=False)
print('Beta1 =', beta1, 'pf1 =', pf1, 'cov_pf1 =', cov_pf1, 'nsimul1 =', nsimul1, 'ttotal1 =', ttotal1)

#

#
# FORM method for gfunction2
#
x0 = [60.00, 60.00, 60.00, 60.00, 60.00, 20.00, 25.00]
test = Reliability(xvar, gfunction2, x0, None)
beta2, x0, alpha[1, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)
#
# MC method
#
test = Reliability(xvar, gfunction2, x0)
beta, pf2, cov_pf2, nsimul2, ttotal2 = test.bucher(100, 5000, 0.05, 1.50, igraph=False)
print('Beta2 =', beta2, 'pf2 =', pf2, 'cov_pf2 =', cov_pf2, 'nsimul2 =', nsimul2, 'ttotal2 =', ttotal2)

#

#
# FORM method for gfunction3
#
x0 = [60.00, 60.00, 60.00, 60.00, 60.00, 20.00, 25.00]
test = Reliability(xvar, gfunction3, x0, None)
beta3, x0, alpha[2, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)
#
# MC method
#
test = Reliability(xvar, gfunction3, x0, None)
beta, pf3, cov_pf3, nsimul3, ttotal3 = test.bucher(100, 5000, 0.05, 1.50, igraph=False)
print('Beta3 =', beta3, 'pf3 =', pf3, 'cov_pf3 =', cov_pf3, 'nsimul3 =', nsimul3, 'ttotal3 =', ttotal3)

#
pff = pf1 + pf2 + pf3
print('Final results:')
print('pf1 =', pf1)
print('cov_pf1 =', cov_pf1)
print('nsimul1 =', nsimul1)
print('pf2 =', pf2)
print('cov_pf2 =', cov_pf2)
print('nsimul2 =', nsimul2)
print('pf3 =', pf3)
print('cov_pf3 =', cov_pf3)
print('nsimul3 =', nsimul3)
print('pff =', pff)

pf = np.array([pf1, pf2, pf3])

pfinf = pf.max()
pfsup = pf.sum()
print('pfinf =', pfinf)
print('pfsup =', pfsup)

b = np.array([beta1, beta2, beta3])
print(b)
print(alpha)
ro = np.ones((3, 3))
for i in range(3):
    for j in range(3):
        ro[i, j] = np.dot(alpha[i], alpha[j])

print(ro)

pa = np.zeros((3, 3))
pb = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        if i != j:
            pa[i, j] = norm.cdf(-b[i]) * norm.cdf(-((b[j]-ro[i, j]*b[i])/np.sqrt(1.-ro[i, j]**2)))
            pb[i, j] = norm.cdf(-b[j]) * norm.cdf(-((b[i]-ro[i, j]*b[j])/np.sqrt(1.-ro[i, j]**2)))

print(pa)
print(pb)

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