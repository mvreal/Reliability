"""
Created on Tue Apr 15 11:13:00 2021
Reliabilty Analysis
Example Beck - 2019
BECK, A.T. Confiabilidade e Segurança das Estruturas.
Elsevier, ISBN 978-85-352-8895-7.
@author: MVREAL
"""
from realpy import *

#
# Step 0 - Limit state funcion
#


def gR1(x):
    r1 = 4.00
    a1 = np.pi * r1 ** 2
    v = 300
    k = 1 / 2
    h = k * v
    le = np.sqrt(v ** 2 + h ** 2)

    g = a1*x[0] / 1000 - le / (2 * v) * (x[2] - x[3] / k)

    return g


def gE2(x):
    r2 = 5.20
    i2 = (np.pi * r2 ** 4) / 4
    v = 300
    k = 1 / 2
    h = k * v
    le = np.sqrt(v ** 2 + h ** 2)

    g = np.pi ** 2 * x[1] * i2 / le ** 2 - le / (2 * v) * (+x[2] + x[3] / k)

    return g

def gE1(x):
    r1 = 4.00
    i1 = (np.pi * r1 ** 4) / 4
    v = 300
    k = 1 / 2
    h = k * v
    le = np.sqrt(v ** 2 + h ** 2)

    g = np.pi ** 2 * x[1] * i1 / le ** 2 - le / (2 * v) * (-x[2] + x[3] / k)

    return g


#
# Data input
#
ng = 3
nvar = 4
alpha = np.zeros((ng, nvar))
beta = np.zeros(ng)
# Random variables: name, probability distribution, mean and coefficient of variation


xvar = [
    {'varname': 'S', 'vardist': 'normal', 'varmean': 24.5643, 'varcov': 0.10},
    {'varname': 'E', 'vardist': 'normal', 'varmean': 70.00, 'varcov': 0.03},
    {'varname': 'H', 'vardist': 'normal', 'varmean': 2.00, 'varcov': 0.20},
    {'varname': 'V', 'vardist': 'normal', 'varmean': 1.00, 'varcov': 0.20},
]
#
# FORM method for gfunction1
#
test = Reliability(xvar, gR1, None, None)
beta[0], x0, alpha[0, :], normgradyk, sigmaxneqk = test.form(iHLRF=False, toler=1.e-3)

#

#
# FORM method for gfunction2
#
test = Reliability(xvar, gE2, None, None)
beta[1], x0, alpha[1, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)

#

#
# FORM method for gfunction3
#
test = Reliability(xvar, gE1, None, None)
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

pa = np.zeros((ng, ng))
pb = np.zeros((ng, ng))
for i in range(ng):
    for j in range(ng):
        if i != j:
            pa[i, j] = norm.cdf(-beta[i]) * norm.cdf(-((beta[j]-ro[i, j]*beta[i])/np.sqrt(1.-ro[i, j]**2)))
            pb[i, j] = norm.cdf(-beta[j]) * norm.cdf(-((beta[i]-ro[i, j]*beta[j])/np.sqrt(1.-ro[i, j]**2)))

print('pa =', pa)
print('pb =', pb)

pfij_inf = np.zeros((ng, ng))
pfij_sup = np.zeros((ng, ng))
for i in range(ng):
    for j in range(ng):
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