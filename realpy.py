import numpy as np
from scipy.stats import norm
from scipy.stats import uniform
from scipy.stats import gumbel_r
from scipy.stats import invweibull
from scipy.stats import weibull_min
from scipy.stats import multivariate_normal
from scipy.stats import beta as beta_dist
import scipy.optimize
from scipy import optimize
from scipy.optimize import fsolve
import scipy.linalg
from scipy.special import gamma
import pandas as pd
import matplotlib.pyplot as plt
import time


class Reliability():

    def __init__(self, xvar, dvar, gx, x0=None, Rz=None):
        self.xvar = xvar
        self.dvar = dvar
        self.nxvar = len(xvar)
        self.ndvar = len(dvar)
        self.fel = gx
        self.x0 = x0
        self.Rz = Rz
        #
        # Initial values of the design variables
        #
        i = -1
        self.d = np.zeros(self.ndvar)
        for var in self.dvar:
            i += 1
            # Mean value of the random variables x
            self.d[i] = float(var['varvalue'])
        
        #
        # Setting variables initial values
        #
        i = -1
        for var in self.xvar:
            i += 1
                        
            #
            # Setting standard variable distribution names
            #
            if var['vardist'].lower() in ['norm', 'normal', 'gauss']:
                var['vardist'] = 'gauss'
            elif var['vardist'].lower() in ['uniform', 'uniforme', 'const']:
                var['vardist'] = 'uniform'
            elif var['vardist'].lower() in ['lognormal', 'lognorm', 'log']:
                var['vardist'] = 'lognorm'
            elif var['vardist'].lower() in ['gumbel', 'extvalue1', 'evt1max']:
                var['vardist'] = 'gumbel'
            elif var['vardist'].lower() in ['frechet', 'extvalue2', 'evt2max']:
                var['vardist'] = 'frechet'
            elif var['vardist'].lower() in ['weibull', 'extvalue3', 'evt3min']:
                var['vardist'] = 'weibull'
            elif var['vardist'].lower() in ['Beta', 'beta', 'beta_dist']:
                var['vardist'] = 'beta'

            # For the beta distribution calculates the mean and standard deviation as a function of parameters a, b, q and r
            if var['vardist'] == 'beta':
                a = var['parameter1']
                b = var['parameter2']
                q = var['parameter3']
                r = var['parameter4']
                var['varmean'] = float(a + q / (q + r) * (b - a))
                var['varstd'] = float(np.sqrt((q * r) / ((q + r) **2 * (q + r + 1)) * (b - a) ** 2))
                var['varhmean'] = float(var['varmean'])

        #
        # Setting variables initial values
        #
        i = -1
        for var in self.xvar:
            i += 1
            # Tests if 'varhmean' key exists in var
            if 'varhmean' in var:
                meanhx = float(var['varhmean'])
            else:
                meanhx = float(var['varmean'])
                var.update({'varhmean': meanhx})
            # Tests if 'varstd' key exists in var
            if 'varstd' in var:
                if var['varmean'] > 0:
                    var['varcov'] = float(var['varstd'] / var['varmean'])
                else:
                    var['varcov'] = 1.00
            else:
                var['varstd'] = float(var['varcov'] * var['varmean'])
        #
        # Initial values of the aleatory variables
        #
        if x0 is None:
            #
            # Original mean of the variables x
            #
            i = -1
            self.x0 = np.zeros(self.nxvar)
            for var in self.xvar:
                i += 1
                # Mean value of the random variables x
                self.x0[i] = float(var['varmean'])
                var['varhmean'] = float(var['varmean'])
        else:
            i = -1
            for var in self.xvar:
                i += 1
                # Mean value of the random variables x
                var['varhmean'] = self.x0[i]

        #
        # Initializes the correlation matrix
        #

        if self.Rz is None:
            self.Rz = np.eye(self.nxvar)
        else:
            self.Rz = self.nataf()

#
# Nataf correction of the correlation matrix
#

    def nataf(self):
        """
        Nataf correction of the correlation matrix
        According to:
        Liu, P.-L. and Kiureghian, A.D. Multivariate distribution models with prescribed marginals and covariances
        Probabilistic Engineering Mechanics, 1986, Vol. 1, No.2, p. 105-112
        """
        Rz1 = np.array(self.Rz)
        for i in range(self.nxvar):
            for j in range(i):

                # Variables parameters
                f = 1.00
                ro = self.Rz[i][j]
                cvi = float(self.xvar[i]['varcov'])
                cvj = float(self.xvar[j]['varcov'])

                # Table 4: Xi is gauss and Xj belongs to group 1 - f is constant

                # 1 Xi = gauss and Xj = gauss

                if self.xvar[i]['vardist'] == 'gauss' and self.xvar[j]['vardist'] == 'gauss':
                    f = 1.000

                # 2 Xi = gauss and Xj = uniform

                elif self.xvar[i]['vardist'] == 'gauss' and self.xvar[j]['vardist'] == 'uniform' \
                        or self.xvar[i]['vardist'] == 'uniform' and self.xvar[j]['vardist'] == 'gauss':
                    f = 1.023

                # 3 Xi = gauss and Xj = gumbel

                elif self.xvar[i]['vardist'] == 'gauss' and self.xvar[j]['vardist'] == 'gumbel' \
                        or self.xvar[i]['vardist'] == 'gumbel' and self.xvar[j]['vardist'] == 'gauss':
                    f = 1.031

                # Table 5: Xi is gauss and Xj belongs to group 2 - f depends on cvj

                # 4 Xi = gauss and Xj = lognorm

                elif self.xvar[i]['vardist'] == 'gauss' and self.xvar[j]['vardist'] == 'lognorm' \
                        or self.xvar[i]['vardist'] == 'lognorm' and self.xvar[j]['vardist'] == 'gauss':
                    if self.xvar[i]['vardist'] == 'lognorm':
                        cv = cvi
                    else:
                        cv = cvj
                    f = cv / (np.sqrt(np.log(1.00 + cv ** 2)))

                # 5 Xi = gauss and Xj = frechet

                elif self.xvar[i]['vardist'] == 'gauss' and self.xvar[j]['vardist'] == 'frechet' \
                        or self.xvar[i]['vardist'] == 'frechet' and self.xvar[j]['vardist'] == 'gauss':
                    if self.xvar[i]['vardist'] == 'frechet':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.030 + 0.238 * cv + 0.364 * cv ** 2

                # 6 Xi = gauss and Xj = weibull

                elif self.xvar[i]['vardist'] == 'gauss' and self.xvar[i]['vardist'] == 'weibull' \
                        or self.xvar[i]['vardist'] == 'weibull' and self.xvar[j]['vardist'] == 'gauss':
                    if self.xvar[i]['vardist'] == 'weibull':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.031 - 0.195 * cv + 0.328 * cv ** 2

                # Table 6: Xi  and Xj belongs to group 2 - f depends on ro

                # 7 Xi = uniform and Xj = uniform

                elif self.xvar[i]['vardist'] == 'uniform' and self.xvar[j]['vardist'] == 'uniform':
                    f = 1.047 - 0.047 * ro ** 2

                # 8 Xi = gumbel and Xj = gumbel

                elif self.xvar[i]['vardist'] == 'gumbel' and self.xvar[j]['vardist'] == 'gumbel':
                    f = 1.064 - 0.069 * ro + 0.005 * ro ** 2

                # 9 Xi = uniform and Xj = gumbel

                elif self.xvar[i]['vardist'] == 'uniform' and self.xvar[j]['vardist'] == 'gumbel' \
                        or self.xvar[i]['vardist'] == 'gumbel' and self.xvar[j]['vardist'] == 'uniform':
                    f = 1.055 + 0.015 * ro ** 2

                # Table 7: Xi belongs to group 1 and Xj belongs to group 2 - f depends on ro and cvj

                # 10 Xi = uniform and Xj = lognorm

                elif self.xvar[i]['vardist'] == 'uniform' and self.xvar[j]['vardist'] == 'lognorm' \
                        or self.xvar[i]['vardist'] == 'lognorm' and self.xvar[j]['vardist'] == 'uniform':
                    if self.xvar[i]['vardist'] == 'lognorm':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.019 - 0.014 * cv + 0.010 * ro ** 2 + 0.249 * cv ** 2

                # 11 Xi = uniform and Xj = frechet

                elif self.xvar[i]['vardist'] == 'uniform' and self.xvar[j]['vardist'] == 'frechet' \
                        or self.xvar[i]['vardist'] == 'frechet' and self.xvar[j]['vardist'] == 'uniform':
                    if self.xvar[i]['vardist'] == 'frechet':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.033 + 0.305 * cv + 0.074 * ro ** 2 + 0.405 * cv ** 2

                # 12 Xi = uniform and Xj = weibull

                elif self.xvar[i]['vardist'] == 'uniform' and self.xvar[j]['vardist'] == 'weibull' \
                        or self.xvar[i]['vardist'] == 'weibull' and self.xvar[j]['vardist'] == 'uniform':
                    if self.xvar[i]['vardist'] == 'weibull':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.061 - 0.237 * cv - 0.005 * ro ** 2 + 0.379 * cv ** 2

                # 13 Xi = gumbel and Xj = lognorm

                elif self.xvar[i]['vardist'] == 'gumbel' and self.xvar[j]['vardist'] == 'lognorm' \
                        or self.xvar[i]['vardist'] == 'lognorm' and self.xvar[j]['vardist'] == 'gumbel':
                    if self.xvar[i]['vardist'] == 'lognorm':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.029 + 0.001 * ro + 0.014 * cv + 0.004 * ro ** 2 + 0.233 * cv ** 2 - 0.197 * ro * cv

                # 14 Xi = gumbel and Xj = frechet

                elif self.xvar[i]['vardist'] == 'gumbel' and self.xvar[j]['vardist'] == 'frechet' \
                        or self.xvar[i]['vardist'] == 'frechet' and self.xvar[j]['vardist'] == 'gumbel':
                    if self.xvar[i]['vardist'] == 'frechet':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.056 - 0.060 * ro + 0.263 * cv + 0.020 * ro ** 2 + 0.383 * cv ** 2 - 0.332 * ro * cv

                # 15 Xi = gumbel and Xj = weibull

                elif self.xvar[i]['vardist'] == 'gumbel' and self.xvar[j]['vardist'] == 'weibull' \
                        or self.xvar[i]['vardist'] == 'weibull' and self.xvar[j]['vardist'] == 'gumbel':
                    if self.xvar[i]['vardist'] == 'weibull':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.064 + 0.065 * ro - 0.210 * cv + 0.003 * ro ** 2 + 0.356 * cv ** 2 - 0.211 * ro * cv

                # Table 8 both Xi and Xj belong to group 2: f depends on ro, cvi e cvj

                # 16 Xi = lognorm and Xj = lognorm

                elif self.xvar[i]['vardist'] == 'lognorm' and self.xvar[j]['vardist'] == 'lognorm':
                    f = np.log(1.00 + ro * cvi * cvj)/(ro * np.sqrt(np.log(1.00 + cvi ** 2) * np.log(1.00 + cvj ** 2)))

                # 17 Xi = lognorm and Xj = frechet

                elif self.xvar[i]['vardist'] == 'lognorm' and self.xvar[j]['vardist'] == 'frechet' \
                        or self.xvar[i]['vardist'] == 'frechet' and self.xvar[j]['vardist'] == 'lognorm':
                    if self.xvar[i]['vardist'] == 'frechet':
                        cvf = cvi
                        cvl = cvj
                    else:
                        cvf = cvj
                        cvl = cvi
                    f = 1.026 + 0.082 * ro - 0.019 * cvl + 0.222 * cvf \
                        + 0.018 * ro ** 2 + 0.288 * cvl ** 2 + 0.379 * cvf ** 2 \
                        - 0.441 * ro * cvl + 0.126 * cvl * cvf - 0.277 * ro * cvf

                # 18 Xi = lognorm and Xj = weibull

                elif self.xvar[i]['vardist'] == 'lognorm' and self.xvar[j]['vardist'] == 'weibull' \
                        or self.xvar[i]['vardist'] == 'weibull' and self.xvar[j]['vardist'] == 'lognorm':
                    if self.xvar[i]['vardist'] == 'weibull':
                        cvw = cvi
                        cvl = cvj
                    else:
                        cvw = cvj
                        cvl = cvi
                    f = 1.031 + 0.052 * ro + 0.011 * cvl - 0.210 * cvw \
                        + 0.002 * ro ** 2 + 0.220 * cvl ** 2 + 0.350 * cvw ** 2 \
                        + 0.005 * ro * cvl + 0.009 * cvl * cvw - 0.174 * ro * cvw

                # 19 Xi = frechet and Xj = frechet

                elif self.xvar[i]['vardist'] == 'frechet' and self.xvar[j]['vardist'] == 'frechet':
                    f = 1.086 + 0.054 * ro + 0.104 * (cvi + cvj) \
                        - 0.055 * ro ** 2 + 0.662 * (cvi ** 2 + cvj ** 2)  \
                        - 0.570 * ro * (cvi + cvj) + 0.203 * cvi * cvj \
                        - 0.020 * ro ** 3 - 0.218 * (cvi ** 3 + cvj ** 3) \
                        - 0.371 * ro * (cvi ** 2 + cvj ** 2) + 0.257 * ro ** 2 * (cvi + cvj) \
                        + 0.141 * cvi * cvj * (cvi + cvj)

                # 20 Xi = frechet and Xj = weibull

                elif self.xvar[i]['vardist'] == 'frechet' and self.xvar[j]['vardist'] == 'weibull' \
                        or self.xvar[i]['vardist'] == 'weibull' and self.xvar[j]['vardist'] == 'frechet':
                    if self.xvar[i]['vardist'] == 'frechet':
                        cvf = cvi
                        cvw = cvj
                    else:
                        cvf = cvj
                        cvw = cvi
                    f = 1.065 + 0.146 * ro + 0.241 * cvf - 0.259 * cvw \
                        + 0.013 * ro ** 2 + 0.372 * cvf ** 2 + 0.435 * cvw ** 2  \
                        + 0.005 * ro * cvf + 0.034 * cvf * cvw - 0.481 * ro * cvw

                # 20 Xi = weibull and Xj = weibull

                elif self.xvar[i]['vardist'] == 'weibull' and self.xvar[j]['vardist'] == 'weibull':
                    f = 1.063 - 0.004 * ro - 0.200 * (cvi + cvj) \
                        - 0.001 * ro ** 2 + 0.337 * (cvi ** 2 + cvj ** 2)  \
                        + 0.007 * ro * (cvi + cvj) - 0.007 * cvi * cvj

                # Application of the correction factor f on the ro coefficient
                ro = f * ro
                Rz1[i, j] = ro
                Rz1[j, i] = ro
#        print('Nataf correlation matrix:')
#        print(Rz1)
        return Rz1

    def form(self, iHLRF, toler=1.e-3, iprint=True):
        """

               Algorithm FORM-iHLRF.

        """
      #
        # FORM - First Order Reliability Method with improved HLRF (iHLRF)
        #
        #
        #
        # Penalty function m(y) for FORM-iHLRF algorithm
        #

        def mfunc(normy, g, c):
            my = 1. / 2. * normy ** 2 + c * np.abs(g)
            return my

        #
        #
        # Evaluation of parameter k for Frechet and Weibull distribution
        #

        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk

        #
        # Equivalent normal distribution parameters
        # xval = value of the variable x (scalar)
        # xpar1,xpar2,xpar3,xpar4 = parameters of the original pdf (scalars)
        # namedist = name of the x probability distribution ('string')
        #

        def normeqv(xval, xpar1, xpar2, xpar3, xpar4, namedist):

            #
            # Normal distribution
            #
            if namedist.lower() == 'gauss':
                mux = xpar1
                sigmax = xpar2
                muxneq = mux
                sigmaxneq = sigmax
            #
            # Uniform or constant distribution
            #
            elif namedist.lower() == 'uniform':
                a = xpar1
                b = xpar2
                c = (b - a)
                pdfx = 1. / c
                cdfx = (xval - a) / c
                zval = norm.ppf(cdfx)
                sigmaxneq = (norm.pdf(zval)) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            # Lognormal distribution
            #
            elif namedist.lower() == 'lognorm':
                mux = xpar1
                sigmax = xpar2
                zetax = np.sqrt(np.log(1. + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.50 * zetax ** 2
                sigmaxneq = zetax * xval
                muxneq = xval * (1. - np.log(xval) + lambdax)
            #
            # Gumbel distribution
            #
            elif namedist.lower() == 'gumbel':
                mux = xpar1
                sigmax = xpar2
                alphan = (np.pi / np.sqrt(6.00)) / (sigmax)
                un = mux - np.euler_gamma / alphan
                cdfx = np.exp(-np.exp(-alphan * (xval - un)))
                pdfx = alphan * np.exp(-alphan * (xval - un)) * cdfx
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Frechet distribution
            #
            elif namedist.lower() == 'frechet':
                mux = xpar1
                sigmax = xpar2
                deltax = sigmax / mux
                kapa0 = 2.50
                gsignal = -1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                vn = mux / gamma(1.00 - 1.00 / kapa)
                cdfx = np.exp(-(vn / xval) ** kapa)
                pdfx = kapa / vn * (vn / xval) ** (kapa + 1) * np.exp(-(vn / xval) ** kapa)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Weibull distribution
            #
            elif namedist.lower() == 'weibull':
                mux = xpar1
                sigmax = xpar2
                epsilon = xpar3
                deltax = sigmax / (mux - epsilon)
                kapa0 = 2.50
                gsignal = 1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsignal))
                w1 = (mux - epsilon) / gamma(1.00 + 1.00 / kapa) + epsilon
                y1 = (xval - epsilon) / (w1 - epsilon)
                pdfx = weibull_min.pdf(y1, kapa) / (w1 - epsilon)
                cdfx = weibull_min.cdf(y1, kapa)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq
            #
            #
            # Beta distribution
            #
            elif namedist.lower() == 'beta':
                a = xpar1
                b = xpar2
                q = xpar3
                r = xpar4
                loc = a
                scale = (b - a)
                pdfx = beta_dist.pdf(xval, q, r, loc, scale)
                cdfx = beta_dist.cdf(xval, q, r, loc, scale)
                zval = norm.ppf(cdfx)
                sigmaxneq = norm.pdf(zval) / pdfx
                muxneq = xval - zval * sigmaxneq

            return muxneq, sigmaxneq

        #
        #
        # Data input
        #
        # Number of variables of the problem

        # Equivalent normal mean and standard deviation of the variables
        muxneqk = np.zeros(self.nxvar)
        sigmaxneqk = np.zeros(self.nxvar)
        namevar = []
        dist = []
        mux0 = []
        sigmax0 = []
        par1 = []
        par2 = []
        par3 = []
        par4 = []
        #
        # Original mean and standard deviation of the variables x
        #

        i = -1
        for var in self.xvar:
            i += 1
            # Names of the random variables x
            namevar.append(str(var['varname']))
            # Names of the probability density functions of the variables x
            dist.append(str(var['vardist']))
            # Mean value of the random variables x
            mux0.append(float(var['varmean']))
            # Standard deviation of the random variables x
            if var['varstd'] == 0.00:
                sigmax0.append(float(var['varcov']) * float(var['varmean']))
            else:
                sigmax0.append(float(var['varstd']))
            # Parameter1
            if 'parameter1' in var:
                par1.append(float(var['parameter1']))
            else:
                par1.append(0.00)
            # Parameter2
            if 'parameter2' in var:
                par2.append(float(var['parameter2']))
            else:
                par2.append(0.00)
            # Parameter3
            if 'parameter3' in var:
                par3.append(float(var['parameter3']))
            else:
                par3.append(0.00)
            # Parameter4
            if 'parameter4' in var:
                par4.append(float(var['parameter4']))
            else:
                par4.append(0.00)
            
           
        #
        # Conversion to array format
        #
        mux0 = np.array(mux0)
        sigmax0 = np.array(sigmax0)
        par1 = np.array(par1)
        par2 = np.array(par2)
        par3 = np.array(par3)
        par4 = np.array(par4)
        #
        #   Algorithm FORM-HLRF: Beck, 2019, pag. 101.
        #
        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrices Jxz and Jzx
        #
        Imatrix = np.eye(self.nxvar)
        #
        # Correlation matrix is self.corrmatrix
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(self.Rz, lower=True)
        Jzy = np.copy(L)
        Jyz = np.linalg.inv(L)
        #
        # Step 2 - Initialize de xk value with mux0
        #
        # Initialization of the variable yk1
        # Jacobian matrices of x==>z and z==>y transformations
        D = sigmax0 * Imatrix
        Jzx = np.linalg.inv(D)
        Jyx = np.dot(Jyz, Jzx)
        Jxz = np.copy(D)
        Jxy = np.dot(Jxz, Jzy)
        yk1 = np.zeros(self.nxvar)
    #    xk1 = mux0 + Jxy.dot(yk1)
        xk1 = np.copy(self.x0)
        #
        # Error tolerance for yk and g(x)
        epsilon = toler
        delta = toler * np.abs(self.fel(xk1, self.d))
        # Initial values for errors and iteration counters
        erro1 = 1000.00
        erro2 = 1000.00
        kiter = 0
        # Value of dx increment for the evaluation of the derivatives
        eps = 1.e-6
        #
        while (erro1 > epsilon or erro2 > delta) and kiter < 100:
            #
            kiter += 1
            xk = np.copy(xk1)
            #
            # Calculation of the equivalent normal distribution parameters for xk
            #
            for i in range(self.nxvar):
                xval = xk[i]
                mux = mux0[i]
                sigmax = sigmax0[i]
                namedist = dist[i]
                xpar1 = mux
                xpar2 = sigmax
                xpar3 = par3[i]
                xpar4 = par4[i]
                if dist[i] == 'beta':
                    xpar1 = par1[i]
                    xpar2 = par2[i]

                muxneqk[i], sigmaxneqk[i] = normeqv(xval, xpar1, xpar2, xpar3, xpar4, namedist)
            #
            # Step 3 - Update of the Jacobian matrices Jyx and Jxy
            #
            Dneq = sigmaxneqk * Imatrix
            Jzx = np.linalg.inv(Dneq)
            Jyx = np.dot(Jyz, Jzx)
            Jxz = np.copy(Dneq)
            Jxy = np.dot(Jxz, Jzy)
            #
            #  Step 4 - Transformation from xk to yk
            #
            yk = Jyx.dot(xk - muxneqk)
            normyk = np.linalg.norm(yk)
            beta = np.linalg.norm(yk)

            #
            #  Step 5 - Evaluation of g(xk)
            #
            gxk = self.fel(xk, self.d)

            #
            # Step 6 - Evaluation of the gradients of g(x) in relation to yk
            #
            #
            # a. Calculation of the partial derivatives of g(x) in relation to xk
            #
            gradxk = optimize.approx_fprime(xk, self.fel, eps, self.d)
            #
            # b. Calculation of the partial derivatives of g(x) in relation to yk
            #
            gradyk = np.transpose(Jxy).dot(gradxk)
            normgradyk = np.linalg.norm(gradyk)
            #
            # c. Calculation of the direction cosines for xk
            #
            # Direction cosines
            alpha = gradyk / normgradyk

            #
            # Step 7. Vector yk updating to yk+1 by HLRF algorithm
            #
            dk = ((np.dot(gradyk, yk) - gxk) / normgradyk ** 2) * gradyk - yk
            lambdak = 1.00
            yk1 = yk + lambdak * dk
            #
            # Parameters of iHLRF method
            #
            if iHLRF:
                gamma0 = 2.0
                a = 0.1
                b = 0.5
                #
                gyk = gxk
                normyk = np.linalg.norm(yk)
                normyk1 = np.linalg.norm(yk1)
                c1 = normyk / normgradyk
                #
                if erro2 > delta:
                    c2 = 0.5 * normyk1 ** 2 / np.abs(gyk)
                    ck = gamma0 * np.max([c1, c2])
                else:
                    ck = gamma0 * c1
                #
                k = -1
                f1 = 1.00
                f2 = 0.00
                while f1 > f2 and k < 10:
                    k += 1
                    lambdak = b ** k
                    yk1 = yk + lambdak * dk
                    xk1 = muxneqk + Jxy.dot(yk1)
                    gyk1 = self.fel(xk1, self.d)
                    normyk1 = np.linalg.norm(yk1)
                    f1 = mfunc(normyk1, gyk1, ck) - mfunc(normyk, gyk, ck)
                    gradm = yk + ck * gradyk * np.sign(gyk)
                    normgradm = np.linalg.norm(gradm)
                    f2 = a * lambdak * np.dot(gradm, dk)
            #        f2=-a*lambdak*normgradm**2 # Beck pg. 85: It does not work!!
            #        res=np.array([k,ck,lambdak,gxk,gyk1,f1,f2])
            #        print(res)
            #
            yk1 = yk + lambdak * dk

            #
            # Step 8. Transformation from yk+1 to xk+1
            #
            xk1 = muxneqk + Jxy.dot(yk1)

            #
            # Step 9. Convergence test for yk and g(x)
            #
            prod = normgradyk * normyk
            # Evaluation of the error in the yk1 vector
            if np.abs(prod) > eps:
                erro1 = 1. - np.abs(np.dot(gradyk, yk) / (normgradyk * normyk))
            else:
                erro1 = 1000.00
            # Evaluation of the error in the limit state function g(x)
            erro2 = np.abs(gxk)
            # Printing of the updated values
            if iprint:
                print('\nIteration number = {0:d} g(x) ={1:0.5e} erro1 ={2:0.5e} Beta ={3:0.4f}'
                      .format(kiter, gxk, erro1, beta))
                datadict = {'xvar': namevar, 'prob_dist': dist, 'mux': muxneqk, 'sigmax': sigmaxneqk,
                            'xk': xk, 'yk': yk, 'alpha': alpha}
                data = pd.DataFrame(datadict)
                print(data)
        #
        pf = norm.cdf(-beta)
        if iprint:
            print('\nProbability of Failure Pf = {0:0.4e}'.format(pf))
        return beta, xk, alpha, normgradyk, sigmaxneqk



    def sorm(self, iHLRF=True, toler=1.e-6, iprint=True):
        """
        Second order reliability method = SORM

        """

        #
        # GRAM-SCHMIDT transformation
        #
        def gramschmidt(A, nxvar):
            rk = np.zeros(nxvar)
            rj = np.zeros(nxvar)
            rk0 = np.zeros(nxvar)
            #
            R = np.zeros((nxvar, nxvar))
            R[nxvar - 1, :] = A[nxvar - 1, :].copy()
            for k in range(nxvar - 2, -1, -1):
                rk0 = A[k, :].copy()
                rk0projection = np.zeros(nxvar)
                for j in range(nxvar - 1, k, -1):
                    rj = R[j, :].copy()
                    projection = (rj.dot(rk0)) / (rj.dot(rj))
                    rk0projection = rk0projection + projection * rj
                rk = rk0 - rk0projection
                R[k, :] = rk.copy()
            for i in range(nxvar):
                R[i, :] = R[i, :] / np.linalg.norm(R[i, :])
            #
            return R

        #
        #
        # Function to calculate the second order derivative: d2g/dxidxj
        #
        def second_order_derivative(x, i, j):
            epsilon = 1.e-4  # tolerance for the increments
            h1 = epsilon  # increments: h1 and h2, when i is not equal to j
            h2 = epsilon  # different increments can be adopted
            h = epsilon  # increment h
            a = x[i]  # reference value for x[i]
            b = x[j]  # reference value for x[j]
            #
            # Code: gmn where m and n are equal to:
            # Index 0 = no increment is applied to the variables i and j
            # Index 1 = a decrement equal to -h is applied to the variable i (or j)
            # Index 2 = an incremente equal to +h is applied to the variable i (or j)
            #
            if i == j:
                x0 = np.copy(x)
                x0[i] = a - h
                g10 = self.fel(x0, self.d)
                x0[i] = a
                g00 = self.fel(x0, self.d)
                x0[i] = a + h
                g20 = self.fel(x0, self.d)
                d2g = (g10 - 2. * g00 + g20) / h ** 2  # second order derivative: d2g/dxi2
            else:
                x0 = np.copy(x)
                x0[i] = a + h1
                x0[j] = b + h2
                g22 = self.fel(x0, self.d)
                x0[i] = a + h1
                x0[j] = b - h2
                g21 = self.fel(x0, self.d)
                x0[i] = a - h1
                x0[j] = b + h2
                g12 = self.fel(x0, self.d)
                x0[i] = a - h1
                x0[j] = b - h2
                g11 = self.fel(x0, self.d)
                d2g = (g22 - g21 - g12 + g11) / (4. * h1 * h2)  # second order derivative: d2g/dxidxj
            #
            return d2g

        #
        # First run FORM-iHLRF algorithm
        #
        n = self.nxvar
        xk = np.zeros(n)
        yk = np.zeros(n)
        gradxk = np.zeros(n)
        alpha = np.zeros(n)
        beta = 0.00
        kiter = 0
        erro1 = 0.00

        beta, xk, alpha, normgradyk, sigmaxneqk = self.form(iHLRF, toler, iprint)
        #
        # Formulation of Second Order Reliability Method - SORM
        #
        #
        # Failure probability calculation
        #
        pfform = norm.cdf(-beta)
        #
        # Calculation of the Hessian Matrix
        #
        bmatrix = np.zeros((n, n))
        dmatrix = np.zeros((n, n))
        amatrix = np.eye(n)
        hmatrix = np.zeros((n, n))

        np.set_printoptions(precision=4)

        #
        # Calculation of the Hessian matrix D: d2g/dyidyj
        #
        for i in range(n):
            for j in range(n):
                dmatrix[i, j] = second_order_derivative(xk, i, j) * sigmaxneqk[i] * sigmaxneqk[j]


        #
        # Calculation of the matrix B
        #
        bmatrix = 1. / normgradyk * dmatrix

        #
        # Calculation of the orthogonal matrix H
        #
        amatrix[n - 1, :] = alpha.copy()
        #
        hmatrix = gramschmidt(amatrix, n)


        #
        # Calculation of the curvature matrix K
        #
        kmatrix = hmatrix.dot(bmatrix.dot(hmatrix.T))

        #
        # Calculation of the failure probability using SORM Breitung equation
        #
        factor = 1.00
        for i in range(n - 1):
            factor = factor * 1. / np.sqrt(1.00 + beta * kmatrix[i, i])
        pfsorm = pfform * factor
        betasorm = -norm.ppf(pfsorm)
        #
        # Print the result
        #
        if iprint:
            print('\nSORM results:')
            print('\nHessian matrix:')
            print(dmatrix)
            print('\nNorm of the gradient of g(y) =', normgradyk)
            print('\nB matrix:')
            print(bmatrix)
            print('\nH matrix:')
            print(hmatrix)
            print('\nK = curvatures matrix:')
            print(kmatrix)
            print('\npfFORM =', pfform)
            print('\nBetaFORM =', beta)
            print('\nfactor =', factor)
            print('\npfSORM =', pfsorm)
            print('\nBetaSORM =', betasorm)

        return

    def var_gen(self, ns, uk_cycle, nsigma=1.00, iprint=True):
        """

           Random variables generator for the Monte Carlo Simulation methods

        """

        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk
        
        def beta_limits(vars, mux, sigmax, q, r):
            a, b = vars
            eq1 = a + q / (q + r) * (b - a) - mux
            eq2 = ((q * r) / ((q + r) ** 2 * (q + r + 1))) ** (0.50) * (b - a) - sigmax
            return [eq1, eq2]

        x = np.zeros((ns, self.nxvar))
        weight = np.ones(ns)
        fx = np.zeros(ns)
        hx = np.zeros(ns)
        fxixj = np.ones(ns)
        yk = np.zeros((ns, self.nxvar))



        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrix Jzy
        #
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(self.Rz, lower=True)
        Jzy = np.copy(L)

        #
        # Generation of Gaussian random numbers
        #

        yk = norm.ppf(uk_cycle)
        zf = np.zeros((ns, self.nxvar))
        zk = np.dot(Jzy, yk.T).T


        #
        i = -1
        for var in self.xvar:
            i += 1
            if var['varstd'] == 0.00:
                var['varstd'] = float(var['varcov']) * float(var['varmean'])
            if iprint:
                print(self.xvar[i])
            #
            #
            # Normal distribution
            #
            namedist = var['vardist']
            if namedist.lower() == 'gauss':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                x[:, i] = muhx + sigmahx * zk[:, i]
                fx = norm.pdf(x[:, i], mufx, sigmafx)
                hx = norm.pdf(x[:, i], muhx, sigmahx)
                zf[:, i] = (x[:, i]-mufx)/sigmafx
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)
            #
            # Uniform or constant distribution
            #
            # *To do* Parameters for the sampling function hx for the uniform distribution
            # c = ?, d = ? ---> Verificar
            #
            elif namedist.lower() == 'uniform':
                a = float(var['varmean'])
                b = float(var['varstd'])
                c = float(var['varmean'])
                d = float(var['varstd'])
                uk = norm.cdf(zk[:, i])
                x[:, i] = c + (d - c) * uk
                zf[:, i] = norm.ppf((x[:, i]-a)/(b-a))
                fx = uniform.pdf(x[:, i], a, b)
                hx = uniform.pdf(x[:, i], c, d)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)
            #
            # Lognormal distribution
            #
            elif namedist.lower() == 'lognorm':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                zetafx = np.sqrt(np.log(1.00 + (sigmafx / mufx) ** 2))
                lambdafx = np.log(mufx) - 0.5 * zetafx ** 2
                zetahx = np.sqrt(np.log(1.00 + (sigmahx / muhx) ** 2))
                lambdahx = np.log(muhx) - 0.5 * zetahx ** 2
                x[:, i] = np.exp(lambdahx + zk[:, i] * zetahx)
                zf[:, i] = (np.log(x[:, i])-lambdafx) / zetafx
                fx = norm.pdf(np.log(x[:, i]), lambdafx, zetafx)
                hx = norm.pdf(np.log(x[:, i]), lambdahx, zetahx)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            # Gumbel distribution
            #
            elif namedist.lower() == 'gumbel':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                alphafn = np.pi / np.sqrt(6.00) / sigmafx
                ufn = mufx - np.euler_gamma / alphafn
                betafn = 1.00 / alphafn
                alphahn = np.pi / np.sqrt(6.00) / sigmahx
                uhn = muhx - np.euler_gamma / alphahn
                betahn = 1.00 / alphahn
                uk = norm.cdf(zk[:, i])
                x[:, i] = uhn - betahn * np.log(np.log(1. / uk))
                cdfx = gumbel_r.cdf(x[:, i], ufn, betafn)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                fx = gumbel_r.pdf(x[:, i], ufn, betafn)
                hx = gumbel_r.pdf(x[:, i], uhn, betahn)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            # Frechet distribution
            #
            elif namedist.lower() == 'frechet':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltafx = sigmafx / mufx
                kapa0 = 2.50
                gsinal = -1.00
                kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
                vfn = mufx / gamma(1.00 - 1.00 / kapaf)
                deltahx = sigmahx / muhx
                kapa0 = 2.50
                gsinal = -1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                vhn = muhx / gamma(1.00 - 1.00 / kapah)
                uk = norm.cdf(zk[:, i])
                x[:, i] = vhn / (np.log(1. / uk)) ** (1. / kapah)
                ynf = x[:, i] / vfn
                ynh = x[:, i] / vhn
                cdfx = invweibull.cdf(ynf, kapaf)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                fx = invweibull.pdf(ynf, kapaf) / vfn
                hx = invweibull.pdf(ynh, kapah) / vhn
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            #
            # Weibull distribution
            #
            elif namedist.lower() == 'weibull':
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                epsilon = float(var['varinf'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                deltafx = sigmafx / (mufx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapaf = scipy.optimize.newton(fkapa, kapa0, args=(deltafx, gsinal))
                w1f = (mufx - epsilon) / gamma(1.00 + 1.00 / kapaf) + epsilon
                deltahx = sigmahx / (muhx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                w1h = (muhx - epsilon) / gamma(1.00 + 1.00 / kapah) + epsilon
                uk = norm.cdf(zk[:, i])
                x[:, i] = (w1h - epsilon) * (np.log(1./(1. - uk))) ** (1. / kapah) + epsilon
                ynf = (x[:, i] - epsilon) / (w1f - epsilon)
                ynh = (x[:, i] - epsilon) / (w1h - epsilon)
                cdfx = weibull_min.cdf(ynf, kapaf)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                fx = weibull_min.pdf(ynf, kapaf) / (w1f - epsilon)
                hx = weibull_min.pdf(ynh, kapah) / (w1h - epsilon)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)

            #
            #
            # Beta distribution
            #
            elif namedist.lower() == 'beta':
                a = float(var['parameter1'])
                b = float(var['parameter2'])
                q = float(var['parameter3'])
                r = float(var['parameter4'])
                mufx = float(var['varmean'])
                sigmafx = float(var['varstd'])
                loc = a
                scale = (b - a)
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmafx
                ah, bh =  fsolve(beta_limits, (1, 1), args= ( muhx, sigmahx, q, r))  
                loch = ah
                scaleh = (bh - ah)        
                uk = norm.cdf(zk[:, i])
                x[:, i] = beta_dist.ppf(uk, q, r, loc, scale)
                fx = beta_dist.pdf(x[:, i], q, r, loc, scale)
                hx = beta_dist.pdf(x[:, i], q, r, loch, scaleh)
                cdfx = beta_dist.cdf(x[:, i], q, r, loc, scale)
                zf[:, i] = norm.ppf(cdfx, 0, 1)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))
                fxixj = fxixj * fx / norm.pdf(zf[:, i], 0, 1)
                

        norm_multivarf = multivariate_normal(mean=None, cov=self.Rz)
        phif = list(map(norm_multivarf.pdf, zf))
        phif = np.array(phif)
        norm_multivarh = multivariate_normal(mean=None, cov=self.Rz)
        phih = list(map(norm_multivarh.pdf, zk))
        phih = np.array(phih)
        weight = weight * phif / phih
        fxixj = fxixj * phif

        return x, weight, fxixj

    def mc(self, nc, ns, delta_lim, nsigma=1.00, igraph=True, iprint=True):
        """
        Monte Carlo Simulation Method
        nc Cycles
        ns Simulations
        Brute force = no adaptive technique

        """
        #
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nc = int(nc)
        ns = int(ns)
        pfc = np.zeros(nc)
        cov_pf = np.zeros(nc)
        pf_mean = np.zeros(nc)
        sum1 = 0.00
        sum2 = 0.00
        fxmax_cycle = np.zeros(nc)
        uk_cycle = np.zeros((ns, self.nxvar))
        #
        # Correlation matrix is self.Rz
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)
        #
        # Standard deviation multiplier for MC-IS
        #
        #
        nsigma = 1.00

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.nxvar) for ns Monte Carlo simulations and self.nxvar random variables
        #
        xp = np.zeros((ns, self.nxvar))
        wp = np.ones(ns)
        fx = np.ones(ns)
        zf = np.zeros((ns, self.nxvar))
        zh = np.zeros((ns, self.nxvar))

        # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables

        dmatrix = np.array([self.d.T] * ns)

        #
        # Adaptive cycles
        #

        for icycle in range(nc):
            kcycle = icycle + 1

            #
            # Monte Carlo Simulations
            #
            #
            # Generation of uniform random numbers
            #
            index = icycle % 2
            uk_new = np.random.rand(ns, self.nxvar)
            if index == 0:
                uk_cycle = uk_new.copy()
            else:
                uk_cycle = 1.00 - uk_cycle

            #
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            xp, wp, fx = self.var_gen(ns, uk_cycle, nsigma, iprint)
            #
            #
            # Step 2 - Evaluation of the limit state function g(x)
            #
            gx = list(map(self.fel, xp, dmatrix))
            gx = np.array(gx)

            #
            #
            # Step 3 - Evaluation of the indicator function I[g(x)]
            #
            igx = np.where(gx <= 0.00, wp, 0)
            nfail = sum(igx)
            pfc[icycle] = nfail/ns
            sum1 += pfc[icycle]
            sum2 += pfc[icycle] ** 2
            fxmax_cycle[icycle] = fx.max()

            #
            #  Step 6 - Evaluation of the error in the estimation of Pf
            #

            pf_mean[icycle] = sum1 / kcycle
            pf = pf_mean[icycle]
            if pf > 0.00 and kcycle > 1:
                cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
            else:
                cov_pf[icycle] = 0.00
            delta_pf = cov_pf[icycle]
            # Probability of failure in this cycle
            if iprint:
                print('Cycle =', kcycle, self.xvar)
                print(f'Probability of failure pf ={pf}')
                print(f'Coefficient of variation of pf ={delta_pf}')
            if delta_pf < delta_lim and kcycle > 3:
                break

        beta = -norm.ppf(pf, 0, 1)
        nsimul = kcycle * ns
        tf = time.time()
        ttotal = tf - ti
        #
        if iprint:
            print('*** Resultados do Método Monte Carlo ***')
            print(f'\nReliability Index Beta = {beta}')
            print(f'Probability of failure pf ={pf}')
            print(f'COV of pf ={delta_pf}')
            print('nimul = {0:0.4f} '.format(nsimul))
            print(f'Function g(x): mean = {gx.mean()}, std = {gx.std()} ')
            print(f'Processing time = {ttotal} s')

        if igraph:
            # Plot results:
            cycle = np.arange(1, kcycle + 1, 1)

            plt.figure(1, figsize=(8.5, 6))
            plt.plot(cycle, pf_mean[:kcycle])
            plt.title("Convergence of Probability of Failure")
            plt.xlabel("Cycle")
            plt.ylabel("Pf")
            plt.show()

            plt.figure(2, figsize=(8.5, 6))
            plt.plot(cycle, cov_pf[:kcycle])
            plt.title("CoV of the Probability of Failure")
            plt.xlabel("Cycle")
            plt.ylabel("CoV Pf")
            plt.show()

        return beta, pf, delta_pf, nsimul, ttotal

    def adaptive(self, nc, ns, delta_lim, nsigma=1.50, igraph=True, iprint=True):
        """
        Monte Carlo Simulations with Importance Sampling (MC-IS)
        Importance sampling with adaptative technique
        Melchers, R.E. Search-based importance sampling.
        Structural Safety, 9 (1990) 117-128

        """
        #
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nfail = 0
        niter = 0
        nc = int(nc)
        ns = int(ns)
        pfc = np.zeros(nc)
        cov_pf = np.zeros(nc)
        pf_mean = np.zeros(nc)
        sum1 = 0.00
        sum2 = 0.00
        fxmax = 0.00
        fxmax_cycle = np.zeros(nc)
        uk_cycle = np.zeros((ns, self.nxvar))

        #
        # Correlation matrix is self.Rz
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.nxvar) for ns Monte Carlo simulations and self.nxvar random variables
        #
        xp = np.zeros((ns, self.nxvar))
        wp = np.ones(ns)
        fx = np.ones(ns)
        zf = np.zeros((ns, self.nxvar))
        zh = np.zeros((ns, self.nxvar))


        # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables

        dmatrix = np.array([self.d.T] * ns)

        #
        # Adaptive cycles
        #

        for icycle in range(nc):
            kcycle = icycle + 1

            #
            # Monte Carlo Simulations
            #
            #
            # Generation of uniform random numbers
            #
            index = icycle % 2
            uk_new = np.random.rand(ns, self.nxvar)
            if index == 0:
                uk_cycle = uk_new.copy()
            else:
                uk_cycle = 1.00 - uk_cycle

            #
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            xp, wp, fx = self.var_gen(ns, uk_cycle, nsigma, iprint)
            #
            #
            # Step 2 - Evaluation of the limit state function g(x)
            #
            gx = list(map(self.fel, xp, dmatrix))
            gx = np.array(gx)

            #
            #
            # Step 3 - Evaluation of the indicator function I[g(x)]
            #
            igx = np.where(gx <= 0.00, wp, 0)
            nfail = sum(igx)
            pfc[icycle] = nfail / ns
            sum1 += pfc[icycle]
            sum2 += pfc[icycle] ** 2
            fxmax_cycle[icycle] = fx.max()

            #
            #  Step 4 - Select adaptative mean
            #
            if nfail == 0:
                #
                # No failures in ns simulations
                #
                imin = np.argmin(gx)
                #
                i = -1
                for var in self.xvar:
                    i += 1
                    var['varhmean'] = xp[imin, i]

            else:
                #
                # Ocurrence of nfail failures in ns simulations
                #
                if fxmax_cycle[icycle] > 1.02 * fxmax:
                    fxmax = fxmax_cycle[icycle]
                    imax = np.argmax(fx)
                    #
                    i = -1
                    for var in self.xvar:
                        i += 1
                        var['varhmean'] = xp[imax, i]

            #
            #  Step 6 - Evaluation of the error in the estimation of Pf
            #

            pf_mean[icycle] = sum1 / kcycle
            pf = pf_mean[icycle]
            if pf > 0.00 and kcycle > 1:
                cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
            else:
                cov_pf[icycle] = 0.00
            delta_pf = cov_pf[icycle]
            # Probability of failure in this cycle
            if iprint:
                print('Cycle =', kcycle)
                print(f'Probability of failure pf ={pf}')
                print(f'Coefficient of variation of pf ={delta_pf}')
            if delta_pf < delta_lim and kcycle > 3:
                break

        beta = -norm.ppf(pf, 0, 1)
        nsimul = kcycle * ns
        tf = time.time()
        ttotal = tf - ti
        #
        if iprint:
            print('*** Resultados do Método Monte Carlo ***')
            print(f'\nReliability Index Beta = {beta}')
            print(f'Probability of failure pf ={pf}')
            print(f'COV of pf ={delta_pf}')
            print('nimul = {0:0.4f} '.format(nsimul))
            print(f'Function g(x): mean = {gx.mean()}, std = {gx.std()} ')
            print(f'Processing time = {ttotal} s')

        if igraph:
            # Plot results:
            cycle = np.arange(1, kcycle+1, 1)

            plt.figure(1, figsize=(8.5, 6))
            plt.plot(cycle, pf_mean[:kcycle])
            plt.title("Convergence of Probability of Failure")
            plt.xlabel("Cycle")
            plt.ylabel("Pf")
            plt.show()

            plt.figure(2, figsize=(8.5, 6))
            plt.plot(cycle, cov_pf[:kcycle])
            plt.title("CoV of the Probability of Failure")
            plt.xlabel("Cycle")
            plt.ylabel("CoV Pf")
            plt.show()

        return beta, pf, delta_pf, nsimul, ttotal

    def bucher(self, nc, ns, delta_lim, nsigma=1.50, igraph=True, iprint=True):
        """
        Monte Carlo Simulations with Importance Sampling (MC-IS)
        Importance sampling with adaptive technique
        BUCHER, C.G. Adaptive sampling – an iterative fast Monte Carlo procedure. Structural
        safety, v. 5, n. 2, p. 119-126, 1988.

        """
        #
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nc = int(nc)
        ns = int(ns)
        xm = np.zeros(self.nxvar)
        sum_xwig = np.zeros(self.nxvar)
        sum_wig = 0.00
        pfc = np.zeros(nc)
        cov_pf = np.zeros(nc)
        pf_mean = np.zeros(nc)
        sum1 = 0.00
        sum2 = 0.00
        uk_cycle = np.zeros((ns, self.nxvar))

        #
        # Correlation matrix is self.Rz
        #
        if iprint:
            print('Correlation Matrix after Nataf correction:')
            print(self.Rz)

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.nxvar) for ns Monte Carlo simulations and self.nxvar random variables
        #
        xp = np.zeros((ns, self.nxvar))
        wp = np.ones(ns)
        fx = np.ones(ns)

        # Matrix dmatrix(ns, self.ndvar) for ns Monte Carlo simulations and self.ndvar design variables

        dmatrix = np.array([self.d.T] * ns)

        #
        # Adaptive cycles
        #

        for icycle in range(nc):
            kcycle = icycle + 1

            #
            # Monte Carlo Simulations
            #
            #
            # Generation of uniform random numbers
            #
            index = icycle % 2
            uk_new = np.random.rand(ns, self.nxvar)
            if index == 0:
                uk_cycle = uk_new.copy()
            else:
                uk_cycle = 1.00 - uk_cycle

            #
            #
            # Step 1 - Generation of the random numbers according to their appropriate distribution
            #

            xp, wp, fx = self.var_gen(ns, uk_cycle, nsigma, iprint)
            #
            #
            # Step 2 - Evaluation of the limit state function g(x)
            #
            gx = list(map(self.fel, xp, dmatrix))
            gx = np.array(gx)

            #
            #
            # Step 3 - Evaluation of the indicator function I[g(x)]
            #
            igx = np.where(gx <= 0.00, wp, 0)
            nfail = sum(igx)
            pfc[icycle] = nfail / ns
            sum1 += pfc[icycle]
            sum2 += pfc[icycle] ** 2
            wig = np.copy(igx)

            #
            #  Step 4 - Select adaptive mean
            #
            if nfail == 0:
                #
                # No failures in ns simulations
                #
                imin = np.argmin(gx)
                #
                i = -1
                for var in self.xvar:
                    i += 1
                    xm[i] = xp[imin, i]
                    var['varhmean'] = xm[i]

            else:
                #
                # Ocurrence of nfail failures in ns simulations
                #
                sum_xwig += np.dot(wig.T, xp)
                sum_wig += sum(wig)
                #
                i = -1
                for var in self.xvar:
                    i += 1
                    xm[i] = sum_xwig[i] / sum_wig
                    var['varhmean'] = xm[i]

            #
            #  Step 6 - Evaluation of the error in the estimation of Pf
            #

            pf_mean[icycle] = sum1 / kcycle
            pf = pf_mean[icycle]
            if pf > 0.00 and kcycle > 1:
                cov_pf[icycle] = 1. / (pf * np.sqrt(kcycle * (kcycle - 1))) * np.sqrt(sum2 - 1. / kcycle * sum1 ** 2)
            else:
                cov_pf[icycle] = 0.00
            delta_pf = cov_pf[icycle]
            nc_final = icycle
            # Probability of failure in this cycle
            if iprint:
                print('Cycle =', kcycle)
                print(f'Probability of failure pf ={pf}')
                print(f'Coefficient of variation of pf ={delta_pf}')
            if delta_pf < delta_lim and kcycle > 3:
                break

        beta = -norm.ppf(pf, 0, 1)
        nsimul = kcycle * ns
        tf = time.time()
        ttotal = tf - ti
        #
        if iprint:
            print('*** Resultados do Método Monte Carlo ***')
            print(f'\nReliability Index Beta = {beta}')
            print(f'Probability of failure pf ={pf}')
            print(f'COV of pf ={delta_pf}')
            print('nimul = {0:0.4f} '.format(nsimul))
            print(f'Function g(x): mean = {gx.mean()}, std = {gx.std()} ')
            print(f'Processing time = {ttotal} s')

        if igraph:
            # Plot results:
            cycle = np.arange(1, kcycle+1, 1)

            plt.figure(1, figsize=(8.5, 6))
            plt.plot(cycle, pf_mean[:kcycle])
            plt.title("Convergence of Probability of Failure")
            plt.xlabel("Cycle")
            plt.ylabel("Pf")
            plt.show()

            plt.figure(2, figsize=(8.5, 6))
            plt.plot(cycle, cov_pf[:kcycle])
            plt.title("CoV of the Probability of Failure")
            plt.xlabel("Cycle")
            plt.ylabel("CoV Pf")
            plt.show()

        return beta, pf, delta_pf, nsimul, ttotal

    def multig(self, xvar, dvar, glist, iprint=True):
        """
        Solution of the problem of the reliability of serial system with multiple limit state functions
        According to:
        BECK, A.T.
        Confiabilidade e Segurança das Estruturas
        Rio de Janeiro, Elsevier, 2019.
        """
        ng = int(len(glist))
        nvar = int(len(xvar))
        beta = np.zeros(ng)
        alpha = np.zeros((ng, nvar))
        pf = np.zeros(ng)
        pa = np.zeros((ng, ng))
        pb = np.zeros((ng, ng))
        pfij_inf = np.zeros((ng, ng))
        pfij_sup = np.zeros((ng, ng))

        i = -1
        for gfunction in glist:
            i += 1
            #
            # FORM method for  the multiple g(x) functions
            #
            test = Reliability(xvar, dvar, gfunction, None, None)
            beta[i], x0, alpha[i, :], normgradyk, sigmaxneqk = test.form(iHLRF=True, toler=1.e-3)
            pf[i] = norm.cdf(-beta[i])
        #

        pfinf = pf.max()
        pfsup = pf.sum()
        alpha_sign = np.sign(alpha)
        alpha2 = alpha_sign * alpha ** 2

        #
        # Sort arrays pf, beta, alpha, alpha2, glist in decrescent order of probability of failure
        #
        glist = np.array(glist)
        ig = (-pf).argsort()
        pf = pf[ig]
        beta = beta[ig]
        alpha = alpha[ig, :]
        alpha2 = alpha2[ig, :]
        glist = glist[ig]


        #
        # Calculation of the correlation coefficients between the limit state functions
        #

        ro = np.dot(alpha, alpha.T)

        #
        # Calculation of the p(Aij) and p(Bij) matrices
        #

        for i in range(ng):
            for j in range(ng):
                if i != j:
                    pa[i, j] = norm.cdf(-beta[i]) * norm.cdf(
                        -((beta[j] - ro[i, j] * beta[i]) / np.sqrt(1. - ro[i, j] ** 2)))
                    pb[i, j] = norm.cdf(-beta[j]) * norm.cdf(
                        -((beta[i] - ro[i, j] * beta[j]) / np.sqrt(1. - ro[i, j] ** 2)))


        #
        # Calculation of the terms P(Fi.Fj)sup and P(Fi.Fj)inf
        #

        for i in range(ng):
            for j in range(ng):
                if i != j:
                    if ro[i, j] >= 0.00:
                        pfij_inf[i, j] = pa[i, j] + pb[i, j]
                        pfij_sup[i, j] = np.max([pa[i, j], pb[i, j]])
                    else:
                        pfij_inf[i, j] = np.min([pa[i, j], pb[i, j]])
                        pfij_sup[i, j] = 0.00

        #
        # Calculation of inferior and superior limits for the probability of failure of the system
        #
        #
        # Inferior limit: pf_inf
        #
        pf_inf = pf[0]
        for i in range(1, ng, 1):
            pfi_inf = pf[i]
            for j in range(0, i, 1):
                pfi_inf -= pfij_inf[i, j]
            pf_inf += np.max([0, pfi_inf])

        #
        # Superior limit: pf_sup
        #
        pf_sup = sum(pf)
        for i in range(1, ng, 1):
            pf_sup -= np.max(pfij_sup[i, 0:i], axis=0)

        beta_sup = -norm.ppf(pf_inf)
        beta_inf = -norm.ppf(pf_sup)
        glist = list(glist)

        if iprint:
            #
            # Print the initial results
            #
            print('Initial results:')
            print('pf =', pf)
            print('beta =', beta)
            print('pfinf =', pfinf)
            print('pfsup =', pfsup)
            print('alpha =', alpha)
            print(('alpha2 =', alpha2))
            print('ro =', ro)
            print('pa =', pa)
            print('pb =', pb)
            #
            # Print final results
            #
            print('Final results:')
            print('g list =', glist)
            print('pf =', pf)
            print('beta =', beta)
            print('alpha =', alpha)
            print(('alpha2 =', alpha2))
            print('pf_inf =', pf_inf)
            print('pf_sup =', pf_sup)
            print('beta_inf =', beta_inf)
            print('beta_sup =', beta_sup)



