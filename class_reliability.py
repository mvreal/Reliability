import numpy as np
from scipy.stats import norm
from scipy.stats import uniform
from scipy.stats import gumbel_r
from scipy.stats import invweibull
from scipy.stats import weibull_min
from scipy.stats import multivariate_normal
import scipy.optimize
from scipy import optimize
import scipy.linalg
from scipy.special import gamma
import pandas as pd
import time

class Reliability():

    def __init__(self, xvar, gx, x0=None, corrmatrix=None):
        self.xvar = xvar
        self.n = len(xvar)
        self.fel = gx
        self.x0 = x0
        self.corrmatrix = corrmatrix
        if x0 is None:
            # Original mean of the variables x
            #
            i = -1
            self.x0 = np.zeros(self.n)
            for var in self.xvar:
                i += 1
                # Mean value of the random variables x
                self.x0[i] = float(var['varmean'])
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

        if self.corrmatrix is None:
            self.corrmatrix = np.eye(self.n)
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
        Rz = np.array(self.corrmatrix)
        for i in range(self.n):
            for j in range(i):

                # Variables parameters
                f = 1.00
                ro = self.corrmatrix[i][j]
                cvi = float(self.xvar[i]['varcov'])
                cvj = float(self.xvar[j]['varcov'])

                # Table 4: Xi is gauss and Xj belongs to group 1 - f is constant

                # 1 Xi = gauss and Xj = gauss

                if self.xvar[i]['vardist'] is 'gauss' and self.xvar[j]['vardist'] is 'gauss':
                    f = 1.000

                # 2 Xi = gauss and Xj = uniform

                elif self.xvar[i]['vardist'] is 'gauss' and self.xvar[j]['vardist'] is 'uniform' \
                        or self.xvar[i]['vardist'] is 'uniform' and self.xvar[j]['vardist'] is 'gauss':
                    f = 1.023

                # 3 Xi = gauss and Xj = gumbel

                elif self.xvar[i]['vardist'] is 'gauss' and self.xvar[j]['vardist'] is 'gumbel' \
                        or self.xvar[i]['vardist'] is 'gumbel' and self.xvar[j]['vardist'] is 'gauss':
                    f = 1.031

                # Table 5: Xi is gauss and Xj belongs to group 2 - f depends on cvj

                # 4 Xi = gauss and Xj = lognorm

                elif self.xvar[i]['vardist'] is 'gauss' and self.xvar[j]['vardist'] is 'lognorm' \
                        or self.xvar[i]['vardist'] is 'lognorm' and self.xvar[j]['vardist'] is 'gauss':
                    if self.xvar[i]['vardist'] is 'lognorm':
                        cv = cvi
                    else:
                        cv = cvj
                    f = cv / (np.sqrt(np.log(1.00 + cv ** 2)))

                # 5 Xi = gauss and Xj = frechet

                elif self.xvar[i]['vardist'] is 'gauss' and self.xvar[j]['vardist'] is 'frechet' \
                        or self.xvar[i]['vardist'] is 'frechet' and self.xvar[j]['vardist'] is 'gauss':
                    if self.xvar[i]['vardist'] is 'frechet':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.030 + 0.238 * cv + 0.364 * cv ** 2

                # 6 Xi = gauss and Xj = weibull

                elif self.xvar[i]['vardist'] is 'gauss' and self.xvar[i]['vardist'] is 'weibull' \
                        or self.xvar[i]['vardist'] is 'weibull' and self.xvar[j]['vardist'] is 'gauss':
                    if self.xvar[i]['vardist'] is 'weibull':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.031 - 0.195 * cv + 0.328 * cv ** 2

                # Table 6: Xi  and Xj belongs to group 2 - f depends on ro

                # 7 Xi = uniform and Xj = uniform

                elif self.xvar[i]['vardist'] is 'uniform' and self.xvar[j]['vardist'] is 'uniform':
                    f = 1.047 - 0.047 * ro ** 2

                # 8 Xi = gumbel and Xj = gumbel

                elif self.xvar[i]['vardist'] is 'gumbel' and self.xvar[j]['vardist'] is 'gumbel':
                    f = 1.064 - 0.069 * ro + 0.005 * ro ** 2

                # 9 Xi = uniform and Xj = gumbel

                elif self.xvar[i]['vardist'] is 'uniform' and self.xvar[j]['vardist'] is 'gumbel' \
                        or self.xvar[i]['vardist'] is 'gumbel' and self.xvar[j]['vardist'] is 'uniform':
                    f = 1.055 + 0.015 * ro ** 2

                # Table 7: Xi belongs to group 1 and Xj belongs to group 2 - f depends on ro and cvj

                # 10 Xi = uniform and Xj = lognorm

                elif self.xvar[i]['vardist'] is 'uniform' and self.xvar[j]['vardist'] is 'lognorm' \
                        or self.xvar[i]['vardist'] is 'lognorm' and self.xvar[j]['vardist'] is 'uniform':
                    if self.xvar[i]['vardist'] is 'lognorm':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.019 - 0.014 * cv + 0.010 * ro ** 2 + 0.249 * cv ** 2

                # 11 Xi = uniform and Xj = frechet

                elif self.xvar[i]['vardist'] is 'uniform' and self.xvar[j]['vardist'] is 'frechet' \
                        or self.xvar[i]['vardist'] is 'frechet' and self.xvar[j]['vardist'] is 'uniform':
                    if self.xvar[i]['vardist'] is 'frechet':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.033 + 0.305 * cv + 0.074 * ro ** 2 + 0.405 * cv ** 2

                # 12 Xi = uniform and Xj = weibull

                elif self.xvar[i]['vardist'] is 'uniform' and self.xvar[j]['vardist'] is 'weibull' \
                        or self.xvar[i]['vardist'] is 'weibull' and self.xvar[j]['vardist'] is 'uniform':
                    if self.xvar[i]['vardist'] is 'weibull':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.061 - 0.237 * cv - 0.005 * ro ** 2 + 0.379 * cv ** 2

                # 13 Xi = gumbel and Xj = lognorm

                elif self.xvar[i]['vardist'] is 'gumbel' and self.xvar[j]['vardist'] is 'lognorm' \
                        or self.xvar[i]['vardist'] is 'lognorm' and self.xvar[j]['vardist'] is 'gumbel':
                    if self.xvar[i]['vardist'] is 'lognorm':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.029 + 0.001 * ro + 0.014 * cv + 0.004 * ro ** 2 + 0.233 * cv ** 2 - 0.197 * ro * cv

                # 14 Xi = gumbel and Xj = frechet

                elif self.xvar[i]['vardist'] is 'gumbel' and self.xvar[j]['vardist'] is 'frechet' \
                        or self.xvar[i]['vardist'] is 'frechet' and self.xvar[j]['vardist'] is 'gumbel':
                    if self.xvar[i]['vardist'] is 'frechet':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.056 - 0.060 * ro + 0.263 * cv + 0.020 * ro ** 2 + 0.383 * cv ** 2 - 0.332 * ro * cv

                # 15 Xi = gumbel and Xj = weibull

                elif self.xvar[i]['vardist'] is 'gumbel' and self.xvar[j]['vardist'] is 'weibull' \
                        or self.xvar[i]['vardist'] is 'weibull' and self.xvar[j]['vardist'] is 'gumbel':
                    if self.xvar[i]['vardist'] is 'weibull':
                        cv = cvi
                    else:
                        cv = cvj
                    f = 1.064 + 0.065 * ro - 0.210 * cv + 0.003 * ro ** 2 + 0.356 * cv ** 2 - 0.211 * ro * cv

                # Table 8 both Xi and Xj belong to group 2: f depends on ro, cvi e cvj

                # 16 Xi = lognorm and Xj = lognorm

                elif self.xvar[i]['vardist'] is 'lognorm' and self.xvar[j]['vardist'] is 'lognorm':
                    f = np.log(1.00 + ro * cvi * cvj)/(ro * np.sqrt(np.log(1.00 + cvi ** 2) * np.log(1.00 + cvj ** 2)))

                # 17 Xi = lognorm and Xj = frechet

                elif self.xvar[i]['vardist'] is 'lognorm' and self.xvar[j]['vardist'] is 'frechet' \
                        or self.xvar[i]['vardist'] is 'frechet' and self.xvar[j]['vardist'] is 'lognorm':
                    if self.xvar[i]['vardist'] is 'frechet':
                        cvf = cvi
                        cvl = cvj
                    else:
                        cvf = cvj
                        cvl = cvi
                    f = 1.026 + 0.082 * ro - 0.019 * cvl + 0.222 * cvf \
                        + 0.018 * ro ** 2 + 0.288 * cvl ** 2 + 0.379 * cvf ** 2 \
                        - 0.441 * ro * cvl + 0.126 * cvl * cvf - 0.277 * ro * cvf

                # 18 Xi = lognorm and Xj = weibull

                elif self.xvar[i]['vardist'] is 'lognorm' and self.xvar[j]['vardist'] is 'weibull' \
                        or self.xvar[i]['vardist'] is 'weibull' and self.xvar[j]['vardist'] is 'lognorm':
                    if self.xvar[i]['vardist'] is 'weibull':
                        cvw = cvi
                        cvl = cvj
                    else:
                        cvw = cvj
                        cvl = cvi
                    f = 1.031 + 0.052 * ro + 0.011 * cvl - 0.210 * cvw \
                        + 0.002 * ro ** 2 + 0.220 * cvl ** 2 + 0.350 * cvw ** 2 \
                        + 0.005 * ro * cvl + 0.009 * cvl * cvw - 0.174 * ro * cvw

                # 19 Xi = frechet and Xj = frechet

                elif self.xvar[i]['vardist'] is 'frechet' and self.xvar[j]['vardist'] is 'frechet':
                    f = 1.086 + 0.054 * ro + 0.104 * (cvi + cvj) \
                        - 0.055 * ro ** 2 + 0.662 * (cvi ** 2 + cvj ** 2)  \
                        - 0.570 * ro * (cvi + cvj) + 0.203 * cvi * cvj \
                        - 0.020 * ro ** 3 - 0.218 * (cvi ** 3 + cvj ** 3) \
                        - 0.371 * ro * (cvi ** 2 + cvj ** 2) + 0.257 * ro ** 2 * (cvi + cvj) \
                        + 0.141 * cvi * cvj * (cvi + cvj)

                # 20 Xi = frechet and Xj = weibull

                elif self.xvar[i]['vardist'] is 'frechet' and self.xvar[j]['vardist'] is 'weibull' \
                        or self.xvar[i]['vardist'] is 'weibull' and self.xvar[j]['vardist'] is 'frechet':
                    if self.xvar[i]['vardist'] is 'frechet':
                        cvf = cvi
                        cvw = cvj
                    else:
                        cvf = cvj
                        cvw = cvi
                    f = 1.065 + 0.146 * ro + 0.241 * cvf - 0.259 * cvw \
                        + 0.013 * ro ** 2 + 0.372 * cvf ** 2 + 0.435 * cvw ** 2  \
                        + 0.005 * ro * cvf + 0.034 * cvf * cvw - 0.481 * ro * cvw

                # 20 Xi = weibull and Xj = weibull

                elif self.xvar[i]['vardist'] is 'weibull' and self.xvar[j]['vardist'] is 'weibull':
                    f = 1.063 - 0.004 * ro - 0.200 * (cvi + cvj) \
                        - 0.001 * ro ** 2 + 0.337 * (cvi ** 2 + cvj ** 2)  \
                        + 0.007 * ro * (cvi + cvj) - 0.007 * cvi * cvj

                # Application of the correction factor f on the ro coefficient
                ro = f * ro
                Rz[i, j] = ro
                Rz[j, i] = ro
        print('Nataf correlation matrix:')
        print(Rz)
        return Rz

    def form(self, iHLRF):
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

            return muxneq, sigmaxneq

        #
        #
        # Data input
        #
        # Number of variables of the problem

        # Equivalent normal mean and standard deviation of the variables
        muxneqk = np.zeros(self.n)
        sigmaxneqk = np.zeros(self.n)
        namevar = []
        dist = []
        mux0 = []
        sigmax0 = []
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
            sigmax0.append(float(var['varcov']) * float(var['varmean']))
        #
        # Conversion to array format
        #
        mux0 = np.array(mux0)
        sigmax0 = np.array(sigmax0)
        #
        #   Algorithm FORM-HLRF: Beck, 2019, pag. 101.
        #
        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrices Jxz and Jzx
        #
        Imatrix = np.eye(self.n)
        #
        # Correlation matrix is self.corrmatrix
        #
        Rz = np.eye(self.n)
        Rz = self.nataf()
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(Rz, lower=True)
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
        yk1 = np.zeros(self.n)
    #    xk1 = mux0 + Jxy.dot(yk1)
        xk1 = np.copy(self.x0)
        #
        # Error tolerance for yk and g(x)
        epsilon = 1e-6
        delta = 1e-6 * np.abs(self.fel(xk1))
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
            for i in range(self.n):
                xval = xk[i]
                mux = mux0[i]
                sigmax = sigmax0[i]
                namedist = dist[i]
                muxneqk[i], sigmaxneqk[i] = normeqv(xval, mux, sigmax, 0, 0, namedist)
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
            gxk = self.fel(xk)

            #
            # Step 6 - Evaluation of the gradients of g(x) in relation to yk
            #
            #
            # a. Calculation of the partial derivatives of g(x) in relation to xk
            #
            gradxk = optimize.approx_fprime(xk, self.fel, eps)
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
                    gyk1 = self.fel(xk1)
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
            print('\nIteration number = {0:d} g(x) ={1:0.5e} erro1 ={2:0.5e} Beta ={3:0.4f}'
                  .format(kiter, gxk, erro1, beta))
            datadict = {'xvar': namevar, 'prob_dist': dist, 'mux': muxneqk, 'sigmax': sigmaxneqk,
                        'xk': xk, 'yk': yk, 'alpha': alpha}
            data = pd.DataFrame(datadict)
            print(data)
        #
        pf = norm.cdf(-beta)
        print('\nProbability of Failure Pf = {0:0.4e}'.format(pf))
        return beta, xk, alpha, normgradyk, sigmaxneqk



    def sorm(self):
        """
        Second order reliability method = SORM

        """

        #
        # GRAM-SCHMIDT transformation
        #
        def gramschmidt(A, n):
            rk = np.zeros(n)
            rj = np.zeros(n)
            rk0 = np.zeros(n)
            #
            R = np.zeros((n, n))
            R[n - 1, :] = A[n - 1, :].copy()
            for k in range(n - 2, -1, -1):
                rk0 = A[k, :].copy()
                rk0projection = np.zeros(n)
                for j in range(n - 1, k, -1):
                    rj = R[j, :].copy()
                    projection = (rj.dot(rk0)) / (rj.dot(rj))
                    rk0projection = rk0projection + projection * rj
                rk = rk0 - rk0projection
                R[k, :] = rk.copy()
            for i in range(n):
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
                g10 = self.fel(x0)
                x0[i] = a
                g00 = self.fel(x0)
                x0[i] = a + h
                g20 = self.fel(x0)
                d2g = (g10 - 2. * g00 + g20) / h ** 2  # second order derivative: d2g/dxi2
            else:
                x0 = np.copy(x)
                x0[i] = a + h1
                x0[j] = b + h2
                g22 = self.fel(x0)
                x0[i] = a + h1
                x0[j] = b - h2
                g21 = self.fel(x0)
                x0[i] = a - h1
                x0[j] = b + h2
                g12 = self.fel(x0)
                x0[i] = a - h1
                x0[j] = b - h2
                g11 = self.fel(x0)
                d2g = (g22 - g21 - g12 + g11) / (4. * h1 * h2)  # second order derivative: d2g/dxidxj
            #
            return d2g

        #
        # First run FORM-iHLRF algorithm
        #
        n = self.n
        xk = np.zeros(n)
        yk = np.zeros(n)
        gradxk = np.zeros(n)
        alpha = np.zeros(n)
        beta = 0.00
        kiter = 0
        erro1 = 0.00

        beta, xk, alpha, normgradyk, sigmaxneqk = self.form(iHLRF=True)
        #
        # Formulation of Second Order Reliability Method - SORM
        #
        print('\nSORM results:')
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

        print('\nHessian matrix:')
        print(dmatrix)

        #
        # Calculation of the matrix B
        #
        bmatrix = 1. / normgradyk * dmatrix
        print('\nNorm of the gradient of g(y) =', normgradyk)
        print('\nB matrix:')
        print(bmatrix)

        #
        # Calculation of the orthogonal matrix H
        #
        amatrix[n - 1, :] = alpha.copy()
        #
        hmatrix = gramschmidt(amatrix, n)

        print('\nH matrix:')

        print(hmatrix)

        #
        # Calculation of the curvature matrix K
        #
        kmatrix = hmatrix.dot(bmatrix.dot(hmatrix.T))
        print('\nK = curvatures matrix:')
        print(kmatrix)

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
        print('\npfFORM =', pfform)
        print('\nfactor =', factor)
        print('\npfSORM =', pfsorm)
        print('\nBetaSORM =', betasorm)

        return

    def var_gen(self, ns, nsigma=1.00):
        """

           Algorithm Monte Carlo Method: importance sampling.

        """

        def fkapa(kapa, deltax, gsignal):
            fk = 1.00 + deltax ** 2 - gamma(1.00 + gsignal * 2.00 / kapa) / gamma(1.00 + gsignal * 1.00 / kapa) ** 2
            return fk

        x = np.zeros((ns, self.n))
        weight = np.ones(ns)
        fx = np.zeros(ns)
        hx = np.zeros(ns)

        #
        # Step 1 - Determination of equivalent correlation coefficients and
        #          Jacobian matrix Jzy
        #
        #
        # Correlation matrix is self.corrmatrix
        #
        Rz = np.eye(self.n)
        Rz = self.nataf()
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(Rz, lower=True)
        Jzy = np.copy(L)
        yk = np.random.normal(0.00, 1.00, [ns, self.n])
        zf = np.zeros((ns, self.n))
        zk = np.dot(Jzy, yk.T).T


        #
        i = -1
        for var in self.xvar:
            i += 1
            var['varstd'] = float(var['varcov']) * float(var['varmean'])
            print(self.xvar[i])
            #

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
            #
            # Uniform or constant distribution
            #
            # *To do* Parameters for the sampling function hx for the uniform distribution
            # c = ?, d = ? ---> Verificar
            #
            elif namedist.lower() == 'uniform':
                a = float(var['a'])
                b = float(var['b'])
                c = float(var['c'])
                d = float(var['d'])
                uk = norm.cdf(zk[:, i])
                x[:, i] = c + (d - c) * uk
                zf[:, i] = norm.ppf((x[:, i]-a)/(b-a))
                fx = uniform.pdf(x[:, i], a, b)
                hx = uniform.pdf(x[:, i], a, b)
                weight = weight * ((fx/norm.pdf(zf[:, i], 0, 1)) / (hx/norm.pdf(zk[:, i], 0, 1)))

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

        norm_multivarf = multivariate_normal(mean=None, cov=Rz)
        phif = list(map(norm_multivarf.pdf, zf))
        phif = np.array(phif)
        norm_multivarh = multivariate_normal(mean=None, cov=Rz)
        phih = list(map(norm_multivarh.pdf, zk))
        phih = np.array(phih)
        weight = weight * phif / phih



        return x, weight

    def mc(self, ns, nsigma=1.00):
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nfail = 0
        niter = 0
        ns = int(ns)
        mean_x = np.zeros(self.n)
        std_x = np.zeros(self.n)
        meanh_x = np.zeros(self.n)
        stdh_x = np.zeros(self.n)
        #
        # Standard deviation multiplier for MC-IS
        #
        #
        i = -1
        for var in self.xvar:
            i += 1
            mean_x[i] = float(var['varmean'])
            std_x[i] = float(var['varcov']) * float(var['varmean'])
            if var.get('varhmean', None) == None:
                var['varhmean'] = var['varmean']
            meanh_x[i] = float(var['varhmean'])
            stdh_x[i] = nsigma * std_x[i]

        #
        #
        # Number of Monte Carlo simulations
        #
        #
        # Matrix xp(ns, self.n) for ns Monte Carlo simulations and self.n random variables
        #
        xp = np.zeros((ns, self.n))
        wp = np.ones(ns)
        zf = np.zeros((ns, self.n))
        zh = np.zeros((ns, self.n))
        Rz = np.array(self.corrmatrix)
        #
        #
        # Step 1 - Generation of the random numbers according to their appropriate distribution
        #

        xp, wp = self.var_gen(ns, nsigma)

        #
        #
        # Step 2 - Evaluation of the limit state function g(x)
        #
        gx = list(map(self.fel, xp))
        gx = np.array(gx)

        #
        #
        # Step 3 - Evaluation of the indicator function I[g(x)]
        #
        igx = np.where(gx <= 0.00, wp, 0)
        nfail = sum(igx)

        #
        #  Step 4 - Evaluation of the failure probability Pf
        #
        pf = nfail / ns
        beta = -norm.ppf(pf)

        #
        #  Step 6 - Evaluation of the error in the estimation of Pf
        #
        if pf > 0.00:
            delta_pf = 1. / (pf * np.sqrt(ns * (ns - 1))) * np.sqrt(
                sum((igx * wp) ** 2) - 1. / ns * (sum(igx * wp)) ** 2)
        else:
            delta_pf = 9999

        tf = time.time()
        ttotal = tf - ti
        #
        print('*** Resultados do Método Monte Carlo ***')
        print(f'\nReliability Index Beta = {beta}')
        print(f'Probability of failure pf ={pf}')
        print(f'COV of pf ={delta_pf}')
        print('nimul = {0:0.4f} '.format(ns))
        print(f'Function g(x): mean = {gx.mean()}, std = {gx.std()} ')
        print(f'Processing time = {ttotal} s')

        return beta, pf, delta_pf,ttotal

