"""
Generation of correlated random variables
Nataf model

"""

import numpy as np
from scipy.stats import norm
from scipy.stats import uniform
from scipy.stats import gumbel_r
from scipy.stats import invweibull
from scipy.stats import weibull_min
import scipy.optimize
from scipy import optimize
import scipy.linalg
from scipy.special import gamma
import pandas as pd
import time


class GenCorr():

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
        if self.corrmatrix is None:
            self.corrmatrix = np.eye(self.n)

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
        Rz = np.array(self.corrmatrix)
        #
        # Cholesky decomposition of the correlation matrix
        #
        L = scipy.linalg.cholesky(Rz, lower=True)
        Jzy = np.copy(L)
        yk = np.random.normal(0.00, 1.00, [ns, self.n])
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
            if namedist.lower() in ['norm', 'normal', 'gauss']:
                mux = float(var['varmean'])
                sigmax = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmax
                x[:, i] = mux + sigmax * zk[:, i]
                fx = norm.pdf(x[:, i], mux, sigmax)
                hx = norm.pdf(x[:, i], muhx, sigmahx)
                weight = weight * (fx / hx)
            #
            # Uniform or constant distribution
            #
            elif namedist.lower() in ['uniform', 'uniforme', 'const']:
                a = float(var['a'])
                b = float(var['b'])
                uk = norm.cdf(zk[:, i])
                x[:, i] = a + (b - a) * uk
                fx = uniform.pdf(x[:, i], a, b)
                hx = uniform.pdf(x[:, i], a, b)
                weight = weight * (fx / hx)

            #
            # Lognormal distribution
            #
            elif namedist.lower() in ['lognormal', 'lognorm', 'log']:
                mux = float(var['varmean'])
                sigmax = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmax
                zetax = np.sqrt(np.log(1.00 + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.5 * zetax ** 2
                zetahx = np.sqrt(np.log(1.00 + (sigmahx / muhx) ** 2))
                lambdahx = np.log(muhx) - 0.5 * zetahx ** 2
                x[:, i] = np.exp(lambdax+zk[:, i]*zetax)
                fx = norm.pdf(np.log(x[:, i]), lambdax, zetax)
                hx = norm.pdf(np.log(x[:, i]), lambdahx, zetahx)
                weight = weight * (fx / hx)

            #
            # Gumbel distribution
            #
            elif namedist.lower() in ['gumbel', 'extvalue1', 'evt1max']:
                mux = float(var['varmean'])
                sigmax = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmax
                alphan = np.pi / np.sqrt(6.00) / sigmax
                un = mux - np.euler_gamma / alphan
                betan = 1.00 / alphan
                alphahn = np.pi / np.sqrt(6.00) / sigmahx
                uhn = muhx - np.euler_gamma / alphahn
                betahn = 1.00 / alphahn
                uk = norm.cdf(zk[:, i])
                x[:, i] = un - betan*np.log(np.log(1./uk))
                fx = gumbel_r.pdf(x[:, i], un, betan)
                hx = gumbel_r.pdf(x[:, i], uhn, betahn)
                weight = weight * (fx / hx)
            #
            # Frechet distribution
            #
            elif namedist.lower() in ['frechet', 'extvalue2', 'evt2max']:
                mux = float(var['varmean'])
                sigmax = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmax
                deltax = sigmax / mux
                kapa0 = 2.50
                gsinal = -1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsinal))
                vn = mux / gamma(1.00 - 1.00 / kapa)
                deltahx = sigmahx / muhx
                kapa0 = 2.50
                gsinal = -1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                vhn = muhx / gamma(1.00 - 1.00 / kapah)
                uk = norm.cdf(zk[:, i])
                x[:, i] = vn / (np.log(1./uk)) ** (1./kapa)
                ynf = x[:, i] / vn
                ynh = x[:, i] / vhn
                fx = invweibull.pdf(ynf, kapa) / vn
                hx = invweibull.pdf(ynh, kapah) / vhn
                weight = weight * (fx / hx)
            #
            #
            # Weibull distribution
            #
            elif namedist.lower() in ['weibull', 'extvalue3', 'evt3min']:
                mux = float(var['varmean'])
                sigmax = float(var['varstd'])
                epsilon = float(var['varinf'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmax
                deltax = sigmax / (mux - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapa = scipy.optimize.newton(fkapa, kapa0, args=(deltax, gsinal))
                w1 = (mux - epsilon) / gamma(1.00 + 1.00 / kapa) + epsilon
                deltahx = sigmahx / (muhx - epsilon)
                kapa0 = 2.50
                gsinal = 1.00
                kapah = scipy.optimize.newton(fkapa, kapa0, args=(deltahx, gsinal))
                w1h = (muhx - epsilon) / gamma(1.00 + 1.00 / kapah) + epsilon
                uk = norm.cdf(zk[:, i])
                x[:, i] = (w1 - epsilon) * (np.log(1. - uk)) ** (1. / kapa)
                ynf = (x[:, i] - epsilon) / (w1 - epsilon)
                ynh = (x[:, i] - epsilon) / (w1h - epsilon)
                fx = weibull_min.pdf(ynf, kapa) / (w1 - epsilon)
                hx = weibull_min.pdf(ynh, kapa) / (w1h - epsilon)
                weight = weight * (fx / hx)
        #
        return x, fx, hx, weight

    def mc(self, ns, nsigma=1.00):
        #
        ti = time.time()
        #
        # Number of variables of the problem
        #
        nfail = 0
        niter = 0
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
        wx = np.ones(ns)
        pdf_fx = np.zeros(ns)
        pdf_hx = np.zeros(ns)
        #
        #
        # Step 1 - Generation of the random numbers according to their appropriate distribution
        #

        xp, pdf_fx, pdf_hx, wp = self.var_gen(ns, nsigma)
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

