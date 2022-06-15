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
        if self.corrmatrix is None:
            self.corrmatrix = np.eye(self.n)

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
            if namedist.lower() in ['norm', 'normal', 'gauss']:
                mux = xpar1
                sigmax = xpar2
                muxneq = mux
                sigmaxneq = sigmax
            #
            # Uniform or constant distribution
            #
            elif namedist.lower() in ['uniform', 'uniforme', 'const']:
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
            elif namedist.lower() in ['lognormal', 'lognorm', 'log']:
                mux = xpar1
                sigmax = xpar2
                zetax = np.sqrt(np.log(1. + (sigmax / mux) ** 2))
                lambdax = np.log(mux) - 0.50 * zetax ** 2
                sigmaxneq = zetax * xval
                muxneq = xval * (1. - np.log(xval) + lambdax)
            #
            # Gumbel distribution
            #
            elif namedist.lower() in ['gumbel', 'extvalue1', 'evt1max']:
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
            elif namedist.lower() in ['frechet', 'extvalue2', 'evt2max']:
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
            elif namedist.lower() in ['weibull', 'extvalue3', 'evt3min']:
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
        Rz = np.array(self.corrmatrix)
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
        epsilon = 1e-3
        delta = 1e-3 * np.abs(self.fel(xk1))
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
        return kiter, gxk, erro1, beta, xk, yk, alpha, gradxk


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
        i = -1
        for var in self.xvar:
            i += 1
            var['varstd'] = float(var['varcov']) * float(var['varmean'])
            print(self.xvar[i])
            #
            namedist = var['vardist']
            if namedist.lower() in ['norm', 'normal', 'gauss']:
                mux = float(var['varmean'])
                sigmax = float(var['varstd'])
                muhx = float(var['varhmean'])
                sigmahx = nsigma * sigmax
                x[:, i] = np.random.normal(muhx, sigmahx, ns)
                fx = norm.pdf(x[:, i], mux, sigmax)
                hx = norm.pdf(x[:, i], muhx, sigmahx)
                weight = weight * (fx / hx)
            #
            # Uniform or constant distribution
            #
            elif namedist.lower() in ['uniform', 'uniforme', 'const']:
                a = float(var['a'])
                b = float(var['b'])
                x[:, i] = np.random.uniform(a, b, ns)
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
                x[:, i] = np.random.lognormal(lambdahx, zetahx, ns)
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
                x[:, i] = np.random.gumbel(uhn, betahn, ns)
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
                x[:, i] = invweibull.rvs(kapah, size=ns) * vhn
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
                x[:, i] = weibull_min.rvs(kapah, size=ns) * (w1h - epsilon) + epsilon
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

