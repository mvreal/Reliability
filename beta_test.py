"""

Test for the Beta distribution

"""
import numpy as np
from scipy.stats import norm
from scipy.stats import beta

# Dada input

xpar1 = 0.45
xpar2 = 1.25
xpar3 = 0.22
xpar4 = 0.36
xval = 0.7534482758620691

#
# Beta distribution
#

lower = xpar1
upper = xpar2
a = xpar3
b = xpar4
loc = lower
scale = (upper - lower)
mux = lower + a / (a + b) * (upper - lower)
sigmax = np.sqrt((a * b) / ((a + b) ** 2 * (a + b + 1)) * (upper - lower) ** 2)
pdfx = beta.pdf(xval, a, b, loc, scale)
cdfx = beta.cdf(xval, a, b, loc, scale)
zval = norm.ppf(cdfx)
sigmaxneq = norm.pdf(zval) / pdfx
muxneq = xval - zval * sigmaxneq

print('mux = ', mux)
print('sigmax =', sigmax)
print('xval = ', xval)
print('pdf x = ', pdfx)
print('cdf x = ', cdfx)
print('zval = ', zval)
print('muxneq = ', muxneq)
print('sigmaxneq = ', sigmaxneq)