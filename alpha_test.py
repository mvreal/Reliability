"""
Test of the correlation matrix from BECK, A.T.
"""
import numpy as np

alpha2 = [[0.132, 0.000, -0.434,  0.434],
          [0.000, 0.014,  0.493, -0.493],
          [0.000, 0.101, -0.450, -0.450]]
print(alpha2)

alpha_sign = np.sign(alpha2)
print(alpha_sign)

alpha = alpha_sign * np.sqrt(np.abs(alpha2))
print(alpha)

ro = np.dot(alpha, alpha.T)
print(ro)