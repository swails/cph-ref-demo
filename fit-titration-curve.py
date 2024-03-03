from scipy.optimize import curve_fit
import numpy as np


def f(ph, pka, n):
   return 1.0/(1.0+10.0**(n*(pka-ph)))

pHs = np.array([-3.0, -1.0, 1.0, 3.0])
frac_deprot = np.array([0.000, 0.024, 0.685, 0.996])

params, covariance = curve_fit(f, pHs, frac_deprot)

print(f"pKa = {params[0]}; n (hill coeff.) = {params[1]}")
