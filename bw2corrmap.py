"""
Reads in Bayesian fit file and calculates corrmap between experimentl and
scaled ensemble avergaed intensity.

Requires freesas to be installed: https://github.com/kif/freesas
"""
import sys
import numpy as np
from freesas.corrmap import gof

fit_file = sys.argv[1]
fit_data = np.genfromtxt(fit_file)
experimental = fit_data[:,1]
simulated = fit_data[:,2]

corr = gof(experimental, simulated)
