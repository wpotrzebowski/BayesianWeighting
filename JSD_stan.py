"""
Solving jensen shannon divergence using mc-stan model
"""
#Read in experimental and simulated curves
#Initialize stan mdoel
#Perform simulations


import sys
import numpy as np
import pystan

stan_code = """
data {
  int n_measures;
  int n_structures;
  vector[n_measures] target_curve;
  vector[n_measures] target_errors;
  matrix[n_measures, n_structures] sim_curves;
  vector[n_structures] energy_priors;
}

parameters {
  simplex[n_structures] weights;
  real<lower=0.0001> scale;
  real boltzman_shift;
}

model {
  vector[n_measures] pred_curve;
  vector[n_measures] alphas;
  alphas = exp(-1.717472947*(boltzman_shift+energy_priors));
  weights ~ dirichlet(alphas);
  pred_curve = sim_curves * weights * scale;
  target_curve ~ normal(pred_curve, target_errors);
}
"""

def InfEntropy(weight):
    S = 0.0
    for wi in weight:
        if wi<0.0000001:
                wi = 0.0000001
        S-=wi*np.log2(wi)
    return S

def JensenShannonDiv(weights_a, weights_b):
    jsd = InfEntropy((weights_a+weights_b)*0.5)-0.5*InfEntropy(weights_a)-0.5*InfEntropy(weights_b)
    return jsd

experimental_file = sys.argv[1]
simulated_file = sys.argv[2]
priors_file = sys.argv[3]

experimental = np.genfromtxt(experimental_file)
simulated = np.genfromtxt(simulated_file)
priors = np.genfromtxt(priors_file)

stan_dat = {"sim_curves": simulated,
            "target_curve": experimental[:,1],
            "target_errors": experimental[:,2],
            "n_measures" : np.shape(experimental)[0],
            "n_structures" : np.shape(simulated)[1],
            "alphas":priors}
sm = pystan.StanModel(model_code=stan_code)
fit = sm.sampling(data=stan_dat, iter=1000, chains=2)

print(fit)

#la = fit.extract(permuted=True)  # return a dictionary of arrays
#mu = la['weights']

## return an array of three dimensions: iterations, chains, parameters
results_array = fit.extract(permuted=False)

nsamples = 0
jsd_sum = 0.0
bayesian_weights = np.zeros(np.shape(simulated)[1])
for iteration in results_array:
    for parameters in iteration:
        current_weights = parameters[:-2]
        bayesian_weights+=current_weights
        nsamples+=1
bayesian_weights=bayesian_weights/nsamples
print(bayesian_weights)

for iteration in results_array:
    for parameters in iteration:
        current_weights = parameters[:-2]
        jsd_sum+=JensenShannonDiv(current_weights, bayesian_weights)
print (np.sqrt(jsd_sum/nsamples))