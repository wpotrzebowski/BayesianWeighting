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
  int m_measures;
  int n_structures;
  vector[n_measures] target_saxs;
  vector[n_measures] target_saxserr;
  matrix[n_measures, n_structures] sim_saxs;
  vector[m_measures] target_cs;
  vector[m_measures] target_cserr;
  vector[m_measures] sim_cserr;
  matrix[m_measures, n_structures] sim_css;
  vector[n_structures] energy_priors;
}

parameters {
  simplex[n_structures] weights;
  real<lower=0.0001> scale;
  real boltzman_shift;
}

model {
  vector[n_measures] pred_saxs;
  vector[m_measures] pred_css;
  vector[n_structures] alphas;
  alphas = exp(-1.717472947*(boltzman_shift+energy_priors));
  weights ~ dirichlet(alphas);
  pred_saxs = sim_saxs * weights * scale;
  pred_css = sim_css * weights;
  target_saxs ~ normal(pred_saxs, target_saxserr);
  target_cs ~ normal(pred_css, sim_cserr+target_cserr);
}
"""

def calculateChiCrysol(weightedIns, expIns, expErr):
        """
        Calculates chis same way as it is done in crysol
        """
        #Calculate scaling factor
        chi2_=0.0
        mixed_term_ = 0.0
        square_calc_ = 0.0

        Sindex = 0
        for ins in expIns:
            Iobs=ins
            square_err_ = expErr[Sindex]*expErr[Sindex]
            mixed_term_ += Iobs*weightedIns[Sindex]/square_err_
            square_calc_ += weightedIns[Sindex]*weightedIns[Sindex]/square_err_
            Sindex+=1
        scale_factor = mixed_term_/square_calc_

        Sindex = 0
        for ins in expIns:
            Iobs=ins
            square_err_ = expErr[Sindex]*expErr[Sindex]
            chi2_+=(Iobs-scale_factor*weightedIns[Sindex])*(Iobs-scale_factor*weightedIns[Sindex])/square_err_
            Sindex+=1
        chi2_=chi2_/Sindex
        return chi2_

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
experimental_cs_file = sys.argv[4]
simulated_cs_file = sys.argv[5]
simulated_cserr_file = sys.argv[6]


experimental = np.genfromtxt(experimental_file)
simulated = np.genfromtxt(simulated_file)
priors = np.genfromtxt(priors_file)
experimental_cs = np.genfromtxt(experimental_cs_file)
simulated_cs = np.genfromtxt(simulated_cs_file)
simulated_cserr = np.genfromtxt(simulated_cserr_file)

stan_dat = {"sim_saxs": simulated,
            "target_saxs": experimental[:,1],
            "target_saxserr": experimental[:,2],
            "sim_css": simulated_cs,
            "sim_cserr": simulated_cserr,
            "target_cs": experimental_cs[:,0],
            "target_cserr": experimental_cs[:,1],
            "n_measures" : np.shape(experimental)[0],
            "m_measures" : np.shape(experimental_cs)[0],
            "n_structures" : np.shape(simulated)[1],
            "energy_priors":priors}
sm = pystan.StanModel(model_code=stan_code)
fit = sm.sampling(data=stan_dat, iter=10000, chains=4)

print(fit)

#la = fit.extract(permuted=True)  # return a dictionary of arrays
#mu = la['weights']

## return an array of three dimensions: iterations, chains, parameters
results_array = fit.extract(permuted=False, inc_warmup=False)

nsamples = 0
jsd_sum = 0.0
bayesian_weights = np.zeros(np.shape(simulated)[1])
for iteration in results_array:
    for parameters in iteration:
        current_weights = parameters[:np.shape(simulated)[1]]
        bayesian_weights+=current_weights
        nsamples+=1
bayesian_weights=bayesian_weights/nsamples
print(bayesian_weights)

for iteration in results_array:
    for parameters in iteration:
        current_weights = parameters[:np.shape(simulated)[1]]
        jsd_sum+=JensenShannonDiv(current_weights, bayesian_weights)
print (np.sqrt(jsd_sum/nsamples))


print "Crysol Chi: ", calculateChiCrysol(np.dot(bayesian_weights,np.transpose(simulated)), experimental[:,1],
                                          experimental[:,2])

fig = fit.plot(pars="weights")
fig.savefig("stan_weights_ep_cs.png")