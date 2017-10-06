"""
Solving jensen shannon divergence using mc-stan model
"""
#Read in experimental and simulated curves
#Initialize stan mdoel
#Perform simulations


import sys
import numpy as np
import pystan
from MCEvidence import MCEvidence

stan_code = """
data {
  int n_measures;
  int n_structures;
  vector[n_measures] target_curve;
  vector[n_measures] target_errors;
  matrix[n_measures, n_structures] sim_curves;
  vector[n_structures] priors;
}

parameters {
  simplex[n_structures] weights;
  real<lower=0.0001> scale;
}

model {
  vector[n_measures] pred_curve;
  vector[n_structures] alphas;
  alphas = priors;
  weights ~ dirichlet(alphas);
  //target += dirichlet_lpdf(weights | alphas);
  pred_curve = sim_curves * weights * scale;
  //target_curve ~ normal(pred_curve, target_errors);
  target += normal_lpdf(target_curve | pred_curve, target_errors);
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
names_file = sys.argv[4]

experimental = np.genfromtxt(experimental_file)
simulated = np.genfromtxt(simulated_file)
priors = np.genfromtxt(priors_file)
#file_names = open(sys.argv[4]).readlines[0].split(" ")
file_names = np.genfromtxt(names_file,dtype='unicode')
#TODO: Otherwise initialze as a array from standard list
sim_curves = simulated
target_curve = experimental[:,1]
target_errors = experimental[:,2]
n_measures = np.shape(experimental)[0]
n_structures = np.shape(simulated)[1]
alphas = priors
log_file = open("full_bayesian.log","w")

threshold = 0.001
#I beleive model can be compiled once only
sm = pystan.StanModel(model_code=stan_code)
for iteration in range(10):
    log_file.write("Starting iteration "+str(iteration)+" with "
                   +str(n_structures)+" models\n")
    stan_dat = {"sim_curves": sim_curves,
            "target_curve": target_curve,
            "target_errors": target_errors,
            "n_measures" : n_measures,
            "n_structures" : n_structures,
            "priors":alphas}

    fit = sm.sampling(data=stan_dat, iter=2000, chains=4, n_jobs=8)

    #Calculating model evidence using external MCEvidence module
    stan_chain=fit.extract(permuted=True)
    if 'weights' in stan_chain.keys(): stan_chain['samples']=stan_chain.pop('weights')
    if 'lp__' in stan_chain.keys(): stan_chain['loglikes']=stan_chain.pop('lp__')
    try:
        mce=MCEvidence(stan_chain,ischain=True).evidence(prewhiten=False)
        log_file.write("Model evidence: ")
        log_file.writelines(np.array_str(mce)+"\n")
    except:
        log_file.write("Failed to estimate model evidence\n")

    current_weights = fit.summary()['summary'][:,0][:n_structures]
    #fit = sm.vb(data=stan_dat, iter=100000, grad_samples=100)
    #current_weights = fit["mean_pars"][:n_structures]
    #lp_ = fit.summary()['summary'][:,0][-1]
    sim_curves = sim_curves[:,current_weights>threshold]
    alphas = alphas[current_weights>threshold]
    n_structures = np.shape(sim_curves)[1]
    file_names = file_names[current_weights>threshold]
    np.savetxt("fit.txt",fit.summary()['summary'][:,0],delimiter=" ")
    #log_file.writelines(current_weights)
    fig = fit.plot()
    fig.savefig("stan_fit_"+str(iteration)+".png")
    print(fit)

print(file_names)
log_file.close()
## return an array of three dimensions: iterations, chains, parameters
results_array = fit.extract(permuted=False, inc_warmup=False)
nsamples = 0
jsd_sum = 0.0
bayesian_weights = np.zeros(np.shape(sim_curves)[1])
for iteration in results_array:
    for parameters in iteration:
        current_weights = parameters[:n_structures]
        bayesian_weights+=current_weights
        nsamples+=1
bayesian_weights=bayesian_weights/nsamples
print(bayesian_weights)

for iteration in results_array:
    for parameters in iteration:
        current_weights = parameters[:n_structures]
        jsd_sum+=JensenShannonDiv(current_weights, bayesian_weights)
print (np.sqrt(jsd_sum/nsamples))

print "Crysol Chi: ", calculateChiCrysol(np.dot(bayesian_weights,
        np.transpose(sim_curves)), experimental[:,1], experimental[:,2])
fig = fit.plot(pars="weights")
fig.savefig("stan_weights.png")

print "Crysol Chi: ", calculateChiCrysol(np.dot(bayesian_weights,
         np.transpose(sim_curves)), experimental[:,1], experimental[:,2])
