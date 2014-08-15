import os
import sys
import time
import cPickle

from numpy        import *
from numpy.random import *
from tssb         import *
from logistic     import *
from util         import *

rand_seed     = 1234
max_data      = 50000
burnin        = 0
num_samples   = 10
checkpoint    = 1000
dp_alpha      = 25.0
dp_gamma      = 1.0
init_drift    = 0.1
alpha_decay   = 0.25
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

seed(rand_seed)

codes = cifar100_codes(max_data)
dims  = codes.shape[1]
root  = Logistic( dims=dims, drift=init_drift )
tssb  = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
              root_node=root, data=codes )

dp_alpha_traces    = zeros((num_samples, 1))
dp_gamma_traces    = zeros((num_samples, 1))
alpha_decay_traces = zeros((num_samples, 1))
drift_traces       = zeros((num_samples, dims))
cd_llh_traces      = zeros((num_samples, 1))

intervals = zeros((7))
print "Starting MCMC run..."
for iter in range(-burnin,num_samples):

    times = [ time.time() ]

    tssb.resample_assignments()
    times.append(time.time())

    tssb.cull_tree()
    times.append(time.time())

    tssb.resample_node_params()
    times.append(time.time())

    root.resample_hypers()
    times.append(time.time())

    tssb.resample_sticks()
    times.append(time.time())

    tssb.resample_stick_orders()
    times.append(time.time())

    if iter > 0:
        tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
    times.append(time.time())
 
    intervals = intervals + diff(array(times)) 

    if iter >= 0:
        dp_alpha_traces[iter]    = tssb.dp_alpha
        dp_gamma_traces[iter]    = tssb.dp_gamma
        alpha_decay_traces[iter] = tssb.alpha_decay
        drift_traces[iter]       = root.drift()
        cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
   
    if iter > 0 and mod(iter, checkpoint) == 0:
        filename = "checkpoints/cifar100-50k-%s-%06d.pkl" % (codename, iter)
        fh = open(filename, 'w')
        cPickle.dump(tssb, fh)
        fh.close()

    if True or mod(iter, 2) == 0:
        (weights, nodes) = tssb.get_mixture()
        print codename, iter, len(nodes), cd_llh_traces[iter], mean(root._drift), tssb.dp_alpha, tssb.dp_gamma, tssb.alpha_decay, " ".join(map(lambda x: "%0.2f" % x, intervals.tolist())), float(root.hmc_accepts)/(root.hmc_accepts+root.hmc_rejects), root.hmc_accepts, root.hmc_rejects
        intervals = zeros((7))
        
    if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
        print "\t%f is best per-data complete data likelihood so far." % (cd_llh_traces[iter]/max_data)
        filename = "bests/cifar100-50k-%s-best.pkl" % (codename)
        fh = open(filename, 'w')
        cPickle.dump(tssb, fh)
        fh.close()


filename = "checkpoints/cifar100-50k-%s-final.pkl" % (codename)
fh = open(filename, 'w')
cPickle.dump({ 'tssb'               : tssb,
               'dp_alpha_traces'    : dp_alpha_traces,
               'dp_gamma_traces'    : dp_gamma_traces,
               'alpha_decay_traces' : alpha_decay_traces,
               'drift_traces'       : drift_traces,
               'cd_llh_traces'      : cd_llh_traces }, fh)
fh.close()

