import os
import sys
import time
import cPickle
import csv
import ipdb
import yaml

from collections import namedtuple

from ordereddict import OrderedDict

from numpy         import *
from numpy.random  import *
from tssb          import *
from pyclone_test  import *
from util          import *
from scipy.stats   import itemfreq
from config        import *
    
error_rate    = 0.001
rand_seed     = 1234
max_data      = 100
burnin        = 0
num_samples   = 1
checkpoint    = 50000
dp_alpha      = 1e-1
dp_gamma      = 1e-1
init_drift    = 0.1
alpha_decay   = 0.5
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

## seed(rand_seed)

PyCloneBinomialData = namedtuple('PyCloneBinomialData',
                                 ['b', 'd', 'cn_n', 'cn_r', 'cn_v', 'mu_n', 'mu_r', 'mu_v', 'log_pi'])

def get_log_pi(weights):
    pi = [x / sum(weights) for x in weights]    
    return tuple([log(x) for x in pi]) 

def get_pyclone_data(mutation, error_rate):
    a = mutation.ref_counts
    b = mutation.var_counts
    
    d = a + b 
    
    cn_n = tuple([x.cn_n for x in mutation.states])
    cn_r = tuple([x.cn_r for x in mutation.states])
    cn_v = tuple([x.cn_v for x in mutation.states])
    
    mu_n = tuple([x.get_mu_n(error_rate) for x in mutation.states])
    mu_r = tuple([x.get_mu_r(error_rate) for x in mutation.states])
    mu_v = tuple([x.get_mu_v(error_rate) for x in mutation.states])
    
    prior_weights = tuple([x.prior_weight for x in mutation.states])

    log_pi = get_log_pi(prior_weights)

    return PyCloneBinomialData(b, d, cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi)

def load_sample_data(file_name, error_rate):
    '''
    Load data from PyClone formatted input file.
    '''
    data = OrderedDict()
    
    fh = open(file_name)
    
    config = yaml.load(fh)
    
    fh.close()

    for mutation_dict in config['mutations']:
        mutation = load_mutation_from_dict(mutation_dict)
        data[mutation.id] = get_pyclone_data(mutation, error_rate)
    
    return data


file_name = '/home/yuan03/Dropbox/dp/pyclone-data/SRR385941.yaml'

data = load_sample_data(file_name, error_rate)
data = data.values()
data = array(data, dtype=tuple)
data = delete(data,139,0)

dims = 1
root = PyClone_test( dims=dims, drift=init_drift )
tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
             root_node=root, data=data )

dp_alpha_traces    = zeros((num_samples, 1))
dp_gamma_traces    = zeros((num_samples, 1))
alpha_decay_traces = zeros((num_samples, 1))
drift_traces       = zeros((num_samples, dims))
cd_llh_traces      = zeros((num_samples, 1))
nodes_traces       = zeros((num_samples, 1))
tssb_traces        = empty((num_samples, 1),dtype = object)


intervals = zeros((7))
print "Starting MCMC run..."
## ipdb.set_trace()
for iter in range(-burnin,num_samples):

    times = [time.time()]

    tssb.resample_assignments()
    times.append(time.time())

    ipdb.set_trace()
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
        tssb_traces[iter]        = cPickle.dumps(tssb)  
        dp_alpha_traces[iter]    = tssb.dp_alpha
        dp_gamma_traces[iter]    = tssb.dp_gamma
        alpha_decay_traces[iter] = tssb.alpha_decay
        drift_traces[iter]       = root.drift()
        cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
        (weights, nodes)         = tssb.get_mixture()
        nodes_traces[iter]       = len(nodes)
        

    if iter > 0 and mod(iter, checkpoint) == 0:
        filename = "checkpoints/norm1d-test-%s-%06d.pkl" % (codename, iter)
        fh = open(filename, 'w')
        cPickle.dump(tssb, fh)
        fh.close()

  
    if mod(iter, 1) == 0:
        (weights, nodes) = tssb.get_mixture()
        print codename, iter, len(nodes), cd_llh_traces[iter], \
           mean(root._drift), tssb.dp_alpha, tssb.dp_gamma, \
           tssb.alpha_decay, \
           " ".join(map(lambda x: "%0.2f" % x, intervals.tolist())), \
           float(root.hmc_accepts)/(root.hmc_accepts+root.hmc_rejects), \
           root.hmc_accepts, root.hmc_rejects
        intervals = zeros((7))
        

    if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
        print "\t%f is best per-data complete data likelihood so far." \
           % (cd_llh_traces[iter]/max_data)
        best_fit = cPickle.dumps(tssb)
        ## filename_best = "bests/norm1d-test-%s-best.pkl" % (codename)
        ## fh = open(filename_best, 'w')
        ## cPickle.dump(tssb, fh)
        ## fh.close()



## filename = "checkpoints/norm1d-test-%s-final.pkl" % (codename)
## fh = open(filename, 'w')
## cPickle.dump({ 'tssb'               : tssb,
##                'dp_alpha_traces'    : dp_alpha_traces,
##                'dp_gamma_traces'    : dp_gamma_traces,
##                'alpha_decay_traces' : alpha_decay_traces,
##                'drift_traces'       : drift_traces,
##                'cd_llh_traces'      : cd_llh_traces
##                'nodes_traces'       : nodes_traces}, fh)
## fh.close()



nodes_tabular = itemfreq(nodes_traces)
best_num_nodes = nodes_tabular[argmax(nodes_tabular[:,1]),0]
best_num_nodes_llh = cd_llh_traces[nodes_traces==best_num_nodes].max()
best_node_fit = cPickle.loads(tssb_traces[cd_llh_traces==best_num_nodes_llh][0])


## (weights_best, nodes_best) = best_node_fit.get_mixture()
## pp = zeros(yy.shape) 

## for kk in range(len(weights_best)):
##     pp = pp + weights_best[kk]*1/sqrt(2*pi*nodes_best[kk].params[1]**2) * \
##       exp(-0.5*(yy-nodes_best[kk].params[0])**2/ nodes_best[kk].params[1]**2)

## fig1 = plt.figure(1)
## plt.plot(cd_llh_traces)
## plt.savefig('figures/test_trace_40.pdf',format='pdf')
## clf()

## fig2 = plt.figure(2)
## plt.plot(yy,pp, color = 'b')
## plt.plot(data,-0.1*ones(data.shape), linestyle = 'none',
##          color = 'g', marker = 'x')
## #plt.ylim((-0.2,0.8))
## plt.savefig('figures/test_40.pdf', format='pdf')
## clf()


## filename = 'testgraph_pyclone.gdl'
## fh2 = open(filename,'w')
## best_node_fit.print_graph_binomial(fh2)
## fh2.close()


cell_freq_traces = zeros((num_samples, data.shape[0]))

for ii in range(num_samples):
    tmp_tssb = cPickle.loads(tssb_traces[ii][0])
    assignments = tmp_tssb.assignments
    for jj in range(data.shape[0]):
        cell_freq_traces[ii,jj] = assignments[jj].params

## numpy.savetxt("pyclone_cell_freq.csv", cell_freq_traces, delimiter=",")
