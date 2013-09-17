import os
import sys
import time
import cPickle
import csv
import pdb
import yaml
##import ipdb


## for python 2.6
#from collections import namedtuple
#from ordereddict import OrderedDict

## for python 2.7
from collections import namedtuple, OrderedDict



from numpy         import *
from numpy.random  import *
from tssb_ke       import *
from beta_binomial import *
from util          import *
from scipy.stats   import itemfreq
from config        import *
    
error_rate    = 0.001
rand_seed     = 1234
max_data      = 100
burnin        = 0
num_samples   = 20000
checkpoint    = 50000
dp_alpha      = 1.0
dp_gamma      = 5e-3
init_drift    = 0.1
alpha_decay   = 0.1
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename
seed(rand_seed)

file_name = './data/syn_reads.csv'
reader = csv.DictReader(open(file_name),
                        delimiter=',')

var_counts = []
tot_counts = []
data = []
for row in reader:
    var_counts.append(row['V2'])
    tot_counts.append(row['V3'])

var_counts = array(var_counts, dtype = 'float64')
tot_counts = array(tot_counts, dtype = 'float64')
var_counts = var_counts[:,newaxis]
tot_counts = tot_counts[:,newaxis]
error_rate = error_rate * ones((tot_counts.shape))
data = hstack([var_counts,tot_counts, error_rate])

dims = 1
root = Beta_Binomial( dims=dims )
tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
             root_node=root, data=data )

dp_alpha_traces    = zeros((num_samples, 1))
dp_gamma_traces    = zeros((num_samples, 1))
alpha_decay_traces = zeros((num_samples, 1))
drift_traces       = zeros((num_samples, dims))
cd_llh_traces      = zeros((num_samples, 1))
cd_mllh_traces     = zeros((num_samples, 1))
nodes_traces       = zeros((num_samples, 1))
tssb_traces        = empty((num_samples, 1),dtype = object)
balpha_traces      = zeros((num_samples, dims))
bbeta_traces       = zeros((num_samples, dims))
preci_traces       = zeros((num_samples, dims))

intervals = zeros((7))
print "Starting MCMC run..."
## ipdb.set_trace()

start_time = time.time()
for iter in range(-burnin,num_samples):

    times = [time.time()]

    tssb.resample_assignments()
    times.append(time.time())

    
    tssb.cull_tree()
    times.append(time.time())

    ##pdb.set_trace()
    tssb.resample_node_params()
    times.append(time.time())

    ##ipdb.set_trace()
    root.resample_hypers()
    times.append(time.time())

    tssb.resample_sticks()
    times.append(time.time())

    tssb.resample_stick_orders()
    times.append(time.time())

    tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
    times.append(time.time())
 
    intervals = intervals + diff(array(times))

    if iter >= 0:
        tssb_traces[iter]        = cPickle.dumps(tssb)  
        dp_alpha_traces[iter]    = tssb.dp_alpha
        dp_gamma_traces[iter]    = tssb.dp_gamma
        alpha_decay_traces[iter] = tssb.alpha_decay
        balpha_traces[iter]      = root.balpha()
        bbeta_traces[iter]       = root.bbeta()
        cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
        ##cd_mllh_traces[iter]     = tssb.complete_data_log_marginal_likelihood()
        (weights, nodes)         = tssb.get_mixture()
        preci_traces[iter]       = root.preci()
        nodes_traces[iter]       = len(nodes)
        

    if iter > 0 and mod(iter, checkpoint) == 0:
        filename = "checkpoints/norm1d-test-%s-%06d.pkl" % (codename, iter)
        fh = open(filename, 'w')
        cPickle.dump(tssb, fh)
        fh.close()

  
    if mod(iter, 100) == 0:
        (weights, nodes) = tssb.get_mixture()
        print codename, iter, len(nodes), cd_llh_traces[iter], \
           tssb.dp_alpha, tssb.dp_gamma, \
           tssb.alpha_decay, root.preci(),\
           " ".join(map(lambda x: "%0.2f" % x, intervals.tolist()))
           
        intervals = zeros((7))
        

    if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
        print "\t%f is best per-data complete data likelihood so far." \
           % (cd_llh_traces[iter]/max_data)
        best_fit = cPickle.dumps(tssb)
        ## filename_best = "bests/norm1d-test-%s-best.pkl" % (codename)
        ## fh = open(filename_best, 'w')
        ## cPickle.dump(tssb, fh)
        ## fh.close()


elapsed_time = time.time() - start_time
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

best_node_fit.remove_empty_nodes()

filename = 'test_tree_beta_binomial_0.gdl'
fh2 = open(filename,'w')
best_node_fit.print_graph_binomial(fh2)
fh2.close()

best = loads(best_fit)
best.remove_empty_nodes()

filename = 'test_tree_beta_binomial_1.gdl'
fh2 = open(filename,'w')
best.print_graph_binomial(fh2)
fh2.close()


cell_freq_traces = zeros((num_samples, data.shape[0]))

for ii in range(num_samples):
    tmp_tssb = cPickle.loads(tssb_traces[ii][0])
    assignments = tmp_tssb.assignments
    for jj in range(data.shape[0]):
        cell_freq_traces[ii,jj] = assignments[jj].params



traces = hstack([dp_alpha_traces, dp_gamma_traces, \
                alpha_decay_traces, cd_llh_traces, \
                bbeta_traces, preci_traces]) 
numpy.savetxt("test_beta_binomial_cell_freq.csv", cell_freq_traces, delimiter=",")
numpy.savetxt('test_beta_binomial_traces.csv', traces, delimiter = ',',
              header = "dp_alpha_traces,dp_gamma_traces,alpha_decay_traces,cd_llh_traces,bbeta_traces,prei_traces", comments='')
