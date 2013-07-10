import os
import sys
import time
import cPickle
import csv

from numpy         import *
from numpy.random  import *
from tssb          import *
from gaussian_test import *
from util          import *
from scipy.stats   import itemfreq



rand_seed     = 1234
max_data      = 100
burnin        = 5000
num_samples   = 10000
checkpoint    = 50000
dp_alpha      = 1.0
dp_gamma      = 1
init_drift    = 0.1
alpha_decay   = 0.25
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

## seed(rand_seed)
## data = concatenate([0.5 + 0.2*randn(70,1), -0.8 + 0.1*randn(35,1), -1.5+0.1*randn(35,1)])
## data = sigmoid(data)

reader = csv.DictReader(open('./pyclone-data/SRR385940.tsv'),
                        delimiter='\t')

data = []
dataid = []
for row in reader:
    dataid.append(row['mutation_id'])
    data.append(row['variant_freq'])
    
data = array(data,dtype = 'float64')
data = data[:,newaxis]

dims = 2
root = Gaussian_test( dims=dims, drift=init_drift )
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
for iter in range(-burnin,num_samples):

    times = [time.time()]

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
        tssb_traces[iter]        = cPickle.dumps(tssb)  
        dp_alpha_traces[iter]    = tssb.dp_alpha
        dp_gamma_traces[iter]    = tssb.dp_gamma
        alpha_decay_traces[iter] = tssb.alpha_decay
        drift_traces[iter]       = root.drift()
        cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
        (weights, nodes)         = tssb.get_mixture()
        nodes_traces[iter]       = len(nodes)
        

    ## if iter > 0 and mod(iter, checkpoint) == 0:
    ##     filename = "checkpoints/norm1d-test-%s-%06d.pkl" % (codename, iter)
    ##     fh = open(filename, 'w')
    ##     cPickle.dump(tssb, fh)
    ##     fh.close()

  
    if mod(iter, 100) == 0:
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


(weights_best, nodes_best) = best_node_fit.get_mixture()
yy = linspace(0,1,1000)
pp = zeros(yy.shape) 

for kk in range(len(weights_best)):
    pp = pp + weights_best[kk]*1/sqrt(2*pi*nodes_best[kk].params[1]**2) * \
      exp(-0.5*(yy-nodes_best[kk].params[0])**2/ nodes_best[kk].params[1]**2)

fig1 = plt.figure(1)
plt.plot(cd_llh_traces)
plt.savefig('figures/test_trace_40.pdf',format='pdf')
clf()

fig2 = plt.figure(2)
plt.plot(yy,pp, color = 'b')
plt.plot(data,-0.1*ones(data.shape), linestyle = 'none',
         color = 'g', marker = 'x')
#plt.ylim((-0.2,0.8))
plt.savefig('figures/test_40.pdf', format='pdf')
clf()


filename = 'testgraph_40.gdl'
fh2 = open(filename,'w')
best_node_fit.print_graph(fh2)
fh2.close()




