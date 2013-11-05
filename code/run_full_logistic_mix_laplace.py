import os
import sys
import time
import cPickle
import csv
import ipdb

from numpy         import *
from numpy.random  import *
from tssb         import *
from logistic_mixlaplace   import *
from util          import *
from scipy.stats   import itemfreq
import numpy

rand_seed     = 1234
max_data      = 100
burnin        = 0
num_samples   = 2000
checkpoint    = 50000
dp_alpha      = 1.0
dp_gamma      = 1.0
init_drift    = 0.5
init_galpha   = 1
init_gbeta    = 0.5
alpha_decay   = 0.1
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

seed(rand_seed)

reader = csv.DictReader(open('./data/fullsyn8_2000_mutmat.csv'),
                        delimiter=',')
data = []
dataid = []
max_snvpos = 1

for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4']),int(row['V5']),int(row['V6']),int(row['V7']),int(row['V8'])])
data = numpy.array(data)
dims = data.shape[1]

#freq = numpy.zeros([max_snvpos,2])
#for read in data:
#    freq[read[0]-1][0] += read[2]
#    freq[read[0]-1][1] += 1
#    freq[read[1]-1][0] += read[3]
#    freq[read[1]-1][1] += 1

#clonal = -1.5*ones(max_snvpos)
#for index, snv in enumerate(freq):
#    if snv[0]/snv[1] > 0.9:
#        clonal[index] = 1.5\


#filename = "test_beta_benomial_depth1.pkl"
#fh = open(filename, 'r')
#depth_traces = cPickle.load(fh)
#fh.close()

#median_depths = median(depth_traces,0)

root = Logistic( dims=dims, drift=init_drift, galpha=1.0, gbeta=0.1 )
tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
             root_node=root, data=data )

dp_alpha_traces    = zeros((num_samples, 1))
dp_gamma_traces    = zeros((num_samples, 1))
alpha_decay_traces = zeros((num_samples, 1))
galpha_traces      = zeros((num_samples, dims))
gbeta_traces       = zeros((num_samples, dims))
cd_llh_traces      = zeros((num_samples, 1))
nodes_traces       = zeros((num_samples, 1))
tssb_traces        = empty((num_samples, 1),dtype = object)

intervals = zeros((7))
print "Starting MCMC run..."
for iter in range(-burnin,num_samples):

    times = [ time.time() ]
    ##ipdb.set_trace()

    tssb.resample_assignments()
    times.append(time.time())

    tssb.cull_tree()
    times.append(time.time())

    ##ipdb.set_trace()
    tssb.resample_node_params()
    times.append(time.time())

    ##root.resample_hypers()
    times.append(time.time())

    tssb.resample_sticks()
    times.append(time.time())

    tssb.resample_stick_orders()
    times.append(time.time())

    if iter > 0:
    ##print "pre-treesamp:", tssb.unnormalized_postertior()
        tssb.resample_tree_topology_root()
    ##print "post-treesamp:", tssb.unnormalized_postertior()

    
    ##if iter > 0:
    tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
    times.append(time.time())
 

    intervals = intervals + diff(array(times))   
    
    if iter >= 0:
        tssb_traces[iter]        = cPickle.dumps(tssb)  
        dp_alpha_traces[iter]    = tssb.dp_alpha
        dp_gamma_traces[iter]    = tssb.dp_gamma
        alpha_decay_traces[iter] = tssb.alpha_decay
        #drift_traces[iter]       = root.drift()
        #galpha_traces[iter]      = root.galpha()
        #gbeta_traces[iter]       = root.gbeta()
        cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
        (weights, nodes)         = tssb.get_mixture()
        nodes_traces[iter]       = len(nodes)

    ## if iter > 0 and mod(iter, checkpoint) == 0:
    ##     filename = "checkpoints/norm1d-test-%s-%06d.pkl" % (codename, iter)
    ##     fh = open(filename, 'w')
    ##     cPickle.dump(tssb, fh)
    ##     fh.close()
  
    if mod(iter, 1) == 0:
        (weights, nodes) = tssb.get_mixture()
        print codename, iter, len(nodes), cd_llh_traces[iter], \
           tssb.dp_alpha, tssb.dp_gamma, \
           tssb.alpha_decay, \
           " ".join(map(lambda x: "%0.2f" % x, intervals.tolist())), \
           float(root.hmc_accepts)/(root.hmc_accepts+root.hmc_rejects), \
           root.hmc_accepts, root.hmc_rejects, time.strftime('%X %x %Z')
        intervals = zeros((7))

    ## if mod(iter, 1) == 0:
    ##     print "root:"
    ##     print " ".join(map(lambda x: "%0.2f" %x, root.params))
    ##     print "children:"
    ##     for child in root._children:
    ##         print " ".join(map(lambda x: "%0.2f" %x, child.params))
    
    if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
        print "\t%f is best per-data complete data likelihood so far." \
           % (cd_llh_traces[iter]/max_data)
        best_fit = cPickle.dumps(tssb)
        filename_best = "bests/pair-test-%s-best.pkl" % (codename)
        fh = open(filename_best, 'w')
        cPickle.dump(tssb, fh)
        fh.close()



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

## filename = 'mcmc-traces/tssb_traces_%s' % (codename)
## fh = open(filename,'w')
## cPickle.dump(tssb_traces, fh)
## fh.close()


nodes_tabular = itemfreq(nodes_traces)
best_num_nodes = nodes_tabular[argmax(nodes_tabular[:,1]),0]
best_num_nodes_llh = cd_llh_traces[nodes_traces==best_num_nodes].max()
best_node_fit = cPickle.loads(tssb_traces[cd_llh_traces==best_num_nodes_llh][0])

## (weights_best, nodes_best) = best_node_fit.get_mixture()
## yy = linspace(0,1,1000)
## pp = zeros(yy.shape) 

## for kk in range(len(weights_best)):
##     pp = pp + weights_best[kk]*1/sqrt(2*pi*nodes_best[kk].params[1]**2) * \
##       exp(-0.5*(yy-nodes_best[kk].params[0])**2/ nodes_best[kk].params[1]**2)

## fig1 = plt.figure(1)
## plt.plot(cd_llh_traces)
## ## plt.savefig('figures/PairLogistic1_trace_39.pdf',format='pdf')
## ## clf()

## filename_best = "bests/PairLogistic1_test.pkl" 
## fh = open(filename_best, 'w')
## cPickle.dump(best_node_fit, fh)
## fh.close()

## fn = "bests/PairLogistic1_test.pkl"
## fh = open(fn, 'r')
## best_node_fit = cPickle.load(fh)
## fh.close()

## fig2 = plt.figure(2)
## plt.plot(yy,pp, color = 'b')
## plt.plot(data,-0.1*ones(data.shape), linestyle = 'none',
##          color = 'g', marker = 'x')
## #plt.ylim((-0.2,0.8))
## plt.savefig('figures/PairLogistic1_39.pdf', format='pdf')
## clf()

## subplot(1,3,1)
## plot(nodes_best[0].params)
## subplot(1,3,2)
## plot(nodes_best[1].params)
## subplot(1,3,3)
## plot(nodes_best[2].params)


#best = loads(best_fit)
best_node_fit.remove_empty_nodes()

filename = './treescripts/test_tree_full_logistic_mix_laplace.gdl'
fh2 = open(filename,'w')
best_node_fit.print_graph_full_logistic(fh2)
fh2.close()
        
## filename = 'treescripts/testgraph_PairLogistic1.gdl'
## fh2 = open(filename,'w')
## best_node_fit.print_graph(fh2)
## fh2.close()
traces = hstack([dp_alpha_traces, dp_gamma_traces, \
                alpha_decay_traces, cd_llh_traces, \
                nodes_traces])
numpy.savetxt('./mcmc-traces/test_full_logistic_mix_laplacetraces.csv', traces, delimiter = ',',
              header = "dp_alpha_traces,dp_gamma_traces,alpha_decay_traces,cd_llh_traces,node_traces", comments='')
