import os
import sys
import time
import cPickle
import csv
#import ipdb

from numpy         import *
from numpy.random  import *
from tssb          import *
from logistic_mixlaplace_methy_test   import *
from util          import *
from scipy.stats   import itemfreq
from sklearn       import metrics 

rand_seed     = 1234
#rand_seed     = 1264
max_data      = 100
burnin        = 0
num_samples   = 20000
checkpoint    = 50000
dp_alpha      = 1
dp_gamma      = 3e-1
init_drift    = 0.5
init_galpha   = 1
init_gbeta    = 0.5
alpha_decay   = 0.1
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

seed(rand_seed)

files = ['CT_IRX2P_R1.csv', 'CT_IRX2P_R4.csv', 'CT_IRX2P_R5.csv', 'CT_IRX2P_R6.csv']

sampno = 0
for seqsamp in files:
    reader = csv.DictReader(open('./data/Sottoriva/'+seqsamp),delimiter=',')

    data = []

    for row in reader:
        data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4']),int(row['V5']),int(row['V6']),int(row['V7']),int(row['V8'])])

    data = numpy.array(data)
    dims = data.shape[1]

    root = Logistic( dims=dims, mu = 5.0)
    tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
                root_node=root, data=data )

    dp_alpha_traces    = zeros((num_samples, 1))
    dp_gamma_traces    = zeros((num_samples, 1))
    alpha_decay_traces = zeros((num_samples, 1))
    cd_llh_traces      = zeros((num_samples, 1))
    nodes_traces       = zeros((num_samples, 1))
    #tssb_traces        = empty((num_samples, 1),dtype = object)
    bignodes_traces    = zeros((num_samples, 1))
    unnormpost_traces  = zeros((num_samples, 1))
    depth_traces       = zeros((num_samples, 1))
    base_value_traces  = zeros((num_samples, 1))
    std_traces         = zeros((num_samples, 1))
    root_bias_traces   = zeros((num_samples, 1))

    intervals = zeros((7))
    print "Starting MCMC run..." +seqsamp
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
   
        tssb.resample_tree_topology_root()

        tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
        times.append(time.time())

        intervals = intervals + diff(array(times))   
    
        if iter >= 0:
            #tssb_traces[iter]        = cPickle.dumps(tssb)  
            dp_alpha_traces[iter]    = tssb.dp_alpha
            dp_gamma_traces[iter]    = tssb.dp_gamma
            alpha_decay_traces[iter] = tssb.alpha_decay
            cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
            (weights, nodes)         = tssb.get_mixture()
            nodes_traces[iter]       = len(nodes)
            bignodes_traces[iter]    = sum(numpy.array(weights)>0.01)
            unnormpost_traces[iter]  = tssb.unnormalized_postertior()
            depth_traces[iter]       = tssb.deepest_node_depth()
            base_value_traces[iter]  = root.mu_caller()
            std_traces[iter]         = root.std_caller()
            root_bias_traces[iter]   = root.mu0_caller()
            
        if mod(iter, 1) == 0:
            (weights, nodes) = tssb.get_mixture()
            print codename, iter, len(nodes), cd_llh_traces[iter], \
              root.mu_caller(), root.std_caller(), \
              root.mu0_caller()+root.mu_caller() , 'big nodes:', int(bignodes_traces[iter]),\
              tssb.dp_alpha, tssb.dp_gamma, \
              tssb.alpha_decay
          
        intervals = zeros((7))

        if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
            print "\t%f is best per-data complete data likelihood so far." \
                % (cd_llh_traces[iter]/max_data)


    #nodes_tabular = itemfreq(nodes_traces)
    #tree_folder = './treescripts/%s-%s/' %(codename,seqsamp)

    #if not os.path.exists(tree_folder):
    #    os.makedirs(tree_folder)

    #for idx, nn in enumerate(nodes_tabular):
    #    node_num = nn[0]
    #    node_freq = nn[1]/num_samples 
    #    node_num_best_llh = cd_llh_traces[nodes_traces==node_num].max()
    #    node_fit = cPickle.loads(tssb_traces[cd_llh_traces==node_num_best_llh][0])
    #    filename = 'nodes-%i-freq-%0.4f.gdl' % (node_num, node_freq)
    #    fn2 = tree_folder + filename
    #    fh2 = open(fn2,'w')
    #    node_fit.print_graph_full_logistic_different_branch_length(fh2)
    #    fh2.close()

    traces = hstack([dp_alpha_traces, dp_gamma_traces,\
                     alpha_decay_traces, cd_llh_traces,\
                     nodes_traces, bignodes_traces, \
                     unnormpost_traces, depth_traces, \
                     base_value_traces, std_traces,\
                     root_bias_traces])
    tracefile = './mcmc-traces/%s_%s_%s_traces.csv' % (codename,seqsamp,str(rand_seed))
    numpy.savetxt(tracefile, traces, delimiter = ',',
                  header = "dp_alpha_traces,dp_gamma_traces,alpha_decay_traces,cd_llh_traces,node_traces,bignodes_traces,unnormpost_traces,depth_traces,base_value_traces,std_traces,root_bias_traces", comments='')

    sampno = sampno + 1

