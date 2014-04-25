import os
import sys
import time
import cPickle
import csv
import ipdb

from numpy         import *
from numpy.random  import *
from tssb          import *
from final_test_methy_averaged_ratemat import *
from util          import *
from scipy.stats   import itemfreq
from sklearn       import metrics 

rand_seed     = int(round(rand(),4)*10000)
rand_seed     = 1234
burnin        = 0
num_samples   = 10000
checkpoint    = 50000
dp_alpha      = 2.0
dp_gamma      = 0.3
alpha_decay   = 0.1
max_depth     = 15

codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

codename = 'final-test-fix-root-hyper-1'
seed(rand_seed)


files = ['noisy_small_clones_8_200_0_mutmat.csv']

trace_folder = './mcmc-traces/full-methy/%s-%i/' %(codename,rand_seed)
if not os.path.exists(trace_folder):
    os.makedirs(trace_folder)

params_folder = trace_folder+'params_traces/' 
if not os.path.exists(params_folder):
    os.makedirs(params_folder)


for seqsamp in files:

    tree_folder = './treescripts/full-methy/%s-%i/%s/' \
                  % (codename,rand_seed,seqsamp)
    if not os.path.exists(tree_folder):
        os.makedirs(tree_folder)
    
    reader = csv.DictReader(open('./data/full_methy/small-clone/'+seqsamp),
                            delimiter=',')
    data = []
    clone = []

    for row in reader:
        data.append([int(row['V1']),int(row['V2']),int(row['V3']),
                     int(row['V4']),int(row['V5']),int(row['V6']),
                     int(row['V7']),int(row['V8'])])
        clone.append(int(row['V9']))

    data = numpy.array(data)
    dims = data.shape[1]
    max_data = data.shape[0]
    root = Logistic(dims = dims)
    tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
                root_node=root, data=data )

    tree_collect_band  = 100
    #dp_alpha_traces    = zeros((num_samples, 1))
    #dp_gamma_traces    = zeros((num_samples, 1))
    #alpha_decay_traces = zeros((num_samples, 1))
    cd_llh_traces      = zeros((num_samples, 1))
    nodes_traces       = zeros((num_samples, 1))
    tssb_traces        = empty((tree_collect_band, 1),dtype = object)
    v_measure_traces   = zeros((num_samples, 3))
    bignodes_traces    = zeros((num_samples, 1))
    unnormpost_traces  = zeros((num_samples, 1))
    depth_traces       = zeros((num_samples, 1))
    base_value_traces  = zeros((num_samples, 1))
    std_traces         = zeros((num_samples, 1))
    root_bias_traces   = zeros((num_samples, 1))
    branch_traces      = zeros((num_samples, data.shape[0]))
    width_dist         = zeros((num_samples, max_depth))
    mass_dist          = zeros((num_samples, max_depth))
    root_dist          = zeros((num_samples, dims))
    label_traces       = zeros((num_samples, data.shape[0]))
    params_traces      = zeros((num_samples, data.shape[0], data.shape[1]))
    node_depth_traces  = zeros((num_samples, data.shape[0]))

    intervals = zeros((7))
    print "Starting MCMC run..." +seqsamp
    for iter in range(-burnin,num_samples):

        times = [ time.time() ]
   

        tssb.resample_assignments()
        times.append(time.time())

        #tssb.resample_tree_topology_root_1()

        tssb.cull_tree()
        times.append(time.time())

        #ipdb.set_trace()
        tssb.resample_node_params()
        times.append(time.time())

        root.resample_hypers()
        times.append(time.time())

        tssb.resample_sticks()
        times.append(time.time())

        tssb.resample_stick_orders()
        times.append(time.time())
   
        tssb.resample_tree_topology_root()

        #tssb.resample_hypers(dp_alpha=False, alpha_decay=True, dp_gamma=True)
        times.append(time.time())

        intervals = intervals + diff(array(times))   
    
        if iter >= 0:            
            #dp_alpha_traces[iter]    = tssb.dp_alpha
            #dp_gamma_traces[iter]    = tssb.dp_gamma
            #alpha_decay_traces[iter] = tssb.alpha_decay
            cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
            (weights, nodes)         = tssb.get_mixture()
            nodes_traces[iter]       = len(nodes)
            bignodes_traces[iter]    = sum(numpy.array(weights)>0.01)
            unnormpost_traces[iter]  = tssb.unnormalized_postertior()
            depth_traces[iter]       = tssb.deepest_node_depth()
            base_value_traces[iter]  = root.mu_caller()
            std_traces[iter]         = root.std_caller()
            branch_traces[iter]      = array([x.branch_length for x in tssb.assignments]).T
            width_dist[iter]         = tssb.get_width_distribution()
            mass_dist[iter]          = tssb.get_weight_distribtuion()
            root_dist[iter]          = tssb.root['node'].params
            data_assign              = object2label(tssb.assignments, nodes)
            label_traces[iter]        = data_assign
            params_traces[iter]       = array([x.params for x in tssb.assignments])
            node_depth_traces[iter]   = array([x.depth for x in tssb.assignments])
            v_measure_traces[iter]   = \
              metrics.homogeneity_completeness_v_measure(clone,data_assign)
            
        if mod(iter, 10) == 0:
            (weights, nodes) = tssb.get_mixture()
            print seqsamp, codename, iter, len(nodes), cd_llh_traces[iter], \
              v_measure_traces[iter],\
              'big nodes:', int(bignodes_traces[iter]), root.mu_caller(),\
              root.std_caller(), diag(root.ratemat_caller()), mean(branch_traces[iter])
          
        intervals = zeros((7))

        if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
            print "\t%f is best per-data complete data likelihood so far." \
                % (cd_llh_traces[iter]/max_data)

        tssb_traces[mod(iter,tree_collect_band)] = cPickle.dumps(tssb)

        # print treescripts
        if iter + 1 == num_samples:
            iter = iter + 1

        if mod(iter,tree_collect_band) == 0 and iter > 0:
            tmp_nodes_tabular = itemfreq(bignodes_traces[iter-tree_collect_band:iter])

            for idx, nn in enumerate(tmp_nodes_tabular):
                node_num = nn[0]
                node_num_best_llh = \
                  unnormpost_traces[bignodes_traces==node_num].max()
                maxidx = where(unnormpost_traces==node_num_best_llh)[0][0]
                if maxidx > iter - tree_collect_band or iter - tree_collect_band == 0:
                    node_fit = cPickle.loads(
                        tssb_traces[mod(maxidx,tree_collect_band)][0])
                    filename = 'nodes-%i.gdl' % (node_num)
                    fn2 = tree_folder + filename
                    fh2 = open(fn2,'w')
                    node_fit.print_graph_full_logistic_different_branch_length(fh2)
                    fh2.close()


    nodes_tabular = itemfreq(bignodes_traces)
    nodes_tabular[:,1] = nodes_tabular[:,1] / num_samples
    treefreqfile = 'tree-freq'
    header = ['unique_big_node_num', 'freq']
    fh4 = open(tree_folder+treefreqfile, 'wb')
    writer = csv.writer(fh4)
    writer.writerow(header)
    [ writer.writerow(x) for x in nodes_tabular ]
    fh4.close()
    

    traces = hstack([cd_llh_traces,\
                     nodes_traces, bignodes_traces, \
                     unnormpost_traces, depth_traces, \
                     base_value_traces, std_traces,\
                     v_measure_traces,\
                     width_dist, mass_dist, root_dist])
                     
    tracefile = 'traces-%s' % (seqsamp)
    header = ['cd_llh_traces',
              'node_traces','bignodes_traces',
              'unnormpost_traces','depth_traces',
              'base_value_traces',
              'std_traces',
              'homogeneity','completeness',
              'v-measure','w1','w2','w3',
              'w4','w5','w6','w7','w8','w9','w10','w11','w12',
              'w13','w14','w15','m1','m2','m3','m4','m5','m6',
              'm7','m8','m9','m10','m11','m12','m13','m14','m15',
              'r1','r2','r3','r4','r5','r6','r7','r8']

    fh3 = open(trace_folder+tracefile, 'wb')
    writer = csv.writer(fh3)
    writer.writerow(header)
    [ writer.writerow(x) for x in traces ]
    fh3.close()

    write_traces2csv(trace_folder+tracefile+'_label_traces.csv',label_traces)
    write_traces2csv(trace_folder+tracefile+'_node_depth_traces.csv',node_depth_traces)
    write_traces2csv(trace_folder+tracefile+'_branch_traces.csv',branch_traces)

    numfiles = 10 # number of small arrays. Tested 50 files for (5e3,2e3,8) arrays. 
             # If numfiles = 50, the function will generate 50 files in the 
             # params_folder each contains a (1e2,2e3,8) array.
    write_params_traces2file(params_traces, numfiles, params_folder)

 

