import os
import time
import cPickle
import csv
import shutil


from numpy         import *
from numpy.random  import *
from scipy.stats   import itemfreq


from bitphylogeny.tssb          import *
from bitphylogeny.model_spec    import *
from bitphylogeny.util          import *
from bitphylogeny.post_process.compute_vmeasures import *


def run_wrapper(args):
    fin = args.fin
    fout = args.fout
    contains_true_label = args.f
    num_samples = args.n
    burnin = args.b
    thin = args.t
    run_analysis(fin, fout, contains_true_label, num_samples, burnin, thin)


def run_analysis(seqsamp , fout, contains_true_label, num_samples, burnin, thin):

    reader = csv.reader(open(seqsamp, 'rb'), delimiter=',')
    data = array( list(reader) )[1:,:].astype('int')

    if contains_true_label:
        data = data[:,:-1]

    dims = data.shape[1]
    max_data = data.shape[0]

    print max_data, dims
   
    seqsamp = seqsamp.split('/')[len(seqsamp.split('/'))-1]
    tree_folder =  fout + '/' + seqsamp  + '/treescripts/'

    if not os.path.exists(tree_folder):
        os.makedirs(tree_folder)
    else:
        shutil.rmtree(tree_folder)
        os.makedirs(tree_folder)


    trace_folder = fout + '/' + seqsamp + '/mcmc-traces/' 

    if not os.path.exists(trace_folder):
        os.makedirs(trace_folder)
    else:
        shutil.rmtree(trace_folder)
        os.makedirs(trace_folder)

    params_folder = trace_folder+'params_traces/' 
    if not os.path.exists(params_folder):
        os.makedirs(params_folder)
    else:
        shutil.rmtree(params_folder)
        os.makedirs(params_folder)

    rand_seed = int(round(rand(),4)*10000)
    seed(rand_seed)
    

    root = BitPhylogeny( dims=dims, mu = 5.0)
    tssb = TSSB( dp_alpha=2.0, dp_gamma=3e-1, alpha_decay=0.1,
            root_node=root, data=data )

    tree_collect_band  = 5 # has to be smaller than num_samples/thin
    max_depth          = 15
    cd_llh_traces      = zeros((num_samples/thin, 1))
    nodes_traces       = zeros((num_samples/thin, 1))
    tssb_traces        = empty((tree_collect_band, 1),dtype = object)
    bignodes_traces    = zeros((num_samples/thin, 1))
    unnormpost_traces  = zeros((num_samples/thin, 1))
    depth_traces       = zeros((num_samples/thin, 1))
    base_value_traces  = zeros((num_samples/thin, 1))
    std_traces         = zeros((num_samples/thin, 1))
    branch_traces      = zeros((num_samples/thin, data.shape[0]))
    width_dist         = zeros((num_samples/thin, max_depth))
    mass_dist          = zeros((num_samples/thin, max_depth))
    root_dist          = zeros((num_samples/thin, dims))
    label_traces       = zeros((num_samples/thin, data.shape[0]))
    params_traces      = zeros((num_samples/thin, data.shape[0], data.shape[1]))
    node_depth_traces  = zeros((num_samples/thin, data.shape[0]))

    intervals = zeros((7))
    print "Starting MCMC run..." +seqsamp
    for iter in range(-int(burnin),int(num_samples)):

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

        #tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
        times.append(time.time())

        intervals = intervals + diff(array(times))   

        idx = iter/thin

        if idx >= 0 and mod(iter, thin) == 0:
            tssb_traces[mod(idx,tree_collect_band)] = cPickle.dumps(tssb)
            cd_llh_traces[idx]      = tssb.complete_data_log_likelihood()
            (weights, nodes)         = tssb.get_mixture()
            nodes_traces[idx]       = len(nodes)
            bignodes_traces[idx]    = sum(numpy.array(weights)>0.01)
            unnormpost_traces[idx]  = tssb.unnormalized_postertior()
            depth_traces[idx]       = tssb.deepest_node_depth()
            base_value_traces[idx]  = root.mu_caller()
            std_traces[idx]         = root.std_caller()
            branch_traces[idx]      = array([x.branch_length for x in tssb.assignments]).T
            width_dist[idx]         = tssb.get_width_distribution()
            mass_dist[idx]          = tssb.get_weight_distribtuion()
            root_dist[idx]          = tssb.root['node'].params
            label_traces[idx]       = object2label(tssb.assignments, nodes)
            params_traces[idx]      = array([x.params for x in tssb.assignments])
            node_depth_traces[idx]  = array([x.depth for x in tssb.assignments])
        
            if mod(idx, 1) == 0:
                (weights, nodes) = tssb.get_mixture()
                print seqsamp, iter, len(nodes), cd_llh_traces[idx], \
                    root.mu_caller(), root.std_caller(), \
                    'big nodes:', int(bignodes_traces[idx]), \
                    diag(root.ratemat_caller()), mean(branch_traces[idx]),\
                    tssb.dp_alpha, tssb.dp_gamma, \
                    tssb.alpha_decay
      
            intervals = zeros((7))

            if idx > 0 and argmax(unnormpost_traces[:idx+1]) == idx:
                print "\t%f is best per-data complete data likelihood so far." \
                    % (unnormpost_traces[idx]/max_data)

            # print treescripts
            if idx + 1 == num_samples/thin:
                idx = idx + 1

            if mod(idx,tree_collect_band) == 0 and idx > 0:
                tmp_nodes_tabular = itemfreq(bignodes_traces[idx-tree_collect_band:idx])

                for ii, nn in enumerate(tmp_nodes_tabular):
                    node_num = nn[0]
                    node_num_best_llh = unnormpost_traces[bignodes_traces==node_num].max()
                    maxidx = where(unnormpost_traces==node_num_best_llh)[0][0]
                    if maxidx > idx - tree_collect_band or idx - tree_collect_band == 0:
                        node_fit = cPickle.loads(tssb_traces[mod(maxidx,tree_collect_band)][0])
                        filename = 'nodes-%i.gdl' % (node_num)
                        fn2 = tree_folder + filename
                        fh2 = open(fn2,'w')
                        node_fit.print_graph_full_logistic_different_branch_length(fh2)
                        fh2.close()

        nodes_tabular = itemfreq(bignodes_traces)
        nodes_tabular[:,1] = nodes_tabular[:,1] / (num_samples/thin)
        treefreqfile = 'tree-freq'
        header = ['unique_node_num', 'freq']
        write_traces2csv(tree_folder+'tree-freq.csv',nodes_tabular,header)

        traces = hstack([cd_llh_traces,\
                         nodes_traces, bignodes_traces, \
                         unnormpost_traces, depth_traces, \
                         base_value_traces, std_traces,\
                         width_dist, mass_dist, root_dist])

        tracefile = 'other_traces.csv'
        header = ['cd_llh_traces',
                  'node_traces','bignodes_traces',
                  'unnormpost_traces','depth_traces',
                  'base_value_traces',
                  'std_traces','w1','w2','w3',
                  'w4','w5','w6','w7','w8','w9','w10','w11','w12',
                  'w13','w14','w15','m1','m2','m3','m4','m5','m6',
                  'm7','m8','m9','m10','m11','m12','m13','m14','m15',
                  'r1','r2','r3','r4','r5','r6','r7','r8']

        write_traces2csv(trace_folder+'other_traces.csv',traces,header)
        write_traces2csv(trace_folder+'label_traces.csv',label_traces)
        write_traces2csv(trace_folder+'node_depth_traces.csv',node_depth_traces)
        write_traces2csv(trace_folder+'branch_traces.csv',branch_traces)

        numfiles = 10 # number of small arrays. Tested 50 files for (5e3,2e3,8) arrays. 
                      # If numfiles = 50, the function will generate 50 files in the 
                      # params_folder each contains a (1e2,2e3,8) array.
        write_params_traces2file(params_traces, numfiles, params_folder)

def post_process_wrapper(args):
    path = args.fout
    f    = args.fin
    write_vmeasures(path, f)
    

def write_vmeasures(path, f):
    traces = [x for x in os.listdir(path) if 'mcmc-traces' == x][0]
    fp = path+'/'+traces
    true_label = get_true_label( f )
    write_vmeasure_traces(fp, true_label)

    mpear_label = load_mpear_label(fp)
    from sklearn import metrics
    vmeasure = metrics.homogeneity_completeness_v_measure(mpear_label, true_label)
    vmeasure = numpy.array(vmeasure)[numpy.newaxis,:]
    header = ['homogeneity','completeness', 'v_measure']
    write_traces2csv(fp + '/'+'mpear_label_vmeasure.csv',vmeasure, header)
    
