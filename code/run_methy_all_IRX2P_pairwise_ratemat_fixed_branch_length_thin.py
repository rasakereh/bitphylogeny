import os
import sys
import time
import cPickle
import csv

from numpy         import *
from numpy.random  import *
from tssb          import *
from logistic_mixlaplace_methy_pairwise_ratemat_fixed_branch_length\
   import *
from util          import *
from scipy.stats   import itemfreq

os.chdir("/home/sathomas/git/phylo-tree/code/")


rand_seed     = 1264
burnin        = 3e4
num_samples   = 5e4
dp_alpha      = 2.0
dp_gamma      = 3e-1
alpha_decay   = 0.1
max_depth     = 15

codename      = os.popen('./random-word').read().rstrip()
#codename      = 'test_thin'

print "Codename: ", codename



rand_seed = int(round(rand(),4)*10000)

#files = [['CT_IRX2P_R1.csv', 'CT_IRX2P_R4.csv',
#          'CT_IRX2P_R5.csv', 'CT_IRX2P_R6.csv'],
#         ['CT_IRX2P_L2.csv', 'CT_IRX2P_L3.csv',
#          'CT_IRX2P_L7.csv', 'CT_IRX2P_L8.csv'],
#         ['CU_IRX2P_R1.csv', 'CU_IRX2P_R2.csv',
#          'CU_IRX2P_R3.csv', 'CU_IRX2P_R5.csv'],
#         ['CU_IRX2P_L4.csv', 'CU_IRX2P_L6.csv',
#          'CU_IRX2P_L7.csv', 'CU_IRX2P_L8.csv'],
#         ['CX_IRX2P_R1.csv', 'CX_IRX2P_R2.csv',
#          'CX_IRX2P_R6.csv'],
#         ['CX_IRX2P_L3.csv', 'CX_IRX2P_L4.csv',
#          'CX_IRX2P_L5.csv'],
#         ['HA_IRX2P_R5.csv', 'HA_IRX2P_R7.csv',
#          'HA_IRX2P_R8.csv', 'HA_IRX2P_R9.csv',
#         'HA_IRX2P_R10.csv'],
#         ['HA_IRX2P_L1.csv', 'HA_IRX2P_L2.csv',
#          'HA_IRX2P_L3.csv', 'HA_IRX2P_L4.csv',
#         'HA_IRX2P_L6.csv']]

files = ['CT_IRX2P_R1.csv', 'CT_IRX2P_R4.csv',
         'CT_IRX2P_R5.csv', 'CT_IRX2P_R6.csv',
         'CT_IRX2P_L2.csv', 'CT_IRX2P_L3.csv',
         'CT_IRX2P_L7.csv', 'CT_IRX2P_L8.csv',
         'CU_IRX2P_R1.csv', 'CU_IRX2P_R2.csv',
         'CU_IRX2P_R3.csv', 'CU_IRX2P_R5.csv',
         'CU_IRX2P_L4.csv', 'CU_IRX2P_L6.csv',
         'CU_IRX2P_L7.csv', 'CU_IRX2P_L8.csv',
         'CX_IRX2P_R1.csv', 'CX_IRX2P_R2.csv',
         'CX_IRX2P_R6.csv',
         'CX_IRX2P_L3.csv', 'CX_IRX2P_L4.csv',
         'CX_IRX2P_L5.csv',
         'HA_IRX2P_R5.csv', 'HA_IRX2P_R7.csv',
         'HA_IRX2P_R8.csv', 'HA_IRX2P_R9.csv',
         'HA_IRX2P_R10.csv',
         'HA_IRX2P_L1.csv', 'HA_IRX2P_L2.csv',
         'HA_IRX2P_L3.csv', 'HA_IRX2P_L4.csv',
         'HA_IRX2P_L6.csv',
         'noisy_full_methy_8_2000_0_mutmat.csv',
         'noisy_full_methy_8_2000_0.01_mutmat.csv',
         'noisy_full_methy_8_2000_0.02_mutmat.csv',
         'noisy_full_methy_8_2000_0.05_mutmat.csv']

#files = [['noisy_full_methy_8_2000_0_mutmat.csv',
#         'noisy_full_methy_8_2000_0.01_mutmat.csv',
#         'noisy_full_methy_8_2000_0.02_mutmat.csv',
#         'noisy_full_methy_8_2000_0.05_mutmat.csv']]

f = int(sys.argv[1])
#f = 0
if isnan(f) or f<0 or f>35:
    exit()

seqsamp = files[f]

tree_folder = './treescripts/Sottoriva/pairwise/%s-%i/%s/' \
  %(codename,rand_seed,seqsamp)

if not os.path.exists(tree_folder):
    os.makedirs(tree_folder)

trace_folder = './mcmc-traces/Sottoriva/pairwise/%s-%i/%s/' \
  %(codename,rand_seed,seqsamp)
if not os.path.exists(trace_folder):
    os.makedirs(trace_folder)
    
reader = csv.DictReader(open('./data/Sottoriva/'+seqsamp),delimiter=',')
#reader = csv.DictReader(open('./data/full_methy/'+seqsamp),delimiter=',')

data = []

for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),
                 int(row['V4']),int(row['V5']),int(row['V6']),
                 int(row['V7']),int(row['V8'])])

data = numpy.array(data)
dims = data.shape[1]
max_data = data.shape[0]

root = Logistic( dims=dims, mu = 5.0)
tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
            root_node=root, data=data )

tree_collect_band  = 5.0
thin               = 5.0
dp_alpha_traces    = zeros((num_samples/thin, 1))
dp_gamma_traces    = zeros((num_samples/thin, 1))
alpha_decay_traces = zeros((num_samples/thin, 1))
cd_llh_traces      = zeros((num_samples/thin, 1))
nodes_traces       = zeros((num_samples/thin, 1))
tssb_traces        = empty((tree_collect_band, 1),dtype = object)
bignodes_traces    = zeros((num_samples/thin, 1))
unnormpost_traces  = zeros((num_samples/thin, 1))
depth_traces       = zeros((num_samples/thin, 1))
base_value_traces  = zeros((num_samples/thin, 1))
std_traces         = zeros((num_samples/thin, 1))
root_bias_traces   = zeros((num_samples/thin, 1))
branch_traces      = zeros((num_samples/thin, 1))
width_dist         = zeros((num_samples/thin, max_depth))
mass_dist          = zeros((num_samples/thin, max_depth))
root_dist          = zeros((num_samples/thin, dims))
label_traces       = zeros((num_samples/thin, data.shape[0]))
params_traces      = zeros((num_samples/thin, data.shape[0], data.shape[1]))
node_depth_traces  = zeros((num_samples/thin, data.shape[0]))

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

    #tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
    times.append(time.time())

    intervals = intervals + diff(array(times))   

    idx = iter/thin

    if idx >= 0 and mod(iter, thin) == 0:
        tssb_traces[mod(idx,tree_collect_band)] = cPickle.dumps(tssb)
        dp_alpha_traces[idx]    = tssb.dp_alpha
        dp_gamma_traces[idx]    = tssb.dp_gamma
        alpha_decay_traces[idx] = tssb.alpha_decay
        cd_llh_traces[idx]      = tssb.complete_data_log_likelihood()
        (weights, nodes)         = tssb.get_mixture()
        nodes_traces[idx]       = len(nodes)
        bignodes_traces[idx]    = sum(numpy.array(weights)>0.01)
        unnormpost_traces[idx]  = tssb.unnormalized_postertior()
        depth_traces[idx]       = tssb.deepest_node_depth()
        base_value_traces[idx]  = root.mu_caller()
        std_traces[idx]         = root.std_caller()
        root_bias_traces[idx]   = root.mu0_caller()
        branch_traces[idx]      = root.branch_caller()
        width_dist[idx]         = tssb.get_width_distribution()
        mass_dist[idx]          = tssb.get_weight_distribtuion()
        root_dist[idx]          = tssb.root['node'].params
        label_traces[idx]       = object2label(tssb.assignments, nodes)
        params_traces[idx]      = array([x.params for x in tssb.assignments])
        node_depth_traces[idx]  = array([x.depth for x in tssb.assignments])
        
    if mod(idx, 1) == 0:
        print seqsamp, codename, iter, len(nodes), cd_llh_traces[idx], \
          root.mu_caller(), root.std_caller(), \
          root.mu0_caller()+root.mu_caller() , \
          'big nodes:', int(bignodes_traces[idx]), root.branch_caller(),\
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
            print maxidx, idx
            if maxidx > idx - tree_collect_band or idx - tree_collect_band == 0:
                node_fit = cPickle.loads(
                    tssb_traces[mod(maxidx,tree_collect_band)][0])
                filename = 'nodes-%i.gdl' % (node_num)
                fn2 = tree_folder + filename
                fh2 = open(fn2,'w')
                node_fit.print_graph_full_logistic(fh2)
                fh2.close()

nodes_tabular = itemfreq(bignodes_traces)
nodes_tabular[:,1] = nodes_tabular[:,1] / (num_samples/thin)
treefreqfile = 'tree-freq'
header = ['unique_node_num', 'freq']
fh4 = open(tree_folder+treefreqfile, 'wb')
writer = csv.writer(fh4)
writer.writerow(header)
[ writer.writerow(x) for x in nodes_tabular ]
fh4.close()

traces = hstack([dp_alpha_traces, dp_gamma_traces,\
                 alpha_decay_traces, cd_llh_traces,\
                 nodes_traces, bignodes_traces, \
                 unnormpost_traces, depth_traces, \
                 base_value_traces, std_traces,\
                 root_bias_traces, branch_traces, width_dist, \
                 mass_dist, root_dist])

tracefile = 'traces-%s' % (seqsamp)
header = ['dp_alpha_traces','dp_gamma_traces', 
          'alpha_decay_traces','cd_llh_traces',
          'node_traces','bignodes_traces',
          'unnormpost_traces','depth_traces',
          'base_value_traces',
          'std_traces','root_bias_traces','branch_traces','w1','w2','w3',
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

fh = open(trace_folder+tracefile+'_params.pkl','w')
cPickle.dump(params_traces,fh)
fh.close()
