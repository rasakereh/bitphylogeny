import os
import sys
import time
import cPickle
import csv

from numpy         import *
from numpy.random  import *
from tssb          import *
from logistic_mixlaplace_methy_fixed_ratemat_nodewise_branch_length import *
from util          import *
from scipy.stats   import itemfreq

#os.chdir("/home/sathomas/git/phylo-tree/code/")


rand_seed     = 1264
burnin        = 3e4
num_samples   = 5e4
dp_alpha      = 2.0
dp_gamma      = 3e-1
alpha_decay   = 0.1
max_depth     = 15

codename      = os.popen('./random-word').read().rstrip()

print "Codename: ", codename



rand_seed = int(round(rand(),4)*10000)


reader = csv.DictReader(open('./data/full_methy/noisy_small_clones_8_2000_genotype.csv'),
                        delimiter=',')
genotypes = []
for row in reader:
    genotypes.append([int(row['clone1']),int(row['clone2']),int(row['clone3']),int(row['clone4']),int(row['clone5']),int(row['clone6']),int(row['clone7'])])
genotypes = numpy.array(genotypes)


reader = csv.DictReader(open('./data/full_methy/noisy_small_clones_8_2000_0_mutmat.csv'),
                        delimiter=',')
data = []
labels = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),
                 int(row['V4']),int(row['V5']),int(row['V6']),
                 int(row['V7']),int(row['V8'])])
    labels.append(int(row['V9']))
    

data = numpy.array(data)
dims = data.shape[1]
max_data = data.shape[0]


seqsamp="ground_truth_small_clones_0.csv"

tree_folder = './treescripts/Sottoriva/small/%s-%i/%s/' \
  %(codename,rand_seed,seqsamp)

if not os.path.exists(tree_folder):
    os.makedirs(tree_folder)

trace_folder = './mcmc-traces/Sottoriva/small/%s-%i/%s/' \
  %(codename,rand_seed,seqsamp)
if not os.path.exists(trace_folder):
    os.makedirs(trace_folder)

params_folder = trace_folder+'params_traces/' 
if not os.path.exists(params_folder):
    os.makedirs(params_folder)
    

root = Logistic( dims=dims )
c1 = Logistic( parent=root,dims=dims )
root.add_child(c1)
c2 = Logistic( parent=c1, dims=dims )
c1.add_child(c2)
c3 = Logistic( parent=c1, dims=dims )
c1.add_child(c3)
c4 = Logistic( parent=c1, dims=dims )
c1.add_child(c4)
c5 = Logistic( parent=c2, dims=dims )
c2.add_child(c5)
c6 = Logistic( parent=c4, dims=dims )
c4.add_child(c6)
c7 = Logistic( parent=c5, dims=dims )
c5.add_child(c7)
c8 = Logistic( parent=c5, dims=dims )
c5.add_child(c8)
c9 = Logistic( parent=c6, dims=dims )
c6.add_child(c9)
c10 = Logistic( parent=c6, dims=dims )
c6.add_child(c10)
c11 = Logistic( parent=c6, dims=dims )
c6.add_child(c11)

tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, alpha_decay=alpha_decay,
            root_node=root, data=data )
c1.tssb = tssb
c2.tssb = tssb
c3.tssb = tssb
c4.tssb = tssb
c5.tssb = tssb
c6.tssb = tssb
c7.tssb = tssb
c8.tssb = tssb
c9.tssb = tssb
c10.tssb = tssb
c11.tssb = tssb

depth = 1
tssb.root['sticks'] = vstack([ tssb.root['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'].append({ 'node'     : c1,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
depth = 2
tssb.root['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'].append({ 'node'     : c2,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
tssb.root['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'].append({ 'node'     : c3,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
tssb.root['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'].append({ 'node'     : c4,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
depth = 3
tssb.root['children'][0]['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][0]['children'].append({ 'node'     : c5,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })

tssb.root['children'][0]['children'][2]['sticks'] = vstack([ tssb.root['children'][0]['children'][2]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][2]['children'].append({ 'node'     : c6,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
depth=4
tssb.root['children'][0]['children'][0]['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['children'][0]['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][0]['children'][0]['children'].append({ 'node'     : c7,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
tssb.root['children'][0]['children'][0]['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['children'][0]['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][0]['children'][0]['children'].append({ 'node'     : c8,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })

tssb.root['children'][0]['children'][2]['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['children'][2]['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][2]['children'][0]['children'].append({ 'node'     : c9,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
tssb.root['children'][0]['children'][2]['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['children'][2]['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][2]['children'][0]['children'].append({ 'node'     : c10,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })
tssb.root['children'][0]['children'][2]['children'][0]['sticks'] = vstack([ tssb.root['children'][0]['children'][2]['children'][0]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][2]['children'][0]['children'].append({ 'node'     : c11,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })


for n in range(tssb.num_data):
    tssb.assignments[n].remove_datum(n)
tssb.assignments = []
    
for i in range(labels.count(1)):
    root.add_datum(i)
    tssb.assignments.append(root)
for i in range(labels.count(2)):
    ii=i+labels.count(1)
    c1.add_datum(ii)
    tssb.assignments.append(c1)
for i in range(labels.count(3)):
    ii=i+labels.count(1)+labels.count(2)
    c2.add_datum(ii)
    tssb.assignments.append(c2)
for i in range(labels.count(4)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)
    c3.add_datum(ii)
    tssb.assignments.append(c3)
for i in range(labels.count(5)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)
    c4.add_datum(ii)
    tssb.assignments.append(c4)
for i in range(labels.count(6)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)
    c5.add_datum(ii)
    tssb.assignments.append(c5)
for i in range(labels.count(7)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)+labels.count(6)
    c6.add_datum(ii)
    tssb.assignments.append(c6)
for i in range(labels.count(8)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)+labels.count(6)+labels.count(7)
    c7.add_datum(ii)
    tssb.assignments.append(c7)
for i in range(labels.count(9)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)+labels.count(6)+labels.count(7)+labels.count(8)
    c8.add_datum(ii)
    tssb.assignments.append(c8)
for i in range(labels.count(10)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)+labels.count(6)+labels.count(7)+labels.count(8)+labels.count(9)
    c9.add_datum(ii)
    tssb.assignments.append(c9)
for i in range(labels.count(11)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)+labels.count(6)+labels.count(7)+labels.count(8)+labels.count(9)+labels.count(10)
    c10.add_datum(ii)
    tssb.assignments.append(c10)
for i in range(labels.count(12)):
    ii=i+labels.count(1)+labels.count(2)+labels.count(3)+labels.count(4)+labels.count(5)+labels.count(6)+labels.count(7)+labels.count(8)+labels.count(9)+labels.count(10)+labels.count(11)
    c11.add_datum(ii)
    tssb.assignments.append(c11)
    

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
branch_traces      = zeros((num_samples, data.shape[0]))
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

    #tssb.resample_assignments()
    times.append(time.time())

    tssb.cull_tree()
    times.append(time.time())

    tssb.resample_node_params()
    times.append(time.time())

    root.resample_hypers()
    times.append(time.time())

    tssb.resample_sticks()
    times.append(time.time())

    #tssb.resample_stick_orders()
    times.append(time.time())

    #tssb.resample_tree_topology_root()

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
        (weights, nodes)        = tssb.get_mixture()
        nodes_traces[idx]       = len(nodes)
        bignodes_traces[idx]    = sum(numpy.array(weights)>0.01)
        unnormpost_traces[idx]  = tssb.unnormalized_postertior()
        depth_traces[idx]       = tssb.deepest_node_depth()
        base_value_traces[idx]  = root.mu_caller()
        std_traces[idx]         = root.std_caller()
        #root_bias_traces[idx]   = root.mu0_caller()
        branch_traces[iter]     = array([x.branch_length for x in tssb.assignments]).T
        width_dist[idx]         = tssb.get_width_distribution()
        mass_dist[idx]          = tssb.get_weight_distribtuion()
        root_dist[idx]          = tssb.root['node'].params
        label_traces[idx]       = object2label(tssb.assignments, nodes)
        params_traces[idx]      = array([x.params for x in tssb.assignments])
        node_depth_traces[idx]  = array([x.depth for x in tssb.assignments])
        
    if mod(idx, 1) == 0:
        (weights, nodes) = tssb.get_mixture()
        print seqsamp, codename, iter, len(nodes), cd_llh_traces[idx], \
          'big nodes:', int(bignodes_traces[idx]), \
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
            #print maxidx, idx
            if maxidx > idx - tree_collect_band or idx - tree_collect_band == 0:
                node_fit = cPickle.loads(
                    tssb_traces[mod(maxidx,tree_collect_band)][0])
                filename = 'nodes-%i.gdl' % (node_num)
                fn2 = tree_folder + filename
                fh2 = open(fn2,'w')
                node_fit.print_graph_full_logistic_different_branch_length(fh2)
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
                 width_dist, \
                 mass_dist, root_dist])

tracefile = 'traces-%s' % (seqsamp)
header = ['dp_alpha_traces','dp_gamma_traces', 
          'alpha_decay_traces','cd_llh_traces',
          'node_traces','bignodes_traces',
          'unnormpost_traces','depth_traces',
          'base_value_traces',
          'std_traces','w1','w2','w3',
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

## Load saved arrays 
# params_traces = load_params_traces(params_folder)
