import os
import sys
import time
import cPickle
import csv
import ipdb

from numpy         import *
from numpy.random  import *
from tssb          import *
from logistic_mixlaplace_methy_test   import *
from util          import *
from scipy.stats   import itemfreq
from sklearn       import metrics 

#rand_seed     = 1234
rand_seed     = 1264
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

reader = csv.DictReader(open('./data/Sottoriva/CT_IRX2P_R1.csv'),
                        delimiter=',')
data = []
dataid = []
max_snvpos = 1

for row in reader:
    ## data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4']),int(row['V5']),int(row['V6']),int(row['V7']),int(row['V8']),int(row['V9']),int(row['V10']),int(row['V11']),int(row['V12']),int(row['V13']),int(row['V14']),int(row['V15']),int(row['V16'])])
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4']),int(row['V5']),int(row['V6']),int(row['V7']),int(row['V8'])])
data = numpy.array(data)

#data = data[2:,:] # remove the first 2 reads

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

root = Logistic( dims=dims, mu = 5.0,
                 ratemat= array([[-1.0/8.0,1.0/8.0],[1.0/8.0,-1.0/8.0]]))
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
silhouette_traces  = zeros((num_samples, 1))
bignodes_traces    = zeros((num_samples, 1))

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

    root.resample_hypers()
    times.append(time.time())

    tssb.resample_sticks()
    times.append(time.time())

    tssb.resample_stick_orders()
    times.append(time.time())
   
    tssb.resample_tree_topology_root()

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
        bignodes_traces[iter]    = sum(numpy.array(weights)>0.01)
        ##data_assign = array(array(tssb.assignments),dtype='str')
        ## silhouette_traces[iter] = metrics.silhouette_score(data,
        ##                                                    data_assign,
        ##                                                    metric='hamming')

    ## if iter > 0 and mod(iter, checkpoint) == 0:
    ##     filename = "checkpoints/norm1d-test-%s-%06d.pkl" % (codename, iter)
    ##     fh = open(filename, 'w')
    ##     cPickle.dump(tssb, fh)
    ##     fh.close()
  
    if mod(iter, 1) == 0:
        (weights, nodes) = tssb.get_mixture()
        print codename, iter, len(nodes), cd_llh_traces[iter], \
           root.base_value, root.std, \
           root.root_bias+root.base_value , 'big nodes:', int(bignodes_traces[iter]),\
           tssb.dp_alpha, tssb.dp_gamma, \
           tssb.alpha_decay
          
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

nodes_tabular = itemfreq(nodes_traces)
tree_folder = './treescripts/%s-CT_IRX2P_R1/' %(codename)

if not os.path.exists(tree_folder):
    os.makedirs(tree_folder)


for idx, nn in enumerate(nodes_tabular):
    node_num = nn[0]
    node_freq = nn[1]/num_samples 
    node_num_best_llh = cd_llh_traces[nodes_traces==node_num].max()
    node_fit = cPickle.loads(tssb_traces[cd_llh_traces==node_num_best_llh][0])
    filename = 'nodes-%i-freq-%0.2f.gdl' % (node_num, node_freq)
    fn2 = tree_folder + filename
    fh2 = open(fn2,'w')
    node_fit.print_graph_full_logistic(fh2)
    fh2.close()

traces = hstack([dp_alpha_traces, dp_gamma_traces, \
                alpha_decay_traces, cd_llh_traces, \
                nodes_traces])
tracefile = './mcmc-traces/%s_CT_IRX2P_R1_traces.csv' % (codename)
numpy.savetxt(tracefile, traces, delimiter = ',',
              header = "dp_alpha_traces,dp_gamma_traces,alpha_decay_traces,cd_llh_traces,node_traces", comments='')


fn = '/home/yuan03/%s-20k-tssb-trace-CT-IRX2P-R1.pkl' % (codename)
fh = open(fn,'w')
cPickle.dump(tssb_traces,fh)
fh.close()
