import os
import sys
import time
import cPickle
import csv

from numpy         import *
from numpy.random  import *
from tssb          import *
from logistic_mixlaplace_methy  import *
from util          import *
from scipy.stats   import itemfreq
from sklearn       import metrics



rand_seed     = 1235
max_data      = 100
burnin        = 0
num_samples   = 3000
checkpoint    = 50000
dp_alpha      = 2.0
dp_gamma      = 3e-1
init_drift    = 0.5
init_galpha   = 1
init_gbeta    = 0.5
alpha_decay   = 0.1
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

seed(rand_seed)


files = ['CT_IRX2P_R1.csv', 'CT_IRX2P_R4.csv', 'CT_IRX2P_R5.csv', 'CT_IRX2P_R6.csv', 'CT_IRX2P_L2.csv', 'CT_IRX2P_L3.csv', 'CT_IRX2P_L7.csv','CT_IRX2P_L8.csv']
#files = ['../cpgmethyl/output_CT_1_R.dat']

output_depth = zeros((num_samples, len(files)))
output_nodenum = zeros((num_samples, len(files)))
output_unnormpost = zeros((num_samples, len(files)))
output_bignodes = zeros((num_samples, len(files)))

sampno = 0
for seqsamp in files:
    reader = csv.DictReader(open('./data/Sottoriva/'+seqsamp),delimiter=',')
    data = []
    dataid = []
    max_snvpos = 1
    clone = []
    
    for row in reader:
        data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4']),int(row['V5']),int(row['V6']),int(row['V7']),int(row['V8'])])
    data = numpy.array(data)
    dims = data.shape[1]

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
    unnorm_post        = zeros((num_samples, 1))
    nodes_traces       = zeros((num_samples, 1))
    bignodes_traces    = zeros((num_samples, 1))
    tssb_traces        = empty((num_samples, 1),dtype = object)
    
    intervals = zeros((7))
    print "Starting MCMC run..."+seqsamp
    for iter in range(-burnin,num_samples):
    
        times = [ time.time() ]
    
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
            #print "pre-treesamp:", tssb.unnormalized_postertior()
            tssb.resample_tree_topology()
            #print "post-treesamp:", tssb.unnormalized_postertior()
    
        
        if iter > 2000:
            tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
        times.append(time.time())
     
    
        intervals = intervals + diff(array(times))   
        
        if iter >= 0:
            tssb_traces[iter]        = cPickle.dumps(tssb)  
            dp_alpha_traces[iter]    = tssb.dp_alpha
            dp_gamma_traces[iter]    = tssb.dp_gamma
            alpha_decay_traces[iter] = tssb.alpha_decay
            cd_llh_traces[iter]      = tssb.complete_data_log_likelihood()
            unnorm_post[iter]        = tssb.unnormalized_postertior()
            (weights, nodes)         = tssb.get_mixture()
            nodes_traces[iter]       = len(nodes)
            bignodes_traces[iter]    = sum(numpy.array(weights)>0.01)
        
        
        if mod(iter, 1) == 0:
            (weights, nodes) = tssb.get_mixture()
            print codename, iter, len(nodes), cd_llh_traces[iter], \
                tssb.dp_alpha, tssb.dp_gamma, \
                tssb.alpha_decay, unnorm_post[iter], 'big nodes:', int(bignodes_traces[iter])
                ## " ".join(map(lambda x: "%0.2f" % x, intervals.tolist())), \
                ## float(root.hmc_accepts)/(root.hmc_accepts+root.hmc_rejects), \
                ## root.hmc_accepts, root.hmc_rejects, time.strftime('%X %x %Z')
        intervals = zeros((7))
           
        if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
            print "\t%f is best per-data complete data likelihood so far." \
                % (cd_llh_traces[iter]/max_data)
        
        if iter >= 0:        
            depth = tssb.deepest_node_depth()
            output_depth[iter,sampno] = depth
            output_nodenum[iter,sampno] = nodes_traces[iter]
            output_unnormpost[iter,sampno] = unnorm_post[iter]
            output_bignodes[iter,sampno] = bignodes_traces[iter]
       
    sampno = sampno + 1


numpy.savetxt('./mcmc-traces/output_CT_depth'+str(rand_seed)+'_CT.csv',
              output_depth,delimiter=',', header = "'CT_R1','CT_R4','CT_R5','CT_R6','CT_L2','CT_L3','CT_L7','CT_L8'", comments = '')

numpy.savetxt('./mcmc-traces/output_CT_nodenum'+str(rand_seed)+'_CT.csv',
              output_nodenum,delimiter=',', header = "'CT_R1','CT_R4','CT_R5','CT_R6','CT_L2','CT_L3','CT_L7','CT_L8'", comments = '')

numpy.savetxt('./mcmc-traces/output_CT_unnormpost'+str(rand_seed)+'_CT.csv',
              output_unnormpost,delimiter=',', header = "'CT_R1','CT_R4','CT_R5','CT_R6','CT_L2','CT_L3','CT_L7','CT_L8'", comments = '')

numpy.savetxt('./mcmc-traces/output_CT_bignodes'+str(rand_seed)+'_CT.csv',
              output_bignodes,delimiter=',', header = "'CT_R1','CT_R4','CT_R5','CT_R6','CT_L2','CT_L3','CT_L7','CT_L8'", comments = '')

