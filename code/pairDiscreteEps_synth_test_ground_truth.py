import os
import sys
import time
import cPickle
import csv

from numpy         import *
from numpy.random  import *
from tssb          import *
from snvdiscreteepsilon  import *
from util          import *
from scipy.stats   import itemfreq
import numpy


rand_seed     = 1234
max_data      = 100
burnin        = 0
num_samples   = 1000
checkpoint    = 50000
dp_alpha      = 1
dp_gamma      = 1
init_epsilon   = 0.02
alpha_decay   = 0.1
codename      = os.popen('./random-word').read().rstrip()
print "Codename: ", codename

seed(rand_seed)
#data = concatenate([0.5 + 0.2*randn(70,1), -0.8 + 0.1*randn(35,1), -1.5+0.1*randn(35,1)])
#data = sigmoid(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat1.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data1 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat2.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data2 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat3.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data3 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat4.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data4 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat5.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data5 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat6.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data6 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_mutmat7.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data7 = numpy.array(data)

reader = csv.DictReader(open('./data/groundtsyn_genotypes.csv'),
                        delimiter=',')
genotypes = []
for row in reader:
    genotypes.append([int(row['p1']),int(row['p2']),int(row['p3']),int(row['p4']),int(row['p5']),int(row['p6']),int(row['p7'])])
genotypes = numpy.array(genotypes)

dims = 50


reader = csv.DictReader(open('./data/groundtsyn_mutmat.csv'),
                        delimiter=',')
data = []
for row in reader:
    data.append([int(row['V1']),int(row['V2']),int(row['V3']),int(row['V4'])])
data = numpy.array(data)

#filename = "test_beta_benomial_depth2.pkl"
#fh = open(filename, 'r')
#depth_traces = cPickle.load(fh)
#fh.close()

#median_depths = median(depth_traces,0)
genotypes = transpose(genotypes)*0.98
genotypes = genotypes + 0.01*ones((genotypes.shape[0], genotypes.shape[1]))

root = SNVDiscreteEpsilon( dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[0] )
c1 = SNVDiscreteEpsilon( parent=root,dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[1] )
root.add_child(c1)
c2 = SNVDiscreteEpsilon( parent=c1, dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[2] )
c1.add_child(c2)
c3 = SNVDiscreteEpsilon( parent=c1, dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[3] )
c1.add_child(c3)
c4 = SNVDiscreteEpsilon( parent=c1, dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[4] )
c1.add_child(c4)
c5 = SNVDiscreteEpsilon( parent=c2, dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[5] )
c2.add_child(c5)
c6 = SNVDiscreteEpsilon( parent=c4, dims=dims, epsilon=init_epsilon, initial_snvs=genotypes[6] )


tssb = TSSB( dp_alpha=dp_alpha, dp_gamma=dp_gamma, max_depth=4, alpha_decay=alpha_decay,
             root_node=root, data=data )
c1.tssb = tssb
c2.tssb = tssb
c3.tssb = tssb
c4.tssb = tssb
c5.tssb = tssb
c6.tssb = tssb

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
c4.add_child(c6)
tssb.root['children'][0]['children'][2]['sticks'] = vstack([ tssb.root['children'][0]['children'][2]['sticks'], boundbeta(1, tssb.dp_gamma) ])
tssb.root['children'][0]['children'][2]['children'].append({ 'node'     : c6,
                          'main'     : boundbeta(1.0, (tssb.alpha_decay**(depth+1))*tssb.dp_alpha) if tssb.min_depth <= (depth+1) else 0.0,
                          'sticks'   : empty((0,1)),
                          'children' : [] })


for n in range(tssb.num_data):
    tssb.assignments[n].remove_datum(n)
tssb.assignments = []
    
for i in range(data1.shape[0]):
    root.add_datum(i)
    tssb.assignments.append(root)
for i in range(data2.shape[0]):
    ii=i+data1.shape[0]
    c1.add_datum(ii)
    tssb.assignments.append(c1)
for i in range(data3.shape[0]):
    ii=i+data1.shape[0]+data2.shape[0]
    c2.add_datum(ii)
    tssb.assignments.append(c2)
for i in range(data4.shape[0]):
    ii=i+data1.shape[0]+data2.shape[0]+data3.shape[0]
    c3.add_datum(ii)
    tssb.assignments.append(c3)
for i in range(data5.shape[0]):
    ii=i+data1.shape[0]+data2.shape[0]+data3.shape[0]+data4.shape[0]
    c4.add_datum(ii)
    tssb.assignments.append(c4)
for i in range(data6.shape[0]):
    ii=i+data1.shape[0]+data2.shape[0]+data3.shape[0]+data4.shape[0]+data5.shape[0]
    c5.add_datum(ii)
    tssb.assignments.append(c5)
for i in range(data7.shape[0]):
    ii=i+data1.shape[0]+data2.shape[0]+data3.shape[0]+data4.shape[0]+data5.shape[0]+data6.shape[0]
    c6.add_datum(ii)
    tssb.assignments.append(c6)
    

dp_alpha_traces    = zeros((num_samples, 1))
dp_gamma_traces    = zeros((num_samples, 1))
alpha_decay_traces = zeros((num_samples, 1))
epsilon_traces       = zeros((num_samples, dims))
cd_llh_traces      = zeros((num_samples, 1))
nodes_traces       = zeros((num_samples, 1))
tssb_traces        = empty((num_samples, 1),dtype = object)


intervals = zeros((7))
print "Starting MCMC run..."
for iter in range(-burnin,num_samples):

    print "ground truth nomix: data tree fit, node-wise fit ",tssb.complete_data_log_likelihood_nomix()

    times = [ time.time() ]
    print tssb.complete_data_log_likelihood()
    if iter>5:
        tssb.resample_assignments()
    times.append(time.time())
    print tssb.complete_data_log_likelihood()
    tssb.cull_tree()
    times.append(time.time())
    #print tssb.complete_data_log_likelihood()
    if iter>5:
        tssb.resample_node_params()
    times.append(time.time())
    print tssb.complete_data_log_likelihood()
    #root.resample_hypers()
    times.append(time.time())
    print tssb.complete_data_log_likelihood()
    tssb.resample_sticks()
    times.append(time.time())
    print tssb.complete_data_log_likelihood()
    tssb.resample_stick_orders()
    times.append(time.time())
    print tssb.complete_data_log_likelihood()

    filename = 'treescripts/testgraph_pairDiscrete1_groundt62_.gdl'
    fh2 = open(filename,'w')
    tssb.print_graph_pairing(fh2)
    fh2.close()

   
    if iter > 0:
        tssb.resample_hypers(dp_alpha=True, alpha_decay=True, dp_gamma=True)
    times.append(time.time())
 

    intervals = intervals + diff(array(times))

    
    if iter >= 0:
        tssb_traces[iter]        = cPickle.dumps(tssb)  
        dp_alpha_traces[iter]    = tssb.dp_alpha
        dp_gamma_traces[iter]    = tssb.dp_gamma
        alpha_decay_traces[iter] = tssb.alpha_decay
        #drift_traces[iter]       = root.drift()
        epsilon_traces[iter]       = root.epsilon()
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
           mean(root._epsilon), tssb.dp_alpha, tssb.dp_gamma, \
           tssb.alpha_decay, \
           " ".join(map(lambda x: "%0.2f" % x, intervals.tolist())), \
           float(root.hmc_accepts)/(root.hmc_accepts+root.hmc_rejects), \
           root.hmc_accepts, root.hmc_rejects, time.strftime('%X %x %Z')
        intervals = zeros((7))

    if mod(iter, 10) == 0:
        print "root:"
        print " ".join(map(lambda x: "%0.2f" %x, root.params))
        print "children:"
        for child in root._children:
            print " ".join(map(lambda x: "%0.2f" %x, child.params))
    
    if iter > 0 and argmax(cd_llh_traces[:iter+1]) == iter:
        print "\t%f is best per-data complete data likelihood so far." \
           % (cd_llh_traces[iter]/max_data)
        best_fit = cPickle.dumps(tssb)
        filename_best = "bests/pairDiscrete-test_groundt_-%s-best.pkl" % (codename)
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
best_num_nodes = nodes_tabular[argmax(nodes_tabular[:,1]),0]
best_num_nodes_llh = cd_llh_traces[nodes_traces==best_num_nodes].max()
best_node_fit = cPickle.loads(tssb_traces[cd_llh_traces==best_num_nodes_llh][0])

(weights_best, nodes_best) = best_node_fit.get_mixture()
yy = linspace(0,1,1000)
pp = zeros(yy.shape) 

#for kk in range(len(weights_best)):
#    pp = pp + weights_best[kk]*1/sqrt(2*pi*nodes_best[kk].params[1]**2) * \
#      exp(-0.5*(yy-nodes_best[kk].params[0])**2/ nodes_best[kk].params[1]**2)

fig1 = plt.figure(1)
plt.plot(cd_llh_traces)
plt.savefig('figures/pairDiscrete_groundt_trace.pdf',format='pdf')
clf()

#fig2 = plt.figure(2)
#plt.plot(yy,pp, color = 'b')
#plt.plot(data,-0.1*ones(data.shape), linestyle = 'none',
#         color = 'g', marker = 'x')
#plt.ylim((-0.2,0.8))
#plt.savefig('figures/pairDiscrete11.pdf', format='pdf')
#clf()


filename = 'treescripts/testgraph_pairDiscrete1_groundt_.gdl'
fh2 = open(filename,'w')
best_node_fit.print_graph_pairing(fh2)
fh2.close()

filename_best = "bests/pairDiscrete1_groundt_test.pkl" 
fh = open(filename_best, 'w')
cPickle.dump(best_node_fit, fh)
fh.close()

fn = "bests/pairDiscrete1_groundt_test.pkl"
fh = open(fn, 'r')
best_node_fit = cPickle.load(fh)
fh.close()
