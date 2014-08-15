import os
import sys
import time
import cPickle
import numpy.random as npr

from time         import *
from numpy        import *
from scipy.sparse import *
from tssb         import *
from topics       import *
from util         import *

dir              = "."
lda_burnin       = 1000
burnin           = 1000
iterations       = 5000
thinning         = 50
num_topic_dists  = 50
num_topics       = int(sys.argv[1])
lda_alpha        = 50.0/num_topics
lda_beta         = 0.01
init_dp_alpha    = 25.0
init_dp_gamma    = 1.0
init_alpha_decay = 0.5
init_topic_conc  = 200.0

def run_nips_tssb(fold):
    npr.seed(fold)

    raw_nips   = nips()
    vocab      = array(map(str, raw_nips['words']))
    vocab_size = len(vocab)

    fold_file    = "nips-fold-%02d.pkl" % (fold+1)
    fh           = open(fold_file)
    fold_indices = cPickle.load(fh)
    fh.close()
    
    all_counts = raw_nips['counts'][:,fold_indices['train_indices']]
    num_docs   = all_counts.shape[1]

    W     = []
    Z     = []                                           # Z[i][j][k] assignment of word_j, token_k in document_i
    NW    = zeros((vocab_size, num_topics), dtype=int32) # NW[i,j] number of word_i assigned to topic_j
    ND    = zeros((num_docs, num_topics), dtype=int32)   # ND[i,j] number of words in doc_i assigned to topic_j
    NWsum = zeros((num_topics), dtype=int32)             # NWsum[i] number of words in topic_i

    # Convert the data and initialize the Z.
    for n in xrange(num_docs):
        counts = all_counts[:,n]
        words  = []
        topics = []
        for word in nonzero(counts)[0]:
            words.extend([word for x in xrange(counts[word,0])])
            topics.extend(npr.randint(num_topics, size=counts[word,0]))

        topics = array(topics, dtype=int32)
        words  = array(words, dtype=int32)

        for topic in xrange(num_topics):
            indices     = nonzero(topics == topic)
            topic_words = words[indices]

            for word in topic_words:
                NW[word,topic] += 1
            ND[n,topic]  += topic_words.shape[0]
            NWsum[topic] += topic_words.shape[0]

        W.append(words)
        Z.append(topics)
    print >>sys.stderr, "LDA Initialized"

    # Run LDA.
    for iter in xrange(lda_burnin):
        if mod(iter,100) == 0:
            print >>sys.stderr, "LDA %d" % (iter)

        lda.lda_gibbs(W, Z, NW, ND, NWsum, lda_beta, lda_alpha)

    # Convert data to TSSB format.
    data = []
    for n in xrange(num_docs):
        data.append({ 'W': W[n], 'Z': Z[n], 'ND': ND[n] })
    data = array(data)

    # Initialize TSSB.
    root  = Topics(NW=NW, conc=init_topic_conc)
    tssb  = TSSB( dp_alpha=init_dp_alpha, dp_gamma=init_dp_gamma, alpha_decay=init_alpha_decay,
                  root_node=root, max_depth=25, data=data )
    print >>sys.stderr, "TSSB Initialized"

    cd_llh_traces = zeros((iterations, 1))

    results = []
    intervals = zeros((9))
    for iter in range(-burnin, iterations):
 
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

        if iter > -burnin//2:
            tssb.resample_hypers()
        times.append(time.time())

        if iter > -burnin//2:            
            root.resample_words()
        times.append(time.time())

        cd_llh_traces[iter] = tssb.complete_data_log_likelihood()
        times.append(time.time())

        intervals = intervals + diff(array(times))
   
        if mod(iter, 10) == 0:
            
            (weights, nodes) = tssb.get_mixture()

            if iter >= 0:
                print >>sys.stderr, iter, len(nodes), cd_llh_traces[iter], tssb.dp_alpha, tssb.dp_gamma, tssb.alpha_decay,  root._conc,
                print >>sys.stderr, " ".join(map(lambda x: "%0.2f" % x, intervals.tolist()))
            else:
                print >>sys.stderr, iter, len(nodes), tssb.dp_alpha, tssb.dp_gamma, tssb.alpha_decay, 
                print >>sys.stderr, " ".join(map(lambda x: "%0.2f" % x, intervals.tolist()))
                
            intervals = zeros((9))

        if iter > 0 and all(cd_llh_traces[iter] > cd_llh_traces[:iter]):            
            filename = "bests/nips-f%02d-t%03d-best.pkl" % (fold, num_topics)
            fh = open(filename, 'w')
            cPickle.dump(tssb, fh)
            fh.close()

        if iter >= 0 and mod(iter, thinning) == 0:
            topic_dists = []
            for i in xrange(num_topic_dists):
                (node, path) = tssb.find_node(npr.rand())
                topic_dist = node.normed()
                topic_dists.append(topic_dist)
            word_dist = (root.NW + root.lda_beta) / (root.NWsum + root.num_words*root.lda_beta)
            results.append({'topic_dists': topic_dists, 'word_dist': word_dist})
    print >>sys.stderr, "Complete."
    return results

name = "nips_tssb_%03d" % (num_topics)
expt = Experiment(name, run_nips_tssb, dir)
expt.run(10)
