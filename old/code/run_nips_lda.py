import os
import sys
import time
import cPickle
import numpy.random as npr
import lda

from time         import *
from numpy        import *
from scipy.sparse import *
from util         import *

dir              = "."
burnin           = 1000
iterations       = 5000
thinning         = 50
num_topics       = int(sys.argv[1])
alpha            = 50.0/num_topics
beta             = 0.01

def run_nips_lda(fold):
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
    print >>sys.stderr, "Initialized..."

    results = []
    for iter in xrange(-burnin, iterations):
        if mod(iter,100) == 0:
            print >>sys.stderr, iter

        lda.lda_gibbs(W, Z, NW, ND, NWsum, beta, alpha)

        if iter >= 0 and mod(iter, thinning) == 0:
            word_dists = (NW + beta) / (NWsum + vocab_size*beta)
            topic_bm   = alpha*ones((num_topics))
            results.append({ 'word_dists': word_dists, 'topic_bm': topic_bm })
    print >>sys.stderr, "Complete."
            
    return results

name = "nips_lda_%03d" % (num_topics)
expt = Experiment(name, run_nips_lda, dir)
expt.run(10)
