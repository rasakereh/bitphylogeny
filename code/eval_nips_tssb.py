import os
import sys
import time
import cPickle
import numpy        as np
import numpy.random as npr
import gpunumpy     as gpn

from scipy.sparse import *
from util         import *

npr.seed(1234)
dir = "."
num_loops  = 1000
num_topics = int(sys.argv[1])

raw_nips   = nips()
vocab      = np.array(map(str, raw_nips['words']))
vocab_size = len(vocab)

def eval_fold(fold):
    file       = "%s/nips_tssb_%03d/result-%05d.pkl" % (dir, num_topics, fold)
    if not os.path.exists(file):
        print >>sys.stderr, "Fold %d incomplete." % (fold)
        return
    fh         = open(file)
    result     = cPickle.load(fh)
    fh.close()
    print >>sys.stderr, "Fold %d loaded." % (fold)

    fold_file    = "nips-fold-%02d.pkl" % (fold+1)
    fh           = open(fold_file)
    fold_indices = cPickle.load(fh)
    test_indices = fold_indices['test_indices']
    fh.close()
    num_docs = test_indices.shape[0]
    counts   = gpn.gpuarray(raw_nips['counts'][:,test_indices].toarray().T)

    num_comps = num_loops*len(result['result'])
    logprobs = gpn.zeros((num_comps, num_docs))
    idx = 0
    for n, res in enumerate(result['result']):
        if np.mod(n+1,10)==0:
            print >>sys.stderr, "\tResult %d" % (n+1)
        topic_dists      = gpn.gpuarray(res['topic_dists'])
        word_dist        = gpn.gpuarray(res['word_dist'])
        for loop in xrange(num_loops):
            topic_dist       = topic_dists[npr.randint(len(topic_dists))]
            multinom         = gpn.sum(topic_dist[np.newaxis,:] * word_dist, axis=1)
            logprobs[idx,:]  = gpn.sum(counts*gpn.log(multinom), axis=1) - gpn.log(num_comps)
            idx += 1
    logprobs = logprobs.T

    maxes = gpn.max(logprobs, axis=1)
    logprobs = gpn.log(gpn.sum(gpn.exp(logprobs - maxes[:,np.newaxis]), axis=1)) + maxes
    perplexity = gpn.exp(-gpn.mean(logprobs/gpn.sum(counts, axis=1)))

    print fold, gpn.mean(logprobs), perplexity

if len(sys.argv) > 2:
    eval_fold(int(sys.argv[2]))
else:
    for fold in xrange(10):
        eval_fold(fold)
