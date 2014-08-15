import os
import sys
import time
import cPickle
import numpy        as np
import numpy.random as npr

from time         import *
from numpy        import *
from scipy.sparse import *
from util         import *

smooth     = 1.0
raw_nips   = nips()
vocab      = array(map(str, raw_nips['words']))
vocab_size = len(vocab)

for fold in xrange(10):
    fold_file     = "nips-fold-%02d.pkl" % (fold+1)
    fh            = open(fold_file)
    fold_indices  = cPickle.load(fh)
    fh.close()
    
    train_indices = fold_indices['train_indices']
    test_indices  = fold_indices['test_indices']
    train_counts  = raw_nips['counts'][:,train_indices].toarray().T
    test_counts   = raw_nips['counts'][:,test_indices].toarray().T

    probs = np.sum(train_counts, axis=0) + 1.0
    probs = probs / np.sum(probs)

    logprobs   = np.sum(test_counts * np.log(probs), axis=1)
    perplexity = np.exp(-np.mean(logprobs/np.sum(test_counts, axis=1)))
    
    print fold, mean(logprobs), perplexity
