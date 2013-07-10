import sys
import lda
import scipy.stats

from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *
from util         import *
from time         import *

class Topics(Node):

    min_conc = 1.0
    max_conc = 100.0
    epsilon  = finfo(float64).eps

    def __init__(self, parent=None, tssb=None, conc=2.0, NW=None, lda_beta=0.01):
        super(Topics, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            if NW is None:
                raise Exception("Must provide initial word-topic assignments.")

            self.NW         = NW
            self.NWsum      = sum(NW, axis=0)
            self.num_words  = NW.shape[0]
            self.num_topics = NW.shape[1]            
            self.lda_beta   = lda_beta
            self._conc      = conc
            self.gammas     = gamma(self.conc()*ones((self.num_topics))) + self.epsilon
            print >>sys.stderr, "Topic Model: %d topics, %d words" % (self.num_topics, self.num_words)

        else:
            self.num_words  = parent.num_words
            self.num_topics = parent.num_topics
            self.gammas     = gamma(self.conc()*self.parent().normed()) + self.epsilon
            
    def normed(self):
        sane = self.gammas + self.epsilon
        return sane / sum(sane)
            
    def conc(self):
        if self.parent() is None:
            return self._conc
        else:
            return self.parent().conc()
        
    def sample(self, args):
        raise Exception("Can't currently fantasise from topic models.")
    
    def resample_params(self):
        
        data   = self.get_data()
        conc   = self.conc()
        counts = zeros((self.num_topics), dtype=int32)
        for doc in data:
            counts = counts + doc['ND']

        def logpost(gammas):
            if any(gammas <= 0.0):
                return -inf
            
            sane   = gammas + self.epsilon
            normed = sane / sum(sane)

            if self.parent() is None:
                llh = sum(gammapdfln( gammas, self.conc()*ones((self.num_topics)), 1 )) + sum(counts*log(normed))
            else:
                llh = sum(gammapdfln( gammas, self.conc()*self.parent().normed(), 1 )) + sum(counts*log(normed))

            for child in self.children():
                llh += sum(gammapdfln( child.gammas, self.conc()*normed, 1 ))

            if not isfinite(llh):
                print gammas
                print normed
                raise Exception("Got a NaN")
            
            return llh

        self.gammas = slice_sample(self.gammas, logpost, step_out=True, compwise=True)
        
    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpost(conc):
            def descend(root):
                llh    = 0.0
                normed = root.normed()
                for child in root.children():
                    llh += sum(gammapdfln( child.gammas, conc*normed, 1))
                return llh
            return sum(gammapdfln(self.gammas, conc*ones((self.num_topics)), 1)) + descend(self)
        
        upper = self.max_conc
        lower = self.min_conc
        llh_s = log(rand()) + logpost(self._conc)
        while True:
            new_conc = (upper-lower)*rand() + lower
            new_llh  = logpost(new_conc)
            if not isfinite(new_llh):
                raise Exception("Log posterior is not finite: %f (%f)" % (new_llh, new_conc))
            if new_llh > llh_s:
                break
            elif new_conc < self._conc:
                lower = new_conc
            elif new_conc > self._conc:
                upper = new_conc
            else:
                raise Exception("Slice sampler shrank to zero!")
        self._conc = new_conc

    def logprob(self, data):
        counts = zeros((self.num_topics), dtype=int32)
        for doc in data:
            counts = counts + doc['ND']
        return sum( counts*log(self.normed()) )

    def complete_logprob(self):

        logprob  = 0.0
        NW       = self.global_param('NW')
        NWsum    = self.global_param('NWsum')
        lda_beta = self.global_param('lda_beta')
        normed   = self.normed()
        thetas   = (NW + lda_beta) / (NWsum + NW.shape[0]*lda_beta)
        
        for doc in self.get_data():            
            W  = doc['W']
            Z  = doc['Z']

            if True:
                logprob += lda.doc_logprob(W, Z, thetas, normed)
            else:
                for i in xrange(len(Z)):
                    logprob += log(thetas[W[i],Z[i]]) + log(normed[Z[i]])

        return logprob


    def global_param(self, key):
        if self.parent() is None:
            return self.__dict__[key]
        else:
            return self.parent().global_param(key)

    def resample_words(self):
        if self.parent() is not None:
            raise Exception("Can only update topics from root!")

        def descend(node):
            normed   = node.normed()
            NW       = node.global_param('NW')
            NWsum    = node.global_param('NWsum')
            lda_beta = node.global_param('lda_beta')

            for n, doc in enumerate(node.get_data()):

                W  = doc['W']
                Z  = doc['Z']
                ND = doc['ND']

                if True:
                    lda.doc_gibbs(W, Z, NW, ND, NWsum, normed, lda_beta)
                else:                
                    for i in xrange(len(Z)):
                            
                        cur_topic           = Z[i]
                        NW[W[i],cur_topic] -= 1
                        ND[cur_topic]      -= 1
                        NWsum[cur_topic]   -= 1

                        probs = ((NW[W[i],:] + lda_beta) / (NWsum + NW.shape[0]*lda_beta)) * normed
                        probs = probs / sum(probs)

                        new_topic = sum(rand() > cumsum(probs))

                        Z[i]                = new_topic
                        NW[W[i],new_topic] += 1
                        ND[new_topic]      += 1
                        NWsum[new_topic]   += 1
            
            for child in node.children():
                descend(child)
        
        descend(self)
