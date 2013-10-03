import sys
import scipy.stats

from util         import *
from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *
from scipy.stats.distributions import beta, binom
from scipy.stats.distributions import bernoulli


class SNVDiscreteEpsilon(Node):

    mutprob      = 0.1
    backmutprob  = 0.01
    hmc_accepts  = 1
    hmc_rejects  = 1
    
    def __init__(self, parent=None, dims=1, tssb=None, initial_snvs=0, epsilon=0.01, prior_depth=0):
        super(SNVDiscreteEpsilon, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self.init_mean = numpy.array(initial_snvs>0.9,dtype='int')
            self._epsilon  = epsilon
            
            self.params = self.init_mean
            self.params = self.params + numpy.array(self.params==0,dtype='int') * bernoulli.rvs(self.mutprob,size=self.dims) - numpy.array(self.params==1,dtype='int') * bernoulli.rvs(self.backmutprob,size=self.dims) 
            
            self.prior_depth = prior_depth
            self.depth = 1
            
        else:
            self.dims   = parent.dims
            self.prior_depth = parent.prior_depth

            self.depth = parent.depth + 1
            
            self.params = parent.params + numpy.array(parent.params==0,dtype='int') * bernoulli.rvs(self.mutprob,size=self.dims) - numpy.array(parent.params==1,dtype='int') * bernoulli.rvs(self.backmutprob,size=self.dims) 
           
            if isinstance( initial_snvs, Iterable):
                self.params = numpy.array(initial_snvs>0.9,dtype='int')
            
            self._epsilon = parent._epsilon
            
    
    def epsilon(self):
        if self.parent() is None:
            return self._epsilon
        else:
            return self.parent().epsilon()    
    
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):

        data       = self.get_data()
        counts     = sum(data[:,[2,3]], axis=0)
        num_data   = data.shape[0]
        epsilon      = self.epsilon()
        
        def logpost(params):
            if self.parent() is None:
                llh = log(binom.pmf(sum(numpy.array(self.params-self.init_mean == 1, dtype='int')),self.dims-sum(self.init_mean),self.mutprob)) + log(binom.pmf(sum(numpy.array(self.params-self.init_mean == -1, dtype='int')),sum(self.init_mean),self.backmutprob))
            else:
                llh = log(binom.pmf(sum(numpy.array(self.params-self.parent().params == 1, dtype='int')),self.dims-sum(self.parent().params),self.mutprob)) + log(binom.pmf(sum(numpy.array(self.params-self.parent().params == -1, dtype='int')),sum(self.parent().params),self.backmutprob))
            for i in range(num_data):
                llh = llh + log(((1-abs(data[i][2]-params[data[i][0]-1]))*(1-self.epsilon())+abs(data[i][2]-params[data[i][0]-1])*self.epsilon()) * ((1-abs(data[i][3]-params[data[i][1]-1]))*(1-self.epsilon()) + abs(data[i][3]-params[data[i][1]-1])*self.epsilon()))
            for child in self.children():
                llh = llh + log(binom.pmf(sum(numpy.array(child.params-params == 1, dtype='int')),self.dims-sum(params),self.mutprob)) + log(binom.pmf(sum(numpy.array(child.params-params == -1, dtype='int')),sum(params),self.backmutprob))
            return llh
        
        for i in range(self.dims):
            llh_s = logpost(self.params) + log(numpy.random.rand()) 
            self.params[i] = int( not (bool(self.params[i])))
            new_llh = logpost(self.params)
            if new_llh < llh_s:
                self.params[i] = int( not (bool(self.params[i])))

    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpostepsilon(epsilons):
            if any(epsilons < self.min_epsilon) or any(epsilons > self.max_epsilon):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    #llh = llh + sum(mixbetapdfln(child.params, self.alpha_base, self.beta_base, root.params, epsilons, self.mix))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) #+ sum(mixbetapdfln(self.params, self.alpha_base, self.beta_base, self.init_mean, epsilons, self.mix))
        
        self._epsilon  = slice_sample(self._epsilon, logpostcut, step_out=True, compwise=True)

    def logprob(self, x):
        x = transpose(x)
        res = log(((1-abs(x[2]-self.params[x[0]-1]))*(1-self.epsilon())+abs(x[2]-self.params[x[0]-1])*self.epsilon()) * ((1-abs(x[3]-self.params[x[1]-1]))*(1-self.epsilon()) + abs(x[3]-self.params[x[1]-1])*self.epsilon()))
        assert all(res<0)
        return sum(res)

    def complete_logprob(self):
        return self.logprob(self.get_data())