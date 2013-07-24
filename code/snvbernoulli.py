import sys
import scipy.stats

from util         import *
from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *

def truncbeta(a, b):
    res = beta(a,b)
    for index, x in enumerate(res):
        if x == 1.0:
            res[index] = 0.9999
        elif isnan(x):
            res[index] = 0.0001
    return res

class SNVBernoulli(Node):

    #bbeta        = 0.1
    min_bbeta    = 0.1
    max_bbeta    = 0.9
    hmc_accepts  = 1
    hmc_rejects  = 1
    
    def __init__(self, parent=None, dims=1, tssb=None, initial_snvs=0, bbeta=0.6):
        super(SNVBernoulli, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self.init_mean = initial_snvs
            self._bbeta  = bbeta*ones(self.dims)
            self.params = truncbeta(10**((self.init_mean-self.bbeta())/abs(self.init_mean-self.bbeta())),1.2**((self.init_mean-self.bbeta())/abs(self.init_mean-self.bbeta())))
            
        else:
            self.dims   = parent.dims
            #self.params = self.drift()*randn(self.dims) + parent.params
            self.params = truncbeta(10**((parent.params-self.bbeta())/abs(parent.params-self.bbeta())), 1.2**((parent.params-self.bbeta())/abs(parent.params-self.bbeta())))
                        
        self._cache_ln()
    
    def _cache_ln(self):
        assert any(self.params[newaxis,:]<1) 
        assert any(self.params[newaxis,:]>0)
        self._ln    = log(self.params[newaxis,:])
        self._negln = log(1-self.params[newaxis,:])

    def bbeta(self):
        if self.parent() is None:
            return self._bbeta
        else:
            return self.parent().bbeta()    
    
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):

        data       = self.get_data()
        counts     = sum(data[:,[2,3]], axis=0)
        num_data   = data.shape[0]
        bbeta      = self.bbeta()
        
        def logpost(params):
            if any(params>=1) or any(params<=0):
                llh = -sys.maxint
                return llh
            if self.parent() is None:
                llh = sum(betapdfln(params, 10**((self.init_mean-self.bbeta())/abs(self.init_mean-self.bbeta())), 1.2**((self.init_mean-self.bbeta())/abs(self.init_mean-self.bbeta()))))
            else:
                llh = sum(betapdfln(params, 10**((self.parent().params-self.bbeta())/abs(self.parent().params-self.bbeta())), 1.2**((self.parent().params-self.bbeta())/abs(self.parent().params-self.bbeta()))))
            for i in range(num_data):
                #llh = llh + data[i][2]*sigmoidln(params[data[i][0]-1]) + (1.0-data[i][2])*sigmoidln(params[data[i][0]-1]) + data[i][3]*sigmoidln(params[data[i][1]-1]) + (1.0-data[i][3])*sigmoidln(params[data[i][1]-1]) + log(1/float(self.dims*(self.dims-1)))
                llh = llh + data[i][2]*log(params[data[i][0]-1]) + (1.0-data[i][2])*log(params[data[i][0]-1]) + data[i][3]*log(params[data[i][1]-1]) + (1.0-data[i][3])*log(params[data[i][1]-1]) + log(1/float(self.dims*(self.dims-1)))
            for child in self.children():
                #llh = llh + normpdfln( child.params, params, drifts)
                llh = llh + sum(betapdfln( child.params, 10**((params-self.bbeta())/abs(params-self.bbeta())), 1.2**((params-self.bbeta())/abs(params-self.bbeta()))))
            return llh

        def logpost_grad(params):
            if any(params>=1) or any(params<=0):
                grad = -sys.maxint
                return grad
            if self.parent() is None:
                grad = (self.init_mean**(1/2) -1)/params + (1/bbeta -1)/(1-params)
            else:
                grad = (self.parent().params**(1/2) -1)/params + (bbeta -1)/(1-params)
            #grad = grad + counts*(1.0-probs) - (num_data-counts)*probs
            #how about the 1/(N*(N-1))? does the gradient remove this?
            for i in range(num_data):
                grad[data[i][0]-1] = grad[data[i][0]-1] + data[i][2]/(params[data[i][0]-1]) + (data[i][2]-1)/(1-params[data[i][0]-1])
                grad[data[i][1]-1] = grad[data[i][1]-1] + data[i][3]/(params[data[i][1]-1]) + (data[i][3]-1)/(1-params[data[i][1]-1])
            for child in self.children():
                #grad = grad + (child.params - params)/drifts**2
                grad = grad + (params**(1/2) - 1)/child.params + (bbeta -1)/(1-child.params)
            return grad

        #func   = logpost(self.params)
        #eps    = 1e-4
        #mygrad = logpost_grad(self.params)
        #fdgrad = zeros(self.params.shape)
        #for d in range(len(self.params)):
        #    mask      = zeros(self.params.shape)
        #    mask[d]   = 1
        #    fdgrad[d] = (logpost(self.params + eps*mask) - logpost(self.params - eps*mask))/(2*eps)       
        #print "MYGRAD: ", mygrad
        #print "FDGRAD: ", fdgrad
        #error = sum(abs(mygrad-fdgrad))
        #print sum(abs(mygrad-fdgrad))

        #if rand() < 0.1:
        for i in range(50):
            self.params = slice_sample(self.params, logpost, step_out=True, compwise=False)
        #else:
            #self.params, accepted = hmc(self.params, logpost, logpost_grad, 25, exponential(0.001))
            #SNVBernoulli.hmc_rejects += 1 - accepted
            #SNVBernoulli.hmc_accepts += accepted
        
        self._cache_ln()

    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpostb(bbetas):
            if any(bbetas < self.min_bbeta) or any(bbetas > self.max_bbeta):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + sum(betapdfln(child.params, root.params, bbetas))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(betapdfln(self.params, self.init_mean, bbetas))
                    
        #self._drift = slice_sample(self._drift, logpost, step_out=True, compwise=True)
        self._bbeta  = slice_sample(self._bbeta, logpostb, step_out=True, compwise=True)

    def logprob(self, x):
        x = transpose(x)
        res = x[2]*self._ln[0][x[0]-1] + (1.0-x[2])*self._negln[0][x[0]-1] + x[3]*self._ln[0][x[1]-1] + (1.0-x[3])*self._negln[0][x[1]-1] + log(1/float(self.dims*(self.dims-1)))
        assert all(res<0)
        return sum(res)

    def complete_logprob(self):
        return self.logprob(self.get_data())