import sys
import scipy.stats

from util         import *
from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *

def sigmoid(x):
    with warnings.catch_warnings(record=True) as w:
        res = 1.0/(1.0+exp(-x))
        if len(w) == 1:
            print x
        assert len(w) <> 1
    return res
def invsigmoid(z):
    return log(z) - log(1.0-z)
def sigmoidln(x):
    return -log(1.0 + exp(-x))
def normpdfln(x, m, std):
    return sum(-0.5*log(2*pi) - log(std) - 0.5*((x-m)/std)**2)

#very dirty:
#otherwise updating parents before children can results in errors
#when the parent parameter gets a bigger parameter than the child parameter 
def mixgammapdfln(x, a, b):
    res = empty(len(x))
    for index, d in enumerate(x):
        if d == 0:
            d += 0.0000001
        if d > 0:
            mix = 0.95
        else:
            mix = 0.05
            d *= -1
        res[index] = numpy.log(mix)+gammapdfln(d, a[index], b[index])
    return res

class SNVLogistic(Node):

    init_mean    = 0.0
    min_drift    = 0.01
    max_drift    = 1.0
    min_galpha   = 0.1
    max_galpha   = 2.5
    min_gbeta    = 0.01
    max_gbeta    = 0.5
    max_param    = 5
    hmc_accepts  = 1
    hmc_rejects  = 1
    
    def __init__(self, parent=None, dims=1, tssb=None, drift=1.0, galpha=1, gbeta=0.05, initial_snvs=0):
        super(SNVLogistic, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self._drift = drift*ones((dims))
            self._galpha= galpha*ones((dims))
            self._gbeta = gbeta*ones((dims))
            #self.params = self._drift*randn(dims) + self.init_mean
            #self.params = zeros((self.dims))
            #inital values already centered at clonal SNVs:
            self.init_mean = initial_snvs
            self.params = normal(self.init_mean,gbeta**2)
            
        else:
            self.dims   = parent.dims
            #self.params = self.drift()*randn(self.dims) + parent.params
            self.params = parent.params + gamma(galpha,gbeta)
            if any(abs(self.params)>self.max_param):
                for index, x in enumerate(self.params):
                    if abs(x)>self.max_param:
                        self.params[index] = parent.params[index]*0.99999    
                        
        self._cache_sigmoidln()
    
    def _cache_sigmoidln(self):
        self._sigln    = sigmoidln(self.params[newaxis,:])
        self._negsigln = sigmoidln(-self.params[newaxis,:])

    def drift(self):
        if self.parent() is None:
            return self._drift
        else:
            return self.parent().drift()
        
    def galpha(self):
        if self.parent() is None:
            return self._galpha
        else:
            return self.parent().galpha()
        
    def gbeta(self):
        if self.parent() is None:
            return self._gbeta
        else:
            return self.parent().gbeta()    
    
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):

        data       = self.get_data()
        counts     = sum(data[:,[2,3]], axis=0)
        num_data   = data.shape[0]
        #drifts     = self.drift()
        galpha     = self.galpha()
        gbeta      = self.gbeta()
        
        def logpost(params):
            #dirty way of truncating the params space
            ## if( any(abs(params)>self.max_param ) ):
            ##      llh = -sys.maxint
            ##      return llh
            if self.parent() is None:
                #TODO: change here to proper root level prior
                #llh = normpdfln(params, self.init_mean, drifts)
                #llh = sum(gammapdfln(params-self.init_mean, galpha, gbeta))

                llh = normpdfln(params, self.init_mean, gbeta**2)
            else:
                #llh = normpdfln(params, self.parent().params, drifts)
                #llh = sum(mixgammapdfln(params-self.parent().params, galpha, gbeta))

                llh = sum(gammapdfln(params-self.parent().params, galpha, gbeta))
            #llh = llh + sum(counts*sigmoidln(params)) + sum((num_data-counts)*sigmoidln(-params))

            for i in range(num_data):
                llh = llh + data[i][2]*sigmoidln(params[data[i][0]-1]) + (1.0-data[i][2])*sigmoidln(-params[data[i][0]-1]) + \
                  data[i][3]*sigmoidln(params[data[i][1]-1]) + (1.0-data[i][3])*sigmoidln(-params[data[i][1]-1]) #+ log(1/float(self.dims*(self.dims-1)))

            ## count1 = sum(data[:,2], axis=0)

            ## count2 = sum(data[:,3], axis=0)

            for child in self.children():
                
                #llh = llh + normpdfln( child.params, params, drifts)
                #llh = llh + sum(mixgammapdfln( child.params - params, galpha, gbeta))
                
                llh = llh + sum(gammapdfln( child.params - params, galpha, gbeta))
            return llh

        def logpost_grad(params):
            ## if( any(abs(params)>self.max_param ) ):
            ##      grad = -sys.maxint*ones(self.dims)
            ##      return grad
            if self.parent() is None:
                #TODO: change here to proper root level prior
                #grad = -(params-self.init_mean)/drifts**2
                #grad = (galpha-1.0)/(params-self.init_mean) - gbeta
                grad = -(params-self.init_mean)/gbeta**4
            else:
                #grad = -(params-self.parent().params)/drifts**2
                grad = (galpha-1.0)/(params-self.parent().params) - gbeta
            probs = sigmoid(params) 
            #grad = grad + counts*(1.0-probs) - (num_data-counts)*probs
            #how about the 1/(N*(N-1))? does the gradient remove this?
            for i in range(num_data):
                grad[data[i][0]-1] = grad[data[i][0]-1] + data[i][2]*(1.0-probs[data[i][0]-1]) - (1-data[i][2])*probs[data[i][0]-1]
                grad[data[i][1]-1] = grad[data[i][1]-1] + data[i][3]*(1.0-probs[data[i][1]-1]) - (1-data[i][3])*probs[data[i][1]-1]
                ## grad[data[i][0]-1] = grad[data[i][0]-1] + (1-data[i][2])*(1.0-probs[data[i][0]-1])*exp(params[data[i][0]-1]) - data[i][2]*probs[data[i][0]-1]*exp(-params[data[i][0]-1])
                ## grad[data[i][1]-1] = grad[data[i][1]-1] + (1-data[i][3])*(1.0-probs[data[i][1]-1])*exp(params[data[i][1]-1]) - data[i][3]*probs[data[i][1]-1]*exp(-params[data[i][1]-1])

            for child in self.children():
                #grad = grad + (child.params - params)/drifts**2
                grad = grad - (galpha-1.0)/(child.params-params) + gbeta
            return grad

        ## func   = logpost(self.params)
        ## eps    = 1e-4
        ## mygrad = logpost_grad(self.params)
        ## fdgrad = zeros(self.params.shape)
        ## for d in range(len(self.params)):
        ##    mask      = zeros(self.params.shape)
        ##    mask[d]   = 1
        ##    fdgrad[d] = (logpost(self.params + eps*mask) - logpost(self.params - eps*mask))/(2*eps)       
        ## print "MYGRAD: ", mygrad
        ## print "FDGRAD: ", fdgrad
        ## error = sum(abs(mygrad-fdgrad))
        ## print "Error: ", sum(abs(mygrad-fdgrad))
        ## print "params", self.params
        ## print "num_data", num_data
        ## print "data", data
        
        ## if rand() < 0.1:
        ##     self.params = slice_sample(self.params, logpost, step_out=True, compwise=True)
        ## else:
        self.params, accepted = hmc(self.params, logpost, logpost_grad, 25, exponential(0.0001))
        assert any(abs(self.params)<self.max_param)
        SNVLogistic.hmc_rejects += 1 - accepted
        SNVLogistic.hmc_accepts += accepted
        
        self._cache_sigmoidln()

    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        #def logpost(drifts):
        #    if any(drifts < self.min_drift) or any(drifts > self.max_drift):
        #        return -inf
        #    def loglh(root):
        #        llh = 0.0
        #        for child in root.children():
        #            llh = llh + normpdfln(child.params, root.params, drifts)
        #            llh = llh + loglh(child)
        #        return llh
        #    return loglh(self) + normpdfln(self.params, self.init_mean, drifts)
        
        def logposta(galphas):
            if any(galphas < self.min_galpha) or any(galphas > self.max_galpha):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + sum(mixgammapdfln(child.params - root.params, galphas, root.gbeta()))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(mixgammapdfln(self.params - self.init_mean, galphas, self.gbeta()))
        
        def logpostb(gbetas):
            if any(gbetas < self.min_gbeta) or any(gbetas > self.max_gbeta):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + sum(mixgammapdfln(child.params - root.params, root.galpha(), gbetas))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(mixgammapdfln(self.params - self.init_mean, self.galpha(), gbetas))
                    
        #self._drift = slice_sample(self._drift, logpost, step_out=True, compwise=True)
        self._galpha = slice_sample(self._galpha, logposta, step_out=True, compwise=True)
        self._gbeta  = slice_sample(self._gbeta, logpostb, step_out=True, compwise=True)

    def logprob(self, x):
        x = transpose(x)
        #both lines result in the same logprob
        #assert (x[2]*self._sigln[0][x[0]-1] - (1.0-x[2])*self._negsigln[0][x[0]-1] + x[3]*self._sigln[0][x[1]-1] - (1.0-x[3])*self._negsigln[0][x[1]-1] + log(1/float(self.dims*(self.dims-1))))==(-x[2]*log(1+exp(-self.params[x[0]-1])) + (1-x[2])*log(1+exp(self.params[x[0]-1])) - x[3]*log(1+exp(-self.params[x[1]-1])) + (1-x[3])*log(1+exp(self.params[x[1]-1])) + log(1/float(self.dims*(self.dims-1))))
        #return x[2]*self._sigln[0][x[0]-1] - (1.0-x[2])*self._negsigln[0][x[0]-1] + x[3]*self._sigln[0][x[1]-1] - (1.0-x[3])*self._negsigln[0][x[1]-1] + log(1/float(self.dims*(self.dims-1)))
        res = -x[2]*log(1+exp(-self.params[x[0]-1])) - (1-x[2])*log(1+exp(self.params[x[0]-1])) - x[3]*log(1+exp(-self.params[x[1]-1])) - (1-x[3])*log(1+exp(self.params[x[1]-1])) + log(1/float(self.dims*(self.dims-1)))
        assert all(res<0)
        return sum(res)

    def complete_logprob(self):
        return self.logprob(self.get_data())
