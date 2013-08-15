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

def exp_trans(x, opt):
    if opt == 'tran_x':
        return log(x);
    if opt == 'orig_x':
        return exp(x);
    if opt == 'det_jacob':
        return x;

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
    
    def __init__(self, parent=None, dims=1, tssb=None, drift=1.0, \
                 galpha=1, gbeta=0.05, initial_snvs=0):
        super(SNVLogistic, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            ## self._drift = drift*ones((dims))
            self._galpha= galpha*ones((dims))
            self._gbeta = gbeta*ones((dims))
            #self.params = self._drift*randn(dims) + self.init_mean
            #self.params = zeros((self.dims))
            #inital values already centered at clonal SNVs:
            self.init_mean = initial_snvs
            self.params = self.init_mean + gamma(galpha,gbeta)
            if any(abs(self.params)>self.max_param):
                for index, x in enumerate(self.params):
                    if abs(x)>self.max_param:
                        self.params[index] = parent.params[index]*0.99999 
            
            
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

        def get_constraint():

            if  self.parent() is None:
                lower = -sys.maxint*ones(self.dims)
            else:
                 lower = self.parent().params

            if  len(self.children()) < 1:
                upper = inf*ones((self.dims))
            else:
                child_params = zeros((len(self.children()), self.dims))
                for index, child in enumerate(self.children()):
                    child_params[index,:] = child.params;
                upper = child_params.max(0);

            return lower, upper

        lower, upper = get_constraint()
        
        def logpost(params):
            xx = transpose(data)
            #dirty way of truncating the params space

            if ( any( (params-lower) < 0.0 ) ) :
               llh = -sys.maxint
               return llh

            if ( any(upper - params) < 0.0 ):
               llh = -sys.maxint
               return llh
            
            if( any(abs(params)>self.max_param ) ):
                  llh = -sys.maxint
                  return llh
            if self.parent() is None:
                #TODO: change here to proper root level prior                
                #llh = sum(gammapdfln(params-self.init_mean, galpha, gbeta))
                llh = normpdfln(params, self.init_mean, gbeta**2)
            else:
                #llh = sum(mixgammapdfln(params-self.parent().params, galpha, gbeta))
                llh = sum(gammapdfln(params-self.parent().params, galpha, gbeta))
          
            ll = sum( xx[2]*sigmoidln(params[xx[0]-1]) + \
                      (1-xx[2])*sigmoidln(-params[xx[0]-1]) + \
                      xx[3]*sigmoidln(params[xx[1]-1]) + \
                      (1-xx[3])*sigmoidln(-params[xx[1]-1]) )
            llh = llh + ll

            for child in self.children():
                #llh = llh + sum(mixgammapdfln( child.params - params, galpha, gbeta))
                llh = llh + sum(gammapdfln( child.params - params, galpha, gbeta))
            return llh

        def logpost_grad(params):

            if ( any( (params-lower) < 0.0 ) ) :
               grad = -sys.maxint*ones(self.dims)
               return grad
    
            if ( any(upper-params) < 0.0 ):
               grad = -sys.maxint*ones(self.dims)
               return grad
            
            if( any(abs(params)>self.max_param ) ):
                 grad = -sys.maxint*ones(self.dims)
                 return grad
            if self.parent() is None:
                #TODO: change here to proper root level prior
                grad = -(params-self.init_mean)/gbeta**4
            else:             
                grad = (galpha-1.0)/(params-self.parent().params) - gbeta
             
           
            #how about the 1/(N*(N-1))? does the gradient remove this? yes,
            #the gradient removes this term
            probs = sigmoid(params)
            for i in range(num_data):
                grad[data[i][0]-1] = grad[data[i][0]-1] + \
                                     data[i][2]*(1.0-probs[data[i][0]-1]) - \
                                     (1-data[i][2])*probs[data[i][0]-1]
                grad[data[i][1]-1] = grad[data[i][1]-1] + \
                                     data[i][3]*(1.0-probs[data[i][1]-1]) - \
                                     (1-data[i][3])*probs[data[i][1]-1]
              

            for child in self.children():
                grad = grad - (galpha-1.0)/(child.params-params) + gbeta
            return grad

        ## func   = logpost(self.params)
        ## eps    = 1e-4
        ## mygrad = logpost_grad(self.params)
        ## fdgrad = zeros(self.params.shape)
        ## for d in range(len(self.params)):
        ##    mask      = zeros(self.params.shape)
        ##    mask[d]   = 1
        ##    fdgrad[d] = (logpost(self.params + eps*mask) - \
        ##                 logpost(self.params - eps*mask))/(2*eps)       
        ## print "MYGRAD: ", mygrad
        ## print "FDGRAD: ", fdgrad
        ## error = sum(abs(mygrad-fdgrad))
        ## print "Error: ", sum(abs(mygrad-fdgrad))
        ## print "params", self.params
        ## print "num_data", num_data
        ## print "data", data
        
        ## if rand() < 0.1:
        ##     self.params = slice_sample(self.params, \
        ##                                logpost, \
        ##                                step_out=True, \
        ##                                compwise=True)
        ## else:


        ## self.params, accepted = hmc(self.params, logpost, \
        ##                             logpost_grad, \
        ##                             25, exponential(0.0005))

        self.params, accepted = bounded_hmc(self.params, logpost, \
                                            logpost_grad, lower, upper, \
                                            25, exponential(0.0005))
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
            ## if any(galphas < self.min_galpha) or any(galphas > self.max_galpha):
            ##     return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    ## llh = llh + sum(mixgammapdfln(child.params - \
                    ##       root.params, galphas, root.gbeta()))
                    llh = llh + sum(gammapdfln(child.params - root.params, \
                                               exp(galphas), root.gbeta()))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(gammapdfln(self.params - self.init_mean, \
                                                exp(galphas), self.gbeta())) + \
                                                sum(galphas)
        ## loglh(self) + sum(mixgammapdfln(self.params - self.init_mean, \
        ##                                 galphas, self.gbeta()))
        
        
        def logpostb(gbetas):
            ## if any(gbetas < self.min_gbeta) or any(gbetas > self.max_gbeta):
            ##     return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    ## llh = llh + sum(mixgammapdfln(child.params - root.params, \
                    ##                               root.galpha(), gbetas))
                    llh = llh + sum(gammapdfln(child.params - root.params, \
                                               root.galpha(), exp(gbetas)))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(gammapdfln(self.params - self.init_mean, \
                                                self.galpha(), exp(gbetas))) + \
                                                sum(gbetas)
        ## loglh(self) + sum(mixgammapdfln(self.params - self.init_mean, \
        ##                                 self.galpha(), gbetas))

        tmp_galpha = slice_sample(log(self._galpha), logposta, \
                                    step_out=True, compwise=True)
        tmp_gbeta  = slice_sample(log(self._gbeta),  logpostb, \
                                    step_out=True, compwise=True)

        self._galpha = exp(tmp_galpha)
        self._gbeta = exp(tmp_gbeta)
        
    def logprob(self, x):
        x = transpose(x)
        llh = -x[2]*log(1.0+exp(-self.params[x[0]-1])) - \
              (1 -x[2])*log(1.0+exp(self.params[x[0]-1])) - \
              x[3]*log(1.0+exp(-self.params[x[1]-1])) - \
              (1-x[3])*log(1.0+exp(self.params[x[1]-1])) 
        assert all(llh<0)

        ## llh = x[2]*self._sigln[:,x[0]-1] + (1.0 - x[2])*self._negsigln[:,x[0]-1] + \
        ##       x[3]*self._sigln[:,x[1]-1] + (1.0 - x[3])*self._negsigln[:,x[1]-1] - \
        ##       log(self.dims) - log(self.dims-1)

        ## print "llh error: ", sum(res) - sum(llh)
        return sum(llh)
              
    def complete_logprob(self):
        return self.logprob(self.get_data())
    













