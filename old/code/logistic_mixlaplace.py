import sys
import scipy.stats

from util         import *
from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *

def sigmoid(x):
    return 1.0/(1.0+exp(-x))
def invsigmoid(z):
    return log(z) - log(1.0-z)
def sigmoidln(x):
    return -log(1.0 + exp(-x))
def normpdfln(x, m, std):
    return -0.5*log(2*pi) - log(std) - 0.5*((x-m)/std)**2

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

def laplacepdfln(x,m,std):
    return -log(2*std) - abs(x - m)/std

def get_weights(parent_params, depth):
    mapper = (parent_params > ones(len(parent_params)))
    m1 = 0.99
    m2 = (depth + 1.0) / len( parent_params ) 
    return m1 * mapper + ~mapper*m2

def mixlaplaceprior(params, base_value, std, depth):

    weights = get_weights(params, depth)    
    mapper = ( rand(len(weights)) < weights )
    m1 = laplace(base_value*ones(len(weights)), std*ones(len(weights)))
    m2 = laplace(-base_value*ones(len(weights)), std*ones(len(weights)))
    
    return m1 * mapper + ~mapper*m2

def mixlaplacepdfln(x, m, std, parent_params, depth ):

    weights = get_weights(parent_params, depth)

    res = empty(len(x))
    for index , d in enumerate(x):
        l1 = log( weights[index] ) + laplacepdfln(d, m[0], std)
        if weights[index] == 1.0:
            res[index] = l1
        else:
            l2 = log(1- weights[index]) + laplacepdfln(d, m[1], std)
            res[index] = logsumexp( (l1, l2) )
    return sum(res)

def rootprior(params, base_value, std):

    weights = 0.01*ones(len(params))

    mapper = ( rand(len(weights)) < weights )

    m1 = laplace(base_value*ones(len(weights)), std*ones(len(weights)))
    m2 = laplace(-base_value*ones(len(weights)), std*ones(len(weights)))

    return m1 * mapper + ~mapper*m2

def rootpriorpdfln(x, m, std):

    weights = 0.01*ones(len(x))

    res = empty(len(x))
    for index , d in enumerate(x):
        l1 = log( weights[index] ) + laplacepdfln(d, m[0], std)  
        l2 = log(1- weights[index]) + laplacepdfln(d, m[1], std)
        res[index] = logsumexp( (l1, l2) )
    return sum(res)

class Logistic(Node):

    std          = 1e-2
    base_value   = 4.5
    init_mean    = -base_value-4.0
    min_drift    = 0.01
    max_drift    = 1.0
    hmc_accepts  = 1
    hmc_rejects  = 1
    min_galpha   = 0.1
    max_galpha   = 2.5
    min_gbeta    = 0.01
    max_gbeta    = 0.5
    max_param    = 20.0

    def __init__(self, parent=None, dims=1, tssb=None, drift=1.0, \
                 galpha=1.0,gbeta=0.1):
        super(Logistic, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self.depth  = 0.0
            self._galpha = galpha*ones((dims))
            self._gbeta = gbeta*ones((dims))
            self._drift = drift*ones((dims))
            self.params = rootprior(self.init_mean*ones(dims),
                                    -self.init_mean,
                                    self.std)
            ## self.params = mixlaplaceprior(self.init_mean*ones(dims),
            ##                               self.base_value, self.std,
            ##                               self.depth)
    
            
        else:
            self.dims   = parent.dims
            self.depth  = parent.depth + 1.0
            self.params = mixlaplaceprior(parent.params, self.base_value,
                                          self.std, self.depth)

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
        counts     = sum(data, axis=0)
        num_data   = data.shape[0]
        drifts     = self.drift()
        galphas    = self.galpha()
        gbetas     = self.gbeta()
        mm         = (self.base_value, -self.base_value)
        std        = self.std
        depth      = self.depth

        
        def logpost(params):

            if any(params >= self.max_param) or any(params <= -self.max_param):
                llh = float('-inf')

            if self.parent() is None:
                llh = rootpriorpdfln(params, (-self.init_mean, self.init_mean), std)
            else:
                llh = mixlaplacepdfln(params, mm, std, self.parent().params,depth)

            llh = llh + sum(counts*sigmoidln(params)) + \
              sum((num_data-counts)*sigmoidln(-params))

            for child in self.children():
                llh = llh + mixlaplacepdfln(child.params, mm, std, params,depth+1.0)

            return llh


        ## def logpost_grad(params):
        ##     if self.parent() is None:
        ##         grad = -(params-self.init_mean)/drifts**2
        ##     else:
        ##         grad = -(params-self.parent().params)/drifts**2
        ##     probs = sigmoid(params) 
        ##     grad = grad + counts*(1.0-probs) - (num_data-counts)*probs
        ##     for child in self.children():
        ##         grad = grad + (child.params - params)/drifts**2
        ##     return grad

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

        ## if rand() < 0.1:
        self.params = slice_sample(self.params, logpost,
                                   sigma = 5.0, step_out=True, compwise=False)
        
        self._cache_sigmoidln()

    def resample_hypers(self):
        ## TODO sample std or weights of the Laplace mixture
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")
        
        
        ## def logpostb(gbetas):
        ##     if any(gbetas < self.min_gbeta) or any(gbetas > self.max_gbeta):
        ##         return float('-inf')
        ##     def loglh(root):
        ##         llh = 0.0
        ##         for child in root.children():
        ##             llh = llh + sum(mixgammapdfln(child.params - root.params, \
        ##                                           root.galpha(), 1/gbetas))
        ##             ## llh = llh + sum(gammapdfln(child.params - root.params, \
        ##             ##                            root.galpha(), 1/gbetas))
        ##             llh = llh + loglh(child)
        ##         return llh
        ##     return loglh(self) + sum(mixgammapdfln(self.params - self.init_mean, \
        ##                                 self.galpha(), 1/gbetas))

        ## ## loglh(self) + sum(gammapdfln(self.params - self.init_mean, \
        ##            ##                              self.galpha(), 1/gbetas))
        

        ## ## self._galpha = slice_sample( self._galpha, logposta, \
        ## ##                             step_out=True, compwise=True)
        ## self._gbeta  = slice_sample(self._gbeta,  logpostb, \
        ##                             step_out=True, compwise=True)



    def logprob(self, x):
        return sum(x*self._sigln + (1.0-x)*self._negsigln)

    def complete_logprob(self):
        return self.logprob(self.get_data())

    def parameter_log_prior(self):

        if self.parent() is None:
            #TODO: change here to proper root level prior                
            llh = rootpriorpdfln(self.params,
                                 (-self.init_mean, self.init_mean),
                                 self.std)
        else:
            llh = mixlaplacepdfln(self.params,(self.base_value,-self.base_value),
                                  self.std, self.parent().params,self.depth)

        for child in self.children():
            llh = llh + mixlaplacepdfln(child.params,
                                        (self.base_value,-self.base_value),
                                        self.std, self.params,
                                        self.depth + 1.0)
       
        return llh
        
    
        
