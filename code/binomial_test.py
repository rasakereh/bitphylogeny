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
    return sum(-0.5*log(2*pi) - log(std) - 0.5*((x-m)/std)**2)

def logit(x):
    return log(x) - log(1.0-x);
    
def sigmoid_trans(x, opt):
    if opt == 'tran_x':
        return logit(x);
    if opt == 'orig_x':
        return sigmoid(x);
    if opt == 'det_jacob':
        return log(x)+ log(1.0-x);
    if opt == 'tran_grad':
        return x*(1.0-x);
    if opt == 'jacob_grad':
        return 1.0-2*x;


def log_factorial(n):
    return gammaln(n + 1.0)

def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(x) - log_factorial(n - x)

def log_binomial_pdf(x, n, p):
    if p == 0:
        return float('-inf')
    
    if p == 1.0:
        return float('-inf')
    
    return sum(log_binomial_coefficient(n, x) + x * log(p) + (n - x) * log(1.0 - p))

class Binomial_test(Node):

    init_mean    = 1.0
    min_drift    = 0.01
    max_drift    = 1.0
    hmc_accepts  = 1
    hmc_rejects  = 1
    rho          = 0.6

    def __init__(self, parent=None, dims=1, tssb=None, drift=0.1):
        super(Binomial_test, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self._drift = drift*ones((dims))
            #self.params = self._drift*randn(dims) + self.rho*self.init_mean
            self.params = 0.9*ones((self.dims))
        else:
            self.dims   = parent.dims
            self.params = self.drift()*randn(self.dims) + \
              self.rho*parent.params
            ## if self.params < 0.0 or self.params > 1.0:
            ##     self.params = parent.params
    

    def drift(self):
        if self.parent() is None:
            return self._drift
        else:
            return self.parent().drift()
        
    ## def sample(self, args):
    ##     num_data = args['num_data'] if args.has_key('num_data') else 1
    ##     return self.param[0] + self.param[1]*randn(num_data, 1)
    
    def resample_params(self):
        data       = self.get_data()
        num_data   = data.shape[0]
        drifts     = self.drift()   
        
        def logpost(params):
            if self.parent() is None:
                llh = normpdfln(params, self.rho*self.init_mean, drifts)
            else:
                llh = normpdfln(params, self.rho*self.parent().params, drifts) 
            ## binomial node
            llh = llh + log_binomial_pdf(data[:,0],data[:,1], params)
            for child in self.children():
                llh = llh + normpdfln( child.params, self.rho*params, drifts)
            return llh

        def logpost_grad(params):
            
            if self.parent() is None:
                grad = -(params-self.rho*self.init_mean)/drifts**2
            else:
                grad = -(params-self.rho*self.parent().params)/drifts**2

            if params == 1 or params == 0:
                g = float('-inf')
            else:
                g = sum(data[:,0]/params - (data[:,1]-data[:,0])/(1 - params)) 
            grad = grad + g
            
            for child in self.children():
                grad = grad + self.rho*(child.params - self.rho*params) \
                  /drifts**2

            return grad
            

        ## tran_param = sigmoid_trans( self.params, 'tran_x')
        ## eps    = 1e-4
        ## mygrad = logpost_grad(self.params) * sigmoid_trans(self.params, 'tran_grad') +\
        ##          sigmoid_trans(self.params, 'jacob_grad')
        ## fdgrad = zeros(self.params.shape)
        ## for d in range(len(self.params)):
        ##     mask      = zeros(self.params.shape)
        ##     mask[d]   = 1
        ##     v1 = sigmoid_trans((tran_param + eps*mask),'orig_x')
        ##     v2 = sigmoid_trans((tran_param - eps*mask),'orig_x')
        ##     fdgrad[d] = (logpost( v1 ) + sigmoid_trans( v1, 'det_jacob') - \
        ##                  logpost( v2 ) - sigmoid_trans( v2, 'det_jacob') ) / (2*eps)
        ## print "MYGRAD: ", mygrad
        ## print "FDGRAD: ", fdgrad
        ## print norm(mygrad-fdgrad)/norm(mygrad+fdgrad)

        ## if rand() < 0.1:
        ##   self.params = slice_sample(self.params, logpost, step_out=True, compwise=True)
        ## else:
        ##   self.params, accepted = hmc(self.params, logpost,
        ##                               logpost_grad, 25,
        ##                               exponential(0.005))
        ## Binomial_test.hmc_rejects += 1 - accepted
        ## Binomial_test.hmc_accepts += accepted

        self.params, accepted = hmc_trans(self.params, logpost,
                                          logpost_grad, sigmoid_trans,
                                          25, exponential(0.05))
        Binomial_test.hmc_rejects += 1 - accepted
        Binomial_test.hmc_accepts += accepted


    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpost(drifts):
            if any(drifts < self.min_drift) or any(drifts > self.max_drift):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + normpdfln(child.params,
                                          self.rho*root.params, drifts)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + normpdfln(self.params,
                                           self.rho*self.init_mean, drifts)

        self._drift = slice_sample(self._drift,
                                   logpost,
                                   step_out=True,
                                   compwise=True)

    def logprob(self, x):
        return log_binomial_pdf(x[:,0], x[:,1], self.params)


    def complete_logprob(self):
        return self.logprob(self.get_data())
        

    ## def _log_binomial_likelihood(self, b, d,
    ##                              cn_n, cn_r, cn_v,
    ##                              mu_n, mu_r, mu_v,
    ##                              cellular_frequency):  
    ##     f = cellular_frequency
    ##     t = self.params.tumour_content
        
    ##     p_n = (1 - t) * cn_n
    ##     p_r = t * (1 - f) * cn_r
    ##     p_v = t * f * cn_v
        
    ##     norm_const = p_n + p_r + p_v
        
    ##     p_n = p_n / norm_const
    ##     p_r = p_r / norm_const
    ##     p_v = p_v / norm_const
        
    ##     mu = p_n * mu_n + p_r * mu_r + p_v * mu_v
        
    ##     return log_binomial_pdf(b, d, mu)
