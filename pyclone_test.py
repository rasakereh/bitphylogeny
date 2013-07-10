import sys
import scipy.stats

from util         import *
from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *
from math         import isinf

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

## def log_sum_exp(log_X):
##     '''
##     Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])
    
##     Numerically safer than naive method.
##     '''
##     max_exp = max(log_X)
    
##     if isinf(max_exp):
##         return max_exp
    
##     total = 0

##     for x in log_X:
##         total += exp(x - max_exp)
    
##     return log(total) + max_exp

def log_factorial(n):
    return gammaln(n + 1.0)

def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(x) - log_factorial(n - x)

def log_binomial_pdf(x, n, p):
    if p == 0:
        return float('-inf')
    
    if p == 1.0:
        return float('-inf')
    
    return log_binomial_coefficient(n, x) + x * log(p) + (n - x) * log(1.0 - p)

def log_binomial_likelihood(b, d, cn_n, cn_r, cn_v,
                             mu_n, mu_r, mu_v, cellular_frequency):  
        f = cellular_frequency
        t = 1
        
        p_n = (1 - t) * cn_n
        p_r = t * (1 - f) * cn_r
        p_v = t * f * cn_v
        
        norm_const = p_n + p_r + p_v
        
        p_n = p_n / norm_const
        p_r = p_r / norm_const
        p_v = p_v / norm_const
        
        mu = p_n * mu_n + p_r * mu_r + p_v * mu_v
        
        return log_binomial_pdf(b, d, mu)


def log_pyclone(data, params):
        ll = 0       
        for cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi  in zip(data[2],
                                                               data[3],
                                                               data[4],
                                                               data[5],
                                                               data[6],
                                                               data[7],
                                                               data[8]):
            ll+= log_pi + log_binomial_likelihood(data[0],
                                                    data[1],
                                                    cn_n,
                                                    cn_r,
                                                    cn_v,
                                                    mu_n,
                                                    mu_r,
                                                    mu_v,
                                                    params)
        
        return sum(ll)


def log_binomial_grad(b,d,cn_n,cn_r,cn_v,mu_n,mu_r,mu_v,params):
    f = params
    t = 1
    
    p_n = (1 - t) * cn_n
    p_r = t * (1 - f) * cn_r
    p_v = t * f * cn_v

    norm_const = p_n + p_r + p_v

    g1 = -t * cn_r + t * cn_v
    g2 = -t * cn_r
    g3 = t * cn_v

    G1 = - p_n / norm_const**2 * g1
    G2 = 1 / norm_const * g2 - p_r / norm_const**2 * g1
    G3 = 1 / norm_const * g3 - p_v / norm_const**2 * g1

    gg = G1 * mu_n + G2 * mu_r + G3 * mu_v

    P_n = p_n / norm_const
    P_r = p_r / norm_const
    P_v = p_v / norm_const
        
    mu = P_n * mu_n + P_r * mu_r + P_v * mu_v

    grad = b / mu - (d - b) / (1 - mu)

    grad = grad * gg
    
    return grad

def log_pyclone_grad(data, params):
    grad = 0
    for cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi  in zip(data[2],
                                                           data[3],
                                                           data[4],
                                                           data[5],
                                                           data[6],
                                                           data[7],
                                                           data[8]):
        grad += log_binomial_grad(data[0],data[1],cn_n,cn_r,
                                  cn_v,mu_n,mu_r,mu_v,params);

    return grad



class PyClone_test(Node):

    init_mean    = 1.0
    min_drift    = 0.01
    max_drift    = 1.0
    hmc_accepts  = 1
    hmc_rejects  = 1
    rho          = 0.6

    def __init__(self, parent=None, dims=1, tssb=None, drift=0.1):
        super(PyClone_test, self).__init__(parent=parent, tssb=tssb)

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
        num_data   = len(data)
        drifts     = self.drift()   
        
        def logpost(params):
            if self.parent() is None:
                llh = normpdfln(params, self.rho*self.init_mean, drifts)
            else:
                llh = normpdfln(params, self.rho*self.parent().params, drifts) 
            ## binomial node

            ll = 0
            for dd in data:
                ll += ( log_pyclone(dd,params) )
                
            llh = llh + ll
            for child in self.children():
                llh = llh + normpdfln( child.params, self.rho*params, drifts)
            return llh

        def logpost_grad(params):
            
            if self.parent() is None:
                grad = -(params-self.rho*self.init_mean)/drifts**2
            else:
                grad = -(params-self.rho*self.parent().params)/drifts**2

            g = 0
            for dd in data:
                g += log_pyclone_grad(dd, params)
                
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
        PyClone_test.hmc_rejects += 1 - accepted
        PyClone_test.hmc_accepts += accepted


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
        ll = 0
        for dd in x:
            ll += log_pyclone(dd, self.params)
        
        return ll


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
