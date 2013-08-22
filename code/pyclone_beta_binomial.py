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

def digamma(x):
    return scipy.special.psi(x)

    
def log_sum_exp(log_X):
    '''
    Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])
    
    Numerically safer than naive method.
    '''
    max_exp = max(log_X)
    
    if isinf(max_exp):
        return max_exp
    
    total = 0

    for x in log_X:
        total += exp(x - max_exp)
    
    return log(total) + max_exp

def log_sum_exp_prod(log_X, Y):

    max_exp = max(log_X)

    if isinf(max_exp):
        return max_exp

    total = 0

    for x, y in zip(log_X, Y):
        total += exp(x - max_exp) * y

    return total * exp( max_exp )
    

def log_factorial(n):
    return gammaln(n + 1.0)

def log_binomial_coefficient(n, x):
    return log_factorial(n) - log_factorial(x) - log_factorial(n - x)

def log_beta(a, b):
    return log_gamma(a) + log_gamma(b) - log_gamma(a + b)

def log_beta_binomial_pdf(x, n, a, b):
    return log_binomial_coefficient(n, x) + log_beta(a + x, b + n - x) - log_beta(a, b)

def log_binomial_likelihood(b, d, cn_n, cn_r, cn_v,
                             mu_n, mu_r, mu_v, cellular_frequency, precision):  
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

        param_a = mu * precision
        
        param_b = (1 - mu) * precision
        
        return log_beta_binomial_pdf(b, d, param_a, param_b)


def log_pyclone(data, params, precision):
        ll = []       
        for cn_n, cn_r, cn_v, mu_n, mu_r, mu_v, log_pi  in zip(data[2],
                                                               data[3],
                                                               data[4],
                                                               data[5],
                                                               data[6],
                                                               data[7],
                                                               data[8]):
            temp = log_pi + log_binomial_likelihood(data[0],
                                                    data[1],
                                                    cn_n,
                                                    cn_r,
                                                    cn_v,
                                                    mu_n,
                                                    mu_r,
                                                    mu_v,
                                                    params, precision)
            ## if isnan(temp):
            ##     return "params ", params
            ## ll.append(temp)
        return log_sum_exp(ll)


class PyClone_Beta_Binomial(Node):

    init_mean    = 1.0
    min_drift    = 0.01
    max_drift    = 1.0
    min_balpha   = 0.1
    max_balpha   = 5
    min_bbeta    = 0.1
    max_bbeta    = 5
    hmc_accepts  = 1
    hmc_rejects  = 1

    def __init__(self, parent=None, dims=1, tssb=None, precision = 1.0
                 balpha = 0.5, bbeta = 0.5):
        super(PyClone_Beta_Binomial, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
           self.dims    = dims
           self._balpha = balpha*ones((1))
           self._bbeta  = bbeta*ones((1))
           self.params  = 0.99*self.init_mean*ones((1))
            
        else:
            self.dims   = parent.dims
            self.params = parent.params * boundbeta(balpha, bbeta)


    def balpha(self):
        if self.parent() is None:
            return self._balpha
        else:
            return self.parent().balpha()
        
    def bbeta(self):
        if self.parent() is None:
            return self._bbeta
        else:
            return self.parent().bbeta()    
        
    ## def sample(self, args):
    ##     num_data = args['num_data'] if args.has_key('num_data') else 1
    ##     return self.param[0] + self.param[1]*randn(num_data, 1)
    
    def resample_params(self):
        data       = self.get_data()
        num_data   = len(data)
        balphas    = self.balpha()
        bbetas     = self.bbeta()

        def get_constraint():
            if self.parent() is None:
                upper = float('inf')*ones(self.dims)
            else:
                upper = self.parent().params

            if  len(self.children()) < 1:
                lower = float('-inf')*ones((self.dims))
            else:
                child_params = zeros((len(self.children()), self.dims))
                for index, child in enumerate(self.children()):
                    child_params[index,:] = child.params;
                lower = child_params.max(0);

            return lower, upper

        lower, upper = get_constraint()

        def logpost_trans(params):

            orig_params = sigmoid_trans(params, 'orig_x')

            if any(orig_params < lower) or any(orig_params > upper):
                return float('-inf')

            if self.parent() is None:
                llh = betapdfln(orig_params/self.init_mean, balphas, bbetas)
            else:
                llh = betapdfln(orig_params/self.parent().params, \
                                balphas, bbetas)

            ll = 0
            for dd in data:
                ll += ( log_pyclone(dd, orig_params) )
                
            llh = llh + ll

            for child in self.children():
                llh = llh + betapdfln(child.params/orig_params, balphas, bbetas)

            det_jacob = sigmoid_trans(orig_params, 'det_jacob')

            llh = llh + det_jacob

            return llh


        trans_params = sigmoid_trans(self.params, 'tran_x')
        temp = slice_sample(trans_params, logpost_trans,
                            step_out=True, compwise=True)
        self.params = sigmoid_trans(temp, 'orig_x')
        
                                
        def logpost(params):
            if self.parent() is None:
                llh = betapdfln(params/self.init_mean, balphas, bbetas)
            else:
                llh = betapdfln(params/self.parent().params, \
                                balphas, bbetas) 
            ## pyclone node

            ll = 0
            for dd in data:
                ll += ( log_pyclone(dd,params) )
                
            llh = llh + ll
        
            for child in self.children():
                llh = llh + betapdfln(child.params/params, balphas, bbetas)

            det_jacobln = sigmoid_trans(params, 'det_jacob')
                
            return llh

        def logpost_grad(params):
            
            if self.parent() is None:
                grad = (balphas-1.0) / params - (bbetas - 1.0) /\
                      (self.init_mean - params)
            else:
                grad = (balphas-1.0) / params - (bbetas - 1.0) /\
                      (self.parent().params - params)

            g = 0
            for dd in data:
                g += log_pyclone_grad(dd, params)
                
            grad = grad + g

                        
            for child in self.children():
                grad = grad + (1.0-balphas) / params + (bbetas - 1.0) * \
                  child.params / (params *(params - child.params))
               ## print child.params
            return grad

        ## func   = logpost(self.params)
        ## eps    = 1e-4
        ## mygrad = logpost_grad(self.params)
        ## fdgrad = zeros(self.params.shape)
        ## for d in range(len(self.params)):
        ##     mask      = zeros(self.params.shape)
        ##     mask[d]   = 1
        ##     fdgrad[d] = (logpost(self.params + eps*mask) - \
        ##                  logpost(self.params - eps*mask))/(2*eps)
        ## print "MYGRAD: ", mygrad
        ## print "FDGRAD: ", fdgrad
        ## print "ERROR: ", norm(mygrad-fdgrad)/norm(mygrad+fdgrad)
        ## print "PARAMS: ", self.params
        ## print "DATA: ", data.shape

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
        ## print "ERROR: ", norm(mygrad-fdgrad)/norm(mygrad+fdgrad)
        ## print "PARAMS: ", self.params
        ## print "DATA: ", data.shape
        ##print "FUNC: ", func
        ## if rand() < 0.1:
        ##     temp = slice_sample(self.params, logpost, step_out=True, compwise=True)
        ## else:
        ##   self.params, accepted = hmc(self.params, logpost,
        ##                               logpost_grad, 25,
        ##                               exponential(0.005))
        ## Binomial_test.hmc_rejects += 1 - accepted
        ## Binomial_test.hmc_accepts += accepted

        ## self.params, accepted = hmc_trans(self.params, logpost,
        ##                                   logpost_grad, sigmoid_trans,
        ##                                   25, exponential(0.005))
        ## PyClone_test.hmc_rejects += 1 - accepted
        ## PyClone_test.hmc_accepts += accepted


    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpost_balpha(balphas):
            if any(balphas < self.min_balpha) or any(balphas > self.max_balpha):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + betapdfln(child.params/root.params, \
                                          balphas, \
                                          root.bbeta())
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + normpdfln(self.params/self.init_mean, \
                                           balphas, \
                                           self.bbeta())

        self._balpha = slice_sample(self._balpha,
                                   logpost_balpha,
                                   step_out=True,
                                   compwise=True)

        def logpost_bbeta(bbetas):
            if any(bbetas < self.min_balpha) or any(bbetas > self.max_balpha):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + betapdfln(child.params/root.params, \
                                          root.balpha(), \
                                          bbetas)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + normpdfln(self.params/self.init_mean, \
                                           self.balpha(), \
                                           bbetas)

        self._bbeta = slice_sample(self._bbeta,
                                   logpost_bbeta,
                                   step_out=True,
                                   compwise=True)

    def logprob(self, x):
        ll = 0
        for dd in x:
            ll += log_pyclone(dd, self.params)
        return ll


    def complete_logprob(self):
        return self.logprob(self.get_data())
        

