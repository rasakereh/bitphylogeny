import sys
import scipy.stats

from util         import *
#from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *
from scipy.linalg import expm, expm2

def sigmoid(x):
    return 1.0/(1.0+exp(-x))
def invsigmoid(z):
    return log(z) - log(1.0-z)
def sigmoidln(x):
    return -log(1.0 + exp(-x))
def normpdfln(x, m, std):
    return -0.5*log(2*pi) - log(std) - 0.5*((x-m)/std)**2

def laplacepdfln(x,m,std):
    return -log(2*std) - abs(x - m)/std

def get_weights(parent_params, depth, ratemat):
    mapper = (parent_params > zeros(len(parent_params)))
    weights = expm2(ratemat)
    m1 = weights[1,1]
    m2 = weights[0,1]
    return m1 * mapper + ~mapper*m2

def mixlaplaceprior(params, base_value, std, depth, ratemat):

    weights = get_weights(params, depth, ratemat)    
    mapper = ( rand(len(weights)) < weights )
    m1 = laplace(base_value*ones(len(weights)), std*ones(len(weights)))
    m2 = laplace(-base_value*ones(len(weights)), std*ones(len(weights)))
    
    return m1 * mapper + ~mapper*m2

def mixlaplacepdfln(x, m, std, parent_params, depth, ratemat):

    weights = get_weights(parent_params, depth, ratemat)

    res = empty(len(x))
    for index , d in enumerate(x):
        l1 = log( weights[index] ) + laplacepdfln(d, m[0][0], std)
        if isnan(log(1-weights[index])):
            res[index] = laplacepdfln(d, m[0][0], std)
        else:
            l2 = log(1- weights[index]) + laplacepdfln(d, m[1][0], std)
            res[index] = logsumexp( (l1, l2) )
    return sum(res)

def rootprior(params, base_value, std, ratemat):
    #weights = dot(array([1,0]),expm2(ratemat))[1] * ones(len(params))
    weights = expm2(ratemat)[0,1]* ones(len(params))
    mapper = ( rand(len(weights)) < weights )
    m1 = laplace(base_value*ones(len(weights)), std*ones(len(weights)))
    m2 = laplace(-base_value*ones(len(weights)), std*ones(len(weights)))

    return m1 * mapper + ~ mapper*m2

def rootpriorpdfln(x, m, std, ratemat):
   
    weights = expm2(ratemat)[0,1] * ones(len(x))
    res = empty(len(x))
    for index , d in enumerate(x):
        l1 = log( weights[index] ) + laplacepdfln(d, m[0], std)  
        l2 = log(1- weights[index]) + laplacepdfln(d, m[1], std)
        res[index] = logsumexp( (l1, l2) )
    return sum(res)

def estimate_ratemat(parent_params, params, branch_length):

    pp = []
    m = len(parent_params)
    pp1 = (parent_params > zeros(m))
    pp.append(pp1)
    pp2 = (params > zeros(m))
    pp.append(pp2)
    
    pp = array(pp)
    rate1 = mean(sum(pp, axis=1).astype(float64)/m)
            

    if rate1 == 0.0:
        rate1 = 0.0 + numpy.finfo(numpy.float64).eps
    if rate1 == 1.0:
        rate1 = 1.0-numpy.finfo(numpy.float64).eps

    rate2 = 1-rate1

    scalar = 1.0/(2*rate1*rate2)
    
    ratemat = branch_length*scalar*array([[-rate1,rate1],[rate2,-rate2]])

    return ratemat


class Logistic(Node):

    min_mu0   = 2.0
    max_mu0   = 5.0
    min_mu    = 2.0
    max_mu    = 7.0
    min_std   = 1e-3
    max_std   = 30.0
    min_branch_length = 0.0
    max_branch_length = 8.0
    max_param    = 20.0

    def __init__(self, parent=None, dims=1, tssb=None, mu = 3.0, mu0 = 3.0,
                 ratemat = (32.0/7.0)*array([[-1.0/8.0,1.0/8.0],[7.0/8.0,-7.0/8.0]]),
                 branch_length = 1.0, std =0.5 ):
        super(Logistic, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self.depth  = 0.0
            self._root_bias = mu0*ones(1)
            self._base_value = mu*ones(1)
            self._std = std*ones(1)
            self.ratemat = ratemat
            self._branch_length = branch_length*ones(1)
            self.params = rootprior(-(self._base_value+self._root_bias) *ones(dims),
                                    self._base_value + self._root_bias,
                                    self._std, self.ratemat)
            
        else:
            self.dims   = parent.dims
            self.depth  = parent.depth + 1.0
            self.ratemat = parent.ratemat
            self.params = mixlaplaceprior(parent.params, self.mu_caller(),
                                          self.std_caller(), self.depth, self.ratemat)

            if any(abs(self.params)>self.max_param):
                for index, x in enumerate(self.params):
                    if abs(x)>self.max_param:
                        self.params[index] = parent.params[index]*0.99999
                                   
        self._cache_sigmoidln()

    def _cache_sigmoidln(self):
        self._sigln    = sigmoidln(self.params[newaxis,:])
        self._negsigln = sigmoidln(-self.params[newaxis,:])

    def mu_caller(self):
        if self.parent() is None:
            return self._base_value
        else:
            return self.parent().mu_caller()

    def mu0_caller(self):
        if self.parent() is None:
            return self._root_bias
        else:
            return self.parent().mu0_caller()

    def std_caller(self):
        if self.parent() is None:
            return self._std
        else:
            return self.parent().std_caller()

    def branch_caller(self):
        if self.parent() is None:
            return self._branch_length
        else:
            return self.parent().branch_caller()
        
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):
        
        data       = self.get_data()
        counts     = sum(data, axis=0)
        num_data   = data.shape[0]
        mm         = (self.mu_caller(), -self.mu_caller())
        std        = self.std_caller()
        depth      = self.depth
        
        if self.parent() is None:
            self.ratemat = estimate_ratemat(zeros(len(self.params)), self.params,
                                            self.branch_caller())
        else:
            self.ratemat = estimate_ratemat(self.parent().params, self.params,
                                            self.branch_caller())
        ratemat    = self.ratemat

        
        def logpost(params):

            if any(params >= self.max_param) or any(params <= -self.max_param):
                llh = float('-inf')

            if self.parent() is None:
                llh = rootpriorpdfln(params, (self.mu_caller() + self.mu0_caller(),
                                              -self.mu_caller() - self.mu0_caller()),
                                     std,
                                     ratemat)
            else:
                llh = mixlaplacepdfln(params, mm, std, self.parent().params, depth,
                                      ratemat)

            llh = llh + sum(counts*sigmoidln(params)) + \
              sum((num_data-counts)*sigmoidln(-params))

            for child in self.children():
                llh = llh + mixlaplacepdfln(child.params, mm, std, params, depth+1.0,
                                            ratemat)

            return llh

        
        self.params = slice_sample(self.params, logpost,
                                   sigma = 5.0, step_out=True,
                                   compwise=False)
        
        self._cache_sigmoidln()

            
    def resample_hypers(self):

        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpost_branch_length(branch_length):

            if any(branch_length < self.min_branch_length) or \
              any(branch_length > self.max_branch_length):
                return float('-inf')

            root_ratemat = estimate_ratemat(zeros(len(self.params)), self.params,
                                            branch_length)

            def loglh(root):
                llh = 0.0
                for child in root.children():
                    ratemat = estimate_ratemat(root.params,child.params,branch_length)
                    llh = llh + mixlaplacepdfln(child.params,
                                                (child.mu_caller(),
                                                 -child.mu_caller()),
                                                child.std_caller(),
                                                root.params, child.depth,
                                                ratemat)
                    llh = llh + loglh(child)
                return llh

            return loglh(self) + rootpriorpdfln(self.params,
                                                (self.mu_caller()+self.mu0_caller(),
                                                 -self.mu_caller()-self.mu0_caller()),
                                                self.std_caller(), root_ratemat) 

        self._branch_length = slice_sample(self._branch_length,
                                          logpost_branch_length,
                                          sigma = 1.0, step_out=True,
                                          compwise=False)

        def logpost_mu(mu):
            if any(mu < self.min_mu) or any(mu > self.max_mu):
                return float('-inf')
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + mixlaplacepdfln(child.params, (mu, -mu),
                                                child.std_caller(),
                                                root.params, child.depth,
                                                child.ratemat)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + rootpriorpdfln(self.params,
                                                (mu+self.mu0_caller(),
                                                 - mu-self.mu0_caller() ),
                                                self.std_caller(), 
                                                self.ratemat) 

        self._base_value = slice_sample(self._base_value, logpost_mu,
                                       sigma = 1.0, step_out=True, compwise=True)

        def logpost_mu0(mu0):
            if any(mu0 < self.min_mu0) or any(mu0 > self.max_mu0):
                return float('-inf')
           
            return rootpriorpdfln(self.params, (self.mu_caller() + mu0,
                                                -self.mu_caller() - mu0),
                                   self.std_caller(), self.ratemat) 

        self._root_bias = slice_sample(self._root_bias, logpost_mu0,
                                      sigma = 1.0, step_out=True, compwise=True)

        def logpost_std(est_std):
            if any(est_std < self.min_std) or any(est_std > self.max_std):
                return float('-inf')
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + mixlaplacepdfln(child.params,
                                                (child.mu_caller(),
                                                 -child.mu_caller()),
                                                est_std, root.params,
                                                child.depth, child.ratemat)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + rootpriorpdfln(self.params,
                                                (self.mu_caller() + self.mu0_caller(),
                                                 -self.mu_caller() - self.mu0_caller()),
                                                est_std, self.ratemat) 

        self._std = slice_sample(self._std, logpost_std,
                                sigma = 1.0, step_out=True, compwise=True)

       


    def logprob(self, x):
        return sum(x*self._sigln + (1.0-x)*self._negsigln)

    def complete_logprob(self):
        return self.logprob(self.get_data())

    def parameter_log_prior(self):

        if self.parent() is None:
            llh = rootpriorpdfln(self.params,
                                 (self.mu_caller()+self.mu0_caller(),
                                  -self.mu_caller()-self.mu0_caller()),
                                 self.std_caller(), self.ratemat)
        else:
            llh = mixlaplacepdfln(self.params,(self.mu_caller(),-self.mu_caller()),
                                  self.std_caller(), self.parent().params, self.depth,
                                  self.ratemat)

        for child in self.children():
            llh = llh + mixlaplacepdfln(child.params,
                                        (self.mu_caller(),-self.mu_caller()),
                                        self.std_caller(), self.params,
                                        self.depth + 1.0, self.ratemat)
       
        return llh
    
