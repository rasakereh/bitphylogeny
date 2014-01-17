import sys
import scipy.stats

from util         import *
from pylab        import *
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
        if weights[index] == 1.0:
            res[index] = l1
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
            

    if rate1 == 0.0 or rate1 == 1.0:
        rate1 = 1.0/8.0

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
    max_std   = 5.0
    min_branch_length = 0.0
    max_branch_length = 8.0
    max_param    = 20.0

    def __init__(self, parent=None, dims=1, tssb=None, mu = 3.0, mu0 = 3.0,
                 ratemat = (32.0/7.0)*array([[-1.0/8.0,1.0/8.0],[7.0/8.0,-7.0/8.0]]),
                 branch_length = 1.0, std = 1e-2 ):
        super(Logistic, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self.depth  = 0.0
            self.root_bias = mu0*ones(1)
            self.base_value = mu*ones(1)
            self.std = std*ones(1)
            self.ratemat = ratemat
            self.branch_length = branch_length*ones(1)
            self.params = rootprior(-(self.base_value+self.root_bias) *ones(dims),
                                    self.base_value + self.root_bias,
                                    self.std, self.ratemat)
            
        else:
            self.dims   = parent.dims
            self.depth  = parent.depth + 1.0
            self.root_bias = parent.root_bias
            self.base_value = parent.base_value
            self.std = parent.std
            self.ratemat = parent.ratemat
            self.branch_length = parent.branch_length
            self.params = mixlaplaceprior(parent.params, self.base_value,
                                          self.std, self.depth, self.ratemat)

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
            return self._mu
        else:
            return self.parent().mu_caller()
        
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):
        
        data       = self.get_data()
        counts     = sum(data, axis=0)
        num_data   = data.shape[0]
        mm         = (self.base_value, -self.base_value)
        std        = self.std
        depth      = self.depth
        ratemat    = self.ratemat

        
        def logpost(params):

            if any(params >= self.max_param) or any(params <= -self.max_param):
                llh = float('-inf')

            if self.parent() is None:
                llh = rootpriorpdfln(params, (self.base_value + self.root_bias,
                                              -self.base_value - self.root_bias),
                                     std,
                                     ratemat)
            else:
                llh = mixlaplacepdfln(params, mm, std, self.parent().params, depth,
                                      ratemat)

            llh = llh + sum(counts*sigmoidln(params)) + \
              sum((num_data-counts)*sigmoidln(-params))

            for child in self.children():
                llh = llh + mixlaplacepdfln(child.params, mm, std, params,depth+1.0,
                                            ratemat)

            return llh

        
        self.params = slice_sample(self.params, logpost,
                                   sigma = 5.0, step_out=True,
                                   compwise=False)
        
        self._cache_sigmoidln()

        def logpost_branch_length(branch_length):

            if any(branch_length < self.min_branch_length) or \
              any(branch_length > self.max_branch_length):
                return float('-inf')

            params = self.params
            
            if self.parent() is None:
                parent_params = zeros(len(self.params))
                ratemat = estimate_ratemat(parent_params, params, branch_length)
                llh = rootpriorpdfln(params, (self.base_value + self.root_bias,
                                              -self.base_value - self.root_bias),
                                     std, ratemat)
            else:
                parent_params = self.parent().params    
                ratemat = estimate_ratemat(parent_params, params, branch_length)
                llh = mixlaplacepdfln(params, mm, std, parent_params, depth,
                                      ratemat)

            return llh

        self.branch_length = slice_sample(self.branch_length,
                                          logpost_branch_length,
                                          sigma = 1.0, step_out=True,
                                          compwise=False)

        if self.parent() is None:
            self.ratemat = estimate_ratemat(zeros(len(self.params)), self.params,
                                            self.branch_length)
        else:
            self.ratemat = estimate_ratemat(self.parent().params, self.params,
                                            self.branch_length)


            
    def resample_hypers(self):

        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        ## def estimate_ratemat(root):
        ##     (w,nodes) = root.tssb.get_mixture()
        ##     params = []
        ##     for idx, nn in enumerate(nodes):
        ##         pp = nn.params
        ##         pp = (pp > ones(len(pp)))
        ##         params.append(pp)

        ##     params = array(params)
        ##     rate1 = mean(sum(params, axis=1).astype(float64)/len(pp))

        ##     if rate1 == 0.0 or rate1 == 1.0:
        ##         rate1 = 1.0/8.0

        ##     rate2 = 1-rate1
        ##     #scalar = 5.0/(2*rate1*rate2)
        ##     scalar = 1.0
        ##     ratemat = scalar*array([[-rate1,rate1],[rate2,-rate2]])
            
        ##     return ratemat

        ## self.ratemat = estimate_ratemat(self)
         
                
        def logpost_mu(mu):
            if any(mu < self.min_mu) or any(mu > self.max_mu):
                return float('-inf')
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + mixlaplacepdfln(child.params, (mu, -mu), self.std, \
                                          root.params, child.depth, self.ratemat  )
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + rootpriorpdfln(self.params,
                                                (mu+self.root_bias,
                                                 - mu-self.root_bias ), self.std, \
                                                self.ratemat) 

        self.base_value = slice_sample(self.base_value, logpost_mu,
                                       sigma = 1.0, step_out=True, compwise=True)


        def logpost_mu0(mu0):
            if any(mu0 < self.min_mu0) or any(mu0 > self.max_mu0):
                return float('-inf')
           
            return rootpriorpdfln(self.params, (self.base_value + mu0,
                                                -self.base_value - mu0),
                                   self.std, self.ratemat) 

        self.root_bias = slice_sample(self.root_bias, logpost_mu0,
                                      sigma = 1.0, step_out=True, compwise=True)

        def logpost_std(est_std):
            if any(est_std < self.min_std) or any(est_std > self.max_std):
                return float('-inf')
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + mixlaplacepdfln(child.params,
                                                (self.base_value, -self.base_value),
                                                est_std, root.params,
                                                child.depth, self.ratemat  )
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + rootpriorpdfln(self.params,
                                                (self.base_value + self.root_bias,
                                                 -self.base_value - self.root_bias),
                                                est_std, self.ratemat) 

        self.std = slice_sample(self.std, logpost_std,
                                sigma = 1.0, step_out=True, compwise=True)


    def logprob(self, x):
        return sum(x*self._sigln + (1.0-x)*self._negsigln)

    def complete_logprob(self):
        return self.logprob(self.get_data())

    def parameter_log_prior(self):

        if self.parent() is None:
            llh = rootpriorpdfln(self.params,
                                 (self.base_value+self.root_bias,
                                  -self.base_value-self.root_bias),
                                 self.std, self.ratemat)
        else:
            llh = mixlaplacepdfln(self.params,(self.base_value,-self.base_value),
                                  self.std, self.parent().params, self.depth,
                                  self.ratemat)

        for child in self.children():
            llh = llh + mixlaplacepdfln(child.params,
                                        (self.base_value,-self.base_value),
                                        self.std, self.params,
                                        self.depth + 1.0, self.ratemat)
       
        return llh
        
    
        
