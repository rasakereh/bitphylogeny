import sys
import scipy.stats

from numpy        import *
from numpy.random import *
from scipy.linalg import expm2

from bitphylogeny.util import *
from bitphylogeny.node import *

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

def boundprob(x):
    return (1.0-finfo(float64).eps)*(x-0.5)+0.5

def get_weights(parent_params,  ratemat):
    mapper = (parent_params > ones(len(parent_params)))
    weights = expm2(ratemat)
    m1 = weights[1,1]
    m2 = weights[0,1]
    if m1 == 0.0:
        m1 = boundprob(m1)
    if m2 == 0.0:
        m2 = boundprob(m2)
    return m1 * mapper + ~mapper*m2

def mixlaplaceprior(params, base_value, std, ratemat):
    weights = get_weights(params, ratemat)    
    mapper = ( rand(len(weights)) < weights )
    m1 = laplace(base_value*ones(len(weights)), std*ones(len(weights)))
    m2 = laplace(-base_value*ones(len(weights)), std*ones(len(weights)))
    return m1 * mapper + ~mapper * m2

def mixlaplacepdfln(x, m, std, parent_params, ratemat):
    weights = get_weights(parent_params, ratemat)
    #print weights, ratemat, parent_params 
    logw1  = log(weights)
    mapper = ( logw1 != 0 )
    lap1 = laplacepdfln(x, m[0], std)
    lap2 = laplacepdfln(x, m[1], std)
    ll1 = logw1[mapper] + lap1[mapper]
    logw2 = log(1-weights[mapper])
    ll2 = logw2 + lap2[mapper]
    ll  = sum(logsumexp(vstack([ll1,ll2]),0))
    ll  = ll + sum(lap1[~mapper])
    return ll

def rootprior(base_value, std, dims):
    m2 = laplace(-5.0*ones(dims), 0.5*ones(dims))
    return m2

def rootpriorpdfln(x, m, std, ratemat):
    ll = laplacepdfln(x, -5.0, 0.5) 
    return sum(ll)

class BitPhylogeny(Node):

    min_mu    = 2.0
    max_mu    = 7.0
    min_std   = 0.3
    max_std   = 20.0
    max_param    = 10.0
    min_branch_length = 0.0
    max_branch_length = 1.0

    def __init__(self, parent=None, dims=1, tssb=None, mu = 5.0,
                 ratemat = 32.0/7.0*array([[-1.0/8.0,1.0/8.0],[7.0/8.0,-7.0/8.0]]),
                 std = 1.0, mode = 'methylation'):
        super(BitPhylogeny, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.mode   = mode
            self.dims   = dims
            self.depth  = 0.0
            self._base_value = mu*ones(1)
            self._std = std*ones(1)
            self._ratemat = ratemat
            self.branch_length = zeros(1) 
            self.params = rootprior(self._base_value,
                                    self._std, dims)
            
        else:
            self.mode = parent.mode
            self.dims   = parent.dims
            self.depth  = parent.depth + 1.0
            self.branch_length = uniform(0,1)*ones(1)
            self.params = mixlaplaceprior(parent.params, self.mu_caller(),
                                          self.std_caller(),
                                          self.ratemat_caller()*self.branch_length)

            if any(abs(self.params)>self.max_param):
                for index, x in enumerate(self.params):
                    if abs(x)>self.max_param:
                        if x >0.0:
                            self.params[index] = self.max_param*0.99999
                        else :
                            self.params[index] = -self.max_param*0.99999
                            
                                   
        self._cache_sigmoidln()

    def _cache_sigmoidln(self):
        self._sigln    = sigmoidln(self.params[newaxis,:])
        self._negsigln = sigmoidln(-self.params[newaxis,:])

    def mu_caller(self):
        if self.parent() is None:
            return self._base_value
        else:
            return self.parent().mu_caller()

    def std_caller(self):
        if self.parent() is None:
            return self._std
        else:
            return self.parent().std_caller()

    def ratemat_caller(self):
        if self.parent() is None:
            return self._ratemat
        else:
            return self.parent().ratemat_caller()

    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):
        
        mode       = self.mode
        data       = self.get_data()
        counts     = sum(data, axis=0)
        num_data   = data.shape[0]
        mm         = (self.mu_caller(), -self.mu_caller())
        est_std    = self.std_caller()
        ratemat    = self.ratemat_caller()
        est_branch = self.branch_length
        
        def logpost(params):

            if any(params >= self.max_param) or any(params <= -self.max_param):
                llh = float('-inf')

            if self.parent() is None:
                llh = rootpriorpdfln(params, self.mu_caller(),
                                     est_std,
                                     ratemat)
            else:
                llh = mixlaplacepdfln(params, mm, est_std, self.parent().params,
                                      ratemat*est_branch)

            if mode == "methylation":
                llh = llh + sum(counts*sigmoidln(params)) + \
                      sum((num_data-counts)*sigmoidln(-params))
            elif num_data > 0:
                llh = llh + nansum(data*sigmoidln(params) + \
                                   (1.0 - data)*sigmoidln(-params))

            for child in self.children():
                llh = llh + mixlaplacepdfln(child.params, mm, est_std, params,
                                            ratemat*child.branch_length)

            return llh

        #for i in range(5):
        self.params = slice_sample(self.params, logpost,
                                   sigma = 5.0, step_out=True, compwise=False)      
        self._cache_sigmoidln()

        def logpost_branch_length(branch_length):

            if any(branch_length < self.min_branch_length) or \
              any(branch_length > self.max_branch_length):
                return float('-inf')

            params = self.params
            parent_params = self.parent().params    
            llh = mixlaplacepdfln(params, mm, est_std, parent_params,
                                  ratemat*branch_length)

            return llh
        if self.parent() is None:
        
            self.branch_length = zeros(1)
        
        else:

            self.branch_length = slice_sample(self.branch_length,
                                              logpost_branch_length,
                                              sigma = 1.0, step_out=True,
                                              compwise=False)

    def resample_hypers(self):

        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def estimate_ratemat(root):
            (w,nodes) = root.tssb.get_mixture()
            params = []
            for idx, nn in enumerate(nodes):
                pp = nn.params
                pp = (pp > ones(len(pp)))
                params.append(pp)

            params = array(params)
            rate1 = mean(sum(params, axis=1).astype(float64)/len(pp))

            if (rate1 == 0.0) or (rate1 == 1.0):
                rate1 = boundprob(rate1)
                
            rate2 = 1-rate1

            if self.mode == "methylation":
                scalar = 1.0/(2*rate1*rate2)
                ratemat = scalar*array([[-rate1,rate1],[rate2,-rate2]])
            else:
                scalar = 1.0/(rate1*rate2)
                ratemat = scalar*array([[-rate1,rate1],[0.0,0.0]])
                
            return ratemat

        
        self._ratemat = estimate_ratemat(self)
                
        def logpost_mu(mu):
            if any(mu < self.min_mu) or any(mu > self.max_mu):
                return float('-inf')
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + mixlaplacepdfln(child.params, (mu, -mu),
                                                child.std_caller(),
                                                root.params, 
                                                child.ratemat_caller()*child.branch_length)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) 

        self._base_value = slice_sample(self._base_value, logpost_mu,
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
                                                child.ratemat_caller()*child.branch_length)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) 

        self._std = slice_sample(self._std, logpost_std,
                                 sigma = 1.0, step_out=True, compwise=True)


    def logprob(self, x):
        if self.mode == "methylation":
            return sum(x*self._sigln + (1.0-x)*self._negsigln)
        else:
            return nansum(x*self._sigln + (1.0-x)*self._negsigln)

    def complete_logprob(self):
        return self.logprob(self.get_data())

    def parameter_log_prior(self):

        if self.parent() is None:
            llh = rootpriorpdfln(self.params,self.mu_caller(),
                                 self.std_caller(), self.ratemat_caller())
        else:
            llh = mixlaplacepdfln(self.params,(self.mu_caller(),-self.mu_caller()),
                                  self.std_caller(), self.parent().params,
                                  self.ratemat_caller()*self.branch_length)
        return llh
