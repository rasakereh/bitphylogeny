import sys
import scipy.stats

from util         import *
from pylab        import *
from numpy        import *
from numpy.random import *
from node         import *
from scipy.stats.distributions import beta

def truncbeta(a, b):
    res = numpy.random.beta(a,b)
    for index, x in enumerate(res):
        if x == 1.0:
            res[index] = 0.9999
        elif isnan(x):
            res[index] = 0.0001
    return res

def mixbetapdfln(params, a, b, parentprobs, mutcut, mix=0.5):
    parentmut = numpy.array(parentprobs>mutcut, dtype='int')
    parentnotmut = numpy.array(parentprobs<=mutcut, dtype='int')
    res = parentmut*betapdfln(params, a, b)
    res = res + parentnotmut*log(mix*beta.pdf(params, a, b) + (1.0-mix)*beta.pdf(params,b,a))
    return res

def mixbeta(a, b, parentprobs, mutcut, mix):
    parentmut = numpy.array(parentprobs>mutcut, dtype='int')
    parentnotmut = numpy.array(parentprobs<=mutcut, dtype='int')
    res = parentmut*truncbeta(a,b)
    u = uniform(size = parentprobs.size)
    highcomponent = numpy.array(u<mix,dtype='int')
    lowcomponent = numpy.array(u>=mix,dtype='int')
    res = res + parentnotmut * highcomponent * truncbeta(a, b)
    res = res + parentnotmut * lowcomponent * truncbeta(b,a)
    return res

class SNVBernoulli(Node):

    #mutcut        = 0.1
    min_mutcut    = 0.1
    max_mutcut    = 0.9
    #alpha_base   = 50 # 50 soapy; 10 irrupting (bernoulli 8)
    #beta_base    = 1.1
    alpha_base   = 80.0 # 300.0 
    beta_base    = 5.0 #20.0
    hmc_accepts  = 1
    hmc_rejects  = 1
    
    def __init__(self, parent=None, dims=1, tssb=None, initial_snvs=0, mutcut=0.6, prior_depth=0):
        super(SNVBernoulli, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims   = dims
            self.init_mean = initial_snvs
            self._mutcut  = mutcut*ones(self.dims)
            #self.params = truncbeta(self.alpha_base**((self.init_mean-self.mutcut())/abs(self.init_mean-self.mutcut())),self.beta_base**((self.init_mean-self.mutcut())/abs(self.init_mean-self.mutcut())))
            #TODO make nicer, i.e., sample from distribution:
            
            self.params = self.init_mean
            
            self.prior_depth = prior_depth
            self.depth = 1
            self.mix = 0.2 #* numpy.array(self.prior_depth <= self.depth, dtype='int') + 0.02 * numpy.array(self.prior_depth > self.depth, dtype='int')
            
        else:
            self.dims   = parent.dims
            self.prior_depth = parent.prior_depth

            #frequency prior:
            self.depth = parent.depth + 1
            #self.alpha_base = self.alpha_base * numpy.array(self.prior_depth <= self.depth, dtype='int') + self.beta_base * numpy.array(self.prior_depth > self.depth, dtype='int')
            
            #TODO remove hard coded mixing parameters:
            #mixing paramter for prior (transition kernel): 0.3 if prior depth is smaller than act depth, 0.1 otherwise 
            #self.mix = 0.2 #* numpy.array(self.prior_depth <= self.depth, dtype='int') + 0.02 * numpy.array(self.prior_depth > self.depth, dtype='int')
            self.mix = self.depth*0.2
            if(self.mix >= 1):
                self.mix = 0.8
            
            self.params = mixbeta( self.alpha_base*ones(self.dims), self.beta_base*ones(self.dims), parent.params, self.mutcut(), self.mix)
            #self.params = truncbeta(self.alpha_base**((parent.params-self.mutcut())/abs(parent.params-self.mutcut())), self.beta_base**((parent.params-self.mutcut())/abs(parent.params-self.mutcut())))
           
            if isinstance( initial_snvs, Iterable):
                self.params = initial_snvs
            
        self.alpha_base = ones(self.dims)*self.alpha_base
        self.beta_base = ones(self.dims)*self.beta_base
        self._cache_ln()
    
    def _cache_ln(self):
        assert any(self.params[newaxis,:]<1) 
        assert any(self.params[newaxis,:]>0)
        self._ln    = log(self.params[newaxis,:])
        self._negln = log(1-self.params[newaxis,:])

    def mutcut(self):
        if self.parent() is None:
            return self._mutcut
        else:
            return self.parent().mutcut()    
    
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return rand(num_data, self.dims) < sigmoid(self.params[newaxis,:])
    
    def resample_params(self):

        data       = self.get_data()
        counts     = sum(data[:,[2,3]], axis=0)
        num_data   = data.shape[0]
        mutcut      = self.mutcut()
        
        def logpost(params):
            if any(params>=1) or any(params<=0):
                llh = -sys.maxint
                return llh
            if self.parent() is None:
                #llh = sum(betapdfln(params, self.alpha_base**((self.init_mean-self.mutcut())/abs(self.init_mean-self.mutcut())), self.beta_base**((self.init_mean-self.mutcut())/abs(self.init_mean-self.mutcut()))))
                llh = sum(mixbetapdfln(params, self.alpha_base, self.beta_base, self.init_mean, self.mutcut(), self.mix))
                #llh = sum(betapdfln(params, 1.0/10000000, 1.0/10000000))
            else:
                #llh = sum(betapdfln(params, self.alpha_base**((self.parent().params-self.mutcut())/abs(self.parent().params-self.mutcut())), self.beta_base**((self.parent().params-self.mutcut())/abs(self.parent().params-self.mutcut()))))
                llh = sum(mixbetapdfln(params, self.alpha_base, self.beta_base, self.parent().params, self.mutcut(), self.mix))
            for i in range(num_data):
                #llh = llh + data[i][2]*sigmoidln(params[data[i][0]-1]) + (1.0-data[i][2])*sigmoidln(params[data[i][0]-1]) + data[i][3]*sigmoidln(params[data[i][1]-1]) + (1.0-data[i][3])*sigmoidln(params[data[i][1]-1]) + log(1/float(self.dims*(self.dims-1)))
                llh = llh + data[i][2]*log(params[data[i][0]-1]) + (1.0-data[i][2])*log(params[data[i][0]-1]) + data[i][3]*log(params[data[i][1]-1]) + (1.0-data[i][3])*log(params[data[i][1]-1]) + log(1/float(self.dims*(self.dims-1)))
            for child in self.children():
                #llh = llh + normpdfln( child.params, params, drifts)
                #llh = llh + sum(betapdfln( child.params, self.alpha_base**((params-self.mutcut())/abs(params-self.mutcut())), self.beta_base**((params-self.mutcut())/abs(params-self.mutcut()))))
                llh = llh + sum(mixbetapdfln(child.params, self.alpha_base, self.beta_base, params, self.mutcut(), self.mix))
                
            #just for TESTING; penalization of parameters between 0.2 and 0.8:
            x = sum(numpy.array([abs(params-0.5)<0.4]))
            llh = llh - x*300  
                
            return llh

        #def logpost_grad(params):
        #    if any(params>=1) or any(params<=0):
        #        grad = -sys.maxint
        #        return grad
        #    if self.parent() is None:
        #        grad = (self.init_mean**(1/2) -1)/params + (1/mutcut -1)/(1-params)
        #    else:
        #        grad = (self.parent().params**(1/2) -1)/params + (mutcut -1)/(1-params)
            #grad = grad + counts*(1.0-probs) - (num_data-counts)*probs
            #how about the 1/(N*(N-1))? does the gradient remove this?
        #    for i in range(num_data):
        #        grad[data[i][0]-1] = grad[data[i][0]-1] + data[i][2]/(params[data[i][0]-1]) + (data[i][2]-1)/(1-params[data[i][0]-1])
        #        grad[data[i][1]-1] = grad[data[i][1]-1] + data[i][3]/(params[data[i][1]-1]) + (data[i][3]-1)/(1-params[data[i][1]-1])
        #    for child in self.children():
        #        #grad = grad + (child.params - params)/drifts**2
        #        grad = grad + (params**(1/2) - 1)/child.params + (mutcut -1)/(1-child.params)
        #    return grad

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
        for i in range(20):
            self.params = slice_sample(self.params, logpost, step_out=True, compwise=False)
        #else:
            #self.params, accepted = hmc(self.params, logpost, logpost_grad, 25, exponential(0.001))
            #SNVBernoulli.hmc_rejects += 1 - accepted
            #SNVBernoulli.hmc_accepts += accepted
        
        self._cache_ln()

    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpostcut(mutcuts):
            if any(mutcuts < self.min_mutcut) or any(mutcuts > self.max_mutcut):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + sum(mixbetapdfln(child.params, self.alpha_base, self.beta_base, root.params, mutcuts, self.mix))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(mixbetapdfln(self.params, self.alpha_base, self.beta_base, self.init_mean, mutcuts, self.mix))
        
        def logpostalpha(alpha):
            if any(alpha <= 0):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + sum(mixbetapdfln(child.params, alpha, self.beta_base, root.params, self.mutcut(), self.mix))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(mixbetapdfln(self.params, alpha, self.beta_base, self.init_mean, self.mutcut(), self.mix))
        
        def logpostbeta(beta):
            if any(beta <= 0):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + sum(mixbetapdfln(child.params, self.alpha_base, beta, root.params, self.mutcut(), self.mix))
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + sum(mixbetapdfln(self.params, self.alpha_base, beta, self.init_mean, self.mutcut(), self.mix))
        
        self._mutcut  = slice_sample(self._mutcut, logpostcut, step_out=True, compwise=True)
        self.alpha_base = slice_sample(self.alpha_base, logpostalpha, step_out=True, compwise=True)
        self.beta_base = slice_sample(self.beta_base, logpostbeta, step_out=True, compwise=True)

    def logprob(self, x):
        x = transpose(x)
        res = x[2]*self._ln[0][x[0]-1] + (1.0-x[2])*self._negln[0][x[0]-1] + x[3]*self._ln[0][x[1]-1] + (1.0-x[3])*self._negln[0][x[1]-1] #+ log(1/float(self.dims*(self.dims-1)))
        assert all(res<0)
        #if any(res>-0.020):
        #    print "stop"
        #    print sum(numpy.array(res>-0.02,dtype='int'))
        return sum(res)

    def complete_logprob(self):
        return self.logprob(self.get_data())