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
    
def trans(x, opt):
    if opt == 'tran_x':
        return numpy.array([logit(x[0]),x[1]]);
    if opt == 'orig_x':
        return numpy.array([sigmoid(x[0]), x[1]]);
    if opt == 'det_jacob':
        return numpy.log(x[0]) + log(1.0-x[0]);
    if opt == 'tran_grad':
        return numpy.array([x[0]*(1.0-x[0]), 1.0]);
    if opt == 'jacob_grad':
        return numpy.array([1.0-2*x[0], 0.0]);

class Gaussian_test(Node):

    init_mean    = 1.0
    init_std     = 0.01
    min_drift    = 0.01
    max_drift    = 1.0
    hmc_accepts  = 1
    hmc_rejects  = 1
    min_balpha   = 0.1
    max_balpha   = 5
    min_bbeta    = 0.1
    max_bbeta    = 5

    def __init__(self, parent=None, dims=1, tssb=None, drift=1.0, balpha = 0.5, bbeta = 0.5 ):
        super(Gaussian_test, self).__init__(parent=parent, tssb=tssb)

        if parent is None:
            self.dims    = dims
            self._drift  = drift*ones((1))
            self._balpha = balpha*ones((1))
            self._bbeta  = bbeta*ones((1))
            self.params  = numpy.array([0.99 * self.init_mean, self.init_std ])

        else:
            self.dims   = parent.dims
            # Gaussian random walk prior for std, beta prior for mean
            p0 = parent.params[0] * boundbeta(balpha, bbeta)
            p1 = parent.params[1] + self.drift()*randn(1)
            self.params = numpy.array([p0,p1])

    def drift(self):
        if self.parent() is None:
            return self._drift
        else:
            return self.parent().drift()

    def balpha(self):
        if self.parent() is None:
            return self._balpha
        else:
            return self.parent().balpha()

    def bbeta(self):
        if self.parent() is None:
            return self._balpha
        else:
            return self.parent().bbeta()
        
    def sample(self, args):
        num_data = args['num_data'] if args.has_key('num_data') else 1
        return self.param[0] + self.param[1]*randn(num_data, 1)
    
    def resample_params(self):
        data       = self.get_data()
        num_data   = data.shape[0]
        drifts     = self.drift()
        balphas    = self.balpha()
        bbetas     = self.bbeta()
        
        def logpost(params):
            if self.parent() is None:
                llh = normpdfln(params[1], self.init_std, drifts) + \
                      betapdfln(params[0]/self.init_mean, balphas, bbetas)
            else:
                llh = normpdfln(params[1], self.parent().params[1], drifts) + \
                      betapdfln(params[0]/self.parent().params[0], \
                                balphas, bbetas)              
            ## Gaussian node
            llh = llh + normpdfln(data, params[0], params[1])
            
            for child in self.children():
                llh = llh + normpdfln(child.params[1], params[1], drifts) + \
                      betapdfln(child.params[0]/params[0], balphas, bbetas)
            return llh

        def logpost_grad(params):

            grad = zeros(self.params.shape)
            if self.parent() is None:
                gp1 = -(params[1] -self.init_std)/drifts**2
                gp0 = (balphas-1.0) / params[0] - (bbetas - 1.0) /\
                      (self.init_mean - params[0])
            else:
                gp1 = -(params[1] - self.parent().params[1])/drifts**2
                gp0 = (balphas-1.0) / params[0] - (bbetas - 1.0) /\
                      (self.parent().params[0] - params[0])
            grad[0] = gp0
            grad[1] = gp1
            
            g1 = sum(data - params[0]) / params[1]**2
            g2 = sum((data - params[0])**2) / params[1]**3 - num_data/params[1] 
            grad = grad + numpy.array([g1, g2])
            
            for child in self.children():
                gg0 = (1.0-balphas) / params[0] + (bbetas - 1.0) * \
                  child.params[0] / (params[0] *(params[0] - child.params[0]))
                gg1 = (child.params[1] - params[1])/drifts**2
                grad[0] = grad[0] + gg0
                grad[1] = grad[1] + gg1
            return grad

        ## ## Uncomment to check normal gradient
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
        ## print "PARAMS: ", self.params
        ## print "DATA: ", data.shape
        ## print norm(mygrad-fdgrad)/norm(mygrad+fdgrad)

        ## ## Uncomment to check transformed gradient
        ## tran_param = trans(self.params, 'tran_x')
        ## eps    = 1e-4
        ## mygrad = logpost_grad(self.params) * trans(self.params, 'tran_grad') +\
        ##          trans(self.params, 'jacob_grad')
        ## fdgrad = zeros(self.params.shape)
        ## for d in range(len(self.params)):
        ##     mask      = zeros(self.params.shape)
        ##     mask[d]   = 1
        ##     v1 = trans((tran_param + eps*mask),'orig_x')
        ##     v2 = trans((tran_param - eps*mask),'orig_x')
        ##     fdgrad[d] = (logpost( v1 ) + trans( v1, 'det_jacob') - \
        ##                  logpost( v2 ) - trans( v2, 'det_jacob') ) / (2*eps)
        ## print "MYGRAD: ", mygrad
        ## print "FDGRAD: ", fdgrad
        ## print "ERROR: ",  norm(mygrad-fdgrad)/norm(mygrad+fdgrad)
        ## print "PARAMS: ", self.params
        ## print "Trans Params ", tran_param
        ## print "DATA: ", data.shape
        
        #if rand() < 0.1:
        #   self.params = slice_sample(self.params, logpost, step_out=True, compwise=True)
        #else:
        self.params, accepted = hmc_trans(self.params, logpost,
                                          logpost_grad, trans,
                                          25, exponential(0.005))
        Gaussian_test.hmc_rejects += 1 - accepted
        Gaussian_test.hmc_accepts += accepted
        
        ## self.params, accepted = hmc(self.params, logpost,
        ##                             logpost_grad, 25,
        ##                             exponential(0.008))
        ## Gaussian_test.hmc_rejects += 1 - accepted
        ## Gaussian_test.hmc_accepts += accepted


    def resample_hypers(self):
        if self.parent() is not None:
            raise Exception("Can only update hypers from root!")

        def logpost(drifts):
            if any(drifts < self.min_drift) or any(drifts > self.max_drift):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + normpdfln(child.params[1],
                                          root.params[1], drifts)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + normpdfln(self.params[1],
                                           self.init_std, drifts)

        self._drift = slice_sample(self._drift,
                                   logpost,
                                   step_out=True,
                                   compwise=True)

        def logpost_balpha(balphas):
            if any(balphas < self.min_balpha) or any(balphas > self.max_balpha):
                return -inf
            def loglh(root):
                llh = 0.0
                for child in root.children():
                    llh = llh + betapdfln(child.params[0]/root.params[0], \
                                          balphas, \
                                          root.bbeta())
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + betapdfln(self.params[0]/self.init_mean, \
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
                    llh = llh + betapdfln(child.params[0]/root.params[0], \
                                          root.balpha(), \
                                          bbetas)
                    llh = llh + loglh(child)
                return llh
            return loglh(self) + betapdfln(self.params[0]/self.init_mean, \
                                           self.balpha(), \
                                           bbetas)

        self._bbeta = slice_sample(self._bbeta,
                                   logpost_bbeta,
                                   step_out=True,
                                   compwise=True)

    def logprob(self, x):
        return  normpdfln(x, self.params[0], self.params[1])


    def complete_logprob(self):
        return self.logprob(self.get_data())
    
    
        
