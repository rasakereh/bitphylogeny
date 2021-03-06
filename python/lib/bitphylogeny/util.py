import os
import sys
import cPickle
import numpy
import numpy.random
import scipy.special
import scipy.stats
import scipy.io
from _abcoll import Iterable
import csv
import re

def object2label(obj,nodes):
    labels = []
    for idx, item in enumerate(obj):
        for ii, nn in enumerate(nodes):
            if item == nn:
                labels.append(ii)
    labels = numpy.array(labels)
    return(labels)

def cluster_with_maxpear(X):
    try:         
        import rpy2
    except:
        raise Exception('''Clustering with maximum PEAR requires rpy2 package. See http://rpy.sourceforge.net/rpy2.html.''')
    
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from rpy2.robjects.packages import importr    
    
    rpy2.robjects.numpy2ri.activate()
    
    importr('mcclust')

    r = robjects.r
    
    robjects.globalenv['X'] = r.matrix(X, nrow=X.shape[0])
    
    clusters = r(
      '''
          class(X) <- 'interger'
          clusters <- maxpear( comp.psm(X+1) )$cl;
       '''
    )

    return clusters


def write_traces2csv(file_name, traces, header=None):
    fh = open(file_name, 'wb')
    writer = csv.writer(fh)
    if header:
        writer.writerow(header)
    [writer.writerow(x) for x in traces]
    fh.close()

def write_traces2h5(file_name, traces):
    try:
        import h5py
    except:
        raise Exception("Need h5py")
    f = h5py.File(file_name, 'w')
    f['dset'] = traces
    f.close()



def summary_clone(best_cluster_labels, params_traces, node_depth_traces):
    
    member = []
    for idx1 , item1 in enumerate(numpy.unique(best_cluster_labels)):
        mm = []
        for idx2, bc in enumerate(best_cluster_labels):
            if bc == item1:
                mm.append(idx2)
                
        member.append(mm)

    genotype = []
    depth = []
    for idx, item in enumerate(numpy.unique(best_cluster_labels)):

        genotype.append(numpy.mean(
            numpy.mean(params_traces[:,member[idx],:], axis = 0),axis = 0))
        depth.append(numpy.mean(
            node_depth_traces[:,member[idx]]) )

    data_genotype = numpy.zeros((params_traces.shape[1], params_traces.shape[2]))
    data_depth = numpy.zeros((params_traces.shape[1]))

    for idx in range(params_traces.shape[1]):
        data_genotype[idx] = genotype[best_cluster_labels[idx]-1]
        data_depth[idx] = depth[best_cluster_labels[idx]-1]
    return([data_genotype , data_depth] ) 


def write_params_traces2file(trace, x, filepath,ftype='.npz'):
    n = trace.shape[0]
    m = n/x
    for i in range(x):
        tmp_trace = trace[i*m:(i+1)*m]
        filename = 'array_{0}'.format(i)
        numpy.savez(filepath+filename,
                                   tmp_trace)
       
    

## Functions for natural sorting
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
        
def load_params_traces(filepath,ftype='.npz'):
    f = [x for x in os.listdir(filepath) if x.endswith(ftype)]
    f.sort(key=natural_keys)
    t = []
    for file in f:
        temp = numpy.load( filepath+file )
        t.append(temp['arr_0'])
    return(numpy.vstack(t))

def log_sum_exp(log_X):
    '''
    Given a list of values in log space, log_X. Compute exp(log_X[0] + log_X[1] + ... log_X[n])
    
    Numerically safer than naive method.
    '''
    max_exp = max(log_X)
    
    if numpy.isinf(max_exp):
        return max_exp
    
    total = 0

    for x in log_X:
        total += numpy.exp(x - max_exp)
    
    return numpy.log(total) + max_exp

def log_sum_exp_prod(log_X, Y):

    max_exp = max(log_X)

    if numpy.isinf(max_exp):
        return max_exp

    total = 0

    for x, y in zip(log_X, Y):
        total += numpy.exp(x - max_exp) * y

    return numpy.log(total) + max_exp

def sigmoid(x):
    res = 1.0/(1.0+numpy.exp(-x))
    return res

def bucket(edges, value):
    return numpy.sum(value > edges)

def sticks_to_edges(sticks):
    return 1.0 - numpy.cumprod(1.0 - sticks)

#def normpdfln(x, m, prec):
#    return numpy.sum(-0.5*numpy.log(2*numpy.pi) + 0.5*numpy.log(prec) - 0.5*prec*(x-m)**2, axis=1)

#def normpdfln_1d(x, m, prec):
#    return numpy.sum(-0.5*numpy.log(2*numpy.pi) + 0.5*numpy.log(prec) - 0.5*prec*(x-m)**2)

def gammaln(x):
    #small  = numpy.nonzero(x < numpy.finfo(numpy.float64).eps)
    result = scipy.special.gammaln(x)
    #result[small] = -numpy.log(x[small])
    return result

def gammapdfln(x, a, b):
    return -gammaln(a) + a*numpy.log(b) + (a-1.0)*numpy.log(x) - b*x

def exp_gammapdfln(y, a, b):
    return a*numpy.log(b) - gammaln(a) + a*y - b*numpy.exp(y)

def betapdfln(x, a, b):
    if not isinstance(x, Iterable):
        x = numpy.array([x])
    if any(x == 0.0) or any(x == 1.0):
            return float('-inf')
    if any(x < 0.0) or any(x > 1.0):
        print x
    return gammaln(a+b) - gammaln(a) - gammaln(b) + (a-1.0)*numpy.log(x) + \
      (b-1.0)*numpy.log(1.0-x) 

def boundbeta(a,b):
    return (1.0-numpy.finfo(numpy.float64).eps)*(numpy.random.beta(a,b)-0.5) + 0.5
    #return numpy.random.beta(a,b)

def lnbetafunc(a):
    return numpy.sum(gammaln(a)) - gammaln(numpy.sum(a))

def dirichletpdfln(p, a):
    return -lnbetafunc(a) + numpy.sum((a-1)*numpy.log(p))

def logsumexp(X, axis=None):
    maxes = numpy.max(X, axis=axis)
    return numpy.log(numpy.sum(numpy.exp(X - maxes), axis=axis)) + maxes

def merge(l):
    return [item for sublist in l for item in sublist]

hotshotProfilers = {}
def hotshotit(func):
    def wrapper(*args, **kw):
        import hotshot
        global hotshotProfilers
        prof_name = func.func_name+".prof"
        profiler = hotshotProfilers.get(prof_name)
        if profiler is None:
            profiler = hotshot.Profile(prof_name)
            hotshotProfilers[prof_name] = profiler
        return profiler.runcall(func, *args, **kw)
    return wrapper

try:
    if bin(0): pass
except NameError, ne:
    def bin(x):
        """
        bin(number) -> string
        
        Stringifies an int or long in base 2.
        """
        if x < 0: return '-' + bin(-x)
        out = []
        if x == 0: out.append('0')
        while x > 0:
            out.append('01'[x & 1])
            x >>= 1
            pass
        try: return '0b' + ''.join(reversed(out))
        except NameError, ne2: out.reverse()
        return '0b' + ''.join(out)

def cifar10_codes(num_data=60000):
    filename = "cifar10-codes.pkl"
    fh       = open(filename)
    raw      = cPickle.load(fh)['codes']
    fh.close()
        
    binarized = []
    for i in range(num_data):
        vector = []
        for j in range(raw.shape[1]):
            boolish = map(lambda x: bool(int(x)), bin(raw[i][j])[2:])
            for k in range(len(boolish),32):
                boolish.insert(0, False)
            vector.extend(boolish)
        binarized.append(vector)
    return numpy.array(binarized)

def cifar100_codes(num_data=50000):
    filename = "cifar100-codes.pkl"
    fh       = open(filename)
    raw      = cPickle.load(fh)['codes']
    fh.close()

    binarized = []
    for i in range(num_data):
        vector = []
        for j in range(raw.shape[1]):
            boolish = map(lambda x: bool(int(x)), bin(raw[i][j])[2:])
            for k in range(len(boolish),32):
                boolish.insert(0, False)
            vector.extend(boolish)
        binarized.append(vector)
    return numpy.array(binarized)

def slice_sample(init_x, logprob, sigma=1.0, step_out=True, max_steps_out=1000, 
                 compwise=False, verbose=False):
    def direction_slice(direction, init_x):
        def dir_logprob(z):
            return logprob(direction*z + init_x)
    
        upper = sigma*numpy.random.rand()
        lower = upper - sigma
        llh_s = numpy.log(numpy.random.rand()) + dir_logprob(0.0)
    
        l_steps_out = 0
        u_steps_out = 0
        if step_out:
            while dir_logprob(lower) > llh_s and l_steps_out < max_steps_out:
                l_steps_out += 1
                lower       -= sigma
            while dir_logprob(upper) > llh_s and u_steps_out < max_steps_out:
                u_steps_out += 1
                upper       += sigma
            
        steps_in = 0
        while True:
            steps_in += 1
            new_z     = (upper - lower)*numpy.random.rand() + lower
            new_llh   = dir_logprob(new_z)
            if numpy.isnan(new_llh):
                print new_z, direction*new_z + init_x, new_llh, llh_s, init_x, logprob(init_x)
                raise Exception("Slice sampler got a NaN")
            if new_llh > llh_s:
                break
            elif new_z < 0:
                lower = new_z
            elif new_z > 0:
                upper = new_z
            else:
                raise Exception("Slice sampler shrank to zero!")

        if verbose:
            print "Steps Out:", l_steps_out, u_steps_out, " Steps In:", steps_in

        return new_z*direction + init_x
    
    dims = init_x.shape[0]
    if compwise:
        ordering = range(dims)
        numpy.random.shuffle(ordering)
        cur_x = init_x.copy()
        for d in ordering:
            direction    = numpy.zeros((dims))
            direction[d] = 1.0
            cur_x = direction_slice(direction, cur_x)
        return cur_x
            
    else:
        direction = numpy.random.randn(dims)
        direction = direction / numpy.sqrt(numpy.sum(direction**2))
        return direction_slice(direction, init_x)

def nips(num_docs=1740):
    file         = 'nips_1-12.mat'
    mat          = scipy.io.loadmat(file, struct_as_record=True)
    authors      = map(lambda x: str(x[0][0]), mat['authors_names'])
    counts       = mat['counts'][:,:num_docs]
    docs_authors = mat['docs_authors'][:num_docs,:]
    words        = map(lambda x: str(x[0][0]), mat['words'])
    years        = [] #mat['years'].T[0]
    return { 'authors'      : authors,
             'counts'       : counts,
             'docs_authors' : docs_authors,
             'words'        : words,
             'years'        : years }

def hmc(init_x, logprob, logprob_grad, num_steps, step_size):
    energy      = lambda x: -logprob(x)
    energy_grad = lambda x: -logprob_grad(x)

    g = energy_grad(init_x)
    E = energy(init_x)
    p = numpy.random.randn(init_x.shape[0])
    H = numpy.sum(p**2)/2.0 + E

    new_x = init_x
    new_g = g
    for step in range(num_steps):
        p     = p - step_size * new_g / 2.0
        new_x = new_x + step_size * p
        new_g = energy_grad(new_x)
        p     = p - step_size * new_g / 2.0;

    new_E = energy(new_x)
    new_H = sum(p**2)/2.0 + new_E

    if numpy.log(numpy.random.rand()) < (H - new_H):
        return new_x, True
    else:
        return init_x, False

def bounded_hmc(init_x, logprob, logprob_grad, lower, upper, num_steps, step_size):
    energy      = lambda x: -logprob(x)
    energy_grad = lambda x: -logprob_grad(x)

    g = energy_grad(init_x)
    E = energy(init_x)
    p = numpy.random.randn(init_x.shape[0])
    H = numpy.sum(p**2)/2.0 + E

    new_x = init_x
    new_g = g
    
    for step in range(num_steps):
        
        p     = p - step_size * new_g / 2.0
        new_x = new_x + step_size * p

        reflections = 0
        for k in range(init_x.shape[0]):
            if new_x[k] < lower[k]:
                new_x[k] = lower[k] + (lower[k] - new_x[k])
                p[k] = -p[k]
                reflections = reflections + 1

            if new_x[k] > upper[k]:
                new_x[k] = upper[k] - (new_x[k] - upper[k])
                p[k] = -p[k]
                reflections = - reflections + 1

        new_g = energy_grad(new_x)
        p     = p - step_size * new_g / 2.0;

    new_E = energy(new_x)
    new_H = sum(p**2)/2.0 + new_E

    if numpy.log(numpy.random.rand()) < (H - new_H):
        return new_x, True
    else:
        return init_x, False
    
def hmc_trans(init_x, logprob, logprob_grad, trans, num_steps, step_size):
    energy      = lambda x: -logprob(x)
    energy_grad = lambda x: -logprob_grad(x)
    transformer = lambda x,y: trans(x,y)
    
    g = energy_grad(init_x)*transformer(init_x,'tran_grad') - \
       transformer(init_x,'jacob_grad')
    E = energy(init_x) - transformer(init_x,'det_jacob')

    if E == float('-inf'):
        return init_x, False
    
    p = numpy.random.randn(init_x.shape[0])
    H = numpy.sum(p**2)/2.0 + E

    new_x = init_x
    tran_x = transformer(new_x, 'tran_x')
    new_g = g
    for step in range(num_steps):
        p     = p - step_size * new_g / 2.0
        tran_x = tran_x + step_size * p
        ## inverse transformation         
        new_x = transformer(tran_x, 'orig_x')
        new_g = energy_grad(new_x) * \
          transformer(new_x, 'tran_grad') - \
          transformer(new_x, 'jacob_grad')
        p     = p - step_size * new_g / 2.0;

        
    new_E = energy(new_x) - transformer(new_x, 'det_jacob')
    new_H = sum(p**2)/2.0 + new_E

    if numpy.log(numpy.random.rand()) < (H - new_H):
        return new_x, True
    else:
        return init_x, False



    

    
import os
import sys
import time
import atexit
import cPickle
import socket

def fail(msg):
    print >>sys.stderr, "ERROR:", msg
    sys.exit(-1)

class Experiment:

    def __init__(self, name, func, exptdir):
        self.name      = name
        self.func      = func
        self.exptdir   = exptdir
        self.dir       = os.path.join(exptdir, name)
        self.locks     = {}

        if not os.path.exists(exptdir):
            fail("Experiment directory '%s' does not exist" % (exptdir))
        elif not os.path.isdir(exptdir):
            fail("Path '%s' exists but is not a directory" % (exptdir))

        if not os.path.exists(self.dir):
            print >>sys.stderr, "Creating directory '%s'" % (self.dir)
            os.mkdir(self.dir)
        elif not os.path.isdir(self.dir):
            fail("Path '%s' exists but is not a directory" % (self.dir))

        atexit.register(self._atexit)

        self._lock("jobs.pkl")
        if not os.path.exists(os.path.join(self.dir, "jobs.pkl")):
            fh = open(os.path.join(self.dir, "jobs.pkl"), "w")
            jobs = {}
            cPickle.dump(jobs, fh)
            fh.close()
        self._unlock("jobs.pkl")

    def run(self, num_runs):
        while not self._runjob(num_runs):
            time.sleep(10)

    def cleanup(self):
        self._lock("jobs.pkl")

        fh = open(os.path.join(self.dir, "jobs.pkl"))
        jobs = cPickle.load(fh)
        fh.close()

        for key in jobs.keys():
            if jobs[key] == 'running':
                print >>sys.stderr, "Setting job %d back to 'available' from 'running'" % key
                jobs[key] = 'available'

        fh = open(os.path.join(self.dir, "jobs.pkl"), "w")
        cPickle.dump(jobs, fh)
        fh.close()
        
        self._unlock("jobs.pkl")

    def reset(self, job):
        self._lock("jobs.pkl")

        fh = open(os.path.join(self.dir, "jobs.pkl"))
        jobs = cPickle.load(fh)
        fh.close()
        
        print >>sys.stderr, "Setting job %d back to 'available'" % (job)
        jobs[job] = 'available'

        fh = open(os.path.join(self.dir, "jobs.pkl"), "w")
        cPickle.dump(jobs, fh)
        fh.close()
        
        self._unlock("jobs.pkl")

    def _runjob(self, num_runs):

        ####################################################################################
        # Phase I: Try to find a job to work on.
        self._lock("jobs.pkl")

        fh = open(os.path.join(self.dir, "jobs.pkl"))
        jobs = cPickle.load(fh)
        fh.close()

        for i in range(num_runs):
            if not jobs.has_key(i):
                jobs[i] = 'available'
        
        available = filter(lambda key: jobs[key] == 'available', jobs.keys())
        
        if len(available) > 0:
            jobs[available[0]] = 'running'

            fh = open(os.path.join(self.dir, "jobs.pkl"), "w")
            cPickle.dump(jobs, fh)
            fh.close()
                    
        self._unlock("jobs.pkl")

        ####################################################################################
        # Phase II: Work on the job, if there is one.
        if len(available) > 0:
            t0 = time.time()
            result = self.func(available[0])
            t1 = time.time()
            print >>sys.stderr, "Completed job %d in %0.2f seconds" % (available[0], (t1-t0))

            filename = "result-%05d.pkl" % (available[0])
            if not self._trylock(filename):
                fail("Someone else has %s locked - that is bad." % (filename))

            pickledata = { 'time_elapsed' : t1-t0,
                           'result'       : result,
                           'host'         : socket.gethostname() }

            outfile = os.path.join(self.dir, filename)
            fh = open(outfile, "w")
            cPickle.dump(pickledata, fh, cPickle.HIGHEST_PROTOCOL)
            fh.close()
            print >>sys.stderr, "Wrote data to '%s' for job %d" % (outfile, available[0])
        
            self._unlock(filename)

            self._lock("jobs.pkl")
            fh = open(os.path.join(self.dir, "jobs.pkl"))
            jobs = cPickle.load(fh)
            fh.close()

            jobs[available[0]] = 'complete'

            fh = open(os.path.join(self.dir, "jobs.pkl"), "w")
            cPickle.dump(jobs, fh)
            fh.close()
                    
            incomplete = filter( lambda key: jobs[key] != 'complete', jobs.keys())

            if len(incomplete) == 0:

                # DO NOT ATTEMPT WITH THESE GIANT F**CKING FILES
                #if not os.path.exists(os.path.join(self.dir, "results.pkl")):
                #    results = []
                #    for i in range(len(jobs)):
                #        filename = "result-%05d.pkl" % (i)
                #        if not self._trylock(filename):
                #            fail("Someone else has %s locked - that is bad." % (filename))
                #        fh = open(os.path.join(self.dir, filename))
                #        result_i = cPickle.load(fh)
                #        results.append(result_i)
                #        fh.close()
                #        os.remove(os.path.join(self.dir, filename))
                #        self._unlock(filename)
                #    filename = "results.pkl"
                #    if not self._trylock(filename):
                #        fail("Someone else has %s locked - that is bad." % (filename))
                #    fh = open(os.path.join(self.dir, filename), "w")
                #    cPickle.dump(results, fh, cPickle.HIGHEST_PROTOCOL)
                #    fh.close()
                #    self._unlock(filename)
                complete = True
            else:
                complete = False

            self._unlock("jobs.pkl")
            return complete
        else:
            return True

    def load(self):
        pass

    def _atexit(self):
        for file in self.locks.keys():
            print >>sys.stderr, "Unlocking hanging lock '%s'" % (file)
            self._unlock(file)

    def _trylock(self, file):
        fullpath = os.path.join(self.dir, "." + file + ".lock")
        if self.locks.has_key(file):
            return True
        command  = "ln -s /dev/null \"%s\" 2> /dev/null" % (fullpath)
        success  = os.system(command) == 0
        if success:
            self.locks[file] = True
        return success

    def _lock(self, file, timeout=60):
        start = time.time()
        while not self._trylock(file):
            time.sleep(1)
            if (time.time() - start) > timeout:
                fail("Timed out waiting for lock on '%s'" % (file))
            else:
                print >>sys.stderr, "Waiting for lock on '%s'" % (file)
        print >>sys.stderr, "Got lock on '%s'" % (file)
        
    def _unlock(self, file):
        fullpath = os.path.join(self.dir, "." + file + ".lock")
        if not self.locks.has_key(file):
            fail("Attempted to unlock a file that was not locked (%s)" % (fullpath))
        fullpath_del = fullpath + ".del"
        command = "mv \"%s\" \"%s\" && rm \"%s\"" % (fullpath, fullpath_del, fullpath_del)
        if os.system(command) != 0:
            fail("Could not perform atomic move of file '%s'" % (fullpath))
        self.locks.pop(file)
        print >>sys.stderr, "Unlocked file '%s'" % (file)

