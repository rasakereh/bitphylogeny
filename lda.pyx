import numpy        as np
import numpy.random as npr

cimport numpy as np
cimport cython

ctypedef np.int32_t dtype_t
        
cdef extern from "math.h":
    double log(double x)

@cython.boundscheck(False)
@cython.wraparound(False)
def lda_gibbs(list W,
              list Z,
              np.ndarray[dtype_t, ndim=2] NW not None,
              np.ndarray[dtype_t, ndim=2] ND not None,
              np.ndarray[dtype_t, ndim=1] NWsum not None,
              double beta, double alpha):
    cdef int num_docs   = ND.shape[0]
    cdef int num_topics = ND.shape[1]
    cdef int vocab_size = NW.shape[0]
    cdef int cur_topic, NDsum, lenZn
    cdef np.ndarray[np.float64_t, ndim=1] probs = np.zeros((num_topics))
    cdef np.ndarray[dtype_t, ndim=1] Wn
    cdef np.ndarray[dtype_t, ndim=1] Zn
    cdef Py_ssize_t n, i, j, Wni, new_topic
    cdef double total
    cdef np.ndarray[np.float64_t, ndim=1] rands

    for n in range(num_docs):
        lenZn = len(Z[n])
        Zn = Z[n]
        Wn = W[n]
        rands = npr.rand(lenZn)
        for i in range(lenZn):

            Wni = Wn[i]
            
            cur_topic              = Zn[i]
            NW[Wni,cur_topic]     -= 1
            ND[n,cur_topic]       -= 1
            NWsum[cur_topic]      -= 1
            NDsum                  = lenZn - 1

            total = 0.0
            for j in range(num_topics):
                probs[j] = ((NW[Wni,j] + beta) / (NWsum[j] + vocab_size*beta)) * ((ND[n,j] + alpha) / (NDsum + num_topics*alpha))
                total += probs[j]
            rands[i] = total*rands[i]
            total = 0.0
            new_topic = num_topics-1
            for j in range(num_topics):
                if rands[i] < total:
                    new_topic = j-1
                    break
                total += probs[j]

            #probs = ((NW[Wni] + beta) / (NWsum + vocab_size*beta)) * ((ND[n] + alpha) / (NDsum + num_topics*alpha))            
            #probs = probs / np.sum(probs)
            #new_topic = np.sum(rands[i] > np.cumsum(probs))

            Zn[i]              = new_topic
            NW[Wni,new_topic] += 1
            ND[n,new_topic]   += 1
            NWsum[new_topic]  += 1

@cython.boundscheck(False)
@cython.wraparound(False)
def doc_gibbs(np.ndarray[dtype_t, ndim=1] W not None,
              np.ndarray[dtype_t, ndim=1] Z not None,
              np.ndarray[dtype_t, ndim=2] NW not None,
              np.ndarray[dtype_t, ndim=1] ND not None,
              np.ndarray[np.int64_t, ndim=1] NWsum not None,
              np.ndarray[np.float64_t, ndim=1] normed not None,
              double beta):
    
    cdef int num_topics = ND.shape[0]
    cdef int vocab_size = NW.shape[0]
    cdef int lenZn
    cdef np.ndarray[np.float64_t, ndim=1] probs = np.zeros((num_topics))
    cdef Py_ssize_t i, j, Wni, new_topic, cur_topic
    cdef double total
    cdef np.ndarray[np.float64_t, ndim=1] rands

    lenZn = Z.shape[0]
    rands = npr.rand(lenZn)
    for i in range(lenZn):
        Wni = W[i]
        
        cur_topic              = Z[i]
        NW[Wni,cur_topic]     -= 1
        ND[cur_topic]         -= 1
        NWsum[cur_topic]      -= 1

        total = 0.0
        for j in range(num_topics):
            probs[j] = ((NW[Wni,j] + beta) / (NWsum[j] + vocab_size*beta)) * normed[j]
            #print NW[Wni,j], beta, NWsum[j], vocab_size*beta, normed[j], probs[j]
            total += probs[j]
        rands[i] = total*rands[i]
        total = 0.0
        new_topic = num_topics-1
        for j in range(num_topics):
            if rands[i] < total:
                new_topic = j-1
                break
            total += probs[j]

        Z[i]              = new_topic
        NW[Wni,new_topic] += 1
        ND[new_topic]     += 1
        NWsum[new_topic]  += 1

@cython.boundscheck(False)
@cython.wraparound(False)
def doc_logprob(np.ndarray[dtype_t, ndim=1] W not None,
                np.ndarray[dtype_t, ndim=1] Z not None,
                np.ndarray[np.float64_t, ndim=2] thetas not None,
                np.ndarray[np.float64_t, ndim=1] normed not None):
    cdef int lenZn
    cdef Py_ssize_t i, Wni, Zni
    cdef double logprob = 0.0

    lenZn = Z.shape[0]
    for i in range(lenZn):
        Wni = W[i]
        Zni = Z[i]
        logprob += log(thetas[Wni,Zni]) + log(normed[Zni])
        
    return logprob
