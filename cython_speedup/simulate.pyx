import cython

import numpy as np
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
def __cy__multivariate_normal_evolving_mixture(
        np.ndarray[double, ndim=3, mode="c"] means not None,
        np.ndarray[double, ndim=4, mode="c"] covs not None,
        np.ndarray[double, ndim=2, mode="c"] p not None,
        np.ndarray[long, ndim=1, mode="c"] size not None):
    cdef int t = means.shape[0]
    cdef int k = means.shape[1]
    cdef int n = means.shape[2]

    cdef long T = np.sum(size)
    cdef np.ndarray[double, ndim=2, mode="c"] result = \
            np.zeros([T, n], dtype=np.float64)

    cdef int i, j, c
    cdef long cur_size, r

    cdef np.ndarray[long, ndim=1, mode="c"] picks
    c = 0
    for i in range(t):
        picks = np.random.multinomial(size[i], p[i])
        for j in range(k):
            result[c:c+picks[j]] = \
                    np.random.multivariate_normal(means[i][j],
                            covs[i][j], size=picks[j])
            c += picks[j]
    return result
