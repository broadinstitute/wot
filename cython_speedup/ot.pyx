import cython

import numpy
cimport numpy

from libc.math cimport exp

@cython.cdivision(True)
@cython.wraparound(False)
@cython.boundscheck(False)
def __cy__transport_stablev3(
        numpy.ndarray[double, ndim=2] C not None,
        numpy.ndarray[double, ndim=1] g not None,
        double lambda1, double lambda2, double epsilon,
        long scaling_iter, long extra_iter):
    cdef int n = C.shape[0]
    cdef int m = C.shape[1]

    cdef double tau = 10e100
    cdef double epsilon0 = 100
    cdef double epsilon_final = epsilon


    cdef double epsilon_i = epsilon
    cdef numpy.ndarray[double, ndim=1] dx = numpy.ones(C.shape[0]) / C.shape[0]
    cdef numpy.ndarray[double, ndim=1] dy = numpy.ones(C.shape[1]) / C.shape[1]
    cdef numpy.ndarray[double, ndim=1] q = numpy.ones(C.shape[1]) * numpy.average(g)

    cdef numpy.ndarray[double, ndim=1] a
    cdef numpy.ndarray[double, ndim=1] u = numpy.zeros(len(g))
    cdef numpy.ndarray[double, ndim=1] v = numpy.zeros(len(q))
    cdef numpy.ndarray[double, ndim=1] b = numpy.ones(len(q))
    cdef numpy.ndarray[double, ndim=2] K = numpy.exp(-C / epsilon_i)

    cdef double alpha1 = lambda1 / (lambda1 + epsilon_i)
    cdef double alpha2 = lambda2 / (lambda2 + epsilon_i)
    cdef int epsilon_index = 0
    cdef int iterations_since_epsilon_adjusted = 0
    cdef int num_inner_iter_max = 100


    for i in range(scaling_iter + extra_iter):
        # scaling iteration
        a = (g /   (K.dot(numpy.multiply(b, dy)))) ** alpha1 * numpy.exp(-u / (lambda1 + epsilon_i))
        b = (q / (K.T.dot(numpy.multiply(a, dx)))) ** alpha2 * numpy.exp(-v / (lambda2 + epsilon_i))

        # stabilization
        iterations_since_epsilon_adjusted += 1
        if i < scaling_iter and (max(max(abs(a)), max(abs(b))) > tau or iterations_since_epsilon_adjusted == num_inner_iter_max):
            if iterations_since_epsilon_adjusted == num_inner_iter_max :
                epsilon_index += 1
                iterations_since_epsilon_adjusted = 0
                u = u + epsilon_i * numpy.log(a)
                v = v + epsilon_i * numpy.log(b)  # absorb
                epsilon_i = (epsilon0 - epsilon_final) * exp(-epsilon_index) \
                        + epsilon_final
                alpha1 = lambda1 / (lambda1 + epsilon_i)
                alpha2 = lambda2 / (lambda2 + epsilon_i)
            else :
                u = u + epsilon_i * numpy.log(a)
                v = v + epsilon_i * numpy.log(b)  # absorb
            K = numpy.exp((numpy.array([u]).T - C + numpy.array([v])) / epsilon_i)
            a = numpy.ones(len(g))
            b = numpy.ones(len(q))

    return (K.T * a).T * b
