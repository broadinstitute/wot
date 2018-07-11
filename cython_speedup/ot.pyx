import cython

import numpy
cimport numpy

from libc.math cimport exp

cdef double get_reg(double epsilon0, double epsilon_final, int n):  # exponential decreasing
    return (epsilon0 - epsilon_final) * exp(-n) + epsilon_final

cdef inline numpy.ndarray[double, ndim=1] update(
        numpy.ndarray[double, ndim=1] g,
        numpy.ndarray[double, ndim=2] K,
        numpy.ndarray[double, ndim=1] b,
        numpy.ndarray[double, ndim=1] d,
        numpy.ndarray[double, ndim=1] u,
        double alpha, double l, double e):
    # equivalent of :
    # (g / (K.dot(numpy.multiply(b, d)))) ** alpha * numpy.exp(-u / (l + e))
    return (g / (K.dot(numpy.multiply(b, d)))) ** alpha * numpy.exp(-u / (l + e))


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


    for i in range(scaling_iter):
        # scaling iteration
        a = update(g, K, b, dy, u, alpha1, lambda1, epsilon_i)
        b = update(q, K.T, a, dx, v, alpha2, lambda2, epsilon_i)

        # stabilization
        iterations_since_epsilon_adjusted += 1
        if (max(max(abs(a)), max(abs(b))) > tau):
            u = u + epsilon_i * numpy.log(a)
            v = v + epsilon_i * numpy.log(b)  # absorb
            K = numpy.exp((numpy.array([u]).T - C + numpy.array([v])) / epsilon_i)
            a = numpy.ones(len(g))
            b = numpy.ones(len(q))

        if iterations_since_epsilon_adjusted == num_inner_iter_max :
            epsilon_index += 1
            iterations_since_epsilon_adjusted = 0
            u = u + epsilon_i * numpy.log(a)
            v = v + epsilon_i * numpy.log(b)  # absorb
            epsilon_i = get_reg(epsilon0, epsilon_final, epsilon_index)
            alpha1 = lambda1 / (lambda1 + epsilon_i)
            alpha2 = lambda2 / (lambda2 + epsilon_i)
            K = numpy.exp((numpy.array([u]).T - C + numpy.array([v])) / epsilon_i)
            a = numpy.ones(len(g))
            b = numpy.ones(len(q))

    for i in range(extra_iter):
        a = update(g, K, b, dy, u, alpha1, lambda1, epsilon_i)
        b = update(q, K.T, a, dx, v, alpha2, lambda2, epsilon_i)

    return (K.T * a).T * b
