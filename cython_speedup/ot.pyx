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
    cdef numpy.ndarray[double, ndim=1, mode="c"] dx = numpy.ones(n) / n
    cdef numpy.ndarray[double, ndim=1, mode="c"] dy = numpy.ones(m) / m
    cdef numpy.ndarray[double, ndim=1, mode="c"] q = numpy.ones(m) * numpy.average(g)

    cdef numpy.ndarray[double, ndim=1, mode="c"] a = numpy.empty(n)
    cdef numpy.ndarray[double, ndim=1, mode="c"] u = numpy.zeros(n)
    cdef numpy.ndarray[double, ndim=1, mode="c"] v = numpy.zeros(m)
    cdef numpy.ndarray[double, ndim=1, mode="c"] b = numpy.ones(m)
    cdef numpy.ndarray[double, ndim=2, mode="c"] K = numpy.ascontiguousarray(numpy.exp(-C / epsilon_i))

    cdef double alpha1 = lambda1 / (lambda1 + epsilon_i)
    cdef double alpha2 = lambda2 / (lambda2 + epsilon_i)
    cdef int epsilon_index = 0
    cdef int iterations_since_epsilon_adjusted = 0
    cdef int num_inner_iter_max = 100

    cdef int i, j, k

    for i in range(scaling_iter + extra_iter):
        # scaling iteration
        # a = (g / (numpy.dot(K, b * dy))) ** alpha1 * numpy.exp(-u / (lambda1 + epsilon_i))
        a = numpy.multiply(numpy.power(numpy.divide(g, (numpy.dot(K, numpy.multiply(b, dy)))), alpha1), numpy.exp(numpy.divide(numpy.negative(u), (lambda1 + epsilon_i))))
        # b = (q / (numpy.dot(K.T, a * dx))) ** alpha2 * numpy.exp(-v / (lambda2 + epsilon_i))
        b = numpy.multiply(numpy.power(numpy.divide(q, (numpy.dot(numpy.transpose(K), numpy.multiply(a, dx)))), alpha2), numpy.exp(numpy.divide(numpy.negative(v), (lambda2 + epsilon_i))))

        # stabilization
        iterations_since_epsilon_adjusted += 1
        if i < scaling_iter and (max(max(numpy.absolute(a)), max(numpy.absolute(b))) > tau or iterations_since_epsilon_adjusted == num_inner_iter_max):
            if iterations_since_epsilon_adjusted == num_inner_iter_max :
                epsilon_index += 1
                iterations_since_epsilon_adjusted = 0
                u = numpy.add(u, numpy.multiply(epsilon_i, numpy.log(a)))
                v = numpy.add(v, numpy.multiply(epsilon_i, numpy.log(b)))  # absorb
                epsilon_i = (epsilon0 - epsilon_final) * exp(-epsilon_index) \
                        + epsilon_final
                alpha1 = lambda1 / (lambda1 + epsilon_i)
                alpha2 = lambda2 / (lambda2 + epsilon_i)
            else :
                u = numpy.add(u, numpy.multiply(epsilon_i, numpy.log(a)))
                v = numpy.add(v, numpy.multiply(epsilon_i, numpy.log(b)))  # absorb
            K = numpy.exp(numpy.divide(numpy.add(numpy.subtract(numpy.transpose(numpy.array([u])),C), numpy.array([v])), epsilon_i))
            a = numpy.ones(n)
            b = numpy.ones(m)

    return numpy.multiply(numpy.transpose(numpy.multiply(numpy.transpose(K), a)), b)
