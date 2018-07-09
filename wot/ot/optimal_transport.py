# -*- coding: utf-8 -*-

import numpy as np
import ot as pot
import scipy.stats


def transport_stable_learnGrowth(C, lambda1, lambda2, epsilon, scaling_iter, g, numInnerItermax=None, tau=None,
                                 epsilon0=None, extra_iter=1000, growth_iters=3):
    """
    Compute the optimal transport with stabilized numerics.
    Args:

        C: cost matrix to transport cell i to cell j
        lambda1: regularization parameter for marginal constraint for p.
        lambda2: regularization parameter for marginal constraint for q.
        epsilon: entropy parameter
        scaling_iter: number of scaling iterations
        g: growth value for input cells
    """

    for i in range(growth_iters):
        if i == 0:
            rowSums = g
        else:
            rowSums = Tmap.sum(axis=1) / Tmap.shape[1]

        Tmap = transport_stablev2(C, lambda1, lambda2, epsilon,
                                  scaling_iter, rowSums, numInnerItermax=numInnerItermax, tau=tau,
                                  epsilon0=epsilon0)
    return Tmap


def transport_stablev2(C, lambda1, lambda2, epsilon, scaling_iter, g, numInnerItermax=None, tau=None,
                       epsilon0=None, extra_iter=1000):
    """
    Compute the optimal transport with stabilized numerics.
    Args:

        C: cost matrix to transport cell i to cell j
        lambda1: regularization parameter for marginal constraint for p.
        lambda2: regularization parameter for marginal constraint for q.
        epsilon: entropy parameter
        scaling_iter: number of scaling iterations
        g: growth value for input cells
    """
    warm_start = tau is not None
    epsilon_final = epsilon

    def get_reg(n):  # exponential decreasing
        return (epsilon0 - epsilon_final) * np.exp(-n) + epsilon_final

    epsilon_i = epsilon0 if warm_start else epsilon
    dx = np.ones(C.shape[0]) / C.shape[0]
    dy = np.ones(C.shape[1]) / C.shape[1]
    p = g
    q = np.ones(C.shape[1]) * np.average(g)

    u = np.zeros(len(p))
    v = np.zeros(len(q))
    b = np.ones(len(q))
    K = np.exp(-C / epsilon_i)

    alpha1 = lambda1 / (lambda1 + epsilon_i)
    alpha2 = lambda2 / (lambda2 + epsilon_i)
    epsilon_index = 0
    iterations_since_epsilon_adjusted = 0

    for i in range(scaling_iter):
        # scaling iteration
        a = (p / (K.dot(np.multiply(b, dy)))) ** alpha1 * np.exp(-u / (lambda1 + epsilon_i))
        b = (q / (K.T.dot(np.multiply(a, dx)))) ** alpha2 * np.exp(-v / (lambda2 + epsilon_i))

        # stabilization
        iterations_since_epsilon_adjusted += 1
        if (max(max(abs(a)), max(abs(b))) > tau):
            u = u + epsilon_i * np.log(a)
            v = v + epsilon_i * np.log(b)  # absorb
            K = np.exp((np.array([u]).T - C + np.array([v])) / epsilon_i)
            a = np.ones(len(p))
            b = np.ones(len(q))

        if (warm_start and iterations_since_epsilon_adjusted == numInnerItermax):
            epsilon_index += 1
            iterations_since_epsilon_adjusted = 0
            u = u + epsilon_i * np.log(a)
            v = v + epsilon_i * np.log(b)  # absorb
            epsilon_i = get_reg(epsilon_index)
            alpha1 = lambda1 / (lambda1 + epsilon_i)
            alpha2 = lambda2 / (lambda2 + epsilon_i)
            K = np.exp((np.array([u]).T - C + np.array([v])) / epsilon_i)
            a = np.ones(len(p))
            b = np.ones(len(q))

    for i in range(extra_iter):
        a = (p / (K.dot(np.multiply(b, dy)))) ** alpha1 * np.exp(-u / (lambda1 + epsilon_i))
        b = (q / (K.T.dot(np.multiply(a, dx)))) ** alpha2 * np.exp(-v / (lambda2 + epsilon_i))

    return (K.T * a).T * b


def transport_stable(p, q, C, lambda1, lambda2, epsilon, scaling_iter, g):
    """
    Compute the optimal transport with stabilized numerics.
    Args:
        p: uniform distribution on input cells
        q: uniform distribution on output cells
        C: cost matrix to transport cell i to cell j
        lambda1: regularization parameter for marginal constraint for p.
        lambda2: regularization parameter for marginal constraint for q.
        epsilon: entropy parameter
        scaling_iter: number of scaling iterations
        g: growth value for input cells
    """
    u = np.zeros(len(p))
    v = np.zeros(len(q))
    b = np.ones(len(q))
    p = p * g
    q = q * np.average(g)
    K0 = np.exp(-C / epsilon)
    K = np.copy(K0)
    alpha1 = lambda1 / (lambda1 + epsilon)
    alpha2 = lambda2 / (lambda2 + epsilon)
    for i in range(scaling_iter):
        # scaling iteration
        a = (p / (K.dot(b))) ** alpha1 * np.exp(-u / (lambda1 + epsilon))
        b = (q / (K.T.dot(a))) ** alpha2 * np.exp(-v / (lambda2 + epsilon))
        # stabilization
        if (max(max(abs(a)), max(abs(b))) > 1e100):
            u = u + epsilon * np.log(a)
            v = v + epsilon * np.log(b)  # absorb
            K = (K0.T * np.exp(u / epsilon)).T * np.exp(v / epsilon)
            a = np.ones(len(p))
            b = np.ones(len(q))
    return (K.T * a).T * b


def optimal_transport(cost_matrix, growth_rate, p=None, q=None, solver=None,
                      delta_days=1, epsilon=0.1, lambda1=1.,
                      lambda2=1., min_transport_fraction=0.05,
                      max_transport_fraction=0.4, min_growth_fit=0.9,
                      l0_max=100, scaling_iter=250, epsilon_adjust=1.1,
                      lambda_adjust=1.5, numItermax=100, epsilon0=100.0, numInnerItermax=10, tau=1000.0, stopThr=1e-06,
                      growth_iters=3):
    if solver == 'unbalanced':

        g = growth_rate ** delta_days
        transport = transport_stable_learnGrowth(C=cost_matrix, lambda1=lambda1, lambda2=lambda2, epsilon=epsilon,
                                                 scaling_iter=scaling_iter, g=g, numInnerItermax=numInnerItermax,
                                                 tau=tau, epsilon0=epsilon0, growth_iters=growth_iters)
        return {'transport': transport}
    elif solver == 'floating_epsilon':
        return optimal_transport_with_entropy(cost_matrix, growth_rate, p=p, q=q,
                                              delta_days=delta_days, epsilon=epsilon, lambda1=lambda1,
                                              lambda2=lambda2, min_transport_fraction=min_transport_fraction,
                                              max_transport_fraction=max_transport_fraction,
                                              min_growth_fit=min_growth_fit,
                                              l0_max=l0_max, scaling_iter=scaling_iter,
                                              epsilon_adjust=epsilon_adjust,
                                              lambda_adjust=lambda_adjust)
    elif solver == 'sinkhorn_epsilon':
        return sinkhorn_epsilon(cost_matrix, growth_rate, p=p, q=q,
                                delta_days=delta_days, epsilon=epsilon, numItermax=numItermax,
                                epsilon0=epsilon0,
                                numInnerItermax=numInnerItermax, tau=tau, stopThr=stopThr)
    elif solver == 'unregularized':
        return unregularized(cost_matrix, growth_rate, p=p, q=q,
                             delta_days=delta_days)
    else:
        raise ValueError('Unknown solver: ' + solver)


def sinkhorn_epsilon(cost_matrix, growth_rate, p=None, q=None,
                     delta_days=1, epsilon=0.1, numItermax=100, epsilon0=100.0,
                     numInnerItermax=10, tau=1000.0, stopThr=1e-06):
    p = np.ones(cost_matrix.shape[0])
    q = np.ones(cost_matrix.shape[1])

    g = growth_rate ** delta_days

    p = p * g
    q = q / q.sum()
    p = p / p.sum()
    val = pot.bregman.sinkhorn_epsilon_scaling(p, q, cost_matrix, reg=epsilon, numItermax=numItermax, epsilon0=epsilon0,
                                               numInnerItermax=numInnerItermax,
                                               tau=tau, stopThr=stopThr, warmstart=None, verbose=True,
                                               print_period=10,
                                               log=False)
    return {'transport': val, 'lambda1': 1, 'lambda2': 1, 'epsilon': 1}


def unregularized(cost_matrix, growth_rate, p=None, q=None, delta_days=1):
    p = np.ones(cost_matrix.shape[0])
    q = np.ones(cost_matrix.shape[1])

    g = growth_rate ** delta_days

    p = p * g
    q = q / q.sum()
    p = p / p.sum()
    val = pot.emd(p, q, cost_matrix, numItermax=max(10000000, cost_matrix.shape[0] * cost_matrix.shape[1]))
    return {'transport': val, 'lambda1': 1, 'lambda2': 1, 'epsilon': 1}


def optimal_transport_with_entropy(cost_matrix, growth_rate, p=None, q=None,
                                   delta_days=1, epsilon=0.1, lambda1=1.,
                                   lambda2=1., min_transport_fraction=0.05,
                                   max_transport_fraction=0.4, min_growth_fit=0.9,
                                   l0_max=100, scaling_iter=250, epsilon_adjust=1.1,
                                   lambda_adjust=1.5):
    """
    Compute the optimal transport.

    Args:
        cost_matrix (ndarray): A 2D matrix that indicates the cost of
        transporting cell i to cell j.
                               Can be generated by
                               sklearn.metrics.pairwise.pairwise_distances
                               for example.
        growth_rate (ndarray): A 1D matrix that indicates the growth rate of
        cells. A growth rate of 2 means that a cell will have 2 descendants
        after 1 day.
        delta_days (float): Elapsed time in days between time points
        epsilon (float): Controls the entropy of the transport map. An
        extremely large entropy parameter will give a
                         maximally entropic transport map, and an extremely
                         small entropy parameter will give a nearly
                         deterministic transport map (but could also lead to
                         numerical instability in the algorithm)
        lambda1 (float): Regularization parameter that controls the fidelity
        of the constraints on p.
                        As lamda1 gets larger, the constraints become more
                        stringent
        lambda2 (float): Regularization parameter that controls the fidelity
        of the constraints on q.
                        As lamda2 gets larger, the constraints become more
                        stringent
        min_transport_fraction (float): The minimum fraction of cells at time
        t that are transported to time t + 1.
        max_transport_fraction (float): The maximum fraction of cells at time
        t that are transported to time t + 1.
        min_growth_fit (float):
        l0_max (float):
        scaling_iter (int): Number of scaling iterations

    Returns:
        ndarray: A dictionary with transport (the transport map), epsilon,
        lambda1, and lambda2
    """

    if p is None:
        p = np.ones(cost_matrix.shape[0]) / cost_matrix.shape[0]
    if q is None:
        q = np.ones(cost_matrix.shape[1]) / cost_matrix.shape[1]

    g = growth_rate ** delta_days
    l0 = 1.
    e0 = 1.
    while True:
        transport = transport_stable(p, q, cost_matrix, lambda1 * l0,
                                     lambda2 * l0, epsilon * e0, scaling_iter,
                                     g)
        avg_transport = np.average(
            [np.exp(scipy.stats.entropy(trans)) for trans in transport])
        growth_fit = 1 - np.linalg.norm(
            transport.sum(1) - g / cost_matrix.shape[0]) ** 2 / np.linalg.norm(
            g / cost_matrix.shape[0]) ** 2
        if avg_transport == 0:
            while avg_transport == 0:
                e0 *= epsilon_adjust
                transport = transport_stable(p, q, cost_matrix, lambda1 * l0,
                                             lambda2 * l0, epsilon * e0,
                                             scaling_iter, g)
                avg_transport = np.average(
                    [np.exp(scipy.stats.entropy(trans)) for trans in transport])
            break
        elif (growth_fit < min_growth_fit) and (l0 < l0_max):
            l0 *= lambda_adjust
        elif avg_transport < transport.shape[1] * min_transport_fraction:
            e0 *= epsilon_adjust
        elif avg_transport < transport.shape[1] * max_transport_fraction:
            break
        else:
            e0 /= epsilon_adjust
    return {'transport': transport, 'lambda1': lambda1 * l0,
            'lambda2': lambda2 * l0, 'epsilon': epsilon * e0}
