# -*- coding: utf-8 -*-

import numpy as np
import ot as pot
import scipy.sparse
import scipy.stats
import sklearn.decomposition
import time

import wot


def transport_stable_learn_growth(C, lambda1, lambda2, epsilon, scaling_iter, g, pp=None, qq=None, tau=None,
                                  epsilon0=None, growth_iters=3, inner_iter_max=None):
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

        Tmap = transport_stablev2(C=C, lambda1=lambda1, lambda2=lambda2, epsilon=epsilon,
                                  scaling_iter=scaling_iter, g=rowSums, tau=tau,
                                  epsilon0=epsilon0, pp=pp, qq=qq, numInnerItermax=inner_iter_max,
                                  extra_iter=1000)
    return Tmap


def transport_stablev_learn_growth_duality_gap(C, g, lambda1, lambda2, epsilon, batch_size, tolerance, tau, epsilon0,
                                               growth_iters, max_iter, pp=None, qq=None):
    """
    Compute the optimal transport with stabilized numerics and duality gap guarantee.

    Parameters
    ----------
    C : 2D array
    g : 1D array
    pp : 1D array
    qq : 1D array
    lambda1 : float
    lambda2 : float
    epsilon : float
    batch_size : int
    tolerance : float
    tau : float
    epsilon0 : float
    growth_iters : int
    max_iter : int

    Returns
    -------
    tmap : 2D array
        Transport map

    Notes
    -----
    It is guaranteed that the duality gap for the result is under the given threshold, or that max_iter have been performed.
    """
    if batch_size <= 0:
        raise ValueError("Batch size must be positive")

    row_sums = g
    for i in range(growth_iters):
        tmap = transport_stablev1(C, row_sums, pp, qq, lambda1, lambda2, epsilon,
                                  batch_size, tolerance, tau, epsilon0, max_iter)
        row_sums = tmap.sum(axis=1) / tmap.shape[1]
    return tmap


# @ Lénaïc Chizat 2015 - optimal transport
def fdiv(l, x, p, dx):
    return l * np.sum(dx * (x * (np.log(x / p)) - x + p))


def fdivstar(l, u, p, dx):
    return l * np.sum((p * dx) * (np.exp(u / l) - 1))


def primal(C, K, R, dx, dy, p, q, a, b, epsilon, lambda1, lambda2):
    I = len(p)
    J = len(q)
    F1 = lambda x, y: fdiv(lambda1, x, p, y)
    F2 = lambda x, y: fdiv(lambda2, x, q, y)
    with np.errstate(divide='ignore'):
        return F1(np.dot(R, dy), dx) + F2(np.dot(R.T, dx), dy) \
               + (epsilon * np.sum(R * np.nan_to_num(np.log(R)) - R + K) \
                  + np.sum(R * C)) / (I * J)


def dual(C, K, R, dx, dy, p, q, a, b, epsilon, lambda1, lambda2):
    I = len(p)
    J = len(q)
    F1c = lambda u, v: fdivstar(lambda1, u, p, v)
    F2c = lambda u, v: fdivstar(lambda2, u, q, v)
    return - F1c(- epsilon * np.log(a), dx) - F2c(- epsilon * np.log(b), dy) \
           - epsilon * np.sum(R - K) / (I * J)


# end @ Lénaïc Chizat

def transport_stablev1(C, g, pp, qq, lambda1, lambda2, epsilon, batch_size, tolerance, tau=10e100, epsilon0=1.,
                       max_iter=1e7):
    """
    Compute the optimal transport with stabilized numerics, with the guarantee that the duality gap is at most `tolerance`

    Parameters
    ----------
    C : 2-D ndarray
        The cost matrix. C[i][j] is the cost to transport cell i to cell j
    g : 1-D array_like
        Growth value for input cells.
    lambda1 : float, optional
        Regularization parameter for the marginal constraint on p
    lambda2 : float, optional
        Regularization parameter for the marginal constraint on q
    epsilon : float, optional
        Entropy regularization parameter.
    batch_size : int, optional
        Number of iterations to perform between each duality gap check
    tolerance : float, optional
        Upper bound on the duality gap that the resulting transport map must guarantee.
    tau : float, optional
        Threshold at which to perform numerical stabilization
    epsilon0 : float, optional
        Starting value for exponentially-decreasing epsilon
    max_iter : int, optional
        Maximum number of iterations. Print a warning and return if it is reached, even without convergence.

    Returns
    -------
    transport_map : 2-D ndarray
        The entropy-regularized unbalanced transport map
    """
    C = np.asarray(C, dtype=np.float64)
    epsilon_scalings = 5
    scale_factor = np.exp(- np.log(epsilon) / epsilon_scalings)

    I, J = C.shape
    dx, dy = np.ones(I) / I, np.ones(J) / J
    p = g
    if pp is not None:
        pp = pp / np.average(pp)
        p *= pp

    if qq is not None:
        qq = qq / np.average(qq)
        q = qq * np.sum(g * pp) / I
    else:
        q = np.ones(J) * np.average(g)
    u, v = np.zeros(I), np.zeros(J)
    a, b = np.ones(I), np.ones(J)

    start_time = time.time()
    duality_time = 0
    epsilon_i = epsilon0 * scale_factor
    current_iter = 0

    for e in range(epsilon_scalings + 1):
        duality_gap = np.inf
        u = u + epsilon_i * np.log(a)
        v = v + epsilon_i * np.log(b)  # absorb
        epsilon_i = epsilon_i / scale_factor
        _K = np.exp(-C / epsilon_i)
        alpha1 = lambda1 / (lambda1 + epsilon_i)
        alpha2 = lambda2 / (lambda2 + epsilon_i)
        K = np.exp((np.array([u]).T - C + np.array([v])) / epsilon_i)
        a, b = np.ones(I), np.ones(J)
        old_a, old_b = a, b
        threshold = tolerance if e == epsilon_scalings else 1e-6

        while duality_gap > threshold:
            for i in range(batch_size if e == epsilon_scalings else 5):
                current_iter += 1
                old_a, old_b = a, b
                a = (p / (K.dot(np.multiply(b, dy)))) ** alpha1 * np.exp(-u / (lambda1 + epsilon_i))
                b = (q / (K.T.dot(np.multiply(a, dx)))) ** alpha2 * np.exp(-v / (lambda2 + epsilon_i))

                # stabilization
                if (max(max(abs(a)), max(abs(b))) > tau):
                    wot.io.verbose("Stabilizing...")
                    u = u + epsilon_i * np.log(a)
                    v = v + epsilon_i * np.log(b)  # absorb
                    K = np.exp((np.array([u]).T - C + np.array([v])) / epsilon_i)
                    a, b = np.ones(I), np.ones(J)

                if current_iter >= max_iter:
                    print("Warning : Reached max_iter with duality gap still above threshold. Returning")
                    return (K.T * a).T * b

            # The real dual variables. a and b are only the stabilized variables
            _a = a * np.exp(u / epsilon_i)
            _b = b * np.exp(v / epsilon_i)

            # Skip duality gap computation for the first epsilon scalings, use dual variables evolution instead
            if e == epsilon_scalings:
                duality_tmp_time = time.time()
                R = (K.T * a).T * b
                pri = primal(C, _K, R, dx, dy, p, q, _a, _b, epsilon_i, lambda1, lambda2)
                dua = dual(C, _K, R, dx, dy, p, q, _a, _b, epsilon_i, lambda1, lambda2)
                duality_gap = (pri - dua) / abs(pri)
                duality_time += time.time() - duality_tmp_time
                # wot.io.verbose("Current (gap, primal, dual) : {:020.18f} {:020.18f} {:020.18f}".format(duality_gap, pri, dua))
            else:
                duality_gap = max(
                    np.linalg.norm(_a - old_a * np.exp(u / epsilon_i)) / (1 + np.linalg.norm(_a)),
                    np.linalg.norm(_b - old_b * np.exp(v / epsilon_i)) / (1 + np.linalg.norm(_b)))

    if np.isnan(duality_gap):
        raise RuntimeError("Overflow encountered in duality gap computation, please report this incident")
    total_time = time.time() - start_time
    wot.io.verbose("Computed tmap in {:.3f}s. Duality gap: {:.3E} ({:.2f}% of computing time)" \
                   .format(total_time, duality_gap, 100 * duality_time / total_time))
    return R


def transport_stablev2(C, lambda1, lambda2, epsilon, scaling_iter, g, pp, qq, numInnerItermax, tau,
                       epsilon0, extra_iter):
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

    # if pp is not None:
    #     pp = pp / np.average(pp)
    #     dx = dx * pp
    #
    # if qq is not None:
    #     qq = qq / np.average(qq)
    #     dy = dy * qq

    # p = g / np.average(g, weights=dx)
    # q = np.ones(C.shape[1])
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


def optimal_transport(cost_matrix, g=None, pp=None, qq=None, p=None, q=None, solver=None,
                      delta_days=1, epsilon=0.1, lambda1=1.,
                      lambda2=1., min_transport_fraction=0.05,
                      max_transport_fraction=0.4, min_growth_fit=0.9,
                      l0_max=100, scaling_iter=250, epsilon_adjust=1.1,
                      lambda_adjust=1.5, numItermax=100, epsilon0=100.0, numInnerItermax=10, tau=1000.0, stopThr=1e-06,
                      growth_iters=3):
    if g is None:
        growth_rate = np.ones(len(cost_matrix))
    else:
        growth_rate = g

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


def compute_pca(m1, m2, n_components):
    matrices = list()
    matrices.append(m1 if not scipy.sparse.isspmatrix(m1) else m1.toarray())
    matrices.append(m2 if not scipy.sparse.isspmatrix(m2) else m2.toarray())
    x = np.vstack(matrices)
    mean_shift = x.mean(axis=0)
    x = x - mean_shift
    pca = sklearn.decomposition.PCA(n_components=n_components, random_state=58951)
    pca.fit(x.T)
    comp = pca.components_.T
    m1_len = m1.shape[0]
    m2_len = m2.shape[0]
    pca_1 = comp[0:m1_len]
    pca_2 = comp[m1_len:(m1_len + m2_len)]
    return pca_1, pca_2, pca, mean_shift


def get_pca(dim, *args):
    """
    Get a PCA projector for the arguments.

    Parameters
    ----------
    dim : int
        The number of components to use for PCA, i.e. number of dimensions of the resulting space.
    *args : ndarray
        The points to compute dimensionality reduction for. Can be several sets of points.

    Returns
    -------
    pca : sklearn.decomposition.PCA
        A PCA projector. Use `pca.transform(x)` to get the PCA-space projection of x.

    Example
    -------
    >>> pca = get_pca(30, p0_x)
    >>> pca.transform(p0_x)
    >>> # -> project p0_x to PCA space
    >>> pca = get_pca(30, p0_x, p1_x)
    >>> p0, p05, p1 = [ pca.transform(x) for x in [p0_x, p05_x, p1_x] ]
    >>> # -> project all three sets of points to PCA space computed from p0_x and p1_x
    """
    args = [a.toarray() if scipy.sparse.isspmatrix(a) else a for a in args]
    x = np.vstack(args)
    mean = x.mean(axis=0)
    x = x - mean
    pca = sklearn.decomposition.PCA(n_components=dim)
    return pca.fit(x), mean


def pca_transform(pca, mean, arr):
    """
    Apply a PCA transformation to argument

    Parameters
    ----------
    pca : sklearn.decomposition.PCA
        A PCA projector. See wot.ot.get_pca
    arr : numpy ndarray or scipy.sparse matrix
        The array to project.

    Returns
    -------
    result : ndarray
    """
    ndarr = arr.toarray() if scipy.sparse.isspmatrix(arr) else arr
    ndarr = ndarr - mean
    return pca.transform(ndarr)
