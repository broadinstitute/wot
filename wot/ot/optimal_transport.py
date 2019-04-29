# -*- coding: utf-8 -*-

import time

import numpy as np
import scipy.sparse
import scipy.stats
import sklearn.decomposition

import wot


def compute_transport_matrix(solver, **params):
    """
    Compute the optimal transport with stabilized numerics.
    Args:
    g: Growth
    solver: transport_stablev2 or optimal_transport_duality_gap
    growth_iters:
  """

    import gc
    g = params['g']
    growth_iters = params['growth_iters']
    for i in range(growth_iters):
        if i == 0:
            row_sums = g
        else:
            row_sums = tmap.sum(axis=1) / tmap.shape[1]
        params['g'] = row_sums
        tmap = solver(**params)
        gc.collect()
    return tmap, row_sums


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

def optimal_transport_duality_gap(C, g, lambda1, lambda2, epsilon, batch_size, tolerance, tau,
                                  epsilon0, max_iter, **ignored):
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
    q = np.ones(C.shape[1]) * np.average(g)

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


def transport_stablev2(C, lambda1, lambda2, epsilon, scaling_iter, g, tau, epsilon0, extra_iter, inner_iter_max,
                       **ignored):
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

        if (warm_start and iterations_since_epsilon_adjusted == inner_iter_max):
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

    R = (K.T * a).T * b

    return R


def compute_pca(m1, m2, n_components):
    matrices = list()
    matrices.append(m1 if not scipy.sparse.isspmatrix(m1) else m1.toarray())
    matrices.append(m2 if not scipy.sparse.isspmatrix(m2) else m2.toarray())
    x = np.vstack(matrices)
    mean_shift = x.mean(axis=0)
    x = x - mean_shift
    n_components = min(n_components, x.shape[0])  # n_components must be <= ncells
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
