# -*- coding: utf-8 -*-
import numpy as np
import ot as pot
import scipy.sparse
import scipy.sparse
import scipy.stats
import sklearn.decomposition
import sklearn.metrics


def compute_growth_scores(proliferation, apoptosis, beta_max=1.7, beta_center=0.25, delta_max=1.7, delta_min=0.3,
                          beta_min=0.3):
    birth = __beta(proliferation, beta_max=beta_max, center=beta_center, beta_min=beta_min)
    death = __delta(apoptosis, delta_max=delta_max, delta_min=delta_min)
    return np.exp(birth - death)


def __gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + __logistic(p, L=beta_max - beta_min, k=4 / width,
                                 x0=center)


def __beta(p, beta_max=1.7, beta_min=0, pmax=1.0, pmin=-0.5, center=0.25):
    # map proliferation p into interval beta_max, beta_min
    return __gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)


def __delta(a, delta_max=1.7, delta_min=0.3, amax=0.5, amin=-0.4, center=0.1):
    return __gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)


def __logistic(x, L=1, k=1, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f


def sample_randomly(exp1, exp2, tm, g, npairs):
    # if args.npairs is None or args.npairs <= 0:
    #     l = tm / tm.sum(axis=0)
    #     l = l.flatten()
    #     qq = np.where(l > args.tm_thresh)
    #     n = len(qq[0])
    # else:
    n = npairs

    # column_sums = np.sum(tm, axis=0)
    # row_sums = np.sum(tm, axis=1)
    # p = np.outer(row_sums, column_sums)

    p = g
    q = np.ones(exp2.shape[0]) * np.average(g)
    p = np.outer(p, q)
    p = p.flatten()
    probabilties = p / p.sum()
    # pairs = np.random.multinomial(n, p, size=1)
    # pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))

    s = np.random.multinomial(n, probabilties, size=1)
    reshaped_s = s.reshape(exp1.shape[0], exp2.shape[0])
    pairs = np.nonzero(reshaped_s)
    weights = reshaped_s[pairs]
    return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]],
            'indices0': pairs[0],
            'indices1': pairs[1],
            'weights': weights}


def sample_uniformly(exp1, exp2, tm, npairs):
    # if args.npairs is None or args.npairs <= 0:
    #     l = tm / tm.sum(axis=0)
    #     l = l.flatten()
    #     qq = np.where(l > args.tm_thresh)
    #     n = len(qq[0])
    # else:
    n = npairs

    p = np.ones(exp1.shape[0] * exp2.shape[0])
    p = p / p.sum()
    pairs = np.random.multinomial(n, p, size=1)
    pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
    return exp1[pairs[0]], exp2[pairs[1]]


def sample_from_transport_map(exp1, exp2, tm, npairs, t_interpolate):
    probabilties = tm / np.power(tm.sum(axis=0), 1.0 - t_interpolate)
    probabilties = probabilties.flatten()
    # if args.npairs is None or args.npairs <= 0:
    #     # l = l / l.sum()
    #     reshaped_tmap = probabilties.reshape(exp1.shape[0], exp2.shape[0])
    #     q = np.where(reshaped_tmap > args.tm_thresh)
    #     qq = np.where(probabilties > args.tm_thresh)
    #     return exp1[q[0]], exp2[q[1]], probabilties[qq]
    # else:
    probabilties = probabilties / probabilties.sum()
    # pairs = np.random.multinomial(args.npairs, probabilties, size=1)
    # pairs = np.nonzero(pairs.reshape(exp1.shape[0], exp2.shape[0]))
    #  return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]], 'indices0': pairs[0], 'indices1': pairs[1], 'weights': None}
    s = np.random.multinomial(npairs, probabilties, size=1)
    reshaped_s = s.reshape(exp1.shape[0], exp2.shape[0])
    pairs = np.nonzero(reshaped_s)
    weights = reshaped_s[pairs]
    return {'pc0': exp1[pairs[0]], 'pc1': exp2[pairs[1]],
            'indices0': pairs[0],
            'indices1': pairs[1],
            'weights': weights}


def interpolate_with_ot(p0, p1, tmap, interp_frac, size):
    """
    Interpolate between p0 and p1 at fraction t_interpolate knowing a transport map from p0 to p1

    Parameters
    ----------
    p0 : 2-D array
        The genes of each cell in the source population
    p1 : 2-D array
        The genes of each cell in the destination population
    tmap : 2-D array
        A transport map from p0 to p1
    t_interpolate : float
        The fraction at which to interpolate
    size : int
        The number of cells in the interpolated population

    Returns
    -------
    p05 : 2-D array
        An interpolated population of 'size' cells
    """
    p0 = p0.toarray() if scipy.sparse.isspmatrix(p0) else p0
    p1 = p1.toarray() if scipy.sparse.isspmatrix(p1) else p1
    p0 = np.asarray(p0, dtype=np.float64)
    p1 = np.asarray(p1, dtype=np.float64)
    tmap = np.asarray(tmap, dtype=np.float64)
    if p0.shape[1] != p1.shape[1]:
        raise ValueError("Unable to interpolate. Number of genes do not match")
    if p0.shape[0] != tmap.shape[0] or p1.shape[0] != tmap.shape[1]:
        raise ValueError("Unable to interpolate. Tmap size is {}, expected {}"
                         .format(tmap.shape, (len(p0), len(p1))))
    I = len(p0);
    J = len(p1)
    # Assume growth is exponential and retrieve growth rate at t_interpolate
    p = tmap / np.power(tmap.sum(axis=0), 1. - interp_frac)
    p = p.flatten(order='C')
    p = p / p.sum()
    choices = np.random.choice(I * J, p=p, size=size)
    return np.asarray([p0[i // J] * (1 - interp_frac) + p1[i % J] * interp_frac for i in choices], dtype=np.float64)


def interpolate_randomly(p0, p1, t_interpolate, size):
    """
    Interpolate between p0 and p1 at fraction t_interpolate

    Parameters
    ----------
    p0 : 2-D array
        The genes of each cell in the source population
    p1 : 2-D array
        The genes of each cell in the destination population
    t_interpolate : float
        The fraction at which to interpolate
    size : int
        The number of cells in the interpolated population

    Returns
    -------
    p05 : 2-D array
        An interpolated population of 'size' cells
    """
    t = t_interpolate
    p0 = p0.toarray() if scipy.sparse.isspmatrix(p0) else p0
    p1 = p1.toarray() if scipy.sparse.isspmatrix(p1) else p1
    p0 = np.asarray(p0, dtype=np.float64)
    p1 = np.asarray(p1, dtype=np.float64)
    if p0.shape[1] != p1.shape[1]:
        raise ValueError("Unable to interpolate. Number of genes do not match")
    I = len(p0)
    J = len(p1)
    choices = np.random.choice(I * J, size=size)
    return np.asarray([p0[i // J] * (1 - t) + p1[i % J] * t for i in choices], dtype=np.float64)


def interpolate_randomly_with_growth(p0, p1, t, size, g):
    p0 = p0.toarray() if scipy.sparse.isspmatrix(p0) else p0
    p1 = p1.toarray() if scipy.sparse.isspmatrix(p1) else p1
    p0 = np.asarray(p0, dtype=np.float64)
    p1 = np.asarray(p1, dtype=np.float64)

    p = g
    q = np.ones(p1.shape[0]) * np.average(g)

    p = np.outer(p, q)
    p = p.flatten(order='C')
    p = p / p.sum()
    I = len(p0)
    J = len(p1)
    choices = np.random.choice(I * J, p=p, size=size)
    return np.asarray([p0[i // J] * (1 - t) + p1[i % J] * t for i in choices], dtype=np.float64)


def earth_mover_distance(cloud1, cloud2, eigenvals=None, weights1=None, weights2=None):
    """
    Returns the earth mover's distance between two point clouds

    Parameters
    ----------
    cloud1 : 2-D array
        First point cloud
    cloud2 : 2-D array
        Second point cloud

    Returns
    -------
    distance : float
        The distance between the two point clouds
    """
    cloud1 = cloud1.toarray() if scipy.sparse.isspmatrix(cloud1) else cloud1
    cloud2 = cloud2.toarray() if scipy.sparse.isspmatrix(cloud2) else cloud2
    if eigenvals is not None:
        cloud1 = cloud1.dot(eigenvals)
        cloud2 = cloud2.dot(eigenvals)
    if weights1 is None:
        p = np.ones(len(cloud1)) / len(cloud1)
    else:
        weights1 = weights1.astype('float64')
        p = weights1 / weights1.sum()

    if weights2 is None:
        q = np.ones(len(cloud2)) / len(cloud2)
    else:
        weights2 = weights2.astype('float64')
        q = weights2 / weights2.sum()

    pairwise_dist = sklearn.metrics.pairwise.pairwise_distances(
        cloud1, Y=cloud2, metric='sqeuclidean')
    return np.sqrt(pot.emd2(p, q, pairwise_dist, numItermax=1e7))


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
