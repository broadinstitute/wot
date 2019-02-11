# -*- coding: utf-8 -*-
import numpy as np
import ot as pot
import scipy.sparse
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


def point_cloud_distance(c1, c2, a=None, b=None, eigenvals=None):
    if eigenvals is not None:
        c1 = c1.dot(eigenvals)
        c2 = c2.dot(eigenvals)
    cloud_distances = sklearn.metrics.pairwise.pairwise_distances(c1, Y=c2, metric='sqeuclidean')

    if a is None:
        a = np.ones((cloud_distances.shape[0]), dtype=np.float64) / cloud_distances.shape[0]
    else:
        a = a / a.sum()
    if b is None:
        b = np.ones((cloud_distances.shape[1]), dtype=np.float64) / cloud_distances.shape[1]
    else:
        b = b / b.sum()
    return np.sqrt(pot.emd2(a, b, cloud_distances, numItermax=10000000))


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


def get_ids(expr):
    import os.path
    if not os.path.isfile(expr):
        set_ids = expr.split(',')
        return list(map(lambda x: x.strip(), set_ids))
    else:
        import pd
        return list(pd.read_csv(expr, index_col=0, header=None).index.values)


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


def split_in_two(n):
    indices = np.random.choice(n, int(n * 0.5))
    indices_c = np.zeros(n, dtype=bool)
    indices_c[indices] = True
    indices_c = np.invert(indices_c)
    return indices, indices_c


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


def earth_mover_distance(cloud1, cloud2, eigenvals):
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
    p = np.ones(len(cloud1)) / len(cloud1)
    q = np.ones(len(cloud2)) / len(cloud2)
    pairwise_dist = sklearn.metrics.pairwise.pairwise_distances(
        cloud1, Y=cloud2, metric='sqeuclidean')
    return np.sqrt(pot.emd2(p, q, pairwise_dist, numItermax=1e7))
