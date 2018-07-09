# -*- coding: utf-8 -*-
import numpy as np
import ot as pot
import sklearn.metrics


def compute_growth_scores(proliferation, apoptosis, beta_max=1.7, beta_center=0.25, delta_max=1.7, delta_min=0.15,
                          beta_min=0):
    birth = __beta(proliferation, beta_max=beta_max, center=beta_center, beta_min=beta_min)
    death = __delta(apoptosis, delta_max=delta_max, delta_min=delta_min)
    return np.exp(birth - death)


def __gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + __logistic(p, L=beta_max - beta_min, k=4 / width,
                                 x0=center)


def __beta(p, beta_max=1.7, beta_min=0, pmax=1.0, pmin=-0.5, center=0.25):
    # map proliferation p into interval beta_max, beta_min
    return __gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)


def __delta(a, delta_max=1.7, delta_min=0.15, amax=0.5, amin=-0.4, center=0.1):
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
        return list(pd.read_table(expr, index_col=0, header=None).index.values)


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
