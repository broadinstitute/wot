# -*- coding: utf-8 -*-

import numpy as np

def interp(t, tp, fp, left=None, right=None, method='linear', smooth=0):
    """
    Multi-dimensional linear interpolation.
    Returns the N-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    t : array_like
        The t-coordinates of the interpolated values.
    tp : 1-D sequence of floats
        The t-coordinates of the data points, must be increasing.
    fp : N-D sequence of floats.
        The x-coordinates of the data points, same length as `tp`.
    left : optional float or complex corresponding to fp
        Value to return for `t < tp[0]`, default is `fp[0]`.
    right : optional float or complex corresponding to fp
        Value to return for `t > tp[-1]`, default is `fp[-1]`.
    method : either 'linear' or 'quadratic', optional
        Spline to use when fitting the points
    smooth : int, optional
        If above 0, number of points to use for running average smoothing

    Returns
    -------
    x : ndarray or array of ndarray
        The interpolated value(s).

    Raises
    ------
    ValueError
        If `tp` and `fp` have different length
    ValueError
        If all items in `fp` do not have the same dimension
    ValueError
        If `tp` is not a 1-D sequence

    Notes
    -----
    Does not check that the x-coordinate sequence `tp` is increasing.
    If `tp` is not increasing, the results are nonsense.
    A simple check for increasing is : np.all(np.diff(tp) > 0)
    """

    return_array = True
    if isinstance(t, (float, int)):
        return_array = False
        t = [t]
    t  = np.asarray(t,  dtype = np.float64)
    tp = np.asarray(tp, dtype = np.float64)
    fp = np.asarray(fp, dtype = np.float64)

    if tp.ndim != 1:
        raise ValueError("Timepoint sequence tp must be a 1-D sequence")
    if tp.shape[0] != fp.shape[0]:
        raise ValueError("tp and fp are not of the same length")

    if left is None:
        left = fp[0]
    if right is None:
        right = fp[-1]

    if return_array:
        return __interp_func(t, tp, fp, left, right, method, smooth=smooth)
    else:
        return __interp_func(t, tp, fp, left, right, method).item()


def __interp_func(t_seq, tp, fp, left, right, method='linear', smooth=0):
    x = []
    for t in t_seq:
        z = np.zeros_like(fp[0])
        i = 0
        if tp[0] > t:
            x.append(left)
            continue
        while i + 1 < len(tp) and tp[i + 1] < t:
            z = - z + 2 * (fp[i+1] - fp[i]) / (tp[i+1] - tp[i])
            i += 1
        if i + 1 == len(tp):
            x.append(right)
            continue
        if method == 'linear':
            r = ( t - tp[i] ) / ( tp[i + 1] - tp[i] )
            x.append(np.asarray(fp[i] * (1 - r) + fp[i+1] * r, dtype=np.float64))
        elif method == 'quadratic':
            z_1 = 2 * (fp[i+1] - fp[i]) / (tp[i+1] - tp[i]) - z
            x.append(np.asarray(
                fp[i] + z * (t - tp[i]) + (z_1 - z) * ((t - tp[i]) ** 2) / (2 * (tp[i+1] - tp[i])),
                dtype=np.float64))
        else:
            raise ValueError("Unkown interpolation method")
    # running mean
    if smooth > 0:
        smooth += 1 - smooth % 2
        i = (smooth - 1)//2
        c = np.cumsum(np.concatenate([x[:i+1], x, x[-i:]]), axis=0)
        x = (c[smooth:] - c[:-smooth]) / float(smooth)
    return x

def multivariate_normal_mixture(means, covs, p=None, size=1):
    """
    Draw random samples from a mixture of multivariate normal distributions

    Parameters
    ----------
    means : (k, N) ndarray
        Means of the k N-dimensional distributions
    covs : (k, N, N) ndarray
        Covariance matrices of the distributions.
        Alternatively, a (k,N) ndarray may be passed to indicate diagonal covariance matrices.
        A (k,) ndarray will be interpreted as multiples of the identity matrix
    p : 1-D array_like, optional
        The probability associated with each distribution. Defaults to uniform
    size : int, optional
        Number of samples to return. Defaults to a single one if unspecified

    Returns
    -------
    x : (N,) ndarray or (size, N) ndarray
        The samples, as an array if size was specified and above 1

    Raises
    ------
    ValueError
        If `means`, `covs`, and `p` have incompatible sizes
        If `p` does not sum to 1
        If `size` is not positive

    """
    means = np.asarray(means, dtype=np.float64)
    if isinstance(covs, (int, float)):
        covs = [ covs * np.identity(means.shape[1]) ] * len(means)
    covs = np.asarray(covs, dtype=np.float64)
    if covs.ndim == 1:
        covs = np.asarray([ c * np.identity(means.shape[1]) for c in covs ])
    elif covs.ndim == 2:
        covs = np.asarray([ np.diag(cov) for cov in covs])

    if p is None:
        p = np.ones(len(means)) / len(means)
    p = np.asarray(p, dtype=np.float64)

    if means.shape[0] != covs.shape[0]:
        raise ValueError("means and covs are not of the same length")
    if means.shape[0] != p.shape[0]:
        raise ValueError("means and p are not of the same length")
    if not np.isclose(sum(p), 1):
        raise ValueError("p does not sum to 1")
    if not size > 0:
        raise ValueError("size is not positive")

    result = np.empty([size, means.shape[1]], dtype=np.float64)
    picks = np.random.multinomial(size, p)
    c = 0
    for i in range(len(p)):
        result[c:c+picks[i]] = \
                np.random.multivariate_normal(means[i], covs[i], size=picks[i])
        c += picks[i]

    if size == 1:
        return result.item()
    else:
        return result
