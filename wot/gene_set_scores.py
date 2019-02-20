import numpy as np
import scipy.sparse
import scipy.stats


def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)


# from http://www.statsmodels.org/dev/_modules/statsmodels/stats/multitest.html
def fdr(pvals, is_sorted=False, method='indep'):
    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    else:
        pvals_sorted = pvals  # alias

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))  # corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
    ##    elif method in ['n', 'negcorr']:
    ##        cm = np.sum(np.arange(len(pvals)))
    ##        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and negcorr implemented')

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected
        return pvals_corrected_
    else:
        return pvals_corrected


def get_p_value_ci(n, n_s, z):
    # smooth
    n = n + 2
    n_s = n_s + 1
    n_f = n - n_s
    ci = (z / n) * np.sqrt((n_s * n_f) / n)
    return ci


def score_gene_sets(ds, gs, method='mean_z_score', permutations=None,
                    random_state=0, smooth_p_values=True, progress=False):
    """Score gene sets.

    Note that datasets and gene sets must be aligned prior to invoking this method. No check is done.

    mean_z_score: Compute the z-score for each gene in the set. Truncate these z-scores at 5 or −5, and define the signature of the cell to be the mean z-score over all genes in the gene set.

    Parameters
    ----------

    random_state : `int`, optional (default: 0)
        The random seed for sampling.

    Returns
    -------
    Observed scores and permuted p-values if permutations > 0

    """

    if permutations is None:
        permutations = 0
    x = ds.X
    gs_1_0 = gs.X
    if not scipy.sparse.issparse(gs.X) and len(gs.X.shape) is 1:
        gs_1_0 = np.array([gs_1_0]).T

    if not scipy.sparse.issparse(gs_1_0):
        gs_1_0 = scipy.sparse.csr_matrix(gs_1_0)

    gs_indices = (gs_1_0 > 0)
    if hasattr(gs_indices, 'toarray'):
        gs_indices = gs_indices.toarray()
    gs_indices = gs_indices.flatten()
    gs_1_0 = gs_1_0[gs_indices]
    ngenes_in_set = gs_1_0.sum(axis=0)

    if len(x.shape) == 1:
        x = np.array([x])
    # preprocess the dataset
    if method == 'mean_z_score':
        x = x[:, gs_indices]  # only include genes in gene set
        if scipy.sparse.issparse(x):
            x = x.toarray()
        mean = x.mean(axis=0)
        var = x.var(axis=0)
        std = np.sqrt(var)
        x = (x - mean) / std
        x[np.isnan(x)] = 0
        x[x < -5] = -5
        x[x > 5] = 5
    elif method == 'mean_rank':  # need all genes for ranking
        ranks = np.zeros(x.shape)
        is_sparse = scipy.sparse.issparse(x)
        for i in range(x.shape[0]):  # rank each cell separately
            row = x[i, :]
            if is_sparse:
                row = row.toarray()
            ranks[i] = scipy.stats.rankdata(row, method='min')
        x = ranks
        x = x[:, gs_indices]
    else:
        x = x[:, gs_indices]
    observed_scores = x.mean(axis=1)
    if hasattr(observed_scores, 'toarray'):
        observed_scores = observed_scores.toarray()

    if progress:
        print('# of genes in gene set ' + str(ngenes_in_set))

    # gene sets has genes on rows, sets on columns
    # ds has cells on rows, genes on columns
    # scores contains cells on rows, gene sets on columns

    if permutations is not None and permutations > 0:
        if random_state:
            np.random.seed(random_state)
        p_values = np.zeros(x.shape[0])
        permuted_X = x.T.copy()  # put genes on rows to shuffle each row indendently
        for i in range(permutations):
            for _x in permuted_X:
                np.random.shuffle(_x)
            # count number of times permuted score is >= than observed score
            p_values += permuted_X.mean(axis=0) >= observed_scores
            if progress:
                print(
                    'permutation ' + str(i) + '/' + str(permutations))

        k = p_values
        if smooth_p_values:
            p_values = (p_values + 1) / (permutations + 2)
        else:
            p_values = p_values / permutations
        return {'score': observed_scores, 'p_value': p_values, 'fdr': fdr(p_values), 'k': k}

    return {'score': observed_scores}
