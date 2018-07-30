import numpy as np
import scipy.sparse
import scipy.stats
import wot


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
    n_f = n - n_s
    ci = (z / n) * np.sqrt((n_s * n_f) / n)
    return ci


def score_gene_sets(dataset_to_score, gs, method='mean_z_score', permutations=None,
                    nbins=25, bin_by='mean', random_state=0, drop_frequency=1000, drop_p_value_threshold=0.05,
                    smooth_p_values=True, progress=False, quantile_bins=False):
    """Score gene sets.

    Note that datasets and gene sets must be aligned prior to invoking this method. No check is done.

    mean_z_score: Compute the z-score for each gene in the set. Truncate these z-scores at 5 or −5, and define the signature of the cell to be the mean z-score over all genes in the gene set.

    Parameters
    ----------
    nbins : `int`, optional
        Number of bins for sampling.
    random_state : `int`, optional (default: 0)
        The random seed for sampling.

    Returns
    -------
    Observed scores and permuted p-values if permutations > 0

    """
    # gene_indices = (gs_x.sum(axis=1) > 0) & (
    #         background_x.std(axis=0) > 0)  # keep genes that are in gene sets and have standard deviation > 0

    x = dataset_to_score.x
    # preprocess the dataset
    if method == 'mean_z_score':
        # gs_x = gs_x[gene_indices]
        # background_x = background_x[:, gene_indices]
        # dataset_to_score.x = dataset_to_score.x[:, gene_indices]
        mean, var = wot.mean_and_variance(x)
        std = np.sqrt(var)
        x = (x - mean) / std
        x[np.isnan(x)] = 0
        x[x < -5] = -5
        x[x > 5] = 5
    elif method == 'mean_rank':
        ranks = np.zeros(x.shape)
        is_sparse = scipy.sparse.issparse(x)
        for i in range(dataset_to_score.x.shape[0]):  # rank each cell separately
            row = x[i, :]
            if is_sparse:
                row = row.toarray()
            ranks[i] = scipy.stats.rankdata(row, method='min')
        x = ranks

    if gs.x.shape[1] > 1:
        raise ValueError('Only one gene set allowed as input')
    gs_1_0 = gs.x
    if permutations is not None and permutations > 0:
        if bin_by == 'mean' and method == 'mean_z_score':
            bin_by = 'std'
        _bin_values = x.mean(axis=0) if bin_by == 'mean' else np.sqrt(wot.mean_and_variance(x)[1])
        bin_values = _bin_values
        if quantile_bins:
            # ranks go from 1 to len(bin_values)
            bin_values = scipy.stats.rankdata(bin_values, method='min')

        bin_min = np.min(bin_values)
        bin_max = np.max(bin_values)
        bin_width = (bin_max - bin_min) / nbins
        bin_assignments = np.floor((bin_values - bin_min) / bin_width)
        bin_assignments[bin_assignments == nbins] = nbins - 1
        bin_index_to_gene_indices = []
        all_gene_indices = []
        missing_bins = False
        for bin_index in range(nbins):
            gene_indices = np.where(bin_assignments == bin_index)[0]
            ngenes_in_set_bin = np.sum(gs_1_0[gene_indices])
            if progress:
                print('bin ' + str(bin_index) + ' ' + str(ngenes_in_set_bin) + '/' + str(
                    len(gene_indices)) + ', bin range: ' + str(np.min(_bin_values[:, gene_indices])) + ' ' + str(
                    np.max(_bin_values[:, gene_indices])))
            if len(gene_indices) > ngenes_in_set_bin:
                bin_index_to_gene_indices.append(gene_indices)
                all_gene_indices.append(gene_indices)
            else:
                if ngenes_in_set_bin > 0:
                    print('Unable to use bin ' + str(bin_index))
                missing_bins = True
        # permute each bin separately
        input_nbins = nbins
        nbins = len(bin_index_to_gene_indices)
        if nbins != input_nbins:
            print('Using ' + str(nbins) + ' out of ' + str(input_nbins) + ' bins')

    # gene sets has genes on rows, sets on columns
    # ds has cells on rows, genes on columns
    # scores contains cells on rows, gene sets on columns
    if not scipy.sparse.issparse(gs_1_0):
        gs_1_0 = scipy.sparse.csr_matrix(gs_1_0)
    observed_scores = x @ gs_1_0
    if hasattr(observed_scores, 'toarray'):
        observed_scores = observed_scores.toarray()
    ngenes_in_set = gs_1_0.sum(axis=0)
    if progress:
        print('# of genes ' + str(ngenes_in_set))
    # ngenes_in_set[ngenes_in_set == 0] = 1  # avoid divide by zero
    p_value_ci = None
    if permutations is not None and permutations > 0:
        if random_state:
            np.random.seed(random_state)
        p_values = np.zeros(x.shape[0])
        z = scipy.stats.norm.ppf(0.99)
        npermutations = np.zeros(x.shape[0])  # number of permutations per cell
        if missing_bins:  # keep genes in bins only
            x = x[:, np.concatenate(all_gene_indices)]
        # keep track of cells that have low p-values
        cells_to_keep = np.ones(x.shape[0], dtype=bool)
        do_drop = drop_frequency > 0
        block_size = drop_frequency if do_drop else min(1000, permutations)
        bin_sizes = []
        for bin_index in range(nbins):
            subset = gs_1_0[bin_index_to_gene_indices[bin_index]]
            bin_sizes.append(subset.shape[0])
        total_permutations = 0
        for i in range(0, permutations, block_size):
            current_block_size = min(i + block_size, permutations) - i
            total_permutations += current_block_size
            npermutations[cells_to_keep] += current_block_size
            permuted_gs = np.zeros((x.shape[1], current_block_size))
            bin_start = 0
            bin_end = 0
            for bin_index in range(nbins):
                subset = gs_1_0[bin_index_to_gene_indices[bin_index]]
                bin_size = subset.shape[0]
                nchoose = subset.sum()
                bin_end += bin_size
                for j in range(current_block_size):
                    # np.random.shuffle(subset)
                    # permuted_gs[:, j] = subset[:, 0]
                    s = np.zeros(bin_size)
                    s[np.random.choice(bin_size, nchoose, replace=False)] = 1
                    permuted_gs[bin_start:bin_end, j] = s
                bin_start = bin_end
                # shuffle each column independently
            # permuted scores has cells on rows, gene set on columns repeated current_block_size times
            permuted_gs = scipy.sparse.csr_matrix(permuted_gs)
            permuted_scores = x[cells_to_keep] @ permuted_gs
            # count number of times permuted score is >= than observed score
            if hasattr(permuted_scores, 'toarray'):
                permuted_scores = permuted_scores.toarray()
            p_values[cells_to_keep] += (permuted_scores >= observed_scores[cells_to_keep]).sum(axis=1)

            if do_drop:
                p_value_ci = get_p_value_ci(npermutations, p_values, z)
                n_s = p_values
                p_low = n_s / npermutations - p_value_ci
                cells_to_keep = cells_to_keep & (p_low < drop_p_value_threshold)
                ncells = np.sum(cells_to_keep)
                if progress:
                    print(
                        'permutation ' + str(total_permutations) + '/' + str(permutations) + ', ' + str(
                            ncells) + '/' + str(
                            x.shape[0]) + ' cells')
                if ncells == 0:
                    break
            elif progress:
                print(
                    'permutation ' + str(total_permutations) + '/' + str(permutations))

        observed_scores = observed_scores / ngenes_in_set
        k = p_values
        if smooth_p_values:
            p_values = (p_values + 1) / (npermutations + 2)
        else:
            p_values /= npermutations
        fdr_low = None
        fdr_high = None
        if p_value_ci is not None:
            p_value_ci[p_value_ci > 1] = 1
            p_low = k / npermutations - p_value_ci
            p_low[p_low < 0] = 0
            fdr_low = fdr(p_low)
            p_high = k / npermutations + p_value_ci
            p_high[p_high > 1] = 1
            fdr_high = fdr(p_high)
        return {'score': observed_scores, 'p_value': p_values, 'fdr': fdr(p_values), 'k': k, 'n': npermutations,
                'p_value_ci': p_value_ci, 'fdr_low': fdr_low, 'fdr_high': fdr_high}
    return {'score': observed_scores}
