import numpy as np
import scipy.sparse
import scipy.stats
import wot


def score_gene_sets(dataset_to_score, gs, background_ds=None, method='mean_z_score', permutations=None,
                    nbins=25, bin_by='mean', random_state=0, drop_frequency=1000, drop_p_value_threshold=0.05,
                    smooth_p_values=True):
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
    if background_ds is None:
        background = x
    else:
        background = background_ds.x
    if gs.x.shape[1] > 1:
        raise ValueError('Only one gene set allowed as input')
    gs_1_0 = gs.x
    if permutations is not None and permutations > 0:
        if bin_by == 'mean' and method == 'mean_z_score':
            bin_by = 'std'
        bin_values = background.mean(axis=0) if bin_by == 'mean' else np.sqrt(wot.mean_and_variance(background)[1])
        # ranks go from 1 to len(bin_values)
        bin_ranks = scipy.stats.rankdata(bin_values, method='ordinal')
        bin_min = 1
        bin_max = 1 + len(bin_ranks)
        bin_width = (bin_max - bin_min) / nbins
        bin_assignments = np.floor((bin_ranks - bin_min) / bin_width)
        bin_index_to_gene_indices = []
        all_gene_indices = []
        missing_bins = False
        for bin_index in range(nbins):
            gene_indices = np.where(bin_assignments == bin_index)[0]
            if len(gene_indices) > 0 and np.sum(gs_1_0[gene_indices]) > 0:
                bin_index_to_gene_indices.append(gene_indices)
                all_gene_indices.append(gene_indices)
            else:
                missing_bins = True
        # permute each bin separately
        input_nbins = nbins
        nbins = len(bin_index_to_gene_indices)
        if nbins != input_nbins:
            print('Using ' + str(nbins) + ' out of ' + str(input_nbins) + ' bins')

    # preprocess the dataset
    if method == 'mean_z_score':
        # gs_x = gs_x[gene_indices]
        # background_x = background_x[:, gene_indices]
        # dataset_to_score.x = dataset_to_score.x[:, gene_indices]
        mean, var = wot.mean_and_variance(background.x)
        std = np.sqrt(var)
        x = (x - mean) / std
        x[np.isnan(x)] = 0
        x[x < -5] = -5
        x[x > 5] = 5
    elif method == 'mean_rank':
        ranks = np.zeros(dataset_to_score.x.shape)
        is_sparse = scipy.sparse.issparse(x)
        for i in range(dataset_to_score.x.shape[0]):
            row = x[i, :]
            if is_sparse:
                row = row.toarray()
            ranks[i] = scipy.stats.rankdata(row)
        x = ranks
    # gene sets has genes on rows, sets on columns
    # ds has cells on rows, genes on columns
    # scores contains cells on rows, gene sets on columns
    is_sparse = scipy.sparse.issparse(x)
    if is_sparse:
        gs_1_0 = scipy.sparse.csr_matrix(gs_1_0)
    observed_scores = x.dot(gs_1_0)
    if is_sparse:
        observed_scores = observed_scores.toarray()
    ngenes_in_set = gs_1_0.sum(axis=0)
    # ngenes_in_set[ngenes_in_set == 0] = 1  # avoid divide by zero
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
            if is_sparse:
                permuted_gs = scipy.sparse.csr_matrix(permuted_gs)
            permuted_scores = x[cells_to_keep].dot(permuted_gs)
            # count number of times permuted score is >= than observed score
            if is_sparse:
                permuted_scores = permuted_scores.toarray()
            p_values[cells_to_keep] += (permuted_scores >= observed_scores[cells_to_keep]).sum(axis=1)

            if do_drop:
                n = npermutations
                n_s = p_values
                n_f = n - n_s
                p_low = n_s / n - (z / n) * np.sqrt((n_s * n_f) / n)
                cells_to_keep = cells_to_keep & (p_low < drop_p_value_threshold)
                ncells = np.sum(cells_to_keep)
                print(
                    'permutation ' + str(total_permutations) + '/' + str(permutations) + ', ' + str(ncells) + '/' + str(
                        x.shape[0]) + ' cells')
                if ncells == 0:
                    break
            else:
                print(
                    'permutation ' + str(total_permutations) + '/' + str(permutations))

        observed_scores = observed_scores / ngenes_in_set
        k = p_values
        if smooth_p_values:
            p_values = (p_values + 1) / (npermutations + 2)
        else:
            p_values /= npermutations
        # FDR
        sorted_p_values = np.sort(p_values)
        n_features = len(p_values)
        fdr = sorted_p_values / n_features * np.arange(1, n_features + 1)
        return {'score': observed_scores, 'p_value': p_values, 'fdr': fdr, 'k': k, 'n': npermutations}
    return {'score': observed_scores}
