import logging

import numpy as np
import pandas as pd
import scipy.sparse
import statsmodels.stats.multitest

import wot.io

logger = logging.getLogger('wot')

import anndata


def diff_exp(adata: anndata.AnnData, fate_datasets: [anndata.AnnData], cell_days_field: str = 'day',
             nperm: int = 0, min_fold_change: float = 0, smooth_p_values: bool = False,
             compare: str = 'all', delta_days: float = None) -> {str: pd.DataFrame}:
    """
    Compute statistics for differential expression

    Parameters
    ----------
    adata :
        The expression matrix.
    fate_datasets :
        List of fate datasets.
    cell_days_field
        Cell days field in adata.
    nperm
        Number of permutations to perform
    min_fold_change
        Exclude genes for permutations that don't have at least min_fold_change
    smooth_p_values
        Whether to smooth p-values
    compare :
        within, match, all, or fate name
    delta_days
        Difference in days for comparisons

    Returns
    -------
    df :
        A dataframe with the columns name1, name2, day1, day2, distance
    """

    if cell_days_field not in adata.obs:
        raise ValueError(cell_days_field + ' not found')
    adata = anndata.AnnData(adata.X, adata.obs.copy(), adata.var)
    fate_names = []
    if not isinstance(fate_datasets, list) and not isinstance(fate_datasets, tuple):
        fate_datasets = [fate_datasets]

    # add fate dataset to obs
    for fate_ds in fate_datasets:
        fate_names += list(fate_ds.var.index)
        adata.obs = adata.obs.join(
            pd.DataFrame(index=fate_ds.obs.index, data=fate_ds.X, columns=fate_ds.var.index))

    unique_days = np.array(sorted(adata.obs[cell_days_field].unique().astype(float)))
    unique_days = unique_days[np.isnan(unique_days) == False]
    logger.info('{} days'.format(len(unique_days)))

    comparisons = wot.tmap.generate_comparisons(comparison_names=fate_names, compare=compare,
                                                days=unique_days, delta_days=delta_days, reference_day='start')

    current_name1 = None
    current_name2 = None
    df = None
    features = adata.var.index
    results = {}
    for comparison in comparisons:
        names = comparison[0]
        days = comparison[1]
        name1 = names[0]
        name2 = names[1]
        day1 = days[0]
        day2 = days[1]
        if current_name1 != name1 or current_name2 != name2:
            if df is not None:
                results['{}_{}'.format(current_name1, current_name2).replace('/', '-')] = df

            df = pd.DataFrame(index=features)
            current_name1 = name1
            current_name2 = name2

        logger.info('{} vs {}, day {}, day {}'.format(name1, name2, day1, day2))
        values1, weights1 = __get_expression_and_weights(adata, cell_days_field, day1, name1)
        values2, weights2 = __get_expression_and_weights(adata, cell_days_field, day2, name2)

        df = __add_stats(values1, weights1, df, '_{}_{}'.format(name1, day1), features)

        df = __add_stats(values2, weights2, df, '_{}_{}'.format(name2, day2), features)

        df = df.join(__do_comparison(values1, weights1, day1, values2, weights2, day2, nperm=nperm,
                                     min_fold_change=min_fold_change,
                                     smooth_p_values=smooth_p_values, features=features))

    if df is not None:
        results['{}_{}'.format(current_name1, current_name2).replace('/', '-')] = df
    return results


def __add_stats(expression_values, weights, df, suffix, features):
    # expression_values = np.expm1(expression_values)
    mean = np.average(expression_values, weights=weights, axis=0)
    fraction_expressed = weights.dot(expression_values > 0)
    # variance = np.average((expression_values - mean) ** 2, weights=weights, axis=0)
    # variance = np.log1p(variance)
    # mean = np.log1p(mean)
    if 'mean{}'.format(suffix) not in df:
        return df.join(pd.DataFrame(index=features,
                                    data={
                                        'mean{}'.format(suffix): mean,
                                        'fraction_expressed{}'.format(suffix): fraction_expressed
                                    }))
    return df


def __get_expression_and_weights(adata, cell_days_field, day, fate_name):
    ds = adata[
        (adata.obs[cell_days_field] == day) & (
                False == adata.obs[fate_name].isna())]
    weights = ds.obs[fate_name].values
    expression_values = ds.X
    if scipy.sparse.isspmatrix(expression_values):
        expression_values = expression_values.toarray()
    weights = weights / weights.sum()
    return expression_values, weights


def __do_comparison(expression_values1, weights1, day1, expression_values2, weights2, day2, nperm, min_fold_change,
                    smooth_p_values, features):
    # expression_values1 = np.expm1(expression_values1)
    # expression_values2 = np.expm1(expression_values2)
    mean1 = np.average(expression_values1, weights=weights1, axis=0)
    mean2 = np.average(expression_values2, weights=weights2, axis=0)
    # variance1 = np.average((expression_values1 - mean1) ** 2, weights=weights1, axis=0)
    # variance2 = np.average((expression_values2 - mean2) ** 2, weights=weights2, axis=0)
    # fold_change = np.log1p(mean1) - np.log1p(mean2)
    observed = (mean1 - mean2)
    suffix = "_{}_{}".format(day1, day2)
    results = pd.DataFrame(index=features, data={'fold_change' + suffix: observed})

    if nperm is not None and nperm > 0:
        genes_use = np.abs(observed) >= min_fold_change
        if genes_use.sum() > 0:
            expression_values1 = expression_values1[:, genes_use]
            expression_values2 = expression_values2[:, genes_use]
            observed_use = observed[genes_use]
            weights1 = weights1.copy()
            weights2 = weights2.copy()
            p = np.zeros(expression_values2.shape[1])

            for i in range(nperm):
                np.random.shuffle(weights1)
                np.random.shuffle(weights2)
                mean1 = np.average(expression_values1, weights=weights1, axis=0)
                mean2 = np.average(expression_values2, weights=weights2, axis=0)

                # permuted_fold_change = np.log1p(mean1) - np.log1p(mean2)
                permuted = (mean1 - mean2)
                p[(permuted >= observed_use)] += 1
            # 2-sided p-value
            k = p
            if smooth_p_values:
                p = (p + 1) / (nperm + 2)
            else:
                p = p / nperm
            one_minus_p = 1.0 - p;
            expr = one_minus_p < p
            p[expr] = 1 - p[expr]
            p *= 2
            fdr = statsmodels.stats.multitest.multipletests(p)[1]
            results = results.join(
                pd.DataFrame(index=features[genes_use], data={'p_value' + suffix: p,
                                                              'fdr' + suffix: fdr,
                                                              'k' + suffix: k}))
    return results
