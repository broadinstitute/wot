import logging

import anndata
import numpy as np
import pandas as pd
import scipy.sparse
import statsmodels.stats.multitest
from scipy import stats

import wot.io

logger = logging.getLogger('wot')


def diff_exp(adata: anndata.AnnData, fate_datasets: [anndata.AnnData], cell_days_field: str = 'day',
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
    compare :
        within, match, all, or fate name
    delta_days
        Difference in days for comparisons. By default uses the closest day to compare to.

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

    df = None
    features = adata.var.index
    for comparison in comparisons:
        names = comparison[0]
        days = comparison[1]
        name1 = names[0]
        name2 = names[1]
        day1 = days[0]
        day2 = days[1]

        logger.info('{} vs {}, day {}, day {}'.format(name1, name2, day1, day2))
        values1, weights1 = __get_expression_and_weights(adata, cell_days_field, day1, name1)
        values2, weights2 = __get_expression_and_weights(adata, cell_days_field, day2, name2)
        result_df = __do_comparison(expression_values1=values1, weights1=weights1, day1=day1,
            expression_values2=values2,
            weights2=weights2,
            day2=day2,
            features=features)
        result_df['name1'] = name1
        result_df['name2'] = name2
        df = pd.concat((df, result_df)) if df is not None else result_df
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


def __do_comparison(expression_values1, weights1, day1, expression_values2, weights2, day2, features,
                    fraction_expressed_ratio_add=0.0001):
    mean1 = np.average(expression_values1, weights=weights1, axis=0)
    mean2 = np.average(expression_values2, weights=weights2, axis=0)
    fraction_expressed1 = weights1.dot(expression_values1 > 0)
    fraction_expressed2 = weights2.dot(expression_values2 > 0)
    fraction_expressed_diff = (fraction_expressed1 + fraction_expressed_ratio_add) / (
            fraction_expressed2 + fraction_expressed_ratio_add)

    fold_change = (mean1 - mean2)
    variance1 = np.average((expression_values1 - mean1) ** 2, weights=weights1, axis=0)
    variance2 = np.average((expression_values2 - mean2) ** 2, weights=weights2, axis=0)
    with np.errstate(invalid="ignore"):
        scores, pvals = stats.ttest_ind_from_stats(
            mean1=mean1, std1=np.sqrt(variance1), nobs1=len(weights1),
            mean2=mean2, std2=np.sqrt(variance2), nobs2=len(weights2), equal_var=False)  # Welch's
    scores[np.isnan(scores)] = 0
    pvals[np.isnan(pvals)] = 1
    fdr = statsmodels.stats.multitest.multipletests(pvals)[1]
    results = pd.DataFrame(index=features,
        data={'fold_change': np.exp(fold_change),
              't_score': scores,
              'p_val': pvals,
              'fdr': fdr,
              'fraction_expressed_ratio': fraction_expressed_diff,
              'day1': day1,
              'day2': day2})

    return results
