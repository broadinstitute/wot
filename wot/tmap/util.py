import itertools

import anndata
import numpy as np
import pandas as pd
import scipy.sparse


def generate_comparisons(comparison_names, compare, days, delta_days=0, reference_day='start'):
    if compare != 'within':  # within, match, all, or trajectory name
        if compare == 'all':
            comparisons = itertools.combinations(comparison_names, 2)
        elif compare == 'match':
            base_name_to_full_names = {}
            comparisons = []
            for i in range(len(comparison_names)):
                full_trajectory_name = comparison_names[i]
                base_trajectory_name = full_trajectory_name[0:full_trajectory_name.rindex('/')]
                names = base_name_to_full_names.get(base_trajectory_name)
                if names is None:
                    names = []
                    base_name_to_full_names[base_trajectory_name] = names
                names.append(full_trajectory_name)
            for base_name in base_name_to_full_names:
                names = base_name_to_full_names[base_name]
                comparisons += list(itertools.combinations(names, 2))
        else:  # compare all to specified trajectory
            comparisons = itertools.product([compare], comparison_names)
        day_pairs = [(day, day) for day in days]
        return itertools.product(comparisons, day_pairs)
    else:
        # within
        filtered_day_pairs = []
        reference_index = 0 if reference_day == 'start' else len(days) - 1
        for day_index in range(len(days)):
            if day_index == reference_index:
                continue
            day2 = days[day_index]
            if delta_days > 0:
                requested_day = day2 - delta_days
                day1 = days[np.abs(days - requested_day).argmin()]
                if day1 == day2 or np.abs(
                        day1 - day2 - delta_days) > 0.1:  # too big or small a gap
                    continue
            else:
                day1 = days[reference_index]
            filtered_day_pairs.append((day1, day2))

    return itertools.product([(name, name) for name in comparison_names], filtered_day_pairs)


def unique_timepoint(*populations):
    """
    Retrieves the unique time point at which all populations live.

    Returns
    -------
    t : float
        The said timepoint

    Raises
    ------
    ValueError
        If all populations do not live in the same timepoint.
    """
    times = set([pop.time for pop in populations])
    if len(times) > 1:
        raise ValueError("Several populations were given, but they are not from the same day")
    if len(times) < 1:
        raise ValueError("No cells found at the given day")
    return list(times)[0]


def glue_transport_maps(tmap_0, tmap_1):
    """
    Glue two transport maps together

    Parameters
    ----------
    tmap_0 : anndata.AnnData
        The first transport map (from t0 to t1)
    tmap_1 : anndata.AnnData
        The second transport map (from t1 to t2)

    Returns
    -------
    result : anndata.AnnData
        The resulting transport map (from t0 to t2)
    """
    # FIXME: Column sum normalization is needed before gluing. Can be skipped only if lambda2 is high enough
    cells_at_intermediate_tpt = tmap_0.var.index
    cait_index = tmap_1.obs.index.get_indexer_for(cells_at_intermediate_tpt)
    result_x = np.dot(tmap_0.X, tmap_1.X[cait_index, :])
    return anndata.AnnData(result_x, tmap_0.obs.copy(), tmap_1.var.copy())


def trajectory_trends_from_trajectory(trajectory_ds, expression_ds, day_field='day'):
    """
    Computes the mean and variance of each gene over time for the given trajectories

    Parameters
    ----------
    trajectory_ds : anndata.AnnData
       anndata.AnnData returned by wot.tmap.TransportModel.trajectories
    expression_ds : anndata.AnnData
        Dataset used to compute mean
    day_field : str
        The day field name in trajectory_ds.obs

    Returns
    -------
    results : list of anndata.AnnData
        One dataset per trajectory. Each dataset has time on the rows and genes on the columns
    """

    # align gene expression matrix with trajectory matrix

    timepoints = []
    mean_list = []

    for j in range(trajectory_ds.shape[1]):
        mean_list.append(None)
        # variance_list.append(None)
    for day, group in trajectory_ds.obs.groupby(day_field):
        timepoints.append(day)
        trajectory_weights = trajectory_ds[group.index].X
        expression_values = expression_ds[group.index].X
        if scipy.sparse.isspmatrix(expression_values):
            expression_values = expression_values.toarray()
        # if inverse_log:
        #     expression_values = np.expm1(expression_values)

        for j in range(trajectory_ds.shape[1]):  # each trajectory
            weights_per_cell = trajectory_weights[:, j] if len(trajectory_weights.shape) > 1 else trajectory_weights
            mean = np.average(expression_values, weights=weights_per_cell, axis=0)
            # if inverse_log:
            #     mean = np.log1p(mean)
            # var = np.average((expression_values - mean) ** 2, weights=weights_per_cell, axis=0)

            if mean_list[j] is None:
                mean_list[j] = mean.T
                # variance_list[j] = var.T
            else:
                mean_list[j] = np.vstack((mean_list[j], mean.T))
                # variance_list[j] = np.vstack((variance_list[j], var.T))

    results = []
    for j in range(len(mean_list)):
        mean_ds = anndata.AnnData(mean_list[j], pd.DataFrame(index=timepoints), expression_ds.var.copy())
        # variance_ds = anndata.AnnData(variance_list[j], pd.DataFrame(index=timepoints), expression_ds.var.copy())
        # results.append((mean_ds, variance_ds))
        results.append(mean_ds)

    return results
