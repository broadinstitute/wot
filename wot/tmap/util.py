import numpy as np

import wot.io


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
        raise ValueError("Several populations were given, but they do not all live in the same timepoint")
    if len(times) < 1:
        raise ValueError("Empty population list")
    return list(times)[0]


def glue_transport_maps(tmap_0, tmap_1):
    """
    Glue two transport maps together

    Parameters
    ----------
    tmap_0 : wot.Dataset
        The first transport map (from t0 to t1)
    tmap_1 : wot.Dataset
        The second transport map (from t1 to t2)

    Returns
    -------
    result : wot.Dataset
        The resulting transport map (from t0 to t2)
    """
    # FIXME: Column sum normalization is needed before gluing. Can be skipped only if lambda2 is high enough
    cells_at_intermediate_tpt = tmap_0.col_meta.index
    cait_index = tmap_1.row_meta.index.get_indexer_for(cells_at_intermediate_tpt)
    result_x = np.dot(tmap_0.x, tmap_1.x[cait_index, :])
    return wot.Dataset(result_x, tmap_0.row_meta.copy(), tmap_1.col_meta.copy())
