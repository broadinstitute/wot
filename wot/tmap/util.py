import anndata
import numpy as np


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
