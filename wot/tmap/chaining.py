import wot.tmap


def chain_transport_maps(tmap_model, pairs_list):
    """
    Chains the transport maps corresponding to the list of pairs for the OTModel.

    Parameters
    ----------
    tmap_model : wot.tmap.TransportMapModel
        The TransportMapModel whose transport maps are to be chained.
    pairs_list : list of (float, float)
        The list of day pairs correspondig to the transport maps to chain.

    Returns
    -------
    tmap : anndata.AnnData
        The final transport map

    Raises
    ------
    ValueError
        If any consecutive pairs (a, b), (c, d) in pairs_list do not verify b = c.
    ValueError
        If any pair in pairs_list is not a valid day pair for the OTModel.
    ValueError
        If any pair (a, b) in pairs_list has a >= b.
    """
    for i in range(len(pairs_list) - 1):
        if pairs_list[i][1] != pairs_list[i + 1][0]:
            raise ValueError("Transport maps cannot be glued together : invalid path")
    for a, b in pairs_list:
        if a >= b:
            raise ValueError("({}, {}) is not a valid transport map : it goes backwards in time".format(a, b))

    tmap_0 = tmap_model.get_coupling(*pairs_list.pop(0))
    while len(pairs_list) > 0:
        tmap_1 = tmap_model.get_coupling(*pairs_list.pop(0))
        tmap_0 = wot.tmap.glue_transport_maps(tmap_0, tmap_1)
    return tmap_0


def find_path(t0, t1, available_pairs, timepoints):
    """
    Finds a path from t0 to t1 using the available day_pairs. Uses the finest resolution possible.

    Parameters
    ----------
    t0 : float
        Starting point for the transport map.
    t1 : float
        Destination point for the transport map.
    available_pairs : list of (float, float) or None
        All day_pairs elligible as edges for the path.
        If None, all consecutives timepoints are allowed.
    timepoints : list of float
        All timepoints of the considered OTModel.

    Returns
    -------
    path : list of (float, float)
        A path using only the available day_pairs to go from t0 to t1.

    Raises
    ------
    ValueError
        If no path exists from t0 to t1 using available day_pairs.
    ValueError
        If t0 or t1 are not in the list of timepoints.

    Notes
    -----
    You cannot use day_pairs backwards :
    [ (0,2), (1, 2), (1, 3) ] does not allow computation from 0 to 3.
    0 -> 2 -> 1 -> 3 is not an acceptable path, because 2 -> 1 goes backwards in time.
    """
    if t0 not in timepoints or t1 not in timepoints:
        raise ValueError("Unable to build path : t0 and t1 must be present in timepoints list")

    if available_pairs is None:
        t0_i = timepoints.index(t0)
        t1_i = timepoints.index(t1)
        return [(timepoints[i], timepoints[i + 1]) for i in range(t0_i, t1_i)]

    # forall a, b in available_pairs, a < b.
    # the increasing order on timepoints is thus a topological ordering of the graph of days.
    reach = {t: [] for t in timepoints}
    for a, b in available_pairs:
        reach[a].append(b)
    reach = {t: sorted(reach[t]) for t in timepoints}

    # FIXME: the longest path is not guaranteed to be the "finest resolution" in all cases
    dist_prev = {t: (0, None) for t in timepoints}
    for t in timepoints[timepoints.index(t0):]:
        d, _ = dist_prev[t]
        for u in reach[t]:
            if d >= dist_prev[u][0]:
                dist_prev[u] = (d + 1, t)

    if dist_prev[t1][0] == 0:
        raise ValueError("{} is not reachable from {} with the chosen day_pairs".format(t1, t0))

    path = []
    cur = t1
    while cur != t0:
        path.insert(0, (dist_prev[cur][1], cur))
        cur = dist_prev[cur][1]
    return path
