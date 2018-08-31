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


