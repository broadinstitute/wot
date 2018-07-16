# -*- coding: utf-8 -*-

def scan_transport_map_directory(core):
    """
    Scan the transport map directory for cached transport maps.

    Parameters
    ----------
    core : wot.Core
        The Core that is scanning for transport maps

    Returns
    -------
    cached_tmaps : dict of (float, float): str
        The path to each transport map that was found in the directory.
    """
    raise ValueError("Not implemented")


def load_transport_map(core, t1, t2):
    """
    Load the given transport map, either from cache or compute it.

    Parameters
    ----------
    core : wot.Core
        The Core that is trying to load a transport map.
    t1 : int or float
        The source timepoint of the transport map.
    t2 : int or float
        The destination timepoint of the transport map.

    Returns
    -------
    tmap : wot.Dataset
        The given transport map
    """
    raise ValueError("Not implemented")
