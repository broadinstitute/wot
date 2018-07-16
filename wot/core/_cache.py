# -*- coding: utf-8 -*-

import os
import glob
import wot.io

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
    cached_tmaps = {}
    pattern = "*" if core.tmap_prefix is None else core.tmap_prefix
    pattern += '_[0-9]*.[0-9]*_[0-9]*.[0-9]*.*'
    files = glob.glob(os.path.join(core.tmap_dir, pattern))
    for path in files:
        if not os.path.isfile(path):
            continue
        basename, ext = wot.io.get_filename_and_extension(path)
        tokens = basename.split('_')
        try :
            t1 = float(tokens[-2])
            t2 = float(tokens[-1])
        except ValueError:
            continue
        cached_tmaps[(t1, t2)] = path
    return cached_tmaps

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
    if core.tmaps.get((t1, t2), None) is None:
        core.compute_transport_map(t1, t2)

    path = core.tmaps.get((t1, t2))
    return wot.io.read_dataset(path)
