# -*- coding: utf-8 -*-

import os
import glob
import wot.io
import hashlib
import yaml

def serialize_day_pairs(day_pairs):
    """Return a YAML-friendly serialized object reprensenting the day_pairs dict"""
    if day_pairs is None:
        return None
    serialize = lambda x: '->'.join(str(z) for z in x)
    serializable = { serialize(x): day_pairs[x] for x in day_pairs }
    return yaml.dump(serializable).rstrip().split('\n')

def deserialize_day_pairs(serialized):
    """Rebuild the day_pairs dict from the YAML-friendly serialized object"""
    if serialized is None or serialized == 'null':
        return None
    deserialize = lambda x: tuple(float(z) for z in x.split('->'))
    deserializable = yaml.load('\n'.join(serialized))
    return { deserialize(x): deserializable[x] for x in deserializable }

def get_full_ot_configuration(ot_model, serialized=True):
    """
    Gets the OT configuration of a given model

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel whose configuration to extract

    Returns
    -------
    config : dict of str
        The configuration for each parameter
    """
    matrix_hash = hashlib.sha256(ot_model.matrix.x.data.tobytes()).hexdigest()
    if serialized:
        day_pairs = serialize_day_pairs(ot_model.day_pairs)
    else:
        day_pairs = ot_model.day_pairs

    return {
            'day_pairs': day_pairs,
            'matrix_hash': matrix_hash,
            ** ot_model.get_ot_config(),
           }

def parse_ot_configuration_from_stream(stream):
    """Parse the OT configuration from the YAML stream"""
    config = yaml.load(stream)
    config['day_pairs'] = deserialize_day_pairs(config.get('day_pairs', 'null'))
    return config

def are_ot_configurations_compatible(ot_model, cache_config):
    """
    Check that the current configuration of an OTModel is compatible with a cache configuration

    A cache configuration is compatible if all scalar variables match,
    and the cached day_pairs is a subset of the OTModel day_pairs.
    """
    current_config = get_full_ot_configuration(ot_model, serialized=False)

    wot.io.verbose("Current configuration :", current_config)
    wot.io.verbose("Cache configuration :", cache_config)

    for key in current_config:
        if key == 'day_pairs':
            # If no collision is present, keep the caches
            if current_config[key] is None and cache_config[key] is None:
                # Both day_pairs are None. No collisions
                continue
            elif current_config[key] is None or cache_config[key] is None:
                wot.io.verbose("One of the two day_pairs is None. Not compatible")
                return False
            # Now both are not None, check if they are compatible
            for x in cache_config[key]:
                if x not in current_config[key]:
                    wot.io.verbose("Tmap", x, "is no longer present")
                    return False
                if current_config[key][x] != cache_config[key][x]:
                    wot.io.verbose("Tmap configuration for", x, "has changed")
                    return False
        elif key in cache_config and current_config[key] != cache_config[key]:
            wot.io.verbose("Configuration value for", key, "not compatible")
            return False
    return True


def write_config_file(ot_model):
    """
    Writes the current configuration to the tmap prefix file

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel whose configurations needs to be written
    """
    config = get_full_ot_configuration(ot_model)

    config_file = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix) + '.yml'
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

def list_cached_transport_maps(tmap_dir, tmap_prefix, with_covariates=False):
    """Get the list of paths to valid transport map names (complex regexp not available in Python)"""
    valid_tmaps = []
    pattern = tmap_prefix
    if with_covariates:
        pattern += '_[0-9]*.[0-9]*_[0-9]*.[0-9]*_cv[0-9]*_cv[0-9]*.*'
    else:
        pattern += '_[0-9]*.[0-9]*_[0-9]*.[0-9]*.*'
    files = glob.glob(os.path.join(tmap_dir, pattern))
    for path in files:
        if not os.path.isfile(path):
            continue
        basename, ext = wot.io.get_filename_and_extension(path)
        tokens = basename.split('_')
        try :
            if with_covariates:
                t1 = float(tokens[-4])
                t2 = float(tokens[-3])
                cv1 = int(tokens[-2][2:])
                cv2 = int(tokens[-1][2:])
            else:
                t1 = float(tokens[-2])
                t2 = float(tokens[-1])
            valid_tmaps.append(path)
        except ValueError:
            continue
    return valid_tmaps

def purge_invalidated_caches(ot_model):
    """
    Removes all cached transport maps that do not match the current setting

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel whose cache needs to be checked
    """
    config_file = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix) + '.yml'
    try:
        with open(config_file, "r") as stream:
            cached_config = parse_ot_configuration_from_stream(stream)
            purge = not are_ot_configurations_compatible(ot_model, cached_config)
            if purge:
                wot.io.verbose("Configuration incompatible with cache. Triggering purge")
            else:
                wot.io.verbose("Cache configuration OK")
    # FIXME: This except is catching all errors, should be restricted to YAML/IO
    except:
        # File not present or not a valid YAML file, purge prefix
        wot.io.verbose("Cache description not present or not a valid YAML file :", config_file)
        purge = True

    if purge:
        wot.io.verbose("Purging ({},{})".format(ot_model.tmap_dir, ot_model.tmap_prefix))
        if os.path.isfile(config_file):
            os.remove(config_file)
        for cov in [False, True]:
            for path in list_cached_transport_maps(ot_model.tmap_dir, ot_model.tmap_prefix, with_covariates=cov):
                wot.io.verbose("Removing tmap", path)
                os.remove(path)
    write_config_file(ot_model)

def scan_transport_map_directory(ot_model):
    """
    Scan the transport map directory for cached transport maps.

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel that is scanning for transport maps

    Returns
    -------
    cached_tmaps : dict of (float, float): str
        The path to each transport map that was found in the directory.
    cached_cov_tmaps : dict of (float, float, float, float): str
        The path to each covariate-restricted transport map that was found in the directory.
    """
    wot.io.verbose("Scanning transport map directory")
    cached_tmaps = {}
    for path in list_cached_transport_maps(ot_model.tmap_dir, ot_model.tmap_prefix):
        basename, ext = wot.io.get_filename_and_extension(path)
        tokens = basename.split('_')
        t1 = float(tokens[-2])
        t2 = float(tokens[-1])
        cached_tmaps[(t1, t2)] = path
    wot.io.verbose("Cached tmaps found:", cached_tmaps)

    cov_cached_tmaps = {}
    for path in list_cached_transport_maps(ot_model.tmap_dir, ot_model.tmap_prefix, with_covariates=True):
        basename, ext = wot.io.get_filename_and_extension(path)
        tokens = basename.split('_')
        t1 = float(tokens[-4])
        t2 = float(tokens[-3])
        cv1 = int(tokens[-2][2:])
        cv2 = int(tokens[-1][2:])
        cov_cached_tmaps[(t1, t2, cv1, cv2)] = path
    wot.io.verbose("Cached cov-restricted tmaps found:", cov_cached_tmaps)
    return cached_tmaps, cov_cached_tmaps

def load_transport_map(ot_model, t1, t2):
    """
    Load the given transport map, either from cache or compute it.

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel that is trying to load a transport map.
    t1 : int or float
        The source timepoint of the transport map.
    t2 : int or float
        The destination timepoint of the transport map.

    Returns
    -------
    tmap : wot.Dataset
        The given transport map
    """
    if ot_model.tmaps.get((t1, t2), None) is None:
        ot_model.compute_transport_map(t1, t2)

    path = ot_model.tmaps.get((t1, t2))

    try:
        tmap = wot.io.read_dataset(path)
    except:
        wot.io.verbose("Warning : Error when reading tmap ({}, {}),"\
                " file might have been corrupted. Recomputing".format(t1, t2))
        ot_model.compute_transport_map(t1, t2)
        path = ot_model.tmaps.get((t1, t2))
        tmap = wot.io.read_dataset(path)

    return tmap

def load_covariate_restricted_transport_map(ot_model, t0, t1, covariate):
    """
    Load the given covariate-restricted transport map, either from cache or compute it

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel that is trying to load a transport map
    t0 : int or float
        The source timepoint of the transport map
    t1 : int or float
        The destination timepoint of the transport map
    covariate : (int, int)
        The covariate restriction on cells from t0 and t1

    Returns
    -------
    tmap : wot.Dataset
        The transport map
    """
    cv0, cv1 = covariate
    if ot_model.cov_tmaps.get((t0, t1, cv0, cv1), None) is None:
        ot_model.compute_transport_map(t0, t1, covariate=(cv0, cv1))

    path = ot_model.cov_tmaps.get((t0, t1, cv0, cv1))

    try:
        tmap = wot.io.read_dataset(path)
    except:
        wot.io.verbose("Warning : Error when reading tmap ({}, {}) cv ({}, {}),"\
                " file might have been corrupted. Recomputing".format(t0, t1, cv0, cv1))
        ot_model.compute_transport_map(t0, t1, covariate=(cv0, cv1))
        path = ot_model.cov_tmaps.get((t0, t1, cv0, cv1))
        tmap = wot.io.read_dataset(path)

    return tmap
