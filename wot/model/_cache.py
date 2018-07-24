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
    if serialized == 'null':
        return None
    deserialize = lambda x: tuple(float(z) for z in x.split('->'))
    deserializable = yaml.load('\n'.join(serialized))
    return { deserialize(x): deserializable[x] for x in deserializable }

def get_ot_configuration(ot_model):
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

    return {
            'matrix_hash': matrix_hash,
            'day_pairs': serialize_day_pairs(ot_model.day_pairs),
            ** ot_model.get_ot_config(),
           }

def parse_ot_configuration_from_stream(stream):
    """Parse the OT configuration from the YAML stream"""
    config = yaml.load(stream)
    config['day_pairs'] = deserialize_day_pairs(config.get('day_pairs', 'null'))
    return config

def write_config_file(ot_model):
    """
    Writes the current configuration to the tmap prefix file

    Parameters
    ----------
    ot_model : wot.OTModel
        The OTModel whose configurations needs to be written
    """
    config = get_ot_configuration(ot_model)

    config_file = os.path.join(ot_model.tmap_dir, ot_model.tmap_prefix) + '.yml'
    with open(config_file, "w") as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

def list_cached_transport_maps(tmap_dir, tmap_prefix):
    """Get the list of paths to valid transport map names (complex regexp not available in Python)"""
    valid_tmaps = []
    pattern = tmap_prefix
    pattern += '_[0-9]*.[0-9]*_[0-9]*.[0-9]*.*'
    files = glob.glob(os.path.join(tmap_dir, pattern))
    for path in files:
        if not os.path.isfile(path):
            continue
        basename, ext = wot.io.get_filename_and_extension(path)
        tokens = basename.split('_')
        try :
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
            # TODO: choose wether or not to purge here
            purge = True
    except:
        # File not present or not a valid YAML file, purge prefix
        purge = True

    if purge:
        if os.path.isfile(config_file):
            os.remove(config_file)
        for path in list_cached_transport_maps(ot_model.tmap_dir, ot_model.tmap_prefix):
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
    """
    purge_invalidated_caches(ot_model)

    cached_tmaps = {}
    for path in list_cached_transport_maps(ot_model.tmap_dir, ot_model.tmap_prefix):
        basename, ext = wot.io.get_filename_and_extension(path)
        tokens = basename.split('_')
        t1 = float(tokens[-2])
        t2 = float(tokens[-1])
        cached_tmaps[(t1, t2)] = path
    return cached_tmaps

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
    return wot.io.read_dataset(path)
