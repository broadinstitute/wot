# -*- coding: utf-8 -*-
import glob
import os

import anndata
import numpy as np
import pandas as pd
import pegasus as pg
import scipy.sparse
import wot


def group_cell_sets(cell_set_paths, group_by_df, group_by_key='day'):
    """
    Return the cell sets by time points given a cell sets file.

    Parameters
    ----------
    cell_set_paths : str or list of str
        The path(s) to the cell sets file. If several are specified, they are
        merged into one list of cell sets after being parsed.
    group_by_df : pandas.DataFrame
        The dataframe containing the considered cell ids as index.
        These may be a subset of cell ids mentionned in the cell sets
        file, in which case only cells in this dataframe appear in the result.
    group_by_key : str (default: 'day')
        The name of the column indicating time information in the dataframe.

    Returns
    -------
    cs_groups : dict of float: list of { 'set': set  of str, 'name': str }
        The different cell sets for each time point.

    Notes
    -----
    cell_set['name'] is a str, the name and time of that cell set.

    For instance, 'cs1' at time 3 would have name 'cs1_3.0'

    Example
    -------
    >>> cs_groups[1.0]
    [ { 'set': { 'cell_1', 'cell_2' }, 'name': 'cell_set_name_1.0' } ]
    """
    group_to_cell_sets = {}
    if isinstance(cell_set_paths, str):
        cell_set_paths = [cell_set_paths]
    for path in cell_set_paths:
        cell_set_ds = wot.io.read_sets(path)
        for i in range(cell_set_ds.shape[1]):
            cell_set_name = cell_set_ds.var.index.values[i]
            cell_ids_in_set = cell_set_ds.obs.index.values[cell_set_ds[:, i].X > 0]

            grouped = group_by_df[group_by_df.index.isin(cell_ids_in_set)].groupby(group_by_key)
            for name, group in grouped:
                cell_sets = group_to_cell_sets.get(name)
                if cell_sets is None:
                    cell_sets = []
                    group_to_cell_sets[name] = cell_sets
                full_name = cell_set_name + '_' + str(name)
                cell_sets.append({'set': set(group.index.values), 'name': full_name})

    return group_to_cell_sets


def list_transport_maps(input_dir):
    transport_maps_inputs = []  # file, start, end
    is_pattern = not os.path.isdir(input_dir)
    files = os.listdir(input_dir) if not is_pattern else glob.glob(input_dir)
    for path in files:
        path = os.path.join(input_dir, path) if not is_pattern else path
        if os.path.isfile(path):
            file_info = wot.io.get_filename_and_extension(os.path.basename(path))
            basename = file_info[0]
            tokens = basename.split('_')
            t1 = tokens[len(tokens) - 2]
            t2 = tokens[len(tokens) - 1]

            try:
                t1 = float(t1)
                t2 = float(t2)

            except ValueError:
                continue

            transport_maps_inputs.append(
                {'path': path, 't1': t1, 't2': t2})

    transport_maps_inputs.sort(key=lambda x: x['t1'])  # sort by t1 (start time)
    return transport_maps_inputs


def read_transport_maps(input_dir, ids=None, time=None):
    """
    Find and parse all transport maps in a directory.
    Returns a list containing the transport maps and start/end timepoints.

    Parameters
    ----------
    input_dir : str
        The directory in which to look for transport maps.
        Alternatively, a pattern may be given, resulting in shell expansion
        before each directory is processed.
    ids : list of str, optional
        Ids to keep the transport maps for.
        If not None, any id not in this list will be filtered out of the maps.
        The order of ids in the resulting transport maps is also guaranteed
        to be the same as this parameter.
    time : int or float, optional
        If ids is not None, specifies the time at which the ids were measured.

    Returns
    -------
    transport_maps : list of { 't1': float, 't2': float, 'transport_map': anndata.AnnData }
        The list of all transport maps

    Raises
    ------
    ValueError
        If exactly one of (ids, time) is None. Must be both or none.
        If no transport map is found in the given directory.
        If several transport maps are found for the same timepoints.

    Notes
    -----
    Time points are determined by the filename.

    Filenames must end in `_{t1}_{t2}.extension`.
    Any transport map not following this convention will be ignored.
    If any other dataset file is present in the listed directories and
    uses this naming convention, it might be interpreted as a transport
    map, yielding unpredictable results.

    All wot commands are guaranteed to enforce this naming convention.
    """
    transport_maps_inputs = []  # file, start, end
    is_pattern = not os.path.isdir(input_dir)

    files = os.listdir(input_dir) if not is_pattern else glob.glob(input_dir)

    if (ids is None) != (time is None):
        raise ValueError("Only one of time and ids is None. Must be both or none")

    tmap_times = set()
    for path in files:
        path = os.path.join(os.path.dirname(input_dir), path) if not is_pattern else path
        if os.path.isfile(path):
            file_info = wot.io.get_filename_and_extension(os.path.basename(path))
            basename = file_info[0]
            tokens = basename.split('_')
            t1 = tokens[len(tokens) - 2]
            t2 = tokens[len(tokens) - 1]

            try:
                t1 = float(t1)
                t2 = float(t2)

            except ValueError:
                continue
            ds = wot.io.read_dataset(path)
            if ids is not None and t1 == time:
                # subset rows
                indices = ds.obs.index.isin(ids)
                ds = ds[indices]
            if ids is not None and t2 == time:
                # subset columns
                indices = ds.var.index.isin(ids)
                ds = ds[:, indices]

            if (t1, t2) in tmap_times:
                raise ValueError("Multiple transport maps found for times ({},{})".format(t1, t2))
            else:
                tmap_times.add((t1, t2))
            transport_maps_inputs.append(
                {'transport_map': ds, 't1': t1, 't2': t2})

    if not transport_maps_inputs:
        raise ValueError("No transport maps found in the given directories")

    transport_maps_inputs.sort(key=lambda x: x['t1'])  # sort by t1 (start time)
    return transport_maps_inputs


def read_sets(path, feature_ids=None, as_dict=False):
    path = str(path)
    hash_index = path.rfind('#')
    set_names = None
    if hash_index != -1:
        set_names = path[hash_index + 1:].split(',')
        path = path[0:hash_index]
    ext = get_filename_and_extension(path)[1]
    if ext == 'gmt':
        gs = read_gmt(path, feature_ids)
    elif ext == 'gmx':
        gs = read_gmx(path, feature_ids)
    elif ext == 'txt' or ext == 'grp':
        gs = read_grp(path, feature_ids)
    else:
        raise ValueError('Unknown file format "{}"'.format(ext))
    if set_names is not None:
        gs_filter = gs.var.index.isin(set_names)
        gs = gs[:, gs_filter]
    if as_dict:
        return wot.io.convert_binary_dataset_to_dict(gs)
    return gs


def read_grp(path, feature_ids=None):
    with open(path) as fp:
        row_id_lc_to_index = {}
        row_id_lc_to_row_id = {}
        if feature_ids is not None:
            for i in range(len(feature_ids)):
                fid = feature_ids[i].lower()
                row_id_lc_to_index[fid] = i
                row_id_lc_to_row_id[fid] = feature_ids[i]

        ids_in_set = set()
        for line in fp:
            if line == '' or line[0] == '#' or line[0] == '>':
                continue
            value = line.strip()
            if value != '':
                value_lc = value.lower()
                row_index = row_id_lc_to_index.get(value_lc)
                if feature_ids is None:
                    if row_index is None:
                        row_id_lc_to_row_id[value_lc] = value
                        row_index = len(row_id_lc_to_index)
                        row_id_lc_to_index[value_lc] = row_index

                if row_index is not None:
                    ids_in_set.add(value)

        if feature_ids is None:
            feature_ids = np.empty(len(row_id_lc_to_index), dtype='object')
            for rid_lc in row_id_lc_to_index:
                feature_ids[row_id_lc_to_index[rid_lc]] = row_id_lc_to_row_id[rid_lc]

        x = np.zeros(shape=(len(feature_ids), 1), dtype=np.int8)
        for id in ids_in_set:
            row_index = row_id_lc_to_index.get(id.lower())
            x[row_index, 0] = 1

        obs = pd.DataFrame(index=feature_ids)
        var = pd.DataFrame(index=[wot.io.get_filename_and_extension(os.path.basename(path))[0]])
        return anndata.AnnData(X=x, obs=obs, var=var)


def read_gmt(path, feature_ids=None):
    with open(path) as fp:
        row_id_lc_to_index = {}
        row_id_lc_to_row_id = {}
        if feature_ids is not None:
            for i in range(len(feature_ids)):
                fid = feature_ids[i].lower()
                row_id_lc_to_index[fid] = i
                row_id_lc_to_row_id[fid] = feature_ids[i]

        members_array = []
        set_descriptions = []
        set_names = []
        for line in fp:
            if line == '' or line[0] == '#':
                continue
            tokens = line.split('\t')
            if len(tokens) < 3:
                continue
            set_names.append(tokens[0].strip())
            description = tokens[1].strip()
            if 'BLANK' == description:
                description = ''
            set_descriptions.append(description)
            ids = tokens[2:]
            ids_in_set = []
            members_array.append(ids_in_set)
            for i in range(len(ids)):
                value = ids[i].strip()
                if value != '':
                    value_lc = value.lower()
                    row_index = row_id_lc_to_index.get(value_lc)
                    if feature_ids is None:
                        if row_index is None:
                            row_id_lc_to_row_id[value_lc] = value
                            row_index = len(row_id_lc_to_index)
                            row_id_lc_to_index[value_lc] = row_index

                    if row_index is not None:
                        ids_in_set.append(value)

        if feature_ids is None:
            feature_ids = np.empty(len(row_id_lc_to_index), dtype='object')
            for rid_lc in row_id_lc_to_index:
                feature_ids[row_id_lc_to_index[rid_lc]] = row_id_lc_to_row_id[rid_lc]

        x = np.zeros(shape=(len(feature_ids), len(set_names)), dtype=np.int8)
        for j in range(len(members_array)):
            ids = members_array[j]
            for id in ids:
                row_index = row_id_lc_to_index.get(id.lower())
                x[row_index, j] = 1

        obs = pd.DataFrame(index=feature_ids)
        var = pd.DataFrame(data={'description': set_descriptions}, index=set_names)
        return anndata.AnnData(X=x, obs=obs, var=var)


def read_gmx(path, feature_ids=None):
    with open(path) as fp:
        set_ids = fp.readline().split('\t')
        descriptions = fp.readline().split('\t')
        nsets = len(set_ids)
        for i in range(len(set_ids)):
            set_ids[i] = set_ids[i].rstrip()

        row_id_lc_to_index = {}
        row_id_lc_to_row_id = {}
        x = None
        array_of_arrays = None
        if feature_ids is not None:
            for i in range(len(feature_ids)):
                fid = feature_ids[i].lower()
                row_id_lc_to_index[fid] = i
                row_id_lc_to_row_id[fid] = feature_ids[i]
            x = np.zeros(shape=(len(feature_ids), nsets), dtype=np.int8)
        else:
            array_of_arrays = []
        for line in fp:
            tokens = line.split('\t')
            for j in range(nsets):
                value = tokens[j].strip()
                if value != '':
                    value_lc = value.lower()
                    row_index = row_id_lc_to_index.get(value_lc)
                    if feature_ids is None:
                        if row_index is None:
                            row_id_lc_to_row_id[value_lc] = value
                            row_index = len(row_id_lc_to_index)
                            row_id_lc_to_index[value_lc] = row_index
                            array_of_arrays.append(np.zeros(shape=(nsets,), dtype=np.int8))
                        array_of_arrays[row_index][j] = 1
                    elif row_index is not None:
                        x[row_index, j] = 1
        if feature_ids is None:
            feature_ids = np.empty(len(row_id_lc_to_index), dtype='object')
            for rid_lc in row_id_lc_to_index:
                feature_ids[row_id_lc_to_index[rid_lc]] = row_id_lc_to_row_id[rid_lc]

        if array_of_arrays is not None:
            x = np.array(array_of_arrays)
        obs = pd.DataFrame(index=feature_ids)
        var = pd.DataFrame(data={'description': descriptions},
            index=set_ids)
        return anndata.AnnData(x, obs=obs, var=var)


def write_sets(sets, path):
    """
    Save gene sets as a gmt file
    :param gene_sets: (str, list) A dict that maps set name to set ids
    :param path: str
        Output file path
    """
    path = str(path)
    path = check_file_extension(path, 'gmt')
    with open(path, 'w') as f:
        for gset in sets:
            f.write('{}\t{}\t{}\n'.format(gset, '-', '\t'.join(sets[gset])))


def convert_binary_dataset_to_dict(ds):
    cell_sets = {}
    for i in range(ds.shape[1]):
        selected = np.where(ds[:, i].X == 1)[0]
        cell_sets[ds.var.index[i]] = list(ds.obs.index[selected])
    return cell_sets


def read_dataset(path, obs=None, var=None, obs_filter=None, var_filter=None, **keywords):
    """
    Read h5ad, loom, mtx, 10X h5, and csv formatted files

    Parameters
    ----------
    path: str
        File name of data file.
    obs: {str, pd.DataFrame}
        Path to obs data file or a data frame
    var: {str, pd.DataFrame}
        Path to var data file or a data frame
    obs_filter {str, pd.DataFrame}
        File with one id per line, name of a boolean field in obs, or a list of ids
    var_filter: {str, pd.DataFrame}
        File with one id per line, name of a boolean field in obs, or a list of ids
    Returns
    -------
    Annotated data matrix.
    """
    if str(path).lower().endswith('.txt'):
        df = pd.read_csv(path, engine='python', header=0, sep=None, index_col=0)
        adata = anndata.AnnData(X=df.values, obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns))
    else:
        adata = pg.read_input(path, **keywords)

    def get_df(meta):
        if not isinstance(meta, pd.DataFrame):
            tmp_path = None
            if meta.startswith('gs://'):
                tmp_path = download_gs_url(meta)
                meta = tmp_path
            meta = pd.read_csv(meta, sep=None, index_col='id', engine='python')
            if tmp_path is not None:
                os.remove(tmp_path)
        return meta

    if obs is not None:
        if not isinstance(obs, list) and not isinstance(obs, tuple):
            obs = [obs]
        for item in obs:
            adata.obs = adata.obs.join(get_df(item))
    if var is not None:
        if not isinstance(var, list) and not isinstance(var, tuple):
            var = [var]
        for item in var:
            adata.var = adata.var.join(get_df(item))

    return filter_adata(adata, obs_filter=obs_filter, var_filter=var_filter)


def write_dataset(ds, path, output_format='txt'):
    path = str(path)
    if not path.lower().endswith('.' + output_format):
        path += '.' + output_format
    if output_format == 'txt':
        x = ds.X.toarray() if scipy.sparse.isspmatrix(ds.X) else ds.X
        pd.DataFrame(x, index=ds.obs.index, columns=ds.var.index).to_csv(path,
            index_label='id', sep='\t', doublequote=False)
    elif output_format == 'h5ad':
        ds.write(path)
    elif output_format == 'loom':
        ds.write_loom(ds, path)
    else:
        raise ValueError('Unknown file format')


def download_gs_url(gs_url):
    from google.cloud import storage
    client = storage.Client()
    path = gs_url[len('gs://'):]
    slash = path.find('/')
    bucket_id = path[0:slash]
    blob_path = path[slash + 1:]
    bucket = client.get_bucket(bucket_id)
    blob = bucket.blob(blob_path)
    dot = path.rfind('.')
    suffix = None
    if dot != -1:
        suffix = path[dot:]
    import tempfile
    tmp = tempfile.mkstemp(suffix=suffix)
    path = tmp[1]
    blob.download_to_filename(path)
    return path


def check_file_extension(name, output_format):
    if not str(name).lower().endswith('.' + output_format):
        name += '.' + output_format
    return name


def get_filename_and_extension(name):
    name = os.path.basename(name)
    dot_index = name.rfind('.')
    ext = ''
    basename = name
    if dot_index != -1:
        ext = name[dot_index + 1:].lower()
        basename = name[0:dot_index]
        if ext == 'txt':  # check for .gmt.txt e.g.
            dot_index2 = basename.rfind('.')
            if dot_index2 != -1:
                ext2 = basename[dot_index2 + 1:].lower()
                if ext2 in set(['gmt', 'grp', 'gmx']):
                    basename = basename[0:dot_index2]
                    return basename, ext2
    return basename, ext


def filter_adata(adata, obs_filter=None, var_filter=None):
    if obs_filter is not None:
        if os.path.exists(obs_filter):
            adata = adata[adata.obs.index.isin(wot.io.read_sets(obs_filter).obs.index)].copy()
        else:
            obs_filter = obs_filter.split(',')
            if len(obs_filter) == 1 and obs_filter[0] in adata.obs:  # boolean field in obs
                adata = adata[adata.obs[obs_filter] == True].copy()
            else:  # list of ids
                adata = adata[adata.obs.index.isin(obs_filter)].copy()
    if var_filter is not None:
        if os.path.exists(var_filter):
            adata = adata[:, adata.var.index.isin(wot.io.read_sets(var_filter).obs.index)].copy()
        else:
            var_filter = var_filter.split(',')
            if len(var_filter) == 1 and var_filter[0] in adata.var:  # boolean field in var
                adata = adata[:, adata.var[var_filter[0]]].copy()
            else:  # list of ids
                adata = adata[:, adata.var.index.isin(var_filter)].copy()
    return adata


def read_days_data_frame(path):
    return pd.read_csv(path, index_col='id',
        engine='python', sep=None, dtype={'day': np.float64})


def add_row_metadata_to_dataset(dataset, days=None, growth_rates=None, covariate=None):
    if days is not None:
        if not os.path.exists(days):
            raise ValueError(days + ' not found')
        dataset.obs = dataset.obs.join(read_days_data_frame(days))
    if growth_rates is not None:
        if not os.path.exists(growth_rates):
            raise ValueError(growth_rates + ' not found')
        dataset.obs = dataset.obs.join(pd.read_csv(growth_rates, index_col='id', engine='python', sep=None))
        # if 'cell_growth_rate' not in dataset.obs:
        #     raise ValueError('Cell growth rates must that the column headers id and cell_growth_rate')
    else:
        dataset.obs['cell_growth_rate'] = 1.0
    # if sampling_bias_path is not None:
    #     dataset.obs = dataset.obs.join(
    #         pd.read_csv(sampling_bias_path, index_col='id', engine='python', sep=None))
    if covariate is not None:
        if not os.path.exists(covariate):
            raise ValueError(covariate + ' not found')
        dataset.obs = dataset.obs.join(pd.read_csv(covariate, index_col='id', engine='python', sep=None))


def read_day_pairs(day_pairs):
    if os.path.isfile(day_pairs):
        target = day_pairs
        args = {'engine': 'python', 'sep': None}
    else:
        import io
        target = io.StringIO(day_pairs)
        args = {'sep': ',', 'lineterminator': ';'}
    return pd.read_csv(target, **args)
