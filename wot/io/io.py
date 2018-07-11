# -*- coding: utf-8 -*-
import glob
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.io
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
    cs_groups : dict of float: list of cell_set
        The different cell sets for each time point.

    Notes
    -----
    The cell_set type is a dictionary with two keys: 'set' and 'name'
    cell_set['set'] : set of str
        The ids of the cells in that cell set.
    cell_set['name'] : str
        The name of that cell set, with the time appended to it.
        For instance, 'cs1' at time 3 would have name 'cs1_3.0'

    Example
    -------
    cs_groups[1.0] : [ { 'set': { 'cell_1', 'cell_2' },
                         'name': 'cell_set_name_1.0' } ]
    """
    group_to_cell_sets = {}
    if isinstance(cell_set_paths, str):
        cell_set_paths = [cell_set_paths]
    for path in cell_set_paths:
        cell_set_ds = wot.io.read_gene_sets(path)
        for i in range(cell_set_ds.x.shape[1]):
            cell_set_name = cell_set_ds.col_meta.index.values[i]
            cell_ids_in_set = cell_set_ds.row_meta.index.values[cell_set_ds.x[:, i] > 0]

            grouped = group_by_df[group_by_df.index.isin(cell_ids_in_set)].groupby(group_by_key)
            for name, group in grouped:
                cell_sets = group_to_cell_sets.get(name)
                if cell_sets is None:
                    cell_sets = []
                    group_to_cell_sets[name] = cell_sets
                full_name = cell_set_name + '_' + str(name)
                cell_sets.append({'set': set(group.index.values), 'name': full_name})

    return group_to_cell_sets


def filter_ds_from_command_line(ds, args):
    params = vars(args)
    if params.get('gene_filter') is not None:
        prior = ds.x.shape[1]
        gene_ids = pd.read_table(args.gene_filter, index_col=0, header=None).index.values
        column_indices = ds.col_meta.index.isin(gene_ids)
        nkeep = np.sum(column_indices)
        if params.get('verbose') and len(gene_ids) > nkeep:
            print(str(len(gene_ids) - nkeep) + ' are in gene filter, but not in matrix')

        ds = wot.Dataset(ds.x[:, column_indices], ds.row_meta, ds.col_meta.iloc[column_indices])
        if params.get('verbose'):
            print('Keeping ' + str(ds.x.shape[1]) + '/' + str(prior) + ' genes')

    if params.get('cell_filter') is not None:
        prior = ds.x.shape[0]
        if not os.path.isfile(args.cell_filter):
            import re
            expr = re.compile(args.cell_filter)
            cell_ids = [elem for elem in ds.row_meta.index.values if expr.match(elem)]
        else:
            cell_ids = pd.read_table(args.cell_filter, index_col=0, header=None).index.values

        # row_indices = np.isin(ds.row_meta.index.values, cell_ids, assume_unique=True)
        row_indices = ds.row_meta.index.isin(cell_ids)
        nkeep = np.sum(row_indices)
        if params.get('verbose') and len(cell_ids) > nkeep:
            print(str(len(cell_ids) - nkeep) + ' are in cell filter, but not in matrix')

        ds = wot.Dataset(ds.x[row_indices], ds.row_meta.iloc[row_indices], ds.col_meta)
        if params.get('verbose'):
            print('Keeping ' + str(ds.x.shape[0]) + '/' + str(prior) + ' cells')
    return ds


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
    transport_maps_inputs = []  # file, start, end
    is_pattern = not os.path.isdir(input_dir)

    files = os.listdir(input_dir) if not is_pattern else glob.glob(input_dir)

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
                indices = ds.row_meta.index.isin(ids)
                ds = wot.Dataset(ds.x[indices], ds.row_meta.iloc[indices], ds.col_meta)
            if ids is not None and t2 == time:
                # subset columns
                indices = ds.col_meta.index.isin(ids)
                ds = wot.Dataset(ds.x[:, indices], ds.row_meta, ds.col_meta.iloc[indices])
            transport_maps_inputs.append(
                {'transport_map': ds, 't1': t1, 't2': t2})

    transport_maps_inputs.sort(key=lambda x: x['t1'])  # sort by t1 (start time)
    return transport_maps_inputs


def read_gene_sets(path, feature_ids=None):
    path = str(path)
    basename_and_extension = get_filename_and_extension(path)
    ext = basename_and_extension[1]
    if ext == 'gmt':
        return read_gmt(path, feature_ids)
    elif ext == 'gmx':
        return read_gmx(path, feature_ids)
    elif ext == 'txt' or ext == 'grp':
        return read_grp(path, feature_ids)
    else:
        raise ValueError('Unknown file format')


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
            if line == '' or line[0] == '#':
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

        row_meta = pd.DataFrame(index=feature_ids)
        col_meta = pd.DataFrame(index=[wot.io.get_filename_and_extension(os.path.basename(path))[0]])
        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)


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

        row_meta = pd.DataFrame(index=feature_ids)
        col_meta = pd.DataFrame(data={'description': set_descriptions}, index=set_names)
        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)


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
        row_meta = pd.DataFrame(index=feature_ids)
        col_meta = pd.DataFrame(data={'description': descriptions},
                                index=set_ids)
        return wot.Dataset(x, row_meta=row_meta, col_meta=col_meta)


def write_gene_sets(gene_sets, path, format=None):
    path = str(path)
    basename, ext = get_filename_and_extension(path)

    if path is None or path in ['STDOUT', 'stdout', '/dev/stdout']:
        f = sys.stdout
    else:
        if format is not None and ext != format:
            path = path + '.' + format
        f = open(path, 'w')

    if format == 'gmt':
        write_gmt(gene_sets, f)
    elif format == 'gmx' or format == 'txt' or format == 'grp':
        raise ValueError("Filetype not supported for writing")
    else:
        raise ValueError("Unkown file format for gene sets")

    if f is not sys.stdout:
        f.close()


def write_gmt(gene_sets, f):
    for gset in gene_sets:
        f.write('{}\t{}\t{}\n'.format(gset, '-', '\t'.join(gene_sets[gset])))


def read_dataset(path, chunks=(500, 500), use_dask=False, genome10x=None, row_filter=None, col_filter=None,
                 force_sparse=False, backed=False):
    path = str(path)
    basename_and_extension = get_filename_and_extension(path)
    ext = basename_and_extension[1]

    if ext == 'mtx':
        # look for .barcodes.txt and .genes.txt
        sp = os.path.split(path)
        row_meta = None
        for f in (
                os.path.join(sp[0],
                             basename_and_extension[0] + '.barcodes.tsv'),
                os.path.join(sp[0],
                             basename_and_extension[0] + '.barcodes.txt'),
                os.path.join(sp[0], 'barcodes.tsv')):
            if os.path.isfile(f) or os.path.isfile(f + '.gz'):
                row_meta = pd.read_table(f if os.path.isfile(f) else f + '.gz', index_col=0, sep='\t',
                                         header=None)
                break
        col_meta = None
        for f in (os.path.join(sp[0], basename_and_extension[0] +
                                      '.genes.tsv'),
                  os.path.join(sp[0], basename_and_extension[0] +
                                      '.genes.txt'),
                  os.path.join(sp[0], 'genes.tsv')):
            if os.path.isfile(f) or os.path.isfile(f + '.gz'):
                col_meta = pd.read_table(f if os.path.isfile(f) else f + '.gz', index_col=0, sep='\t',
                                         header=None)
                break

        x = scipy.io.mmread(path)
        x = scipy.sparse.csr_matrix(x.T)
        if col_meta is None:
            print(basename_and_extension[0] + '.genes.txt not found')
            col_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=x.shape[1], step=1))
        if row_meta is None:
            print(basename_and_extension[0] + '.barcodes.txt not found')
            row_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=x.shape[0], step=1))

        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
    elif ext == 'hdf5' or ext == 'h5' or ext == 'loom' or ext == 'h5ad':
        f = h5py.File(path, 'r')
        if ext == 'h5ad':
            h5_x = '/X'
            h5_row_meta = '/obs'
            h5_col_meta = '/var'
        elif ext == 'loom':
            h5_x = '/matrix'
            h5_row_meta = '/row_attrs'
            h5_col_meta = '/col_attrs'
        else:
            if genome10x is None:
                keys = list(f.keys())
                if len(keys) > 0:
                    genome10x = keys[0]
            group = f['/' + genome10x]

            M, N = group['shape'][()]
            data = group['data'][()]
            x = scipy.sparse.csr_matrix((data, group['indices'][()], group['indptr'][()]), shape=(N, M))
            col_meta = pd.DataFrame(index=group['gene_names'][()].astype(str),
                                    data={'ensembl': group['genes'][()].astype(str)})
            row_meta = pd.DataFrame(index=group['barcodes'][()].astype(str))
            f.close()

            return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
        if ext == 'h5ad':
            row_meta = pd.DataFrame.from_records(f[h5_row_meta][()], index='index')
            col_meta = pd.DataFrame.from_records(f[h5_col_meta][()], index='index')
            row_meta.index = row_meta.index.values.astype(str)
            col_meta.index = col_meta.index.values.astype(str)

        else:
            row_attrs = read_h5_attrs(f, h5_row_meta, row_filter)
            nrows = len(row_attrs['indices']) if row_attrs['indices'] is not None else f[h5_x].shape[0]
            row_meta = pd.DataFrame(row_attrs['attrs'], index=pd.RangeIndex(start=0, stop=nrows, step=1))
            if row_meta.get('id') is not None:
                row_meta.set_index('id', inplace=True)

            col_attrs = read_h5_attrs(f, h5_col_meta, col_filter)
            ncols = len(col_attrs['indices']) if col_attrs['indices'] is not None else f[h5_x].shape[1]

            col_meta = pd.DataFrame(col_attrs['attrs'], index=pd.RangeIndex(start=0, stop=ncols, step=1))
            if col_meta.get('id') is not None:
                col_meta.set_index('id', inplace=True)

        if not use_dask:
            x = f[h5_x]
            is_x_sparse = False
            if type(x) == h5py.Group:
                data = x['data'][()]
                x = scipy.sparse.csr_matrix((data, x['indices'][()], x['indptr'][()]), shape=x.attrs['h5sparse_shape'])
                backed = False
            else:
                is_x_sparse = x.attrs.get('sparse')
                if not backed and (is_x_sparse or force_sparse) and (row_filter is None and col_filter is None):
                    # read in blocks of 1000
                    chunk_start = 0
                    chunk_step = min(nrows, 1000)
                    chunk_stop = chunk_step
                    nchunks = int(np.ceil(max(1, nrows / chunk_step)))
                    sparse_arrays = []
                    for chunk in range(nchunks):
                        chunk_stop = min(nrows, chunk_stop)
                        subset = scipy.sparse.csr_matrix(x[chunk_start:chunk_stop])
                        sparse_arrays.append(subset)
                        chunk_start += chunk_step
                        chunk_stop += chunk_step

                    x = scipy.sparse.vstack(sparse_arrays)
                else:
                    if row_filter is None and col_filter is None and not backed:
                        x = x[()]
                    elif row_filter is not None and col_filter is not None:
                        x = x[row_attrs['indices']]
                        x = x[:, col_attrs['indices']]
                    elif row_filter is not None:
                        x = x[row_attrs['indices']]
                    elif col_filter is not None:
                        x = x[:, col_attrs['indices']]

                    if not backed and (is_x_sparse or force_sparse):
                        x = scipy.sparse.csr_matrix(x)
            if not backed:
                f.close()

            return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
        else:

            import dask.array as da
            x = da.from_array(f[h5_x], chunks=chunks)
            # TODO load in chunks
            # row_meta = dd.from_pandas(row_meta, npartitions=4, sort=False)
            # col_meta = dd.from_pandas(col_meta, npartitions=4, sort=False)
        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
    elif ext == 'gct':
        return wot.io.read_gct(path)
    else:
        with open(path) as fp:
            row_ids = []
            header = fp.readline()

            for s in ['\t', ' ', ',']:
                test = header.split(s)
                if len(test) > 1:
                    sep = s
                    column_ids = test
                    break

            column_ids = column_ids[1:]
            column_ids[len(column_ids) - 1] = column_ids[
                len(column_ids) - 1].rstrip()

            i = 0
            np_arrays = []
            for line in fp:
                line = line.rstrip()
                if line != '':
                    tokens = line.split(sep)
                    row_ids.append(tokens[0])
                    np_arrays.append(np.array(tokens[1:], dtype=np.float64))
                    i += 1

            return wot.Dataset(x=np.array(np_arrays),
                               row_meta=pd.DataFrame(index=row_ids),
                               col_meta=pd.DataFrame(index=column_ids))


def read_h5_attrs(f, path, filter):
    g = f[path]
    indices = None
    if filter is not None:
        for key in filter.keys():
            values = g[key][()]
            f = filter[key]
            if values.dtype.kind == 'S':
                values = values.astype(str)
            _indices = []
            for i in range(len(values)):
                if f(values[i]):
                    _indices.append(i)
            if indices is not None:
                indices = np.intersect1d(_indices, indices, assume_unique=True)
            else:
                indices = _indices

        if len(indices) is 0:
            raise ValueError('No indices passed filter')

    data = {}
    for key in g:
        values = g[key]
        if indices is None:
            values = values[()]
        else:
            values = values[indices]
        if values.dtype.kind == 'S':
            values = values.astype(str)
        data[key] = values

    return {'attrs': data, 'indices': indices}


def check_file_extension(name, output_format):
    expected = None
    if output_format == 'csv':
        expected = '.csv'
    elif output_format == 'txt':
        expected = '.txt'
    elif output_format == 'txt.gz':
        expected = '.txt.gz'
    elif output_format == 'loom':
        expected = '.loom'
    elif output_format == 'gct':
        expected = '.gct'
    if expected is not None:
        if not str(name).lower().endswith(expected):
            name += expected
    return name


def get_filename_and_extension(name):
    dot_index = name.rfind('.')
    ext = ''
    basename = name
    if dot_index != -1:
        ext = name[dot_index + 1:]
        if ext == 'gz':
            return get_filename_and_extension(name[0:dot_index])

    if dot_index != -1:
        basename = name[0:dot_index]
    return basename, ext


def write_dataset(ds, path, output_format='txt', txt_full=True):
    path = check_file_extension(path, output_format)
    if output_format == 'txt' or output_format == 'txt.gz' or output_format == 'gct':
        if txt_full or output_format == 'gct':
            f = open(path, 'w')
            # write columns ids

            if output_format == 'gct':
                f.write('#1.3\n')
                f.write(str(ds.x.shape[0]) + '\t' + str(ds.x.shape[1]) + '\t' + str(len(ds.row_meta.columns)) +
                        '\t' + str(len(ds.col_meta.columns)) + '\n')
            f.write('id\t')
            f.write('\t'.join(ds.row_meta.columns))
            if len(ds.row_meta.columns) > 0:
                f.write('\t')
            f.write('\t'.join(ds.col_meta.index.values))
            f.write('\n')
            spacer = ''.join(np.full(len(ds.row_meta.columns), '\t', dtype=object))
            # column metadata fields + values
            for field in ds.col_meta.columns:
                f.write(field)
                f.write(spacer)
                for val in ds.col_meta[field].values:
                    f.write('\t')
                    f.write(str(val))

                f.write('\n')
            # TODO write as sparse array
            pd.DataFrame(index=ds.row_meta.index, data=np.hstack(
                (ds.row_meta.values, ds.x.toarray() if scipy.sparse.isspmatrix(ds.x) else ds.x))).to_csv(f, sep='\t',
                                                                                                         header=False
                                                                                                         )
            f.close()
        else:

            pd.DataFrame(ds.x.toarray() if scipy.sparse.isspmatrix(ds.x) else ds.x, index=ds.row_meta.index,
                         columns=ds.col_meta.index).to_csv(path,
                                                           index_label='id',
                                                           sep='\t',
                                                           doublequote=False,
                                                           compression='gzip' if output_format == 'txt.gz'
                                                           else None)
    elif output_format == 'loom':
        f = h5py.File(path, 'w')
        x = ds.x
        is_sparse = scipy.sparse.isspmatrix(x)
        is_dask = str(type(x)) == "<class 'dask.array.core.Array'>"
        save_in_chunks = is_sparse or is_dask

        dset = f.create_dataset('/matrix', shape=x.shape, chunks=(1000, 1000) if
        x.shape[0] >= 1000 and x.shape[1] >= 1000 else None,
                                maxshape=(None, x.shape[1]),
                                compression='gzip', compression_opts=9,
                                data=None if save_in_chunks else x)
        if is_dask:
            chunks = tuple((c if i == 0 else (sum(c),))
                           for i, c in enumerate(x.chunks))

            x = x.rechunk(chunks)
            xstart = 0
            xend = 0
            for xchunk in x.chunks[0]:
                xend += xchunk
                dset[slice(xstart, xend)] = x[slice(xstart, xend)].compute()
                xstart = xend


        elif is_sparse:
            dset.attrs['sparse'] = True
            # write in chunks of 1000
            start = 0
            step = min(x.shape[0], 1000)
            stop = step
            nchunks = int(np.ceil(max(1, x.shape[0] / step)))
            for i in range(nchunks):
                stop = min(x.shape[0], stop)
                dset[start:stop] = x[start:stop].toarray()
                start += step
                stop += step

        f.create_group('/layers')
        f.create_group('/row_graphs')
        f.create_group('/col_graphs')
        # for key in ds.layers:
        #     x = ds.layers[key]
        #     f.create_dataset('/layers/' + key, shape=x, chunks=(1000, 1000),
        #                      maxshape=(None, x.shape[1]),
        #                      compression='gzip', compression_opts=9,
        #                      data=x)

        wot.io.save_loom_attrs(f, False, ds.row_meta, ds.x.shape[0])
        wot.io.save_loom_attrs(f, True, ds.col_meta, ds.x.shape[1])

        f.close()

    else:
        raise Exception('Unknown file output_format')


def write_dataset_metadata(ds, path, metadata_name):
    if metadata_name not in ds.row_meta.columns:
        raise ValueError("Metadata not present in dataset: \"{}\"".format(metadata_name))
    ds.row_meta[[metadata_name]].to_csv(path, index_label='id', sep='\t', doublequote=False)


def save_loom_attrs(f, is_columns, metadata, length):
    attrs_path = '/col_attrs' if is_columns else '/row_attrs'
    f.create_group(attrs_path)

    def save_metadata_array(path, array):
        # convert object or unicode to string
        if array.dtype.kind == 'U' or array.dtype.kind == 'O':
            array = array.astype('S')
        f[path] = array

    if metadata is not None:
        save_metadata_array(attrs_path + '/' + ('id' if metadata.index.name is
                                                        None or metadata.index.name is 0 else
                                                str(metadata.index.name)), metadata.index.values)
        for name in metadata.columns:
            save_metadata_array(attrs_path + '/' + str(name), metadata[name].values)
    else:
        save_metadata_array(attrs_path + '/id', np.array(range(1, length + 1)).astype('S'))


def read_days_data_frame(path):
    return pd.read_table(path, index_col='id',
                         engine='python', sep=None, dtype={'day': np.float64})


def incorporate_days_information_in_dataset(dataset, path):
    days_data_frame = read_days_data_frame(path)
    dataset.row_meta = dataset.row_meta.join(days_data_frame)
    return dataset
