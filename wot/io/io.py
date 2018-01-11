# -*- coding: utf-8 -*-
import h5py
import dask.array as da
import dask.dataframe as dd
import pandas as pd
import wot
import numpy as np
import os


def read_gene_sets(path, feature_ids=None):
    path = str(path)
    basename_and_extension = get_file_basename_and_extension(path)
    ext = basename_and_extension[1]
    if ext == 'gmt':
        return read_gmt(path, feature_ids)
    elif ext == 'gmx':
        return read_gmx(path, feature_ids)
    else:
        raise ValueError('Unknown file format')


def read_gmt(path, feature_ids=None):
    with open(path) as fp:
        row_id_to_index = {}
        if feature_ids is not None:
            for i in range(len(feature_ids)):
                row_id_to_index[feature_ids[i]] = i
            row_ids = feature_ids
        else:
            row_ids = []

        members = []
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
            set_ids = []
            for i in range(len(ids)):
                id = ids[i].strip()
                if id != '':
                    set_ids.append(id)
                    if feature_ids is None:
                        index = row_id_to_index.get(id)
                        if index is None:
                            index = len(row_id_to_index)
                            row_id_to_index[id] = index
                            row_ids.append(id)
            members.append(set_ids)
        x = np.zeros(shape=(len(row_ids), len(set_names)), dtype=np.int8)
        for j in range(len(members)):
            ids = members[j]
            for id in ids:
                row_index = row_id_to_index.get(id)
                if row_index is not None:
                    x[row_index, j] = 1

        return wot.Dataset(x=x,
                           row_meta=pd.DataFrame(index=row_ids),
                           col_meta=pd.DataFrame(data={'description':
                                                           set_descriptions},
                                                 index=set_names))


def read_gmx(path, feature_ids=None):
    with open(path) as fp:
        ids = fp.readline().split('\t')
        descriptions = fp.readline().split('\t')
        ncols = len(ids)
        ids[len(ids) - 1] = ids[len(ids) - 1].rstrip()
        descriptions[len(ids) - 1] = descriptions[len(ids) - 1].rstrip()
        row_id_to_index = {}
        matrix = None
        array_of_arrays = None
        if feature_ids is not None:
            for i in range(len(feature_ids)):
                row_id_to_index[feature_ids[i]] = i
            matrix = np.zeros(shape=(len(feature_ids), ncols), dtype=np.int8)
        else:
            array_of_arrays = []
        for line in fp:
            tokens = line.split('\t')
            for j in range(ncols):
                value = tokens[j].strip()
                if value != '':
                    i = row_id_to_index.get(value)
                    if feature_ids is None:
                        if i is None:
                            i = len(row_id_to_index)
                            row_id_to_index[value] = i
                            array_of_arrays.append(np.zeros(shape=(ncols,),
                                                            dtype=np.int8))
                    elif i is not None:
                        matrix[i, j] = 1
        if feature_ids is None:
            feature_ids = np.empty(len(row_id_to_index), dtype='object')
            for rid in row_id_to_index:
                feature_ids[row_id_to_index[rid]] = rid
        return wot.Dataset(
            x=np.array(array_of_arrays),
            row_meta=pd.DataFrame(index=feature_ids),
            col_meta=pd.DataFrame(data={'description': descriptions},
                                  index=ids))


def read_dataset(path, chunks=(200, 200), h5_x=None, h5_row_meta=None,
                 h5_col_meta=None, use_dask=True):
    path = str(path)
    basename_and_extension = get_file_basename_and_extension(path)
    ext = basename_and_extension[1]

    if ext == 'mtx':
        # look for .barcodes.txt and .genes.txt
        sp = os.path.split(path)
        row_meta = None
        for f in (
                os.path.join(sp[0],
                             basename_and_extension[0] + ".barcodes.tsv"),
                os.path.join(sp[0],
                             basename_and_extension[0] + ".barcodes.txt"),
                os.path.join(sp[0], "barcodes.tsv")):
            if os.path.isfile(f) or os.path.isfile(f + '.gz'):
                data = np.genfromtxt(f
                                     if os.path.isfile(
                    f) else f + '.gz', dtype=str)
                if len(data.shape) > 1:
                    data = data[:, 0]  # TODO add other columns to df
                row_meta = pd.DataFrame(index=data)
                break
        col_meta = None
        for f in (os.path.join(sp[0], basename_and_extension[0] +
                                      ".genes.tsv"),
                  os.path.join(sp[0], basename_and_extension[0] +
                                      ".genes.txt"),
                  os.path.join(sp[0], "genes.tsv")):
            if os.path.isfile(f) or os.path.isfile(f + '.gz'):
                data = np.genfromtxt(f
                                     if os.path.isfile(
                    f) else f + '.gz', dtype=str)
                index = data
                df_data = {}
                if len(data.shape) > 1:
                    index = data[:, 0]  # TODO add other columns to df
                    df_data = {'Symbol': data[:, 1]}

                col_meta = pd.DataFrame(index=index, data=df_data)
                break

        if path.endswith('.gz'):
            import gzip
            stream = gzip.open(path, 'rb')
        elif path.endswith('.bz2'):
            import bz2
            stream = bz2.BZ2File(path, 'rb')
        else:
            stream = open(path, 'rb')
        from numpy.compat import asstr
        # % % MatrixMarket matrix coordinate real general
        # %
        # 32738 2700 2286884
        line = stream.readline()
        mmid, matrix, format, field, symmetry = \
            [asstr(part.strip()) for part in line.split()]
        if not mmid.startswith('%%MatrixMarket'):
            raise ValueError('source is not in Matrix Market format')
        if not matrix.lower() == 'matrix':
            raise ValueError("Problem reading file header: " + line)
        # skip comments
        while line.startswith(b'%'):
            line = stream.readline()

        line = line.split()
        if not len(line) == 3:
            raise ValueError("Header line not of length 3: " + line)
        rows, cols, entries = map(int, line)

        x = np.zeros(shape=(cols, rows), dtype=np.float32)

        while line:
            line = stream.readline()
            if not line or line.startswith(b'%'):
                continue

            l = line.split()
            i, j = map(int, l[:2])
            i, j = i - 1, j - 1
            aij = float(l[2])
            x[j, i] = aij
        stream.close()
        if col_meta is None:
            col_meta = pd.DataFrame(index=pd.RangeIndex(start=0,
                                                        stop=x.shape[1],

                                                        step=1))
        if row_meta is None:
            row_meta = pd.DataFrame(
                index=pd.RangeIndex(start=0, stop=x.shape[0], step=1))

        if use_dask:
            return wot.Dataset(x=da.from_array(x, chunks=chunks),
                               row_meta=dd.from_pandas(row_meta,
                                                       npartitions=4,
                                                       sort=False),
                               col_meta=dd.from_pandas(col_meta,
                                                       npartitions=4,
                                                       sort=False))
        else:
            return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
    elif ext == 'hdf5' or ext == 'h5' or ext == 'loom':
        if ext == 'loom':
            h5_x = '/matrix'
            h5_row_meta = '/row_attrs'
            h5_col_meta = '/col_attrs'
        f = h5py.File(path, 'r')
        g = f[h5_row_meta]
        data = {}
        for key in g:
            values = g[key][()]
            if values.dtype.kind == 'S':
                values = values.astype(str)
            data[key] = values

        row_meta = pd.DataFrame(data,
                                index=pd.RangeIndex(start=0, stop=f[
                                    h5_x].shape[0],
                                                    step=1))
        g = f[h5_col_meta]
        data = {}
        for key in g:
            values = g[key][()]
            if values.dtype.kind == 'S':
                values = values.astype(str)
            data[key] = values
        col_meta = pd.DataFrame(data,
                                index=pd.RangeIndex(start=0, stop=f[
                                    h5_x].shape[1],
                                                    step=1))
        if not use_dask:
            x = f[h5_x][()]
            f.close()
            return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
        else:
            f = h5py.File(path, 'r')
            x = da.from_array(f[h5_x], chunks=chunks)
            # TODO load in chunks
            row_meta = dd.from_pandas(row_meta, npartitions=4, sort=False)
            col_meta = dd.from_pandas(col_meta, npartitions=4, sort=False)
        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
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
            dask_arrays = []
            np_arrays = []
            for line in fp:
                line = line.rstrip()
                if line != '':
                    tokens = line.split(sep)
                    row_ids.append(tokens[0])
                    np_arrays.append(np.array(tokens[1:], dtype=np.float32))
                    if use_dask and len(np_arrays) == 100:
                        dask_arrays.append(
                            da.from_array(np.array(np_arrays), chunks=chunks))
                        np_arrays = []
                    i += 1
            if use_dask:
                if len(np_arrays) > 0:
                    dask_arrays.append(
                        da.from_array(np.array(np_arrays), chunks=chunks))
                x = da.concatenate(dask_arrays, axis=0)
                return wot.Dataset(x=x, row_meta=dd.from_pandas(pd.DataFrame(
                    index=row_ids), npartitions=4, sort=False),
                                   col_meta=dd.from_pandas(pd.DataFrame(
                                       index=column_ids), npartitions=4,
                                       sort=False))
            else:
                return wot.Dataset(x=np.array(np_arrays),
                                   row_meta=pd.DataFrame(index=row_ids),
                                   col_meta=pd.DataFrame(index=column_ids))


def check_file_extension(name, format):
    expected = None
    if format == 'txt':
        expected = '.txt'
    elif format == 'txt.gz':
        expected = '.txt.gz'
    elif format == 'loom':
        expected = '.loom'
    if expected is not None:
        if not str(name).lower().endswith(expected):
            name += expected
    return name


def get_file_basename_and_extension(name):
    dot_index = name.rfind('.')
    ext = ''
    basename = name
    if dot_index != -1:
        ext = name[dot_index + 1:]
        if ext == 'gz':
            return get_file_basename_and_extension(name[0:dot_index])

    if dot_index != -1:
        basename = name[0:dot_index]
    return basename, ext


def write_dataset(ds, path, output_format='txt'):
    if output_format == 'txt' or output_format == 'txt.gz':
        f = pd.DataFrame(data=ds.x, index=ds.row_meta.index,
                         columns=ds.col_meta.index)
        f.to_csv(path,
                 index_label="id",
                 sep='\t',
                 compression='gzip' if output_format == 'txt.gz'
                 else None)
    elif output_format == 'loom':
        f = h5py.File(path, 'w')
        x = ds.x
        f.create_dataset("/matrix", shape=x.shape, chunks=(1000, 1000) if
        x.shape[0] >= 1000 and x.shape[1] >= 1000 else None,
                         maxshape=(None, x.shape[1]),
                         compression="gzip", compression_opts=9,
                         data=x)
        # f.create_group('/layers')
        # for key in ds.layers:
        #     x = ds.layers[key]
        #     f.create_dataset("/layers/" + key, shape=x, chunks=(1000, 1000),
        #                      maxshape=(None, x.shape[1]),
        #                      compression="gzip", compression_opts=9,
        #                      data=x)

        f.create_group('/row_attrs')
        f.create_group('/col_attrs')

        def save_metadata_array(path, array):
            # convert object or unicode to string
            if array.dtype.kind == 'U' or array.dtype.kind == 'O':
                array = array.astype('S')
            f[path] = array

        save_metadata_array('/row_attrs/' + ('id' if ds.row_meta.index.name is
                                                     None else
        ds.row_meta.index.name),
                            ds.row_meta.index.values)
        save_metadata_array('/col_attrs/' + ('id' if ds.col_meta.index.name is
                                                     None else
        ds.col_meta.index.name), ds.col_meta.index.values)

        for name in ds.row_meta.columns:
            save_metadata_array('/row_attrs/' + name,
                                ds.row_meta[name].values())
        for name in ds.col_meta.columns:
            save_metadata_array('/col_attrs/' + name,
                                ds.col_meta[name].values())
        f.close()

    else:
        raise Exception('Unknown file output_format')
