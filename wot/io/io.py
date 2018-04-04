# -*- coding: utf-8 -*-
import h5py
import pandas as pd
import wot
import numpy as np
import os
import scipy.sparse
import glob


def list_transport_maps(input_dir):
    transport_maps_inputs = []  # file, start, end
    for path in glob.glob(input_dir):
        if os.path.isfile(path):
            file_info = wot.io.get_file_basename_and_extension(os.path.basename(path))
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
    for path in glob.glob(input_dir):
        if os.path.isfile(path):
            file_info = wot.io.get_file_basename_and_extension(os.path.basename(path))
            basename = file_info[0]
            tokens = basename.split('_')
            t1 = tokens[len(tokens) - 2]
            t2 = tokens[len(tokens) - 1]

            try:
                t1 = float(t1)
                t2 = float(t2)

            except ValueError:
                continue
            transport_map = pd.read_table(path, index_col=0)
            if ids is not None and t1 == time:
                # subset rows
                transport_map = transport_map[transport_map.index.isin(ids)]
            if ids is not None and t2 == time:
                # subset columns
                transport_map = transport_map[ids]
            transport_maps_inputs.append(
                {'transport_map': transport_map, 't1': t1, 't2': t2})

    transport_maps_inputs.sort(key=lambda x: x['t1'])  # sort by t1 (start time)
    return transport_maps_inputs


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
            members.append(set_ids)
            for i in range(len(ids)):
                id = ids[i].strip()
                if id != '':
                    row_index = row_id_to_index.get(id)
                    if feature_ids is None:
                        set_ids.append(id)
                        if row_index is None:
                            row_index = len(row_id_to_index)
                            row_id_to_index[id] = row_index
                            row_ids.append(id)
                    elif row_index is not None:
                        set_ids.append(id)

        x = np.zeros(shape=(len(row_ids), len(set_names)), dtype=np.int8)
        for j in range(len(members)):
            ids = members[j]
            for id in ids:
                row_index = row_id_to_index.get(id)
                x[row_index, j] = 1

        row_meta = pd.DataFrame(index=row_ids)
        col_meta = pd.DataFrame(data={'description': set_descriptions}, index=set_names)
        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)


def read_gmx(path, feature_ids=None):
    with open(path) as fp:
        ids = fp.readline().split('\t')
        descriptions = fp.readline().split('\t')
        ncols = len(ids)
        ids[len(ids) - 1] = ids[len(ids) - 1].rstrip()
        descriptions[len(ids) - 1] = descriptions[len(ids) - 1].rstrip()
        row_id_to_index = {}
        x = None
        array_of_arrays = None
        if feature_ids is not None:
            for i in range(len(feature_ids)):
                row_id_to_index[feature_ids[i]] = i
            x = np.zeros(shape=(len(feature_ids), ncols), dtype=np.int8)
        else:
            array_of_arrays = []
        for line in fp:
            tokens = line.split('\t')
            for j in range(ncols):
                value = tokens[j].strip()
                if value != '':
                    row_index = row_id_to_index.get(value)
                    if feature_ids is None:
                        if row_index is None:
                            row_index = len(row_id_to_index)
                            row_id_to_index[value] = row_index
                            array_of_arrays.append(np.zeros(shape=(ncols,), dtype=np.int8))
                        array_of_arrays[row_index][j] = 1
                    elif row_index is not None:
                        x[row_index, j] = 1
        if feature_ids is None:
            feature_ids = np.empty(len(row_id_to_index), dtype='object')
            for rid in row_id_to_index:
                feature_ids[row_id_to_index[rid]] = rid
        if array_of_arrays is not None:
            x = np.array(array_of_arrays)
        row_meta = pd.DataFrame(index=feature_ids)
        col_meta = pd.DataFrame(data={'description': descriptions},
                                index=ids)
        return wot.Dataset(x, row_meta=row_meta, col_meta=col_meta)


def read_dataset(path, chunks=(500, 500), h5_x=None, h5_row_meta=None,
                 h5_col_meta=None, use_dask=False, genome10x=None, row_filter=None, col_filter=None):
    path = str(path)
    basename_and_extension = get_file_basename_and_extension(path)
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

        import scipy.io
        x = scipy.io.mmread(path)
        if col_meta is None:
            col_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=x.shape[1], step=1))
        if row_meta is None:
            row_meta = pd.DataFrame(index=pd.RangeIndex(start=0, stop=x.shape[0], step=1))
        return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
    elif ext == 'hdf5' or ext == 'h5' or ext == 'loom':
        f = h5py.File(path, 'r')
        if genome10x is not None or f.get('/mm10') is not None:
            genome10x = '/mm10' if genome10x is None else genome10x
            group = f['/' + genome10x]

            from scipy.sparse import csr_matrix
            M, N = group['shape'][()]
            data = group['data'][()]
            # if data.dtype == np.dtype('int32'):
            #     data = data.astype('float32')
            x = csr_matrix((data, group['indices'][()], group['indptr'][()]), shape=(N, M))
            col_meta = pd.DataFrame(index=group['genes'][()].astype(str),
                                    data={'gene_names': group['gene_names'][()].astype(str)})
            row_meta = pd.DataFrame(index=group['barcodes'][()].astype(str))
            f.close()
            return wot.Dataset(x=x, row_meta=row_meta, col_meta=col_meta)
        if ext == 'loom':
            h5_x = '/matrix'
            h5_row_meta = '/row_attrs'
            h5_col_meta = '/col_attrs'

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
            import scipy
            if x.attrs.get('sparse') and row_filter is None and col_filter is None:
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
                if row_filter is None and col_filter is None:
                    x = x[()]
                elif row_filter is not None and col_filter is not None:
                    x = x[row_attrs['indices']]
                    x = x[:, col_attrs['indices']]
                elif row_filter is not None:
                    x = x[row_attrs['indices']]
                else:
                    x = x[:, col_attrs['indices']]
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


def check_file_extension(name, format):
    expected = None
    if format == 'txt':
        expected = '.txt'
    elif format == 'txt.gz':
        expected = '.txt.gz'
    elif format == 'loom':
        expected = '.loom'
    elif format == 'gct':
        expected = '.gct'
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


def write_dataset(ds, path, output_format='txt', txt_full=False, progress=None):
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
                                                                                                         header=False)
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
                if progress:
                    progress(xstart)

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
                if progress:
                    progress(start)

        # f.create_group('/layers')
        # for key in ds.layers:
        #     x = ds.layers[key]
        #     f.create_dataset('/layers/' + key, shape=x, chunks=(1000, 1000),
        #                      maxshape=(None, x.shape[1]),
        #                      compression='gzip', compression_opts=9,
        #                      data=x)

        wot.io.save_loom_attrs(f, False, ds.row_meta)
        wot.io.save_loom_attrs(f, True, ds.col_meta)

        f.close()

    else:
        raise Exception('Unknown file output_format')


def save_loom_attrs(f, is_columns, metadata):
    attrs_path = '/col_attrs' if is_columns else '/row_attrs'
    f.create_group(attrs_path)

    def save_metadata_array(path, array):
        # convert object or unicode to string
        if array.dtype.kind == 'U' or array.dtype.kind == 'O':
            array = array.astype('S')
        f[path] = array

    if metadata is not None:
        save_metadata_array(attrs_path + '/' + ('id' if metadata.index.name is
                                                        None else
                                                metadata.index.name), metadata.index.values)
        for name in metadata.columns:
            save_metadata_array(attrs_path + '/' + name, metadata[name].values)
