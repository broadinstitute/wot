import pandas as pd
import numpy as np
import scipy.sparse
import wot
import os


def read_gmx(path):
    with open(path) as fp:
        ids = fp.readline().split('\t')

        descriptions = fp.readline().split('\t')
        data = []
        rows = []
        cols = []
        row_id_to_index = {}

        ncols = len(ids)
        nrows = 0
        ids[len(ids) - 1] = ids[len(ids) - 1].rstrip()
        descriptions[len(ids) - 1] = descriptions[len(ids) - 1].rstrip()
        for line in fp:
            tokens = line.split('\t')
            for j in range(ncols):
                value = tokens[j].strip()
                if value != '':
                    cols.append(j)
                    data.append(1)
                    index = row_id_to_index.get(value)
                    if index is None:
                        index = nrows
                        row_id_to_index[value] = index
                        nrows += 1
                    rows.append(index)
        a = scipy.sparse.coo_matrix((data, (rows, cols)), shape=(nrows, ncols),
                                    dtype=np.uint8).tocsr()

        row_ids = np.empty(len(row_id_to_index), dtype='object')
        for rid in row_id_to_index:
            row_ids[row_id_to_index[rid]] = rid

        return wot.Dataset(x=a, row_meta=pd.DataFrame(index=row_ids),
                           col_meta=pd.DataFrame(data={'description':
                                                           descriptions},
                                                 index=ids))


def get_file_basename_and_extension(name):
    dot_index = name.rfind('.')
    ext = ''
    basename = name
    if dot_index != -1:
        ext = name[dot_index:]
        if ext == '.gz':
            return get_file_basename_and_extension(name[0:dot_index])

    if dot_index != -1:
        basename = name[0:dot_index]
    return (basename, ext)


def write_dataset(ds, path, output_format='txt'):
    if output_format == 'txt' or output_format == 'txt.gz':
        f = pd.DataFrame(data=ds.x, index=ds.row_meta.index,
                         columns=ds.col_meta.index)
        f.to_csv(path,
                 index_label="id",
                 sep='\t',
                 compression='gzip' if output_format == 'txt.gz'
                 else None)
    if output_format == 'loom':
        import h5py
        f = h5py.File(path, 'w')

        x = ds.x
        x = x.todense() if scipy.sparse.issparse(x) else x
        f.create_dataset("matrix", shape=x, chunks=(100, 100),
                         maxshape=(None, x.shape[1]),
                         compression="gzip", compression_opts=9,
                         data=x)
        f.create_group('/layers')
        for key in ds.layers:
            x = ds.layers[key]
            x = x.todense() if scipy.sparse.issparse(x) else x
            f.create_dataset("/layers/" + key, shape=x, chunks=(100, 100),
                             maxshape=(None, x.shape[1]),
                             compression="gzip", compression_opts=9,
                             data=x)

        f.create_group('/row_attrs')
        f.create_group('/col_attrs')
        rids = ds.row_meta.index.values
        if rids.dtype.kind == 'U' or rids.dtype.kind == 'O':
            rids = rids.astype(np.string_)
        f['/row_attrs/id'] = rids
        cids = ds.col_meta.index.values
        if cids.dtype.kind == 'U' or cids.dtype.kind == 'O':
            cids = cids.astype(np.string_)
        f['/col_attrs/id'] = cids
        for name in ds.row_meta.columns:
            v = ds.row_meta[name].values()
            if v.dtype.kind == 'U' or v.dtype.kind == 'O':
                v = v.astype(np.string_)
            f['/row_attrs/' + name] = v
        for name in ds.col_meta.columns:
            v = ds.col_meta[name].values()
            if v.dtype.kind == 'U' or v.dtype.kind == 'O':
                v = v.astype(np.string_)
            f['/col_attrs/' + name] = v
        f.close()

    else:
        raise Exception('Unknown file output_format')


def read_dataset(path, sep=None, dtype=np.float32, is_sparse=True):
    path = str(path)
    basename_and_extension = get_file_basename_and_extension(path)
    ext = basename_and_extension[1]
    if ext == 'mtx':
        x = scipy.sparse.csr_matrix(scipy.io.mmread(path).astype(dtype))
        # look for .barcodes.txt and .genes.txt
        sp = os.path.split(path)
        barcodes_files = (
        os.path.join(sp, basename_and_extension[0] + ".barcodes.tsv"),
        os.path.join(sp, basename_and_extension[0] + ".barcodes.txt"),
        os.path.join(sp, "barcodes.tsv"))
        row_meta = None
        for f in barcodes_files:
            if os.path.isfile(f):
                row_meta = pd.DataFrame(index=np.genfromtxt(f, dtype=str))
                break
        if row_meta is None:
            row_meta = pd.DataFrame(
                index=pd.RangeIndex(start=0, stop=x.shape[0], step=1))
        genes_files = (os.path.join(sp, basename_and_extension[0] +
                                    ".genes.tsv"),
                       os.path.join(sp, basename_and_extension[0] +
                                    ".genes.txt"),
                       os.path.join(sp, "genes.tsv"))
        col_meta = None
        for f in genes_files:
            if os.path.isfile(f):
                data = np.genfromtxt(f, dtype=str)
                if len(data.shape) > 1:
                    data = data[:, 0]  # TODO add other columns to df
                col_meta = pd.DataFrame(index=data)
                break
        if col_meta is None:
            col_meta = pd.DataFrame(index=pd.RangeIndex(start=0,
                                                        stop=x.shape[1],
                                                        step=1))
        return wot.Dataset(x=x, row_meta=row_meta,
                           col_meta=col_meta)

    with open(path) as fp:
        data = []
        row_ids = []
        if is_sparse:
            rows = []
            cols = []

        header = fp.readline()
        if sep is None:
            for s in ['\t', ' ', ',']:
                test = header.split(s)
                if len(test) > 1:
                    sep = s
                    column_ids = test
                    break

        column_ids = column_ids[1:]
        column_ids[len(column_ids) - 1] = column_ids[
            len(column_ids) - 1].rstrip()
        ncols = len(column_ids)
        i = 0
        for line in fp:
            line = line.rstrip()
            if line != '':
                tokens = line.split(sep)
                row_ids.append(tokens[0])
                if is_sparse:
                    for j in range(ncols):
                        value = float(tokens[j + 1])
                        if value != 0:
                            cols.append(j)
                            data.append(value)
                            rows.append(i)
                else:
                    data.append(np.array(tokens[1:], dtype=dtype))
                i += 1

        if is_sparse:
            data = np.array(data, copy=False)
            rows = np.array(rows, dtype=np.uint32, copy=False)
            cols = np.array(cols, dtype=np.uint64, copy=False)
            x = scipy.sparse.coo_matrix((data, (rows, cols)), shape=(i, ncols),
                                        dtype=dtype).tocsr()
        else:
            x = np.array(data, dtype=dtype, copy=False)
        return wot.Dataset(x=x, row_meta=pd.DataFrame(index=row_ids),
                           col_meta=pd.DataFrame(index=column_ids))
