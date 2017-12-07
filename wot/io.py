import pandas
import numpy as np
import scipy.sparse
import wot


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

        return wot.Dataset(x=a, row_meta=pandas.DataFrame(index=row_ids),
                           col_meta=pandas.DataFrame(data={'description':
                                                               descriptions},
                                                     index=ids))


def write_dataset(ds, path, output_format='txt'):
    if output_format == 'txt' or output_format == 'txt.gz':
        f = pandas.DataFrame(data=ds.x, index=ds.row_meta.index,
                             columns=ds.col_meta.index)
        f.to_csv(path + '.txt' + ('.gz' if output_format == 'txt.gz' else
        ''),
                 index_label="id",
                 sep='\t',
                 compression='gzip' if output_format == 'txt.gz'
                 else None)
    if output_format == 'loom':
        import h5py
        f = h5py.File(path + '.loom', 'w')
        f.create_dataset("matrix", shape=ds.x.shape, chunks=(100, 100),
                         maxshape=(None, ds.x.shape[1]),
                         compression="gzip", compression_opts=9,
                         data=ds.x.todense() if
                         scipy.sparse.issparse(ds.x) else ds.x)
        f.create_group('/row_attrs')
        f.create_group('/col_attrs')
        for name in ds.row_meta.columns:
            f['/row_attrs/' + name] = ds.row_meta[name].values()
        for name in ds.col_meta.columns:
            f['/col_attrs/' + name] = ds.col_meta[name].values()
        f.close()

    else:
        raise Exception('Unknown file output_format')


def read_dataset(path, sep=None):
    with open(path) as fp:
        data = []
        rows = []
        cols = []
        row_ids = []
        i = 0
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
        for line in fp:
            line = line.rstrip()
            if line != '':
                tokens = line.split(sep)
                row_ids.append(tokens[0])
                for j in range(ncols):
                    value = float(tokens[j + 1])
                    if value != 0:
                        cols.append(j)
                        data.append(value)
                        rows.append(i)
                i += 1

        data = np.array(data, copy=False)
        rows = np.array(rows, dtype=np.uint32, copy=False)
        cols = np.array(cols, dtype=np.uint64, copy=False)
        a = scipy.sparse.coo_matrix((data, (rows, cols)), shape=(i, ncols),
                                    dtype=np.float32).tocsr()
        return wot.Dataset(x=a, row_meta=pandas.DataFrame(index=row_ids),
                           col_meta=pandas.DataFrame(index=column_ids))
