import pandas
import numpy as np
import scipy.sparse


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

        return Dataset(x=a, row_meta=pandas.DataFrame(index=row_ids),
                       col_meta=pandas.DataFrame(data={'description':
                                                           descriptions},
                                                 index=ids))


def read_sparse(path, sep=None):
    with open(path) as fp:
        data = []
        rows = []
        cols = []
        row_ids = []
        i = 0
        header = fp.readline()
        if sep is None:
            for s in ['\t', ' ']:
                test = header.split(s)
                if len(test) > 1:
                    sep = s
                    column_ids = test
                    break

        print(sep)
        column_ids = column_ids[1:]
        ncols = len(column_ids)
        for line in fp:
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
        cols = np.array(cols, dtype=np.uint32, copy=False)

        a = scipy.sparse.coo_matrix((data, (rows, cols)), shape=(i, ncols),
                                    dtype=np.float32).tocsr()
        return Dataset(x=a, row_meta=pandas.DataFrame(index=row_ids),
                       col_meta=pandas.DataFrame(index=column_ids))
