import pandas


class Dataset:
    """
       A Dataset consists of a 2-d matrix x, row metadata, and column metadata.
       Args:
           x (ndarray|scipy.spmatrix): 2-d matrix
           row_meta (pandas.DataFrame) Row metadata
           col_meta (pandas.DataFrame) Column metadata
       """

    def __init__(self, x, row_meta, col_meta):
        self.x = x
        self.row_meta = row_meta
        self.col_meta = col_meta
        if x.shape[0] != row_meta.shape[0]:
            raise Exception('Row dimensions do not match: ' + str(x.shape[0]) +
                            '!=' + str(row_meta.shape[0]))
        if x.shape[1] != col_meta.shape[0]:
            raise Exception(
                'Column dimensions do not match: ' + str(x.shape[1]) +
                '!=' + str(col_meta.shape[0]))

    def transpose(self):
        self.x = self.x.transpose()
        tmp = self.row_meta
        self.row_meta = self.col_meta
        self.col_meta = tmp

    def write(self, path, output_format='txt'):
        if output_format == 'txt' or output_format == 'txt.gz':
            f = pandas.DataFrame(data=self.x, index=self.row_meta.index,
                                 columns=self.col_meta.index)
            f.to_csv(path + '.txt' + ('.gz' if output_format == 'txt.gz' else
            ''),
                     index_label="id",
                     sep='\t',
                     compression='gzip' if output_format == 'txt.gz'
                     else None)
        else:
            raise Exception('Unknown file output_format')
