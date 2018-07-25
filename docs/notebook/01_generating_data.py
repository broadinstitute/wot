# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
covariate_file = 'covariate.txt'
gene_sets_file = 'gene_sets.gmt'

number_of_timepoints = 51
covariates_count = 5

gene_sets = {
        'Stem cells': [ 'Stem_gene' ],
        'Myeloid stem cells': [ 'Myeloid_gene' ],
        'Red blood cells': [ 'RBC_gene' ],
        'Granulocytes': [ 'Granulo_gene' ],
        'Lymphocytes': [ 'Lympho_gene' ],
        }
# --------------------------------------

import numpy
import wot.simulate
from numpy.random import randint
from numpy.random import random

gene_names = [ 'X_gene', 'Y_gene',
        'RBC_gene', 'Granulo_gene', 'Lympho_gene',
        'Myeloid_gene', 'Stem_gene' ]

tips = [
        [ [ -.3, 0, 0, 0, 0, 0, 10 ], [ 0, 2.5, 0, 0, 0, 2, 3 ], [ -1, 5, 0, 0, 0, 4, 1 ], [ -2, 7.5, 5, 0, 0, 0, 0 ], [ -2.5, 10, 10, 0, 0, 0, 0 ] ],
        [ [ -.3, 0, 0, 0, 0, 0, 10 ], [ 0, 2.5, 0, 0, 0, 2, 3 ], [ -1, 5, 0, 0, 0, 4, 1 ], [-.5, 7.5, 0, 5, 0, 0, 0 ], [  -.5, 10, 0, 10, 0, 0, 0 ] ],
        [ [ -.3, 0, 0, 0, 0, 0, 10 ], [ 0, 2.5, 0, 0, 0, 2, 3 ], [ .3, 5, 0, 0, 5, 1, 1 ], [  1, 7.5, 0, 0, 9, 0, 0 ], [    2, 10, 0, 0, 10, 0, 0 ] ]
       ]
times = [ [ i / (len(k) - 1) for i in range(len(k)) ] for k in tips ]

N = number_of_timepoints
timepoints = [ i / ( N - 1) for i in range(N) ]

means = numpy.array(
        [ wot.simulate.interp(timepoints, times[k], tips[k],
              method='linear', smooth=(N // 10))
            for k in range(len(tips)) ])
means = numpy.asarray([ means[:,t] for t in range(N) ])

covs = [ .08, .1, .04, .04, .04, .03, .05 ]
covs = [ [ c * (random() + .5) for c in covs ] for t in timepoints ]

sizes =  [ 5000 + randint(-100, 100) for _ in range(N) ]
splits = [ 0 ]
for i in range(N):
    splits.append(sizes[i] + splits[-1])

data = wot.simulate.multivariate_normal_evolving_mixture(
        means, [ [ c ] * means.shape[1] for c in covs ], size=sizes)
data_to_dataset = lambda i : \
        wot.dataset_from_x(data[splits[i]:splits[i+1],:],
                row_prefix="cell_g{:02}_".format(i),
                columns=gene_names)
dataset_list = [ data_to_dataset(i) for i in range(N) ]


for i in range(N):
    wot.set_cell_metadata(dataset_list[i], 'day', i)
    covariates = randint(0, covariates_count - 1, size=splits[i+1] - splits[i])
    wot.set_cell_metadata(dataset_list[i], 'covariate', covariates)
ds = wot.merge_datasets(*dataset_list)

wot.io.write_gene_sets(gene_sets, gene_sets_file, "gmt")
wot.io.write_dataset(ds, matrix_file, txt_full = False)
wot.io.write_dataset_metadata(ds, days_file, 'day')
wot.io.write_dataset_metadata(ds, covariate_file, 'covariate')
