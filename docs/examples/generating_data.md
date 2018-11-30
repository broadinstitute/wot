---
noheader: true
layout: documentation
location: Examples
---

# Generating simulated data
---------------------------

## Gaussians around a path ##

One way to generate simulated data is to pick a few "paths" in gene-expression
space, defined by piecewise linear interpolation between a set of fixed points,
and then generate cells for those paths by simply adding Gaussian noise to the
paths.

You can generate this kind of simulation, along with days and covariates data
for future use, with the following script :

```python
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
times = [ numpy.linspace(0, 1, num=len(k)) for k in tips ]

N = number_of_timepoints
timepoints = numpy.linspace(0, 1, num=N)

means = numpy.array(
        [ wot.simulate.interp(timepoints, times[k], tips[k],
              method='linear', smooth=(N // 10))
            for k in range(len(tips)) ])
means = numpy.asarray([ means[:,t] for t in range(N) ])

covs = [ .08, .1, .04, .04, .04, .03, .05 ]
covs = [[ c * (random() + .5) for c in covs ]] * len(tips)

sizes =  [ 5000 + randint(-100, 100) for _ in range(N) ]

data = [ wot.simulate.multivariate_normal_mixture(means[i],
    covs, size=sizes[i]) for i in range(N) ]
data_to_dataset = lambda i : \
        wot.dataset_from_x(data[i],
                row_prefix="cell_g{:02}_".format(i),
                columns=gene_names)
dataset_list = [ data_to_dataset(i) for i in range(N) ]


for i in range(N):
    wot.set_cell_metadata(dataset_list[i], 'day', i)
    covariates = randint(0, covariates_count, size=sizes[i])
    wot.set_cell_metadata(dataset_list[i], 'covariate', covariates)
ds = wot.merge_datasets(*dataset_list)

wot.io.write_gene_sets(gene_sets, gene_sets_file, "gmt")
wot.io.write_dataset(ds, matrix_file)
wot.io.write_dataset_metadata(ds.obs, days_file, 'day')
wot.io.write_dataset_metadata(ds.obs, covariate_file, 'covariate')
```


## Plotting two features ##

One of the best ways to visualize your dataset is to use Force-directed Layout Embedding
to project your data into 2 dimensions.

However, you might also want to just plot 2 specific features of you dataset.
The following snippet does that :

```python
import wot.io
import wot.graphics
from matplotlib import pyplot

ds = wot.io.read_dataset('matrix.txt')

pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds, x=0, y=1)
pyplot.show()
```

## Coloring the graph ##

You can easily add colors to this kind of plot by using the 'color' metadata
in your dataset. For instance, here is how to color according to time :

```sh
import numpy
import wot.io
import wot.graphics
from matplotlib import pyplot

ds = wot.io.read_dataset('matrix.txt')

color1 = [ .08, .34, .59 ] # first color
color2 = [ .08, .59, .34 ] # final color
wot.io.incorporate_days_information_in_dataset(ds, 'days.txt')

# you can use any of the columns here, or metadata information :
cell_colors = numpy.asarray(ds.obs['day'])
cell_colors = cell_colors / max(cell_colors)
cell_colors = [ wot.graphics.color_mix(color1, color2, d)
        for d in cell_colors ]
wot.set_cell_metadata(ds, 'color', cell_colors)

pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds, x=0, y=1)
pyplot.show()
```

## Result ##

![generated data plot]({{site.baseurl}}/images/generated_pop.png)
