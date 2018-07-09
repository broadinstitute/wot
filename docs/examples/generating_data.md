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
for future use, with the following snippet :

```python
import numpy
import wot.simulate
from math import sin
from numpy.random import randint

covariates_count = 10

tips = [
        [ [ -3, 0, 0, 0, 0 ], [ 0, 25, 0, 0, 0 ], [ -10, 50, 0, 0, 0 ], [ -20, 75, 10, 0, 0], [ -25, 100, 50, 0, 0] ],
        [ [ -3, 0, 0, 0, 0 ], [ 0, 25, 0, 0, 0 ], [ -10, 50, 0, 0, 0 ], [  -5, 75, 0, 10, 0], [  -5, 100, 0, 50, 0] ],
        [ [ -3, 0, 0, 0, 0 ], [ 0, 25, 0, 0, 0 ], [   3, 50, 0, 0, 0 ], [  10, 75, 0, 0, 10], [  20, 100, 0, 0, 50] ]
       ]
times = [ [ i / (len(k) - 1) for i in range(len(k)) ] for k in tips ]

N = 101
timepoints = [ i / ( N - 1) for i in range(N) ]

means = numpy.array(
        [ wot.simulate.interp(timepoints, times[k], tips[k], method='linear')
            for k in range(len(tips)) ])

covs = [ 15, 8, 5, 5, 5 ]
covs = [ [ c / 2 * (2 + sin(80*t)) for c in covs ] for t in timepoints ]

random_population = lambda t : \
    wot.simulate.multivariate_normal_mixture(
             means[:,t], [covs[t]] * len(means), size=1000)
data = [ random_population(t) for t in range(N) ]

data_to_dataset = lambda i : \
    wot.dataset_from_x(data[i], row_prefix="cell_g{:02}_".format(i))
dataset_list = [ data_to_dataset(i) for i in range(N) ]

for i in range(N):
    wot.set_cell_metadata(dataset_list[i], 'day', i)
    covariates = randint(0, covariates_count - 1, size=data[i].shape[0])
    wot.set_cell_metadata(dataset_list[i], 'covariate', covariates)
ds = wot.merge_datasets(*dataset_list)

wot.io.write_dataset(ds, 'matrix.txt', txt_full = False)
wot.io.write_dataset_metadata(ds, 'days.txt', 'day')
wot.io.write_dataset_metadata(ds, 'covariate.txt', 'covariate')
```


## Plotting two features ##

The best way to visualize you data is using the Force-directed Layout Embedding
to project it into 2 dimensions.

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
cell_colors = numpy.asarray(ds.row_meta['day'])
cell_colors = cell_colors / max(cell_colors)
cell_colors = [ wot.graphics.color_mix(color1, color2, d)
        for d in cell_colors ]
wot.set_cell_metadata(ds, 'color', cell_colors)

pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds, x=0, y=1)
pyplot.show()
```
