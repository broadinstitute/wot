---
noheader: true
layout: documentation
location: Examples
---

# Plotting shared ancestry
-------------------------

## Compute trajectories ##

The first thing you need to do is compute the trajectories of the cell sets.

```sh
wot trajectory --matrix matrix.txt --cell_days days.txt \
 --cell_set cell_sets.gmt --tmap tmaps --time 50
```

Note : You need to have a transport map configuration called `tmaps.yml` for this to work. Please refer to the [CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps) for instructions on how to compute transport maps.

## Compute shared ancestry ##

```python
import wot
import numpy
import pandas
from matplotlib import pyplot

trajectory_filename = 'wot_trajectory.txt'
matrix_filename = 'matrix.txt'

cell_set_1 = 'Red blood cells'
cell_set_2 = 'Granulocytes'

red = (1, 0, 0)
blue = (0, 0, 1)

traj = pandas.read_table(trajectory_filename, index_col='id')
ds = wot.io.read_dataset(matrix_filename)
wot.set_cell_metadata(ds, 'color', '#80808050')

pyplot.figure(figsize=(5, 5))
pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds)

probabilities = traj.loc[:, [cell_set_1, cell_set_2]].values
alphas = numpy.amax(probabilities, axis=1)
t = numpy.log((probabilities.T[0] + 1e-8) / (probabilities.T[1] + 1e-18)) / 8
t = numpy.clip(t + .5, 0, 1)
alphas = numpy.clip(alphas / numpy.mean(alphas), 0, 1)
colors = [ wot.graphics.hexstring_of_rgba([t[i], 0, 1-t[i], alphas[i]])
           for i in range(len(t)) ]
ds.row_meta.loc[:, 'color']  = colors

wot.graphics.plot_2d_dataset(pyplot, ds)
pyplot.show()
```


## Result ##

![Shared Ancestry]({{site.baseurl}}/images/shared_ancestry.png)

