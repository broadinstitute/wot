---
noheader: true
layout: documentation
location: Examples
---

# Plotting cell sets
--------------------

## Generating the cell sets ##

If you don't already have a [cell sets file]({{ site.baseurl }}/cli_documentation#cellset_file),
you may want to create one from a [gene sets file]({{ site.baseurl }}/cli_documentation#geneset_file).

To do so, simply use the **cells_by_gene_set** wot tool :

```sh
wot cells_by_gene_set --matrix matrix.txt --gene_sets gene_sets.gmt \
 --out cell_sets.gmt --format gmt --quantile 0.01
```

## Coloring the cell sets ##

The next step is to assign a color to each cell set. The following function may
be used for that purpose :

```python
def color_cell_set(ds, cell_sets_path, cset_name, color):
    cs_groups = (wot.io.group_cell_sets(cell_sets_path, ds.row_meta))
    for t in cs_groups:
        for cset in cs_groups[t]:
            if cset['name'].startswith(cset_name + '_'):
                wot.set_cell_metadata(ds, 'color',
                        color, indices=cset['set'])

```

You can then quickly visualize your colored gene sets

```python
import numpy
import wot.io
import wot.graphics
from matplotlib import pyplot

ds = wot.io.read_dataset('matrix.txt')

bgcolor = "#80808020"
wot.io.incorporate_days_information_in_dataset(ds, 'days.txt')
wot.set_cell_metadata(ds, 'color', bgcolor)

color_cell_set(ds, 'cell_sets.gmt', 'tip1', '#800000')
color_cell_set(ds, 'cell_sets.gmt', 'tip2', '#008000')
color_cell_set(ds, 'cell_sets.gmt', 'tip3', '#000080')

pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds, x=0, y=1)
pyplot.show()
```
