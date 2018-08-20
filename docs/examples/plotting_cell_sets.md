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
 --out cell_sets.gmt --format gmt --quantile 0.99
```

## Coloring the cell sets ##

The next step is to assign a color to each cell set. The following function may
be used for that purpose :

```python
# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
cell_sets_file = 'cell_sets.gmt'
bg_color = "#80808020"
cell_sets_to_color = [
        [ 'red',    'Red blood cells' ],
        [ 'blue',   'Granulocytes' ],
        [ 'green',  'Lymphocytes' ],
        [ 'purple', 'Myeloid stem cells' ],
        [ 'black',  'Stem cells' ],
        ]
gene_x_plot = 0
gene_y_plot = 1
destination_file = "cell_sets.png"
# --------------------------------------

import wot
from matplotlib import pyplot

ds = wot.io.read_dataset(matrix_file)
wot.io.incorporate_days_information_in_dataset(ds, days_file)

cell_sets = wot.io.read_cell_sets(cell_sets_file)

wot.set_cell_metadata(ds, 'color', bg_color)

for color, cset_name in cell_sets_to_color:
    wot.set_cell_metadata(ds, 'color', color,
            indices=cell_sets[cset_name])

pyplot.figure(figsize=(5,5))
pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds,
        x=gene_x_plot, y=gene_y_plot)
wot.graphics.legend_figure(pyplot, cell_sets_to_color)
pyplot.autoscale(enable=True, tight=True)
pyplot.tight_layout(pad=0)
pyplot.savefig(destination_file)
```

## Result ##

![cell sets plot]({{site.baseurl}}/images/notebook_cell_sets.png)
