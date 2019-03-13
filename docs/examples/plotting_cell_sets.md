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

The next step is to assign a color to each cell set. The following code may
be used for that purpose :

```python
{% include_relative code/03_plotting_cell_sets.py %}

```

## Result ##

![cell sets plot]({{site.baseurl}}/images/cell_sets.png)
