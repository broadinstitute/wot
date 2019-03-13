---
noheader: true
layout: documentation
location: Examples
---

# Plotting ancestors
--------------------


## Generating the cell sets ##

If you don't already have a [cell sets file]({{ site.baseurl }}/cli_documentation#cellset_file),
you may want to create one from a [gene sets file]({{ site.baseurl }}/cli_documentation#geneset_file).

To do so, simply use the **cells_by_gene_set** wot tool :

```sh
wot cells_by_gene_set --matrix matrix.txt --gene_sets gene_sets.gmt \
 --out cell_sets.gmt --format gmt --quantile 0.99
```

## Computing and plotting ancestors ##

Note : You need to have a transport maps computed for this to work. Please refer to the [CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps) for instructions on how to compute transport maps.


```python
{% include_relative code/04_plotting_ancestors.py %}

```

## Result ##

![ancestors plot]({{site.baseurl}}/images/ancestors.png)
