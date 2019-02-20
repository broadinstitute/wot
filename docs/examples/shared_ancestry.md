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
 --cell_set cell_sets.gmt --tmap tmaps --time 7
```

Note : You need to have a transport maps computed for this to work. Please refer to the [CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps) for instructions on how to compute transport maps.

## Compute shared ancestry ##

```python
{% include_relative code/05_plotting_shared_ancestry.py %}

```


## Result ##

![Shared Ancestry]({{site.baseurl}}/images/shared_ancestry.png)

