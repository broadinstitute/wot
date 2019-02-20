---
noheader: true
layout: documentation
location: Examples
---

# Plotting an ancestor census
-----------------------------

### Computing ###

You can compute an ancestor census for your cell sets with:

```sh
wot census --tmap tmaps \
 --cell_set cell_sets.gmt \
 --out census --time 7
```

See the [CLI documentation]({{site.baseurl}}/cli_documentation#ancestor-census)
for more information about this command.

Note : You need to have a transport maps computed for this to work. Please refer to the [CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps) for instructions on how to compute transport maps.


### Plotting ###

This will create several census files, that you can then plot as follows :

```python
{% include_relative code/07_plotting_ancestor_census.py %}
```

### Result ###

![ancestor census plot]({{site.baseurl}}/images/ancestor_census.png)
