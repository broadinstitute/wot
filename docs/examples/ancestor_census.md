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
 --out census --time 15
```

See the [CLI documentation]({{site.baseurl}}/cli_documentation#ancestor-census)
for more information about this command.

Note : You need to have a transport maps computed for this to work. Please refer to the [CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps) for instructions on how to compute transport maps.


### Plotting ###

This will create several census files, that you can then plot as follows :

```python
import wot.graphics
from matplotlib import pyplot

picture_title = "Ancestor census of cells from Tip 1 at day 10"

# Change this to the name of census file you want to plot
census_file = 'census_tip1.txt'

# Choose a color and a name to display for each cell set in the census
legend = [
        [ "#e00000", "Tip 1" ],
        [ "#00e000", "Tip 2" ],
        [ "#0000e0", "Tip 3" ]
        ]

ds = wot.io.read_dataset(census_file)

wot.graphics.legend_figure(pyplot, legend)
for i in range(ds.X.shape[1]):
    pyplot.plot(ds.X[:,i], color=legend[i][0])

pyplot.title(picture_title)
pyplot.xlabel("time")
pyplot.ylabel("Proportion of ancestors per cell set")
pyplot.show()
```

### Result ###

![ancestor census plot]({{site.baseurl}}/images/ancestor_census.png)
