---
noheader: true
layout: documentation
location: Examples
---

# Plotting an ancestor census
-----------------------------

```python
import wot.io
import wot.graphics
from matplotlib import pyplot

ds = wot.io.read_dataset('census_tip1_100.0.txt')

legend = [
        [ "#e00000", "Tip 1" ],
        [ "#00e000", "Tip 2" ],
        [ "#0000e0", "Tip 3" ]
        ]
wot.graphics.legend_figure(pyplot, legend)
for i in range(ds.x.shape[1]):
    pyplot.plot(ds.x[:,i], color=legend[i][0])

pyplot.title("Ancestor census of cells from Tip 1 at day 100")
pyplot.xlabel("time")
pyplot.ylabel("Proportion of ancestors per cell set")
pyplot.show()
```
