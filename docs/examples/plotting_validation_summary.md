---
noheader: true
layout: documentation
location: Examples
---

# Plotting a validation summary
-------------------------------

```python
import numpy import pandas
import wot.graphics
from matplotlib import pyplot

t_interp = .5

vs = pandas.read_table("val_tmaps_validation_summary.txt", engine='python', sep='\t')
vs['time'] = vs['interval_start'] * t_interp + vs['interval_end'] * (1 - t_interp)
vs['type'] = (
        vs['pair0'].astype(str).str[0] +
        vs['pair1'].astype(str).str[0]
        ).apply(lambda s: ''.join(sorted(s)))
# The type is one of PP, IP, IR, PR
res = vs.groupby(['time', 'type'])['distance'].agg([numpy.mean, numpy.std])

legend = {
        'PP': [ "#f000f0", "between real batches" ],
        'PR': [ "#f0f000", "between random and real" ],
        'IP': [ "#00f0f0", "between interpolated and real"]
        }

pyplot.title("Average distance and standard deviation over time")
pyplot.xlabel("time")
pyplot.ylabel("distance")
wot.graphics.legend_figure(pyplot, legend.values())
for p, d in res.groupby('type'):
    if p not in legend.keys():
        continue
    t = numpy.asarray(d.index.get_level_values('time'))
    m = numpy.asarray(d['mean'])
    s = numpy.asarray(d['std'])
    pyplot.plot(t, m, '-o', color=legend[p][0])
    pyplot.fill_between(t, m - s, m + s, color=legend[p][0] + "80")
pyplot.show()
```
