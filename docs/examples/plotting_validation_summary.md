---
noheader: true
layout: documentation
location: Examples
---

# Plotting a validation summary
-------------------------------

## Computing validation summary ##

You can compute the validation summary for a set of constant parameters
from the command line :

```sh
wot optimal_transport_validation --matrix matrix.txt --cell_days days.txt \
 --covariate covariate.txt --out validation_summary.txt
```

This will create a file called 'validation_summary.txt', containing
a tab-separated table with all the validation information.


## Reading DataFrame from a file ##

You can then create a `pandas.DataFrame` from the obtained validation summary file :

```python
import pandas

path_to_file = 'validation_summary.txt'
vs = pandas.read_table(path_to_file)

# Call the function defined below to plot
plot_validation_summary(vs, 'validation_summary.png')
```

## Plotting a DataFrame ##

If you have the validation summary available as `pandas.DataFrame`,
it is very easy to plot the mean and standard deviation of the distances
between all considered populations. The following function does that :

```python
import wot
import numpy
from matplotlib import pyplot

def plot_validation_summary(vs, filename):
    vs['time'] = (vs['interval_start'] + vs['interval_end']) / 2
    vs['type'] = vs['pair0'].astype(str).str[0]
    res = vs.groupby(['time', 'type'])['distance']\
        .agg([numpy.mean, numpy.std])

    legend = {
            'P': [ "#f000f0", "between real batches" ],
            'R': [ "#00f000", "between random and real" ],
            'I': [ "#f00000", "between interpolated and real"],
            'F': [ "#00f0f0", "between first and real"],
            'L': [ "#f0f000", "between last and real"],
            }

    pyplot.figure(figsize=(10, 10))
    pyplot.title("Validation of the OT model")
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
        pyplot.fill_between(t, m - s, m + s, color=legend[p][0] + "50")
    pyplot.savefig(filename)
```

## Getting the DataFrame from the OTModel ##

You can run the command-line tool to create a validation summary file, and then
load it into a DataFrame to plot it, but you could also go even further by
getting this validation summary directly from the OTModel.

This could for instance be used to test several configurations, generating a validation
plot for each of the configurations :

```python
import wot
import pandas

matrix_file = 'matrix.txt'
days_file = 'days.txt'
covariate_file = 'covariate.txt'

for l in [ 1, 10, 50, 80, 120 ]:
    ot_model = wot.initialize_ot_model(matrix_file, days_file,
        covariate=covariate_file, lambda1=l, growth_iters=1, fast=True)
    vs = wot.commands.compute_validation_summary(ot_model)
    c = ot_model.get_ot_config()
    figure_name = "validation_summary_e{}_l{}_l{}_g{}.png"\
            .format(c['epsilon'], c['lambda1'], c['lambda2'],
                    c['growth_iters'])
    plot_validation_summary(vs, figure_name)
```

This would take a very long time to compute though, watch out.

## Result ##

![Validation summary plot]({{site.baseurl}}/images/validation_summary.png)
