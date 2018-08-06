# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
covariate_file = 'covariate.txt'
destination_file = 'validation_summary.png'
# --------------------------------------

import wot
import numpy
from matplotlib import pyplot

ot_model = wot.initialize_ot_model(matrix_file, days_file,
    covariate=covariate_file, growth_iters=1, tmap_prefix='val')
vs = wot.commands.compute_validation_summary(ot_model)
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
pyplot.savefig(destination_file)
