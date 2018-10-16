import numpy as np
from matplotlib import pyplot

import wot.commands
import wot.graphics
# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
covariate_file = 'covariate.txt'
destination_file = 'validation_summary.png'
# --------------------------------------

ot_model = wot.ot.initialize_ot_model(matrix_file, days_file,
                                      covariate=covariate_file, growth_iters=1, tmap_prefix='val')
vs = wot.commands.compute_validation_summary(ot_model)
vs['time'] = (vs['interval_start'] + vs['interval_end']) / 2
vs['type'] = vs['pair0'].astype(str).str[0]
res = vs.groupby(['time', 'type'])['distance'] \
    .agg([np.mean, np.std])

legend = {
    'P': ["#f000f0", "between real batches"],
    'R': ["#00f000", "between random and real"],
    'I': ["#f00000", "between interpolated and real"],
    'F': ["#00f0f0", "between first and real"],
    'L': ["#f0f000", "between last and real"],
}

pyplot.figure(figsize=(10, 10))
pyplot.title("Validation of the OT model")
pyplot.xlabel("time")
pyplot.ylabel("distance")
wot.graphics.legend_figure(pyplot, legend.values())
for p, d in res.groupby('type'):
    if p not in legend.keys():
        continue
    t = np.asarray(d.index.get_level_values('time'))
    m = np.asarray(d['mean'])
    s = np.asarray(d['std'])
    pyplot.plot(t, m, '-o', color=legend[p][0])
    pyplot.fill_between(t, m - s, m + s, color=legend[p][0] + "50")
pyplot.savefig(destination_file)
