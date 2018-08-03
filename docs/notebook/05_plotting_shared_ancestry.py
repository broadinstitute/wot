# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
bg_color = "#80808080"
gene_x_plot = 0
gene_y_plot = 1
cell_set_1 = "Red blood cells"
cell_set_2 = "Granulocytes"
cell_sets_file = 'cell_sets.gmt'
target_timepoint = 50
destination_file = "shared_ancestry.png"
# --------------------------------------

import wot
import numpy
import pandas
from matplotlib import pyplot

ot_model = wot.load_ot_model(matrix_file, days_file, 'tmaps')
cell_sets = wot.io.read_cell_sets(cell_sets_file)
populations = ot_model.population_from_cell_sets(cell_sets,
        at_time=target_timepoint)

traj = wot.commands.compute_trajectories(ot_model, *populations.values())
traj.columns = populations.keys()
ds = wot.io.read_dataset(matrix_file)
wot.set_cell_metadata(ds, 'color', '#80808050')

pyplot.figure(figsize=(5, 5))
pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds)

probabilities = traj.loc[:, [cell_set_1, cell_set_2]].values
alphas = numpy.amax(probabilities, axis=1)
t = numpy.log((probabilities.T[0] + 1e-9) / (probabilities.T[1] + 1e-9))
t = numpy.clip(t / 8 + .5, 0, 1)
alphas = alphas / max(alphas)
colors = [ wot.graphics.hexstring_of_rgba([t[i], 0, 1-t[i], alphas[i]])
           for i in range(len(t)) ]
ds.row_meta.loc[:, 'color']  = colors

wot.graphics.plot_2d_dataset(pyplot, ds)
wot.graphics.legend_figure(pyplot,
        [["#A00000", "Ancestors of " + cell_set_1],
         ["#0000A0", "Ancestors of " + cell_set_2]],
        loc=3 )
pyplot.autoscale(enable=True, tight=True)
pyplot.tight_layout(pad=0)
pyplot.savefig(destination_file)
