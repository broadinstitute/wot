# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
cell_sets_file = 'cell_sets.gmt'
target_cell_set = 'Red blood cells'
target_timepoint = 50
skip_first_n_genes = 2
destination_file = "trajectory_trends.png"
# --------------------------------------

import wot
import numpy
from matplotlib import pyplot

ot_model = wot.load_ot_model(matrix_file, days_file, 'tmaps')
cell_sets = wot.io.read_cell_sets(cell_sets_file)
all_populations = ot_model.population_from_cell_sets(cell_sets,
        at_time=target_timepoint)
population = all_populations[target_cell_set]

timepoints, trends, variances = \
        wot.commands.compute_trajectory_trends(ot_model, population)

pyplot.figure(figsize=(10, 10))

stds = numpy.sqrt(variances)
genes = ot_model.matrix.col_meta.index
for i in range(skip_first_n_genes, trends.shape[1]):
    pyplot.plot(timepoints, trends[:,i], label=genes[i])
    pyplot.fill_between(timepoints, trends[:,i] - stds[:,i],
            trends[:,i] + stds[:,i], alpha=.5)

pyplot.xlabel("Time")
pyplot.ylabel("Gene expression")
pyplot.title("Trajectory trend of {} from time {}"\
        .format(target_cell_set, target_timepoint))
pyplot.legend()
pyplot.savefig(destination_file)
