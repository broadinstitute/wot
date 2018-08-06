# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
cell_sets_file = 'cell_sets.gmt'
target_cell_set = 'Red blood cells'
target_timepoint = 50
destination_file = "ancestor_census.png"
# --------------------------------------

import wot
import numpy
from matplotlib import pyplot

ot_model = wot.load_ot_model(matrix_file, days_file, 'tmaps')
cell_sets = wot.io.read_cell_sets(cell_sets_file)
cs_matrix = wot.io.read_gene_sets(cell_sets_file)
all_populations = ot_model.population_from_cell_sets(cell_sets,
        at_time=target_timepoint)
population = all_populations[target_cell_set]

timepoints, census = \
        wot.commands.compute_ancestor_census(ot_model, cs_matrix, population)

pyplot.figure(figsize=(5, 5))

cset_names = list(cell_sets.keys())
for i in range(census.shape[1]):
    pyplot.plot(timepoints, census[:,i], label=cset_names[i])

pyplot.xlabel("Time")
pyplot.ylabel("Proportion of ancestors")
pyplot.title("Ancestor census of {} from time {}"\
        .format(target_cell_set, target_timepoint))
pyplot.legend()
pyplot.savefig(destination_file)
