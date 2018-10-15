import numpy
from matplotlib import pyplot

import wot

# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
cell_sets_file = 'cell_sets.gmt'
target_cell_set = 'Red blood cells'
target_timepoint = 50
skip_first_n_genes = 2
destination_file = "trajectory_trends.png"
# --------------------------------------

ds = wot.io.read_dataset(matrix_file)
tmap_model = wot.tmap.TransportMapModel.from_directory('tmaps')
cell_sets = wot.io.read_sets(cell_sets_file, as_dict=True)
all_populations = tmap_model.population_from_cell_sets(cell_sets,
                                                       at_time=target_timepoint)
population = all_populations[target_cell_set]

# timepoints, means, variances = tmap_model.compute_trajectory_trends(ds, population)
trajectory_ds = tmap_model.compute_trajectories({target_cell_set: all_populations[target_cell_set]})
results = wot.ot.compute_trajectory_trends_from_trajectory(trajectory_ds, ds)
means, variances = results[0]
timepoints = means.row_meta.index.values
pyplot.figure(figsize=(5, 5))
means = means.x
stds = numpy.sqrt(variances.x)
genes = ds.col_meta.index
for i in range(skip_first_n_genes, means.shape[1]):
    pyplot.plot(timepoints, means[:, i], label=genes[i])
    pyplot.fill_between(timepoints, means[:, i] - stds[:, i],
                        means[:, i] + stds[:, i], alpha=.5)

pyplot.xlabel("Time")
pyplot.ylabel("Gene expression")
pyplot.title("Trajectory trend of {} from time {}" \
             .format(target_cell_set, target_timepoint))
pyplot.legend()
pyplot.savefig(destination_file)
