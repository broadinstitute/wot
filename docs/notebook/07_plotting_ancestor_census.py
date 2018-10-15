from matplotlib import pyplot

import wot

# ------ Configuration variables -------
matrix_file = 'matrix.txt'
cell_sets_file = 'cell_sets.gmt'
target_cell_set = 'Red blood cells'
target_timepoint = 50
destination_file = "ancestor_census.png"
# --------------------------------------


tmap_model = wot.tmap.TransportMapModel.from_directory('tmaps')

cs_matrix = wot.io.read_sets(cell_sets_file)
cell_sets = wot.io.convert_binary_dataset_to_dict(cs_matrix)
all_populations = tmap_model.population_from_cell_sets(cell_sets,
                                                       at_time=target_timepoint)
population = all_populations[target_cell_set]

timepoints, census = tmap_model.compute_ancestor_census(cs_matrix, population)

pyplot.figure(figsize=(5, 5))

cset_names = list(cell_sets.keys())
for i in range(census.shape[2]):
    pyplot.plot(timepoints, census[:, :, i].T, label=cset_names[i])

pyplot.xlabel("Time")
pyplot.ylabel("Proportion of ancestors")
pyplot.title("Ancestor census of {} from time {}" \
             .format(target_cell_set, target_timepoint))
pyplot.legend()
pyplot.savefig(destination_file)
