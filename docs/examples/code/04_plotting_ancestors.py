import numpy as np
import pandas as pd
import wot.graphics
from matplotlib import pyplot

# ------ Configuration variables -------
matrix_file = 'matrix.txt'
bg_color = "#80808080"
cell_sets_file = 'cell_sets.gmt'
target_cell_set = "Red blood cells"
target_timepoint = 7
destination_file = "ancestors.png"
embedding_file = 'embedding.csv'
# --------------------------------------

ds = wot.io.read_dataset(matrix_file)
ds.obs = ds.obs.join(pd.read_csv(embedding_file, index_col='id'))
tmap_model = wot.tmap.TransportMapModel.from_directory('tmaps')

transparent = lambda x: wot.graphics.hexstring_of_rgba((.08, .34, .59, x))


def color_cells(population):
    p = population.p
    if not np.isclose(max(p), 0):
        p = p / max(p)
    color = [transparent(x) for x in p]
    wot.set_cell_metadata(ds, 'color', color,
                          indices=tmap_model.cell_ids(population))


pyplot.figure(figsize=(5, 5))
pyplot.axis('off')
pyplot.scatter(ds.obs['x'], ds.obs['y'], c=bg_color,
               s=0.8, marker=',', edgecolors='none')

cell_sets = wot.io.read_sets(cell_sets_file, as_dict=True)
population = tmap_model.population_from_ids(
    cell_sets[target_cell_set],
    at_time=target_timepoint)[0]
color_cells(population)

while tmap_model.can_pull_back(population):
    population = tmap_model.pull_back(population)
    color_cells(population)

pyplot.scatter(ds.obs['x'], ds.obs['y'], c=ds.obs['color'].values,
               s=1, marker=',', edgecolors='none')

wot.graphics.legend_figure(pyplot,
                           [["#316DA2", "Ancestors of {}".format(target_cell_set)]])
pyplot.autoscale(enable=True, tight=True)
pyplot.tight_layout(pad=0)
pyplot.savefig(destination_file)
