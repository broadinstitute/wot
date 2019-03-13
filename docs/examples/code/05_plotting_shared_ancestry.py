import numpy as np
import pandas as pd
import wot.graphics
from matplotlib import pyplot

# ------ Configuration variables -------
matrix_file = 'matrix.txt'
bg_color = "#80808050"
gene_x_plot = 0
gene_y_plot = 1
cell_set_1 = "Red blood cells"
cell_set_2 = "Granulocytes"
cell_sets_file = 'cell_sets.gmt'
target_timepoint = 7
destination_file = "shared_ancestry.png"
embedding_file = 'embedding.csv'
# --------------------------------------


tmap_model = wot.tmap.TransportMapModel.from_directory('tmaps')
cell_sets = wot.io.read_sets(cell_sets_file, as_dict=True)
populations = tmap_model.population_from_cell_sets(cell_sets,
                                                   at_time=target_timepoint)

trajectories = tmap_model.compute_trajectories(populations)
ds = wot.io.read_dataset(matrix_file)
ds.obs = ds.obs.join(pd.read_csv(embedding_file, index_col='id'))

pyplot.figure(figsize=(5, 5))
pyplot.axis('off')
pyplot.scatter(ds.obs['x'], ds.obs['y'], c=bg_color,
               s=0.8, marker=',', edgecolors='none')

probabilities = trajectories.X[:, trajectories.var.index.get_indexer_for([cell_set_1, cell_set_2])]
alphas = np.amax(probabilities, axis=1)
t = np.log((probabilities.T[0] + 1e-9) / (probabilities.T[1] + 1e-9))
t = np.clip(t / 8 + .5, 0, 1)
alphas = alphas / max(alphas)
colors = [wot.graphics.hexstring_of_rgba([t[i], 0, 1 - t[i], alphas[i]])
          for i in range(len(t))]
ds.obs.loc[:, 'color'] = colors
pyplot.scatter(ds.obs['x'], ds.obs['y'], c=colors,
               s=1, marker=',', edgecolors='none')
wot.graphics.legend_figure(pyplot,
                           [["#A00000", "Ancestors of " + cell_set_1],
                            ["#0000A0", "Ancestors of " + cell_set_2]],
                           loc=3)
pyplot.autoscale(enable=True, tight=True)
pyplot.tight_layout(pad=0)
pyplot.savefig(destination_file)
