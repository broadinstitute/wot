# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
bg_color = "#80808080"
gene_x_plot = 0
gene_y_plot = 1
cell_sets_file = 'cell_sets.gmt'
target_cell_set = "Red blood cells"
target_timepoint = 50
destination_file = "ancestors.png"
# --------------------------------------

import numpy
import wot
from matplotlib import pyplot

ds = wot.io.read_dataset(matrix_file)

ot_model = wot.load_ot_model(matrix_file, days_file, 'tmaps')

transparent = lambda x : wot.graphics.hexstring_of_rgba((.08, .34, .59, x))
def color_cells(population):
    p = population.p
    if not numpy.isclose(max(p), 0):
        p = p / max(p)
    color = [ transparent(x) for x in p ]
    wot.set_cell_metadata(ds, 'color', color,
            indices = ot_model.cell_ids(population))


pyplot.figure(figsize=(5,5))
pyplot.axis('off')
wot.set_cell_metadata(ds, 'color', bg_color)
wot.graphics.plot_2d_dataset(pyplot, ds, x=gene_x_plot, y=gene_y_plot)

cell_sets = wot.io.read_cell_sets(cell_sets_file)
population = ot_model.population_from_ids(
        cell_sets[target_cell_set],
        at_time=target_timepoint)
color_cells(population)

while ot_model.can_pull_back(population):
    population = ot_model.pull_back(population)
    color_cells(population)

wot.graphics.plot_2d_dataset(pyplot, ds,
        x=gene_x_plot, y=gene_y_plot)
wot.graphics.legend_figure(pyplot,
        [["#316DA2", "Ancestors of {}".format(target_cell_set)]])
pyplot.autoscale(enable=True, tight=True)
pyplot.tight_layout(pad=0)
pyplot.savefig(destination_file)
