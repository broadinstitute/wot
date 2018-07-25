# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
gene_sets_file = 'gene_sets.gmt'
quantile_for_cell_sets = .88
cell_sets_file = 'cell_sets.gmt'
bg_color = "#80808020"
cell_sets_to_color = [
        [ 'red',    'Red blood cells' ],
        [ 'blue',   'Granulocytes' ],
        [ 'green',  'Lymphocytes' ],
        [ 'purple', 'Myeloid stem cells' ],
        [ 'black',  'Stem cells' ],
        ]
gene_x_plot = 0
gene_y_plot = 1
destination_file = "cell_sets.png"
# --------------------------------------

import wot
from matplotlib import pyplot

ds = wot.io.read_dataset(matrix_file)
wot.io.incorporate_days_information_in_dataset(ds, days_file)

# Compute the cell sets for the given quantile

gene_sets = wot.io.read_gene_sets(gene_sets_file, wot.cell_names(ds))
cell_sets = wot.commands.get_cells_in_gene_sets(gene_sets, ds,
        quantile=quantile_for_cell_sets)
wot.io.write_gene_sets(cell_sets, cell_sets_file, "gmt")

# Plot the cell sets

wot.set_cell_metadata(ds, 'color', bg_color)

for color, cset_name in cell_sets_to_color:
    wot.set_cell_metadata(ds, 'color', color,
            indices=cell_sets[cset_name])

pyplot.figure(figsize=(5,5))
pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds,
        x=gene_x_plot, y=gene_y_plot)
wot.graphics.legend_figure(pyplot, cell_sets_to_color)
pyplot.autoscale(enable=True, tight=True)
pyplot.tight_layout(pad=0)
pyplot.savefig(destination_file)
