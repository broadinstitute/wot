# -*- coding: utf-8 -*-

from matplotlib import patches
from matplotlib import pyplot


def __make_figure(y = 1, x = 1, projection=None):
    pyplot.clf()
    return pyplot.subplots(y, x, figsize=(8 * x, 6 * y), projection=None)

def plot_2d_dataset(figure, dataset, x = 0, y = 1, title=None):
    colors = "#808080"
    if 'color' in dataset.row_meta.columns:
        colors = dataset.row_meta['color'].values
    figure.scatter(dataset.x[:,x], dataset.x[:,y], c=colors,
            s=.2, marker=',', edgecolors='none')
    if title is not None:
        figure.title.set_text(title)

def legend_figure(figure, legend_list, loc=0):
    patch_list = [ patches.Patch(color=c, label=l) for c, l in legend_list ]
    figure.legend(handles = patch_list, loc = loc)
