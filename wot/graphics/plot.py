# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import patches
from matplotlib import pyplot

import wot.graphics


def __make_figure(y=1, x=1, projection=None):
    pyplot.clf()
    return pyplot.subplots(y, x, figsize=(8 * x, 6 * y), projection=None)


def plot_2d_dataset(figure, dataset, x=0, y=1, s=0.2, colors="#808080", title=None):
    figure.scatter(dataset.X[:, x], dataset.X[:, y], c=colors,
                   s=s, marker=',', edgecolors='none')
    if title is not None:
        figure.title.set_text(title)


def legend_figure(figure, legend_list, loc=0):
    patch_list = [patches.Patch(color=c, label=l) for c, l in legend_list]
    figure.legend(handles=patch_list, loc=loc)


def interpolate(x, xi, yi, sigma):
    val = x - xi
    val *= -val
    diff = val
    sigma2 = 2 * sigma ** 2
    w = np.exp(diff / sigma2)
    fx = (yi * w).sum()
    return fx / w.sum()


def kernel_smooth(xi, yi, start, stop, steps, sigma):
    xlist = np.linspace(start, stop, steps)
    fhat = np.zeros(len(xlist))
    for i in range(len(xlist)):
        fhat[i] = interpolate(xlist[i], xi, yi, sigma)
    return xlist, fhat


ot_validation_legend = {
    'P': ["#e41a1c", "between real batches"],
    'I': ["#377eb8", "between interpolated and real"],
    'F': ["#4daf4a", "between first and real"],
    'L': ["#984ea3", "between last and real"],
    'R': ["#ff7f00", "between random (no growth) and real"],
    'Rg': ["#ffff33", "between random (with growth) and real"],
    'A': ["#bdbdbd", "between first and last"]
}


def plot_ot_validation_summary(df, filename, bandwidth=None):
    df = df.reset_index()
    pyplot.figure(figsize=(10, 10))
    pyplot.title("OT Validation")
    pyplot.xlabel("time")
    pyplot.ylabel("distance")
    legend = {}

    for p, d in df.groupby('name'):
        if p not in ot_validation_legend.keys():
            continue
        t = np.asarray(d['interval_mid'])
        m = np.asarray(d['mean'])
        s = np.asarray(d['std'])
        legend[p] = ot_validation_legend[p]
        if bandwidth is not None:
            x, m = kernel_smooth(t, m, 0, t[len(t) - 1], 1000, bandwidth)
            x, s = kernel_smooth(t, s, 0, t[len(t) - 1], 1000, bandwidth)
            t = x
        pyplot.plot(t, m, '-o', color=ot_validation_legend[p][0])
        pyplot.fill_between(t, m - s, m + s, color=ot_validation_legend[p][0], alpha=0.2)
    wot.graphics.legend_figure(pyplot, legend.values())
    pyplot.savefig(filename)
