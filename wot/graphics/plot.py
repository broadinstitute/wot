# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib import pyplot

import wot.graphics


def __make_figure(y=1, x=1, projection=None):
    pyplot.clf()
    return pyplot.subplots(y, x, figsize=(8 * x, 6 * y), projection=None)


def plot_2d_dataset(figure, dataset, x=0, y=1, title=None):
    colors = "#808080"
    if 'color' in dataset.row_meta.columns:
        colors = dataset.row_meta['color'].values
    figure.scatter(dataset.x[:, x], dataset.x[:, y], c=colors,
                   s=.2, marker=',', edgecolors='none')
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


def group_ot_validation_summary(df, filename):
    df['time'] = (df['interval_start'] + df['interval_end']) / 2
    df['type'] = df['pair0'].astype(str).str[0]
    res = df.groupby(['time', 'type'])['distance'].agg([np.mean, np.std])
    is_first = True
    legend = {
        'P': ["#984ea3", "between real batches"],
        'R': ["#4daf4a", "between random and real"],
        'I': ["#e41a1c", "between interpolated and real"],
        'F': ["#377eb8", "between first and real"],
        'L': ["#ff7f00", "between last and real"],
    }
    with open(filename, 'w') as f:
        for p, d in res.groupby('type'):
            if p not in legend.keys():
                continue
            t = np.asarray(d.index.get_level_values('time'))
            m = np.asarray(d['mean'])
            s = np.asarray(d['std'])
            pd.DataFrame(data={"time": t, "mean": m, "std": s, "type": legend[p][1]}).to_csv(f, sep="\t",
                                                                                             header=is_first,
                                                                                             index=False)
    return res


def plot_ot_validation_summary(res, filename, bandwidth=None):
    legend = {
        'P': ["#984ea3", "between real batches"],
        'R': ["#4daf4a", "between random and real"],
        'I': ["#e41a1c", "between interpolated and real"],
        'F': ["#377eb8", "between first and real"],
        'L': ["#ff7f00", "between last and real"],
    }

    pyplot.figure(figsize=(10, 10))
    pyplot.title("OT Validation")
    pyplot.xlabel("time")
    pyplot.ylabel("distance")
    wot.graphics.legend_figure(pyplot, legend.values())

    for p, d in res.groupby('type'):
        if p not in legend.keys():
            continue
        t = np.asarray(d.index.get_level_values('time'))
        m = np.asarray(d['mean'])
        s = np.asarray(d['std'])
        if bandwidth is not None:
            x, m = kernel_smooth(t, m, 0, t[len(t) - 1], 1000, bandwidth)
            x, s = kernel_smooth(t, s, 0, t[len(t) - 1], 1000, bandwidth)
            t = x
        pyplot.plot(t, m, '-o', color=legend[p][0])
        pyplot.fill_between(t, m - s, m + s, color=legend[p][0] + "50")
    pyplot.savefig(filename)
