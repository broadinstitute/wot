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


ot_validation_legend = {
    'P': ["#a6cee3", "between real batches"],
    'I': ["#1f78b4", "between interpolated and real"],
    'F': ["#b2df8a", "between first and real"],
    'L': ["#33a02c", "between last and real"],
    'R': ["#fb9a99", "between random (no growth) and real"],
    'Rg': ["#e31a1c", "between random (with growth) and real"]
}


def group_ot_validation_summary(df):
    df['time'] = (df['interval_start'] + df['interval_end']) / 2
    full_df = df[df['cv0'] == 'full']
    cv_df = df[df['cv0'] != 'full']
    cv_df = cv_df.copy()
    cv_df['type'] = cv_df['pair0'].astype(str).str.split('_').str.get(0)
    cv_agg = cv_df.groupby(['time', 'type'])['distance'].agg([np.mean, np.std])
    result = []
    result_columns = ['time', 'mean', 'std', 'type']
    # value from full batches, std from mean, except for P vs P
    for type, d in cv_agg.groupby('type'):
        t = np.asarray(d.index.get_level_values('time'))
        mean = d['mean'] if type is 'P' else full_df[full_df['pair0'] == type]['distance']
        std = d['std']
        result.append([t[0], np.asarray(mean)[0], np.asarray(std)[0], type])
    grouped_df = pd.DataFrame(result, columns=result_columns)
    return grouped_df


def plot_ot_validation_summary(grouped_df, filename, bandwidth=None):
    pyplot.figure(figsize=(10, 10))
    pyplot.title("OT Validation")
    pyplot.xlabel("time")
    pyplot.ylabel("distance")
    wot.graphics.legend_figure(pyplot, ot_validation_legend.values())

    for p, d in grouped_df.groupby('type'):
        if p not in ot_validation_legend.keys():
            print('skipping ' + str(p))
            continue
        t = np.asarray(d['time'])
        m = np.asarray(d['mean'])
        s = np.asarray(d['std'])
        if bandwidth is not None:
            x, m = kernel_smooth(t, m, 0, t[len(t) - 1], 1000, bandwidth)
            x, s = kernel_smooth(t, s, 0, t[len(t) - 1], 1000, bandwidth)
            t = x
        pyplot.plot(t, m, '-o', color=ot_validation_legend[p][0])
        pyplot.fill_between(t, m - s, m + s, color=ot_validation_legend[p][0] + "50")
    pyplot.savefig(filename)
