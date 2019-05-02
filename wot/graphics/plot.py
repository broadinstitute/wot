# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import patches
from matplotlib import pyplot as plt

import wot.graphics


def __make_figure(y=1, x=1, projection=None):
    plt.clf()
    return plt.subplots(y, x, figsize=(8 * x, 6 * y), projection=None)


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
    'A': ["#bdbdbd", "between first and last"],
    'I1': ["#a6cee3", "between first and interpolated"],
    'I2': ["#fb9a99", "between last and interpolated"]
}


def plot_ot_validation_ratio(df, filename):
    # (interpolated - real) / (null - real)
    df = df.reset_index()
    df = df.sort_values('interval_mid')

    interpolated_df = df[df['name'] == 'I']
    null_growth = df[df['name'] == 'Rg']
    null_no_growth = df[df['name'] == 'R']
    if (interpolated_df['interval_mid'].values - null_growth['interval_mid'].values).sum() != 0:
        raise ValueError('Timepoints are not aligned')
    if (interpolated_df['interval_mid'].values - null_no_growth['interval_mid'].values).sum() != 0:
        raise ValueError('Timepoints are not aligned')

    plt.figure(figsize=(10, 10))
    with_growth_score = (interpolated_df['mean'].values / null_growth['mean'].values).sum()
    no_growth_score = (interpolated_df['mean'].values / null_no_growth['mean'].values).sum()
    plt.title(
        "OT Validation: \u03A3(interpolated - real)/(null - real), with growth={:.2f}, no growth={:.2f}".format(
            with_growth_score, no_growth_score))
    plt.xlabel("time")
    plt.ylabel("ratio")

    plt.plot(interpolated_df['interval_mid'], interpolated_df['mean'].values / null_growth['mean'].values,
             label='with growth')
    plt.plot(interpolated_df['interval_mid'], interpolated_df['mean'].values / null_no_growth['mean'].values,
             label='no growth')
    plt.legend()
    plt.savefig(filename)


def plot_ot_validation_summary_stats(df, bandwidth=None):
    df = df.reset_index()
    plt.figure(figsize=(10, 10))
    plt.title("OT Validation")
    plt.xlabel("time")
    plt.ylabel("distance")
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
        plt.plot(t, m, '-o', color=ot_validation_legend[p][0])
        plt.fill_between(t, m - s, m + s, color=ot_validation_legend[p][0], alpha=0.2)
    wot.graphics.legend_figure(plt, legend.values())
