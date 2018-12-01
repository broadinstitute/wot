import numexpr as ne
import numpy as np
import scipy
import wot.grn
from scipy.stats import entropy
from sklearn.cluster import SpectralClustering


def dum(x):
    return x


def dumg(x):
    return 1


def compose_transports(Lineage, TP, lag):
    ComposedLineage = []
    for i in range(len(Lineage)):
        if len(Lineage[i]) > 0:
            j = i + 1
            while (j < len(TP)) and (TP[j] > TP[j - 1]) and (TP[j] < TP[i] + lag):
                j += 1
            if j == len(TP):
                break
            if TP[j] < TP[i] + lag:
                ComposedLineage.append([])
            else:
                composed_lineage = Lineage[i]
                for k in range(i + 1, j):
                    composed_lineage = composed_lineage.dot(Lineage[k])
                ComposedLineage.append(composed_lineage)
        else:
            ComposedLineage.append([])
    return ComposedLineage


try:
    from gslrandom import PyRNG, multinomial


    def coupling_random_sample(n, l, s, threads):
        P = np.empty((s, len(l)), dtype=np.uint32)
        l_tile = np.tile(l, (s, 1))
        rngs = [PyRNG(np.random.randint(2 ** 16)) for _ in range(threads)]
        return multinomial(rngs, n, l_tile, P)
except:
    # print('gslrandom not available. Please install for faster random number generation.')
    def coupling_random_sample(n, l, s, threads):
        return np.random.multinomial(n[0], l, size=s)


def coupling_sampler(Lineage, nf=1e-3, s=1, threads=1, nmin=10):
    Pairs = [[] for _ in range(s)]
    for lineage in Lineage:
        if len(lineage) > 0:
            l = lineage / lineage.sum(0)
            l = l / l.sum()
            l = l.flatten()
            sd = np.exp(entropy(l))
            n = max(nmin, int(sd * nf)) * np.ones(s, dtype=np.uint32)
            # P = np.ones((s,len(l)),dtype=np.uint32)
            P = coupling_random_sample(n, l, s, threads)
            for i in range(s):
                pairs = np.nonzero(P[i].reshape(lineage.shape))
                Pairs[i].append(pairs)
        else:
            for i in range(s):
                Pairs[i].append([])
    return Pairs


def get_expression_pairs(Pairs, Lineage, Xg, Xr, TP, lag, differences=True):
    Xgs = [[]]
    Xrs = []
    for i in range(len(Lineage)):
        if len(Lineage[i]) > 0:
            j = i + 1
            while (j < len(TP)) and (TP[j] > TP[j - 1]) and (TP[j] < TP[i] + lag):
                j += 1
            if j == len(TP):
                break
            if TP[j] < TP[i] + lag:
                xgp = []
                xr = []
            else:
                pairs = Pairs[i]
                if differences:
                    a = Xg[j][pairs[1]]
                    b = Xg[i][pairs[0]]
                    xgp = ne.evaluate("a-b")
                else:
                    xgp = Xg[j][pairs[1]]
                # interpolate
                delta_t = TP[j] - TP[i]
                x1 = Xr[i][pairs[0]]
                x2 = Xr[j][pairs[1]]
                xr = x1 * lag / delta_t + x2 * (1. - lag / delta_t)
        else:
            xgp = []
            xr = []
        Xgs.append(xgp)
        Xrs.append(xr)
    Xrs.append([])
    return Xgs, Xrs


def initialize_modules(Xg, N, threads=1, nn=150):
    XG = np.vstack(Xg)
    XG = XG[np.random.choice(range(XG.shape[0]), int(XG.shape[0] / 20), replace=False)]
    SC = SpectralClustering(n_clusters=N, n_neighbors=nn, affinity='nearest_neighbors', n_jobs=threads)
    labels = SC.fit_predict(XG.T)
    U = np.zeros((N, XG.shape[1]))
    for i, c in enumerate(set(labels)):
        U[i, np.where(labels == c)[0]] = 1. / (labels == c).sum() ** .5
    return U


def update_regulation(Lineage, Xg, Xr, TP, lag, Z=[], U=[], num_modules=50, lda_z1=100., lda_z2=100., lda_u=10.,
                      epochs=25, sample_fraction=0.01, inner_iters=1, threads=1, k=1, b=6, y0=.01, x0=0,
                      differences=True, frequent_fa=False, fa_update_itrs=10, epoch_block_size=100, savepath=None):
    # lineage from optimal transport
    # Xg is a list of observed gene expression matrices at each time point
    # Xr is a list of regulatory levels at each time point (lagged by 1 from Xg)
    model = wot.grn.SparseOptimization(threads)
    if len(U) == 0:
        model.U = np.random.random((num_modules, Xg[0].shape[1]))
        model.U = (model.U.T / np.linalg.norm(model.U, axis=1)).T
    else:
        model.U = U
    if len(Z) > 0:
        model.Z = Z
    if differences:
        model.fa = dum
        model.fa_grad = dumg
        model.k, model.b, model.y0, model.x0 = k, b, y0, x0
    else:
        model.set_fa(k, b, y0, x0)
    for epoch_block in range(0, epochs, epoch_block_size):
        nepochs = min(epoch_block_size, epochs - epoch_block)
        Pairs = coupling_sampler(Lineage, nf=sample_fraction, s=nepochs, threads=threads)
        for epoch in range(nepochs):
            Xgs, Xrs = get_expression_pairs(Pairs[epoch], Lineage, Xg, Xr, TP, lag, differences=differences)
            model.withCoupling = False
            model.Lineage, model.Xg, model.Xr = Lineage, Xgs, Xrs
            for innerItr in range(inner_iters):
                model.lda1, model.lda2 = lda_z1, lda_z2
                model.z_shape = (Xr[0].shape[1], num_modules)
                model.update_Z(maxItr=1, with_prints=False, fa_update_freq=10 ** 6, forceNorm=True)
                model.lda1, model.lda2 = lda_u, 0
                model.update_U(maxItr=1, with_prints=False, fa_update_freq=10 ** 6)
            if ((epoch + epoch_block) % 5 == 0) or ((epoch + epoch_block) == epochs - 1):
                if frequent_fa:
                    model.update_fa(fmin_itrs=fa_update_itrs)
                model.print_performance(model.U, epoch + epoch_block)
        if savepath != None:
            np.save('%s/Z.%d.npy' % (savepath, epoch + epoch_block), model.Z)
            np.save('%s/U.%d.npy' % (savepath, epoch + epoch_block), model.U)
            np.save('%s/kbyx.%d.npy' % (savepath, epoch + epoch_block), (model.k, model.b, model.y0, model.x0))
    model.Xg, model.Xr = Xg, Xr
    return model.Z, model.U, model.get_all_Xhat(model.U), model.k, model.b, model.y0, model.x0


def main(argsv):
    import argparse
    import wot.io
    import pandas as pd
    import os
    parser = argparse.ArgumentParser(
        description='Gene Regulatory Networks')

    parser.add_argument('--tmap', help=wot.commands.TMAP_HELP, required=True, action='append')
    parser.add_argument('--tf',
                        help='File with one gene id per line to assign transcription factors', required=True)
    parser.add_argument('--gene_filter',
                        help='File with one gene id per line to include')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--time_lag', help='Time lag', type=float, required=True)
    parser.add_argument('--nmodules', help='Number of gene expression modules', type=int, default=50)
    parser.add_argument('--z1', help='Z L1 regularlization parameter', type=float, default=2)
    parser.add_argument('--z2', help='Z L2 regularlization parameter', type=float, default=0.5)
    parser.add_argument('--u2', help='U L2 regularlization parameter', type=float, default=1.5)

    parser.add_argument('--threads',
                        help='Number of threads to use', type=int)
    parser.add_argument('--epochs',
                        help='Number of epochs', type=int, default=10000)

    parser.add_argument('--sample_fraction',
                        help='For each pair of consecutive time points, we compute the Shannon diversity S of the transport map, then randomly sample max(S × sample_fraction, 10) pairs of points to add to the batch',
                        type=float, default=5e-6)
    parser.add_argument('--epoch_block_size',
                        help='Epoch block size', type=int, default=500)

    parser.add_argument('--exp2',
                        help='Element-wise 2 to the power x minus 1', action='store_true')
    parser.add_argument('--percentile',
                        help='Threshold matrix values at specified percentile (99.5 recommended)', type=float)
    parser.add_argument('--out',
                        help='Prefix for ouput file names')
    args = parser.parse_args(argsv)
    if args.out is None:
        args.out = wot.io.get_filename_and_extension(os.path.basename(args.matrix))[0] + '_grn'
    N = args.nmodules
    epochs = args.epochs
    TimeLag = args.time_lag

    lda_z1 = args.z1
    lda_z2 = args.z2
    lda_u = args.u2

    ds = wot.io.read_dataset(args.matrix)
    if scipy.sparse.isspmatrix(ds.X):
        ds.X = ds.X.toarray()
    if args.exp2:
        ds.X = 2 ** ds.X - 1

    if args.percentile is not None:
        thresh = np.percentile(ds.X, args.percentile)
        ds.X[(ds.X > thresh)] = thresh

    tf_ids = pd.read_table(args.tf, index_col=0, header=None).index.values
    tf_column_indices = ds.var.index.isin(tf_ids)
    ntf_sum = tf_column_indices.sum()
    if ntf_sum == 0:
        print('No transcription factors found')
        exit(1)
    print(str(ntf_sum) + ' transcription factors')

    non_tf_column_indices = ~tf_column_indices
    non_tf_sum = non_tf_column_indices.sum()
    if non_tf_sum == 0:
        print('No non-transcription factors found')
        exit(1)
    print(str(non_tf_sum) + ' non-transcription factors')

    TP = []
    Lineage = []  # list of transport maps
    if args.threads is None:
        threads = os.cpu_count()
    else:
        threads = args.threads

    differences = False

    Xg = []  # list of non-tf expression
    Xr = []  # list of tf expression

    for transport_map_dir in args.tmap:
        transport_maps = wot.io.list_transport_maps(transport_map_dir)
        if len(transport_maps) == 0:
            print('No transport maps found in ' + transport_map_dir)
            exit(1)
        for i in range(len(transport_maps)):
            tmap_dict = transport_maps[i]
            tmap = wot.io.read_dataset(tmap_dict['path'])
            if i == 0:  # first timepoint, align dataset with tmap rows
                aligned_order = ds.obs.index.get_indexer_for(tmap.obs.index.values.astype(str))
                if (aligned_order == -1).sum() > 0:
                    print(tmap.obs.index.values[aligned_order == -1])
                    nmissing = (aligned_order == -1).sum()
                    raise ValueError(str(nmissing) + ' missing ids for ' + tmap_dict['path'])

                ds_t = anndata.AnnData(ds.X[aligned_order], ds.obs.iloc[aligned_order], ds.var)
                Xg.append(ds_t.X[:, non_tf_column_indices])
                Xr.append(ds_t.X[:, tf_column_indices])
                TP.append(tmap_dict['t1'])
                if len(TP) > 0 and TP[len(TP) - 1] < TP[len(TP) - 2]:
                    Lineage.append([])

            # align dataset with tmap columns
            aligned_order = ds.obs.index.get_indexer_for(tmap.var.index.values.astype(str))
            if (aligned_order == -1).sum() > 0:
                print(tmap.var.index.values[aligned_order == -1])
                nmissing = (aligned_order == -1).sum()
                raise ValueError(str(nmissing) + ' missing ids for ' + tmap_dict['path'])

            ds_t = anndata.AnnData(ds.X[aligned_order], ds.obs.iloc[aligned_order], ds.var)
            Xg.append(ds_t.X[:, non_tf_column_indices])
            Xr.append(ds_t.X[:, tf_column_indices])
            TP.append(tmap_dict['t2'])
            Lineage.append(tmap.X)

    # if args.U is None:
    # rows are modules, columns are genes
    U = initialize_modules(Xg, N, threads=threads)

    # else:
    #     U = wot.io.read_dataset(args.U).X
    # TODO ensure in same order as ds

    Uinv = np.linalg.pinv(U)
    Z = []
    XU = []
    for xg, xr in zip(Xg, Xr):
        XU.append(xg.dot(Uinv))

    XU = np.vstack(XU)
    k, b, x0 = XU.max(0), 1.5, -0.2
    y0 = np.array([ki / 10 for ki in k])
    # for differences model:
    # lda_z1,lda_z2,lda_u = 0.1

    # optionally add constant to each day
    for i, x in enumerate(Xr):
        xconst = np.zeros((x.shape[0], len(TP)))
        xconst[:, i] = 1
        Xr[i] = np.hstack([x, xconst])

    # put all features on the same scale
    # (globally)
    Xr_avg = np.average(np.vstack(Xr), axis=0)
    Xr = [x / (Xr_avg + 0.001) for x in Xr]
    # (all TFs have mean 1 on every day)
    # Xr = [x/(np.average(x,axis=0) + 0.01) for x in Xr]

    ComposedLineage = compose_transports(Lineage, TP, TimeLag)

    # for original model: lda_z1=1.5,lda_z2=0.25,lda_u=1.5
    # for sparse model: lda_z1=3,lda_z2=1.5,lda_u=0.25
    # currently using: lda_z1=2,lda_z2=0.5,lda_u=1.5
    Z, U, Xh, k, b, y0, x0 = update_regulation(ComposedLineage, Xg, Xr, TP, TimeLag, Z=Z, U=U, lda_z1=lda_z1,
                                               lda_z2=lda_z2,
                                               lda_u=lda_u, epochs=epochs, sample_fraction=args.sample_fraction,
                                               threads=threads,
                                               inner_iters=1,
                                               k=k, b=b, y0=y0, x0=x0, differences=differences, frequent_fa=False,
                                               num_modules=N, epoch_block_size=args.epoch_block_size,
                                               savepath=None)

    # U has modules on rows, non-TFs on columns
    wot.io.write_dataset(
        anndata.AnnData(U, obs=pd.DataFrame(index=pd.RangeIndex(start=0, stop=U.shape[0], step=1)),
                    var=ds.var.iloc[non_tf_column_indices]),
        args.out + '_U', output_format='loom')

    # # Z has TFs on rows, modules on columns
    wot.io.write_dataset(
        anndata.AnnData(Z, obs=ds.var.iloc[tf_column_indices],
                    var=pd.DataFrame(index=pd.RangeIndex(start=0, stop=Z.shape[1], step=1))),
        args.out + '_Z', output_format='loom')

    with open(args.out + '_kbyx.txt', 'w') as writer:
        writer.write('k' + '\t' + 'b' + '\t' + 'y0' + '\t' + 'x0' + '\n')
        writer.write(str(k))
        writer.write('\t')
        writer.write(str(b))
        writer.write('\t')
        writer.write(str(y0))
        writer.write('\t')
        writer.write(str(x0))
        writer.write('\n')
    # np.save(args.out + '_kbyx.npy', (k, b, y0, x0))

    # for i, tp in enumerate(TP):
    #     if len(Xh[i]) > 0:
    #         # genes on columns, modules on rows
    #         name = args.out + '_regulators_deltaX.day-%d' % tp
    #         wot.io.write_dataset(
    #             anndata.AnnData(Xh[i], obs=pd.DataFrame(index=pd.RangeIndex(start=0, stop=Xh[i].shape[0], step=1)),
    #                         var=ds.var.iloc[non_tf_column_indices]),
    #             name, output_format='loom')
