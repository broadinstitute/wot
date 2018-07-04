#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io
import numpy as np
import os
import pandas as pd
import subprocess
import pkg_resources
import tempfile


def compute_force_layout(ds, n_neighbors=100, n_comps=100, neighbors_diff=20, n_steps=10000):
    import anndata
    import scanpy.api as sc
    if type(ds) == wot.Dataset:
        adata = anndata.AnnData(ds.x, ds.row_meta, ds.col_meta)
    else:
        adata = ds
    if n_neighbors > 0:
        sc.pp.neighbors(adata, use_rep='X', n_neighbors=n_neighbors)
    if n_comps > 0:
        sc.tl.diffmap(adata, n_comps=n_comps)
        sc.pp.neighbors(adata, use_rep='X_diffmap', n_neighbors=neighbors_diff)
    input_graph_file = tempfile.mkstemp(prefix='gephi', suffix='.net')[1]
    output_coord_file = tempfile.mkstemp(prefix='coords', suffix='.txt')[1]
    W = adata.uns['neighbors']['connectivities']
    n_obs = W.shape[0]
    cids = adata.obs.index.values
    with open(input_graph_file, 'w') as writer:
        writer.write("*Vertices {n_obs}\n".format(n_obs=n_obs))
        for i in range(n_obs):
            writer.write("{node} \"{node}\"\n".format(node=i + 1))
        writer.write("*Edges\n")
        rows, cols = W.nonzero()
        for i, j in zip(rows, cols):
            if i < j:
                writer.write("{u} {v} {w:.6g}\n".format(u=i + 1, v=j + 1, w=W[i, j]))
    layout = 'fa'
    import psutil
    memory = int(0.5 * psutil.virtual_memory()[0] * 1e-9)
    classpath = os.path.dirname(
        pkg_resources.resource_filename('wot', 'commands/resources/graph_layout/GraphLayout.class')) + ':' + \
                pkg_resources.resource_filename('wot', 'commands/resources/graph_layout/gephi-toolkit-0.9.2-all.jar')
    subprocess.check_call(['java', '-Djava.awt.headless=true', '-Xmx{memory}g'.format(memory=memory), '-cp', classpath, \
                           'GraphLayout', input_graph_file, output_coord_file, layout, str(n_steps),
                           str(os.cpu_count())])
    # replace numbers with cids
    df = pd.read_table(output_coord_file, header=0, index_col='id')
    os.remove(output_coord_file)
    os.remove(input_graph_file)
    df.index = np.array(cids)[df.index.values - 1]
    return df, adata


def main(argv):
    import anndata
    import scanpy.api as sc
    parser = argparse.ArgumentParser(description='Force-directed layout embedding')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP, required=True)
    parser.add_argument('--neighbors', help='Number of nearest neighbors', type=int, default=100)
    parser.add_argument('--neighbors_diff', help='Number of nearest neighbors to use in diffusion component space',
                        type=int, default=20)
    parser.add_argument('--n_comps', help='Number of diffusion components', type=int, default=20)
    parser.add_argument('--n_steps', help='Number of force layout iterations', type=int, default=1000)
    parser.add_argument('--out', help='Output file name')
    args = parser.parse_args(argv)
    if args.out is None:
        args.out = 'wot'

    if os.path.isfile(args.matrix):
        ds = wot.io.read_dataset(args.matrix)
        compute_force_layout(ds, neighbors=args.neighbors, neighbors_diff=args.neighbors_diff, n_comps=args.n_comps,
                     n_steps=args.n_steps)

    else:
        print('Input matrix is not a file')
        exit(1)
        cids = []
        transport_maps = wot.io.list_transport_maps(args.input)
        if len(transport_maps) == 0:
            print('No transport maps found in ' + args.input)
            exit(1)
        datasets = []
        for tm in transport_maps:
            datasets.append(wot.io.read_dataset(tm['path']))

        with open(input_graph_file, 'w') as writer:
            n_obs = 0
            name_to_id = {}
            for ds_index in range(len(datasets)):
                ds = datasets[ds_index]
                n_obs += ds.x.shape[0]
                if ds_index == len(datasets) - 1:
                    n_obs += datasets[len(datasets) - 1].x.shape[1]
            writer.write("*Vertices {n_obs}\n".format(n_obs=n_obs))
            node_counter = 0
            for ds_index in range(len(datasets)):
                ds = datasets[ds_index]
                for cid in ds.row_meta.index.values:
                    node_counter += 1
                    name_to_id[cid] = node_counter
                    cids.append(cid)
                    writer.write("{node} \"{node}\"\n".format(node=node_counter))
                if ds_index == len(datasets) - 1:
                    for cid in ds.col_meta.index.values:
                        node_counter += 1
                        name_to_id[cid] = node_counter
                        cids.append(cid)
                        writer.write("{node} \"{node}\"\n".format(node=node_counter))

            writer.write("*Edges\n")
            for ds_index in range(len(datasets)):
                ds = datasets[ds_index]

                for i in range(ds.x.shape[0]):
                    start_id = name_to_id[ds.row_meta.index.values[i]]
                    order = np.argsort(ds.x[i])

                    for j in range(n_neighbors):
                        writer.write(
                            "{u} {v} {w:.6g}\n".format(u=start_id, v=name_to_id[ds.col_meta.index.values[order[j]]],
                                                       w=ds.x[i, order[j]]))
                # internal nearest neighbors
                adata = anndata.AnnData(ds.x, ds.row_meta, ds.col_meta)
                sc.pp.neighbors(adata, use_rep='X', n_neighbors=n_neighbors)
                W = adata.uns['neighbors']['connectivities']
                rows, cols = W.nonzero()
                for i, j in zip(rows, cols):
                    if i < j:
                        writer.write("{u} {v} {w:.6g}\n".format(u=i + 1, v=j + 1, w=W[i, j]))

    df, adata = compute_force_layout(ds)
    adata.write(args.out + '.h5ad')
    csv_file = args.out if args.out.lower().endswith('.csv') else args.out + '.csv'
    df.to_csv(csv_file, index_label='id')
