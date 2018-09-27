#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import subprocess
import tempfile

import numpy as np
import pandas as pd
import pkg_resources
import wot.io


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

    run_gephi(input_graph_file, output_coord_file, n_steps)
    # replace numbers with cids
    df = pd.read_table(output_coord_file, header=0, index_col='id')
    os.remove(output_coord_file)
    os.remove(input_graph_file)
    df.index = np.array(cids)[df.index.values - 1]
    return df, adata


def run_gephi(input_graph_file, output_coord_file, n_steps):
    layout = 'fa'
    import psutil
    memory = int(0.5 * psutil.virtual_memory()[0] * 1e-9)
    classpath = os.path.dirname(
        pkg_resources.resource_filename('wot', 'commands/resources/graph_layout/GraphLayout.class')) + ':' + \
                pkg_resources.resource_filename('wot', 'commands/resources/graph_layout/gephi-toolkit-0.9.2-all.jar')
    subprocess.check_call(['java', '-Djava.awt.headless=true', '-Xmx{memory}g'.format(memory=memory), '-cp', classpath, \
                           'GraphLayout', '--input', input_graph_file, '--output', output_coord_file, '--layout',
                           layout, '--nsteps', str(n_steps)])


def main(argv):
    parser = argparse.ArgumentParser(description='Force-directed layout embedding')
    parser.add_argument('--matrix', help=wot.commands.MATRIX_HELP)
    parser.add_argument('--n_comps', help='Number of diffusion components', type=int, default=20)
    parser.add_argument('--neighbors', help='Number of nearest neighbors', type=int, default=100)
    parser.add_argument('--neighbors_diff', help='Number of nearest neighbors to use in diffusion component space',
                        type=int, default=20)

    parser.add_argument('--graph',
                        help='Precomputed graph to use instead of matrix. See https://gephi.org/users/supported-graph-formats/')

    parser.add_argument('--n_steps', help='Number of force layout iterations', type=int, default=1000)
    parser.add_argument('--out', help='Output file name')
    args = parser.parse_args(argv)
    if args.out is None:
        args.out = 'wot'
    import os
    jar_file = pkg_resources.resource_filename('wot', 'commands/resources/graph_layout/gephi-toolkit-0.9.2-all.jar')
    if not os.path.isfile(jar_file):
        import urllib.request
        import shutil
        print('Downloading gephi jar file')
        with urllib.request.urlopen(
                'https://github.com/gephi/gephi-toolkit/releases/download/v0.9.2/gephi-toolkit-0.9.2-all.jar') as response, open(
            jar_file, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
    if args.matrix is not None:
        ds = wot.io.read_dataset(args.matrix)
        df, adata = compute_force_layout(ds, n_neighbors=args.neighbors,
                                         neighbors_diff=args.neighbors_diff, n_comps=args.n_comps,
                                         n_steps=args.n_steps)
        adata.write(args.out + '.h5ad')
        csv_file = args.out if args.out.lower().endswith('.csv') else args.out + '.csv'
        df.to_csv(csv_file, index_label='id')
    elif args.graph is not None:
        run_gephi(args.graph, args.out, args.n_steps)
