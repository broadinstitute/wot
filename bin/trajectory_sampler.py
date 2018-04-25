import wot.ot
import pandas as pd
import argparse
import seaborn as sns
import numpy as np

parser = argparse.ArgumentParser(
    description='Compute cell ancestors/descendants')
parser.add_argument('--dir',
                    help='Directory of transport maps as produced by ot',
                    required=True)
parser.add_argument('--time',
                    help='The time',
                    required=True, type=float)
parser.add_argument('--prefix',
                    help='Prefix for ouput file names.',
                    required=True)
parser.add_argument('--matrix', help='Gene expression matrix')
parser.add_argument('--verbose', action='store_true',
                    help='Print progress information')
parser.add_argument('--gene', help='List of genes', action='append')
parser.add_argument('--gene_sets', help='Gene sets')

parser.add_argument('--cell_sets',
                    help='Grouping of cells into cell sets in gmt or gmx format',
                    required=True)

args = parser.parse_args()
time = args.time
cell_set_ds = wot.io.read_gene_sets(args.cell_sets)
transport_maps = wot.io.list_transport_maps(args.dir)
datasets = []
summaries = []
if args.gene is not None:
    genes = set(np.char.lower(np.array(args.gene)))
    ds = wot.io.read_dataset(args.matrix, col_filter={'id': lambda x: x.lower() in genes})
    datasets.append(ds)
    summaries.append('mean')

if args.gene_sets is not None:
    datasets.append(
        wot.io.read_dataset(args.gene_sets))
    summaries.append(None)

for t in range(len(transport_maps)):
    if transport_maps[t]['t2'] == time:
        time_index = t
        break

if time_index is None:
    raise RuntimeError(
        'Transport transport_map for time ' + str(time) + ' not found.')

result = wot.ot.TrajectorySampler.compute(cell_set_ds=cell_set_ds,
                                          transport_maps=transport_maps, unaligned_datasets=datasets,
                                          summaries=summaries)

sns.set_style("whitegrid")
sns.set_style("ticks", {"xtick.major.size": 2})
g = sns.factorplot(x="t", y="value", row="cell_set", col='name', data=pd.DataFrame(data=result),
                   kind='violin')
g.set_xticklabels(rotation=45)
