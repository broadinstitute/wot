import wot.ot
import pandas as pd

r = wot.ot.Ancestors.from_cmd_line()
result = wot.ot.Ancestors.compute(cell_set_ds=r['cell_set_ds'], transport_maps=r['transport_maps'],
                                  time_index=r['time_index'], t2_index=0, full_ds=r['full_ds'],
                                  gene_set_scores=r['gene_set_scores'], genes=r['genes'], verbose=r['args'].verbose)
import seaborn as sns

sns.set_style("whitegrid")
sns.set_style("ticks", {"xtick.major.size": 2})
g = sns.factorplot(x="t", y="value", row="cell_set", col='name', data=pd.DataFrame(data=result),
                   kind='violin')
g.set_xticklabels(rotation=45)
g.savefig(r['args'].prefix + 'png', dpi=200)
