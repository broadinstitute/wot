---
noheader: true
layout: documentation
location: Examples
---

# Computing Transport maps
--------------------------

**wot** uses the Optimal Transport approach to infer relationships between cells.
It thus relies on the computation of transport maps between consecutive timepoints.

These transport maps can be computed from the command line with :

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt \
 --out tmaps --local_pca -1
```

Several options are available for this command, to choose the Optimal Transport parameters
that best suit your needs. Please refer to [the CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps)
for more information on those.

Here is another example, specifying the default values explicitely :

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt \
 --out tmaps --local_pca -1 --scaling_iter 3000 \
 --epsilon 0.05 --lambda1 1 --lambda2 50
```
