#!/usr/bin/env bash

# generate example data
python 00_generating_data.py

# generate cell sets
wot cells_by_gene_set --matrix matrix.txt --gene_sets gene_sets.gmt --out cell_sets.gmt --format gmt --quantile 0.99

# compute transport maps
wot optimal_transport --matrix matrix.txt --cell_days days.txt

# compute trajectories
wot trajectory --tmap tmaps --cell_set cell_sets.gmt --day 7

# compute trajectory trends
wot trajectory_trends --trajectory wot_trajectory.txt --cell_days days.txt --matrix matrix.txt

# compute validation summary
wot optimal_transport_validation --matrix matrix.txt --cell_days days.txt --covariate covariate.txt
