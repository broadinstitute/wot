#!/usr/bin/env bash

# generate example data
python 00_generating_data.py

# compute gene set scores
wot gene_set_scores --matrix matrix.txt --gene_sets gene_sets.gmt --out gene_set_scores.txt --method mean_z_score

# generate cell sets
wot cells_by_gene_set --score gene_set_scores.txt_Myeloid_stem_cells.txt --score gene_set_scores.txt_Red_blood_cells.txt \
--score gene_set_scores.txt_Granulocytes.txt --score gene_set_scores.txt_Lymphocytes.txt --out cell_sets.gmt

# compute transport maps
wot optimal_transport --matrix matrix.txt --cell_days days.txt

# compute and plot trajectories
wot trajectory --tmap tmaps --cell_set cell_sets.gmt --day 7 --embedding embedding.csv

# compute trajectory divergence
wot trajectory_divergence --trajectory wot_trajectory.txt --cell_days days.txt --matrix matrix.txt --compare all --covariate covariate.txt --plot

# compute and plot trajectory trends
wot trajectory_trends --trajectory wot_trajectory.txt --cell_days days.txt --matrix matrix.txt --plot

# compute differentially expressed genes along trajectory
wot diff_exp --matrix matrix.txt --cell_days days.txt --trajectory wot_trajectory.txt

# compute and plot validation summary
wot optimal_transport_validation --matrix matrix.txt --cell_days days.txt --covariate covariate.txt
