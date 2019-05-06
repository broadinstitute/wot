#!/usr/bin/env bash

# Compute gene set scores, see notebook for an example of converting gene set scores to growth rates
wot gene_set_scores \
--matrix data/ExprMatrix.loom \
--method mean_z_score \
--gene_sets data/gene_sets.gmx

# Compute transport maps
wot optimal_transport \
--matrix data/ExprMatrix.var.genes.loom \
--cell_days data/cell_days.txt \
--cell_filter data/serum_cell_ids.txt \
--growth_iters 3 \
--cell_growth_rates data/growth_gs_init.txt \
--out tmaps/serum \
--verbose

# Compute and plot trajectories
wot trajectory \
--tmap tmaps/serum \
--cell_set data/major_cell_sets.gmt \
--day 18 \
--embedding data/FLE.txt \
--verbose

# Compute trajectory divergence
wot trajectory_divergence \
--trajectory wot_trajectory.txt \
--cell_days data/cell_days.txt \
--matrix data/ExprMatrix.var.genes.loom \
--compare all \
--covariate data/batches.txt \
--plot

# Compute and plot trajectory trends
wot trajectory_trends \
--trajectory wot_trajectory.txt \
--cell_days data/cell_days.txt \
--matrix data/ExprMatrix.loom \
--gene_filter Nanog,Obox6,Shisa8,Zfp42 \
--plot

# Compute differentially expressed genes along trajectory
wot diff_exp \
--matrix data/ExprMatrix.loom \
--cell_days data/cell_days.txt \
--trajectory wot_trajectory.txt \
--nperm 1000 \
--verbose

# Compute and plot validation summary
wot optimal_transport_validation \
--matrix data/ExprMatrix.var.genes.loom \
--cell_days data/cell_days.txt \
--cell_filter data/serum_cell_ids.txt \
--covariate data/batches.txt \
--cell_growth_rates tmaps/serum_g.txt \
--cell_growth_rates_field g2 \
--verbose

