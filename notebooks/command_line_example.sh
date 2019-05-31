#!/usr/bin/env bash

# Compute gene set scores, see notebook for an example of converting gene set scores to growth rates
wot gene_set_scores \
--matrix data/ExprMatrix.h5ad \
--method mean_z_score \
--gene_sets data/gene_sets.gmx

# neighborhood graph
wot neighborhood_graph \
--matrix data/ExprMatrix.var.genes.h5ad \
--space dmap \
--neighbors 50 \
--pca_comps 100 \
--diff_comps 20 \
--out fle-input.gexf

# FLE https://github.com/klarman-cell-observatory/forceatlas2
java -Djava.awt.headless=true -Xmx8g -cp forceatlas2.jar:gephi-toolkit-0.9.2-all.jar kco.forceatlas2.Main \
--input fle-input.gexf \
--output fle \
--nsteps 1000 \
--2d

# Compute transport maps
#wot optimal_transport \
#--matrix data/ExprMatrix.var.genes.h5ad \
#--cell_days data/cell_days.txt \
#--cell_filter data/serum_cell_ids.txt \
#--growth_iters 3 \
#--cell_growth_rates data/growth_gs_init.txt \
#--out tmaps/serum \
#--verbose

# Compute and plot trajectories
wot trajectory \
--tmap tmaps/serum \
--cell_set data/major_cell_sets.gmt \
--day 18 \
--embedding data/fle_coords.txt \
--verbose

# Compute trajectory divergence
wot trajectory_divergence \
--trajectory wot_trajectory.txt \
--cell_days data/cell_days.txt \
--matrix data/ExprMatrix.var.genes.h5ad \
--compare within \
--verbose \
--plot

# Compute and plot trajectory trends
wot trajectory_trends \
--trajectory wot_trajectory.txt \
--cell_days data/cell_days.txt \
--matrix data/ExprMatrix.h5ad \
--gene_filter Nanog,Obox6,Zfp42 \
--plot

wot fates \
--tmap tmaps/serum \
--cell_set data/major_cell_sets.gmt \
--day 17 \
--cell_set_filter IPS \
--out IPS_d17 \
--verbose

wot transition_table \
--tmap tmaps/serum \
--cell_set data/major_cell_sets.gmt \
--start_time 12 \
--end_time 18

# Compute differentially expressed genes at day 14 that are predictive of IPS fate at day 17
wot diff_exp \
--matrix data/ExprMatrix.h5ad \
--cell_days data/cell_days.txt \
--fate IPS_d17_fates.txt \
--fold_change 0 \
--gene_filter data/TFs.txt \
--cell_day_filter 14 \
--verbose

# Compute and plot validation summary
#wot optimal_transport_validation \
#--matrix data/ExprMatrix.var.genes.h5ad \
#--cell_days data/cell_days.txt \
#--cell_filter data/serum_cell_ids.txt \
#--covariate data/batches.txt \
#--cell_growth_rates tmaps/serum_g.txt \
#--cell_growth_rates_field g2 \
#--verbose

