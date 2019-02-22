import numpy as np

import wot.commands
import wot.graphics

# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
covariate_file = 'covariate.txt'
destination_file = 'validation_summary.png'
# --------------------------------------

ot_model = wot.ot.initialize_ot_model(matrix_file, days_file,
                                      covariate=covariate_file, growth_iters=1, tmap_out='val')
summary = wot.commands.compute_validation_summary(ot_model)
summary_stats = summary.groupby(['interval_mid', 'name'])['distance'].agg([np.mean, np.std])
summary_stats.to_csv('cv_validation_summary_stats.txt')
wot.graphics.plot_ot_validation_summary(summary_stats, 'cv_validation_summary.png')
