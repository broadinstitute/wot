# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
# --------------------------------------

import wot

wot.initialize_ot_model(matrix_file, days_file, growth_iters=1,
        epsilon=.02, lambda1=10, lambda2=80, local_pca=0)
