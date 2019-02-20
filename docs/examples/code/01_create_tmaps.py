import wot.ot

# ------ Configuration variables -------
matrix_file = 'matrix.txt'
days_file = 'days.txt'
# --------------------------------------

ot_model = wot.ot.initialize_ot_model(matrix_file, days_file, growth_iters=1,
                                      epsilon=.02, lambda1=10, lambda2=80, local_pca=0)
ot_model.compute_all_transport_maps()
