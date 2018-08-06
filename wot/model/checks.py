# -*- coding: utf-8 -*-

import wot

def check_model_consistency(ot_model):
    check_local_pca(ot_model)
    check_days_information(ot_model)
    pass

def check_local_pca(ot_model):
    local_pca = ot_model.get_ot_config()['local_pca']
    if local_pca > ot_model.matrix.x.shape[1]:
        print("Warning : local_pca set to {}, above gene count of {}. Disabling PCA"\
                .format(local_pca, ot_model.matrix.x.shape[1]))
        ot_model.ot_config['local_pca'] = 0

def check_days_information(ot_model):
    if 'day' not in ot_model.matrix.row_meta.columns:
        raise ValueError("Days information not available for matrix")
    if any(ot_model.matrix.row_meta['day'].isnull()):
        query = ot_model.matrix.row_meta['day'].isnull()
        faulty = list(ot_model.matrix.row_meta.index[query])
        raise ValueError("Days information missing for cells : {}".format(faulty))

def check_covariate_information(ot_model):
    pass
