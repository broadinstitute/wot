import wot.ot

import unittest

import pandas as pd
import numpy as np
import anndata

class TestIO(unittest.TestCase):
    
    def test_costs_single_tmap(self):
        #Construct an aritifical adata with 6 cells (3 per day)
        #and 1000 random genes
        obs_df = pd.DataFrame(data={'day': [1, 1, 1, 2, 2, 2]})
        var_df = pd.DataFrame(data=None, index=np.arange(1000))
        expr_matr = np.random.rand(6, 1000)
        adata = anndata.AnnData(X=expr_matr, obs=obs_df, var=var_df)
        
        #Cells should transport along the diagonal 
        cost_matrix = np.array([[0, 100, 100],
                                [100, 0, 100],
                                [100, 100, 0]])

        #Make an ot_model and transport
        ot_model = wot.ot.OTModel(adata, epsilon=0.01, lambda1=1, lambda2=50)
        tmap = ot_model.compute_transport_map(1, 2, cost_matrix=cost_matrix)

        #Check if the transport map is within tolerance
        expected_tmap = np.array([[1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]])
        are_similar = np.allclose(tmap.X, expected_tmap, atol=0.01, rtol=0)
        
        self.assertTrue(are_similar)
        
    def test_costs_multiple_tmaps(self):
        #Construct an aritifical adata with 9 cells (3 per day)
        #and 1000 random genes
        obs_df = pd.DataFrame(data={'day': [1, 1, 1, 2, 2, 2, 3, 3, 3]})
        var_df = pd.DataFrame(data=None, index=np.arange(1000))
        expr_matr = np.random.rand(9, 1000)
        adata = anndata.AnnData(X=expr_matr, obs=obs_df, var=var_df)
        
        #Cells should transport along the diagonal for each timepoint
        cost_matrix = np.array([[0, 100, 100],
                                [100, 0, 100],
                                [100, 100, 0]])
        cost_matrices = [cost_matrix, cost_matrix]
        
        #Make an ot_model and transport
        ot_model = wot.ot.OTModel(adata, epsilon=0.01, lambda1=1, lambda2=50)
        ot_model.compute_all_transport_maps(tmap_out='tmaps_out/tmap', cost_matrices=cost_matrices)
        
        #Read the transport maps and compare them to expected
        tmap_1 = anndata.read_h5ad('tmaps_out/tmap_1_2.h5ad')
        tmap_2 = anndata.read_h5ad('tmaps_out/tmap_2_3.h5ad')
        
        #Check if the transport maps are within tolerance
        expected_tmap = np.array([[1, 0, 0],
                                [0, 1, 0],
                                [0, 0, 1]])
        are_similar_tmap_1 = np.allclose(tmap_1.X, expected_tmap, atol=0.01, rtol=0)
        are_similar_tmap_2 = np.allclose(tmap_2.X, expected_tmap, atol=0.01, rtol=0)
        
        self.assertTrue(are_similar_tmap_1 and are_similar_tmap_2)

        