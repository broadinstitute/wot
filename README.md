# WOT: Waddington-OT

Uses time-course data to infer how the probability distribution of cells in gene-expression space evolves over time,
by using the mathematical approach of Optimal Transport (OT)

* [Install](#install)
* [Prepare Gene Expression Matrix](#prepare-expression-matrix)
* [Optimal Transport](#optimal_transport)
* [Visualization](#visualization)
* [Optimal Transport Trajectory](#optimal_transport_trajectory)
* [Optimal Transport Transition Table](#optimal_transport_transition_table)
* [Optimal Transport Validation](#optimal_transport_validation)
* [Gene Set Scores](#gene_set_score)
* [File Formats](#file_formats)



## <a name="install"></a> Install
Python 3 is required.

Optionally install dependencies using conda:
```
conda install cython h5py flask gunicorn numexpr numpy pandas scikit-learn scipy simplejson psutil
conda install -c conda-forge pot
```

Install WOT:

```
pip install wot
```

## <a name="prepare-expression-matrix"></a> Prepare Expression Matrix
Apply a pre-processing workflow to normalize and scale your data, and detect variable genes.
Suggested tools include [Seurat](https://satijalab.org/seurat/) or [Scanpy](http://scanpy.readthedocs.io/en/latest/).

    
## WOT Tools
WOT tools are run using the syntax *wot tool*. To see all tool options, type *wot tool -h* (e.g. wot optimal_transport -h)


## <a name="optimal_transport"></a> Optimal Transport
*wot optimal_transport* calculates transport maps between consecutive time points and automatically learns cellular growth and death rates.

Common Options (required in **bold**)

Flag | Description
--- | --- |
**matrix** | [Normalized gene expression matrix](#matrix).
**cell_days** | [Assigns days to cells](#cell_days)
gene_set_scores | Apoptosis and cell cycle scores used to compute growth rates. If not specified, a constant growth rate is used. The wot tool [gene_set_scores](#gene_set_scores) can be used to compute gene set scores.
local_pca | Use principal component analysis to reduce the dimensionality of the expression matrix locally in the space of consecutive days. Thirty components are used by default.
day_pairs | [Pairs of days to compute transport maps for](#day_pairs)
gene_filter | File with one gene id per line to use for computing cost matrices (e.g. variable genes)
out | Base name for output files 


## <a name="validation">Optimal Transport Validation</a>
*wot optimal_transport_validation* tests the performance of optimal transport by comparing interpolated distributions to held-out time points.

optimal_transport_validation has the same options as [optimal_transport](#optimal_transport) and

Flag | Description
--- | --- |
**t_interpolate** | Interpolation fraction between two time points
covariate | Two column file with headers "id" and "covariate" indicating cell ids and covariate value



## <a name="visualization">Visualization</a>
*wot force_layout* generates a force-directed layout using the ForceAtlas2 algorithm. The tool can optionally perform dimensionality reduction using diffusion component embedding of the dataset. 
After you have generated the force layout coordinates, you can use *wot wot_server* to view trajectories and gene expression in the force layout.

**Trajectory Trends** plots the expression of a gene over time based on the transport maps

## <a name="optimal_transport_trajectory">Trajectory<a>
*wot trajectory* generate ancestors and descendants given a starting cell cet and transport maps.
Options

Flag | Description
--- | --- |
**tmap** |Directory of transport maps as produced by [optimal_transport](#optimal_transport)
**cell_set** | [Assigns cells to cell sets](#cell_sets) 
**cell_days** | [Assigns days to cells](#cell_days)
 


## <a name="optimal_transport_transition_table">Transition Table</a>
*wot transition_table* generate a transition table from one cell set to another cell set.
Options

Flag | Description
--- | --- |
**tmap** |Directory of transport maps as produced by [optimal_transport](#optimal_transport)
**cell_set** | [Assigns cells to cell sets](#cell_sets) 
**cell_days** | [Assigns days to cells](#cell_days)
**start_time** | The start time for the cell sets to compute the transitions to cell sets at end_time
**end_time** | The end time. 

    
## <a name="gene_set_scores">Gene Set Scores</a>
*wot gene_set_scores* computes gene set scores for each cell given a gene expression matrix and gene sets.

Options

Flag | Description
--- | --- |
**matrix** | Normalized gene expression matrix to compute gene set scores (e.g. apoptosis and cell cycle scores) for each cell.

## <a name="file_formats"></a> File Formats

#### <a name="matrix">Gene Expression matrix</a> 
Cells on rows and genes (features) on columns. Accepted formats are [Market Exchange Format (MEX)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices), [HDF5 Gene-Barcode Matrix](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices), [loom](http://linnarssonlab.org/loompy/format/index.html), [h5ad](http://scanpy.readthedocs.io/en/latest/), and text
    
Example Text File:
    
<table>
<tr><td>id</td><td>gene_1</td><td>gene_2</td><td>gene_3</td></tr>
<tr><td>cell_1</td><td>1.2</td><td>12.2</td><td>5.4</td></tr>
<tr><td>cell_2</td><td>2.3</td><td>4.1</td><td>5.0</td></tr>
</table>
   
       

#### <a name="cell_days">Cell Days</a>
Two column file with header "id" and "day".

Example:

<table>
<tr><td>id</td><td>day</td></tr>
<tr><td>cell_1</td><td>1</td></tr>
<tr><td>cell_2</td><td>2.5</td></tr>
</table>
  
#### <a name="day_pairs">Days Pairs</a> 
Two column file without header with pairs of days to compute transport maps for.

Example:

<table>
<tr><td>0</td><td>2</td></tr>
<tr><td>2</td><td>4</td></tr>
<tr><td>4</td><td>6</td></tr>
</table>


#### <a name="cell_sets">Cell Sets</a>
Cell sets can be provided in [gmt](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) or [gmx](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29). 
