# WOT

Uses time-course data to infer how the probability distribution of cells in gene-expression space evolves over time,
by using the mathematical approach of Optimal Transport (OT)

* [Install](#install)
* [Prepare Gene Expression Matrix](#prepare-expression-matrix)
* [File Formats](#file_formats)
* [Apoptosis and Cell Cycle Scores](#gene_set_scores)
* [Optimal Transport](#optimal_transport)
* [Visualization](#visualization)



## <a name="install"></a> Install

```
git clone https://github.com/broadinstitute/wot.git
cd wot
git checkout develop
pip install -e .
```

## <a name="prepare-expression-matrix"></a> Prepare Expression Matrix
Apply a pre-processing workflow to normalize and scale your data, and detect variable genes.
Suggested tools include [Seurat](https://satijalab.org/seurat/) or [Scanpy](http://scanpy.readthedocs.io/en/latest/).


## <a name="file_formats"></a> File Formats
* Expression matrix
    * [Market Exchange Format (MEX)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)
    * [HDF5 Gene-Barcode Matrix](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices)
    * [Loom](http://linnarssonlab.org/loompy/format/index.html)
    * Tab-delimited Text
    
    Example:
    ```
    id TAB gene_1 TAB gene_2 TAB gene 3
    cell_1 TAB 1.2 TAB 12.2 TAB 5.4
    cell_2 TAB 1.2 TAB 2 TAB 3.0
    cell_n TAB 12.2 TAB 2 TAB 3
    ```    

* Cell days
    * Two column tab delimited file without header with cell ids and days.

    Example:
    ```
    cell_1 TAB 1
    cell_2 TAB 1
    cell_n TAB 2
    ```
* Day pairs
    * Two column tab delimited file without header with pairs of days to compute transport maps for.

    Example:
    ```
    0 TAB 2
    2 TAB 4
    4 TAB 6
    ```
## <a name="gene_set_scores"></a> Compute apoptosis and cell cycle scores to estimate growth rates (optional)
Required Options

Description | Flag
--- | --- |
**Normalized gene expression matrix with ALL GENES to compute an apoptosis and cell cycle score for each cell.** | --matrix
**Base name for output file** | --prefix


```
wot gene_set_scores --matrix my_expression_matrix.loom
```

## <a name="optimal_transport"></a> Optimal Transport
Required Options

Description | Flag
--- | --- |
**Normalized gene expression matrix of VARIABLE GENES.** | --matrix
**Assigns days to cells** | --cell_days
**Pairs of days to compute transport maps for** | --day_pairs
**Base name for output files** | --prefix

Optional Common Options

Description | Flag
--- | --- |
Apoptosis and cell cycle scores used to compute growth rates. If not specified, a constant growth rate is used. | --gene_set_scores 
Use principal component analysis to reduce the dimensionality of the expression matrix locally in the space of consecutive days to n_components| --local_pca 

To see all options, type:
```
wot ot -h 
```

