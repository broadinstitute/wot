# WOT: Waddington-OT

Uses time-course data to infer how the probability distribution of cells in gene-expression space evolves over time,
by using the mathematical approach of Optimal Transport (OT)

* [Install](#install)
* [Prepare Gene Expression Matrix](#prepare-expression-matrix)
* [Optimal Transport](#optimal_transport)
* [Visualization](#visualization)
* [File Formats](#file_formats)



## <a name="install"></a> Install
Python 3 is required.

```
git clone https://github.com/broadinstitute/wot.git
cd wot
git checkout develop
pip install -e .
```

## <a name="prepare-expression-matrix"></a> Prepare Expression Matrix
Apply a pre-processing workflow to normalize and scale your data, and detect variable genes.
Suggested tools include [Seurat](https://satijalab.org/seurat/) or [Scanpy](http://scanpy.readthedocs.io/en/latest/).

    

## <a name="optimal_transport"></a> Optimal Transport
Required Options

Description | Flag
--- | --- |
**Normalized gene expression matrix.** | --matrix
**Assigns days to cells** | --cell_days

Optional Common Options

Description | Flag
--- | --- |
Apoptosis and cell cycle scores used to compute growth rates. If not specified, a constant growth rate is used. The tool gene_set_scores can be used to compute gene set scores. | --gene_set_scores
Use principal component analysis to reduce the dimensionality of the expression matrix locally in the space of consecutive days. Thirty components are used by default. | --local_pca 
Pairs of days to compute transport maps for | --day_pairs
File with one gene id per line to use for computing cost matrices | --gene_filter
Base name for output files | --out
To see all options, type:
```
wot optimal_transport -h 
```

    
## <a name="visualization">Visualization</a>

## <a name="gene_set_scores">Gene Set Scores</a>
Required Option

Description | Flag
--- | --- |
**Normalized gene expression matrix to compute gene set scores (e.g. apoptosis and cell cycle scores) for each cell.** | --matrix

## <a name="file_formats"></a> File Formats
* Expression matrix
    * [Market Exchange Format (MEX)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices)
    * [HDF5 Gene-Barcode Matrix](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices)
    * [Loom](http://linnarssonlab.org/loompy/format/index.html)
    * Tab-delimited Text
    
    Example:
    <table>
    <tr><td>id</td><td>gene_1</td><td>gene_2</td><td>gene_3</td></tr>
    <tr><td>cell_1</td><td>1.2</td><td>12.2</td><td>5.4</td></tr>
    <tr><td>cell_2</td><td>2.3</td><td>4.1</td><td>5.0</td></tr>
    </table>
   
       

* Cell days
    * Two column tab delimited file header "id" and "day".

    Example:
    
    <table>
        <tr><td>id</td><td>day</td></tr>
        <tr><td>cell_1</td><td>1</td></tr>
        <tr><td>cell_2</td><td>2.5</td></tr>
        </table>
  
    
* Day pairs
    * Two column tab delimited file without header with pairs of days to compute transport maps for.

    Example:
    
    <table>
            <tr><td>0</td><td>2</td></tr>
            <tr><td>2</td><td>4</td></tr>
            <tr><td>4</td><td>6</td></tr>
            </table>
