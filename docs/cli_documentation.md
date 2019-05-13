---
noheader: true
title: Command Line
permalink: cli_documentation/
layout: documentation
location: Documentation
---


## Command Line ##
-----------

**wot** consists of several tools. Each tool can be used with the syntax `wot <tool>`.

Help is available for each tool with `wot <tool> -h`. For instance:

```
wot optimal_transport -h
```

In the following sections, each command is described with an example and
a table containing the core options. Required options are in bold font.

<hr />




### Transport maps ###

```sh
wot optimal_transport --matrix matrix.txt \
 --cell_days days.txt --out tmaps
```

This command will create a file `tmaps_{A}_{B}.h5ad` for each pair `{A}`, `{B}`
of consecutive timepoints. These transport maps can then be converted to any format you find
convenient with the [convert_matrix tool](#matrix_file).


<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--matrix</b></td>
      <td>Normalized gene expression matrix. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--cell_days</b></td>
      <td>Timestamps for each cell. See <a href="#days_file">formats</a></td>
    </tr>
    <tr>
	  <td>--tranpose</td>
	  <td>Swap the rows and column of the input matrix</td>
	</tr>
    <tr>
      <td>--epsilon</td>
      <td>Regularization parameter that controls the entropy of the transport map<br/>default : 0.05</td>
    </tr>
    <tr>
      <td>--lambda1</td>
      <td>Regularization parameter that controls the fidelity of the constraint on p<br/>default : 1</td>
    </tr>
    <tr>
      <td>--lambda2</td>
      <td>Regularization parameter that controls the fidelity of the constraint on q<br/>default : 50</td>
    </tr>
    <tr>
      <td>--local_pca</td>
      <td>Number of dimensions to use when doing PCA<br/>default : 30</td>
    </tr>
    <tr>
          <td>--growth_iters</td>
          <td>Number of growth iterations for learning the growth rate.<br/>default : 1</td>
        </tr>
<tr>
<td>--cell_growth_rates</td>
<td>File with "id" and "cell_growth_rate" headers corresponding to cell id and growth rate per day.</td>
</tr>
    <tr>
    <td>--gene_filter</td>
    <td>File with one gene id per line to use (e.g. variable genes)</td>
    </tr>
  </tbody>
</table>




##### Local PCA #####

The default transport cost uses Principal Component Analysis to reduce the
dimension of the data before computing distances pairs of timepoints.
By default, **wot** uses 30 dimensions.

Dimensionality reduction can be disabled with `--local_pca 0`

If you specify more dimensions for the PCA than your dataset has genes,
**wot** will skip PCA and print a warning.


### Trajectories ###

Ancestors and descendants in **wot** are computed through the `trajectory` tool.

You can select a **cell set** by specifying a [cell set file](#cellset_file).
You can either manually edit this type of file, or generate it from a gene set file
using the [cells_by_gene_set](#cells_by_gene_set) tool. Please note that the order of the cell ids in the output trajectory file
may differ from that of the expression matrix used to generate transport maps.

```sh
wot trajectory --tmap tmaps \
 --cell_set cell_sets.gmt --day 10 --out traj.txt
```

<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--tmap</b></td>
      <td>Prefix of the transport maps configuration file</td>
    </tr>
    <tr>
      <td><b>--cell_set</b></td>
      <td>Target cell set. See <a href="#cellset_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--day</b></td>
      <td>The target timepoint to consider if some cell sets span across multiple timepoints</td>
    </tr>
    <tr>
      <td>--out</td>
      <td>Output filename<br/>default : wot_trajectory.txt. Please note that the order of the cell ids in the output trajectory file
		may differ from that of the expression matrix used to generate transport maps.</td>
    </tr>
  </tbody>
</table>


### Ancestor census ###

The ancestor census command computes the amount of mass of an
ancestor distribution falls into each cell set.


```sh
wot census --tmap tmaps  \
 --cell_set cell_sets.gmt \
 --out census --day 10
```

This command creates several census files named `<prefix>_<cellset>.txt`,
for instance `census_tip1.txt`. See <a href="#census_file">formats</a>
for more information.

<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--tmap</b></td>
      <td>Prefix of the transport maps configuration file</td>
    </tr>
    <tr>
      <td><b>--cell_set</b></td>
      <td>Target cell sets. See <a href="#cellset_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--day</b></td>
      <td>The target timepoint to consider if some cell sets span across multiple timepoints</td>
    </tr>
    <tr>
      <td>--out</td>
      <td>Output filenames prefix.<br/>default : 'census'</td>
    </tr>
  </tbody>
</table>


### Trajectory trends ###

The trajectory trends command computes the mean and variance of each gene in the specified matrix for the given trajectories.
Please use the trajectory tool to compute trajectories.


```sh
wot trajectory_trends --trajectory tmaps \
 --matrix matrix.txt \
 --out trends 
```

This will create a file `trends_<trajectory name>.txt` with the mean expression
profile among ancestors/descendants and `trends_<trajectory name>.variance.txt` with the
variance for each feature at each timepoint

<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--trajectory</b></td>
      <td>The trajectory produced by the trajectory tool</td>
    </tr>
    <tr>
      <td><b>--matrix</b></td>
      <td>Normalized gene expression matrix. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td>--out</td>
      <td>Output filenames prefix.<br/>default : 'trends'</td>
    </tr>
  </tbody>
</table>


### Differential expression ###

The diff_exp command finds the genes that are differentially expressed between two sets of cells. You input one or more fates created using
the fates tool, an expression matrix, and the tool outputs a table with statistics about differentially expressed genes.

### Validation ###

You can validate the transport maps that have been computed above.

Say, for instance, that you have data at time points 0, 1, and 2.
You could compute the transport map from 0 to 2, and then interpolate with
optimal transport to see what distribution is expected at time 1 by this model.
It is then easy to compare with the actual distributions at those time points,
and thus check that the transport maps properly predict the cell's trajectory.

```sh
wot optimal_transport_validation --matrix matrix.txt \
 --cell_days days.txt --covariate covariate.txt
```

This would create a validation summary, as a text file, containing
all the information needed to evaluate the accuracy of the predictions.

<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--matrix</b></td>
      <td>Normalized gene expression matrix. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--cell_days</b></td>
      <td>Timestamps for each cell. See <a href="#days_file">formats</a></td>
    </tr>
    <tr>
      <td>--covariate</td>
      <td>Covariate value for each cell. See <a href="#covariate_file">formats</a></td>
    </tr>
    <tr>
      <td>--out</td>
      <td>The filename for the validation summary<br/>default : validation_summary.txt</td>
    </tr>
  </tbody>
</table>

The validation tool also accepts all options of the
[optimal_transport tool](#transport-maps). Mandatory options are
reproduced here for convenience.



##### Covariate #####

To measure the quality of interpolation, we compare the distance between
batches of cells at the same time point.

The covariate values may be specified in a tab-separated text file.
It must have exactly two headers : "id" and "covariate".
Each subsequent line must consist of a cell name, a tab, and a covariate value.

##### Null hypothesises #####

**wot** tests the results of the Optimal Transport predictions
against three different null hypothesises.

- **First point** :
  The state of the population at time t0 is a good approximation for its
  state between times t0 and t1.
  This is equivalent to saying cells stay where they were at time t0,
  and then instantly move to their next state exactly at time t1.
- **Last point** :
  The state of the population at time t1 is a good approximation for its
  state between times t0 and t1.
  This is equivalent to saying cells are in a given state at time t0, and
  then instantly move to their next state, and stay there until time t1.
- **Randomized interpolation** :
  Linear interpolation between randomly-picked pairs of points from t0 and t1
  is a good approximation for the state of the population between t0 and t1.


#### Validation summary ####

The validation summary, generated as text file, is organized as follows :

Each line contains information about the relation between two cell sets :
 - **interval_start** and **interval_end** indicate which pair of timepoints is being considered.
 - **pair0** and **pair1** indicate which sets are being considered:
   - **P** is the real population at the intermediate timepoint
   - **F** is the estimated population in *first point* hypothesis
   - **L** is the estimated population in *last point* hypothesis
   - **R** is the estimated population in *randomized interpolation* hypothesis
   - **I** is the estimated population in *optimal transport* hypothesis
 - **distance** is the Wasserstein distance between the two sets considered

**pair0** and **pair1** will have a suffix of the form **_cvX_cvY** or **_cvZ**
when covariates are used, to indicate which batch they were extracted from.


### Embedding ###

One method of visualizing data in two or three dimensions is
[Force-directed Layout Embedding (FLE)](https://en.wikipedia.org/wiki/Force-directed_graph_drawing).
We compute this nearest neighbor graph in diffusion component space.
We recommend using [forceatlas2](https://github.com/klarman-cell-observatory/forceatlas2) to run the force directed layout.
You can generate a graph to use as input to the force directed layout using the command neighborhood_graph

```sh
wot neighborhood_graph --matrix matrix.txt
```

<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--matrix</b></td>
      <td>Normalized gene expression matrix. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
		<td>--gene_filter</td>
        <td>File with one gene id per line to use (e.g. variable genes)</td>
    </tr>
    <tr>
		<td>--cell_filter</td>
		<td>File with one cell id per line to include from the matrix</td>
    </tr>
    <tr>
      <td>--transpose</td>
      <td>Transpose the matrix</td>
    </tr>
	<tr>
  		<td>--pca_comps</td>
        <td>Number of PCA components<br/>default : 50</td>
     </tr>
    <tr>
      <td>--diff_comps</td>
      <td>Number of diffusion components<br/>default : 15</td>
    </tr>
    <tr>
      <td>--neighbors</td>
      <td>Number of nearest neighbors<br/>default : 15</td>
    </tr>
      <tr>
		<td>--space</td>
        <td>Space to compute the neighborhood graph in', choices=['dmap', 'pca', 'input']</td>
        </tr>
    <tr>
      <td>--out</td>
      <td>Output file name. The file is saved in gexf format (https://gephi.org/gexf/format/)</td>
    </tr>
  </tbody>
</table>


## Supported file formats ##
----------------------------

### <a name="matrix_file">Gene expression matrices</a> ###

The *matrix* file specifies the gene expression matrix to use.

The following formats are accepted by all tools: *mtx*, *txt*, *h5ad*, *loom*, and *gct*. Please note that *wot* expects 
cells on the rows and genes on the columns, except for the *mtx* format.

##### Text #####

The text format consists of tab or comma separated columns with genes on the columns and cells on the rows.

The first row, the header, must consist of an "id" field, and then the list of genes to be considered.

Each subsequent row will give the expression level of each gene for a given cell.

The first field must be a unique identifier for the cell, and then the tab or comma separated list
of expression levels for each gene/feature.

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>gene_1</td><td>gene_2</td><td>gene_3</td></tr>
<tr><td>cell_1</td><td>1.2</td><td>12.2</td><td>5.4</td></tr>
<tr><td>cell_2</td><td>2.3</td><td>4.1</td><td>5.0</td></tr>
</table>

##### MTX #####

The MTX format is a sparse matrix format with genes on the rows and cells on the columns as output by [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices). 
You should also have TSV files with genes and barcode sequences corresponding to row and column indices, respectively. 
These files must be located in the same folder as the MTX file with the same base file name. For example if the MTX file is my_data.mtx, you should
also have a my_data.genes.txt file and a my_data.barcodes.txt file.


##### GCT #####

A GCT file is a tab-delimited text file. Please see description [here](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT)


##### H5AD #####

A HDF5 file that provides a scalable way of keeping track of data together with learned annotations.. Please see description at [https://anndata.readthedocs.io](https://anndata.readthedocs.io/en/latest/)


##### Loom #####

A HDF5 file for efficient storage and access of large datases. Please see description at [http://loompy.org/](http://loompy.org/)

---

> <button class="btn-info rounded border-0 px-3 py-1" disabled>Note</button>
>
> You can convert input matrices from one format to another :
>
> ```sh
> wot convert_matrix --format txt *.loom
> ```
>
> This command converts all loom files in the current directory to text.
>
> Supported input formats are *mtx*, *hdf5*, *h5*, *h5ad*, *loom*, and *gct*.
>
> Supported output formats are *txt*, *loom*, and *gct*. (default: *loom*)



### <a name="days_file">Cell timestamps (cell day file)</a> ###

The timestamp associated with each cell of the matrix file is specified in the *days* file.
This file must be a tab-separated plain text file, with two header fields: "id" and "day".

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>day</td></tr>
<tr><td>cell_1</td><td>1</td></tr>
<tr><td>cell_2</td><td>2.5</td></tr>
</table>

### OT Configuration file ###

There are several options to specify Optimal Transport parameters in **wot**.

The easiest is to just use constant parameters and specify them when
computing transport maps with the `--epsilon` or `--lambda1` options.

If you want more control over what parameters are used, you can use a
detailed configuration file. There are two kinds of configuration files
accepted by **wot**.

#### Configuration per timepoint ####

You can specify each parameter at each timepoint.
When computing a transport map between two timepoints, the average
of the two parameters for the considered timepoints will be taken into account.

For instance, if you have prior knowledge of the amount of entropy
at each timepoint, you could specify a different value of epsilon for each
timepoint, and those would be used to compute more accurate transport maps.

The configuration file is a tab-separated text file that starts with a header
that must contain a column named `t`, for the timepoint, and then the name
of any parameter you want to set. Any parameter not specified in this
file can be specified as being constant as previously, with the command-line
arguments `--epsilon`, `--lambda1`, `--tolerance`, etc. .

Example:

<table class="table" style="display: table">
<tr><td>t</td><td>epsilon</td></tr>
<tr><td>0</td><td>0.001</td></tr>
<tr><td>1</td><td>0.002</td></tr>
<tr><td>2</td><td>0.005</td></tr>
<tr><td>3</td><td>0.008</td></tr>
<tr><td>3.5</td><td>0.01</td></tr>
<tr><td>4</td><td>0.005</td></tr>
<tr><td>5</td><td>0.001</td></tr>
</table>

#### Configuration per pair of timepoints ####

If you want to be even more explicit about what parameters to use for each
transport map computation, you can specify parameters for pairs of timepoints.

As previously, the configuration is specified in a tab-separated text file.
Its header must have columns `t0` and `t1`, for source and destination timepoints.

Bear in mind though, that any pair of timepoints not specified in this file
will not be computable. That means you should at least put all pairs
of consecutive timepoints if you want to be able to compute full trajectories.

Example:

<table class="table" style="display: table">
<tr><td>t0</td><td>t1</td><td>lambda1</td></tr>
<tr><td>0</td><td>1</td><td>50</td></tr>
<tr><td>1</td><td>2</td><td>80</td></tr>
<tr><td>2</td><td>4</td><td>30</td></tr>
<tr><td>4</td><td>5</td><td>10</td></tr>
</table>

This can for instance be used if you want to skip a timepoint (note how
timepoints 3 or 3.5 are not present here). If a timepoint is present in the
dataset but not in this configuration file, it will be ignored.

You can use as many parameter columns as you want, even none.
All parameters not specified here can be specified as being constant as previously,
with the command-line arguments `--epsilon`, `--lambda1`, `--tolerance`, etc. .

### <a name="geneset_file">Gene/Cell sets</a> ###

Gene or cell sets can be in **gmx** (Gene MatriX), **gmt** (Gene Matrix Transposed), or **grp** format.

The **gmt** format is convenient to store large databases of sets.
However, for a handful of sets, the **gmx** format might offer better
excel-editablity.

More information on these formats can be found
in the [Broad Institute Software Documentation](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats)

##### GMT #####

The **gmt** format consists of one set per line. Each line is a
tab-separated list composed as follows :

- The set name (can contain spaces)
- A commentary / description of the set (may be empty or contain spaces)
- A tab-separated list of set members

Example:

<table class="table" style="display: table">
<tr><td>Set1</td><td>The first set</td><td>gene_2</td><td>gene_1</td></tr>
<tr><td>Set2</td><td>The second set</td><td>gene_3</td></tr>
<tr><td>Set3</td><td>The third set</td><td>gene_4</td><td>gene_1</td></tr>
</table>

##### GMX #####

The **gmx** format is the transposed of the **gmx** format.
Each column represents a set. It is also tab-separated.

Example:

<table class="table" style="display: table">
<tr><td>Set1</td><td>Set2</td><td>Set3</td></tr>
<tr><td>The first set</td><td>The second set</td><td>The third set</td></tr>
<tr><td>gene_2</td><td>gene_3</td><td>gene_4</td></tr>
<tr><td>gene_1</td><td></td><td>gene_1</td></tr>
</table>


##### GRP #####

The **grp** format contains a single set in a simple newline-delimited text format. 

Example:

<table class="table" style="display: table">
<tr><td>gene_1</td></tr>
<tr><td>gene_2</td></tr>
<tr><td>gene_3</td></tr>
</table>


##### <a name="cells_by_gene_set">Cell selecting tool</a> #####

If you want to select a cell sets corresponding to a list of gene sets,
you may use the **cells_by_gene_set** command-line tool provided byt **wot**.

```sh
wot cells_by_gene_set --score gene_set_scores.txt \
 --out cell_sets.gmt --format gmt --quantile 99
```

You can select which proportion of the cells having each gene to select
with the `--quantile` option. The default value is 99, which would
select the top 1% of each gene. Choosing 50 for instance would
select every cell that has all genes above the median in the population.

### <a name="census_file">Census file</a> ###

Census files are datasets files : tab-separated text files with a header.
The header consists of an "id" field, and then the list of cell sets
for the census.

Each subsequent row will give the proportion of ancestors that
pertained in each of the mentionned cell sets.

The id is the time at which the ancestors lived.

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>tip1</td><td>tip2</td><td>tip3</td></tr>
<tr><td>0.0</td><td>0.15</td><td>0.05</td><td>0.05</td></tr>
<tr><td>1.0</td><td>0.28</td><td>0.05</td><td>0.03</td></tr>
<tr><td>2.0</td><td>0.42</td><td>0.03</td><td>0.02</td></tr>
<tr><td>3.0</td><td>0.72</td><td>0.02</td><td>0.01</td></tr>
<tr><td>4.0</td><td>0.89</td><td>0.00</td><td>0.00</td></tr>
<tr><td>5.0</td><td>0.99</td><td>0.00</td><td>0.00</td></tr>
</table>


### <a name="covariate_file">Covariate file</a> ###

The batch associated with each cell of the matrix file is specified in the *covariate* file.
This file must be a tab or comma separated plain text file, with two header fields: "id" and "covariate".

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>covariate</td></tr>
<tr><td>cell_1</td><td>0</td></tr>
<tr><td>cell_2</td><td>1</td></tr>
</table>


