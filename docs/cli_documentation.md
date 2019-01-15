---
---
layout: documentation
location: Documentation
---

Waddington Optimal Transport uses time-course data to infer how the
probability distribution of cells in gene-expression space evolves
over time, by using the mathematical approach of Optimal Transport (OT).

> <button class="btn-success rounded border-0 px-3 py-1" disabled>Interactive</button>
>
> Alternatively, **wot** features an interactive web interface to visualize
> your data and perform many of the tasks described below.
>
> <div class="center-block text-center py-2">
>   <a class="nounderline btn-outline-secondary btn-md rounded border px-3 py-2"
>      role="button" href="{{site.baseurl}}/interactive_documentation">
>      View documentation for web interface &raquo;
>   </a>
> </div>
>


## Installation ##
------------------

### Dependencies ###

This packages depends on [Python 3](https://www.python.org/downloads/).

Several other python packages are required, but they can easily be installed through [pip][pip-install] or [conda](https://www.anaconda.com/download/)

### Install from pip ###

The recommended and easiest way to install **wot** is via [pip][pip-install]

```sh
pip install --user wot
```

**wot** is then installed and ready to use.



## Usage ##
-----------

**wot** consists of several tools. Each tool can be used with the syntax `wot <tool>`.

Help is available for each tool with `wot <tool> -h`. For instance:

```
wot optimal_transport -h
```

In the following sections, each command is described with an example and
a table containing the core options. Required options are in bold font.

<hr />

> <button class="btn-info rounded border-0 px-3 py-1" disabled>Note</button>
>
> To help you get started with **wot**, we provide an archive containing all
> the simulated data that was used to compute all the pictures shown throughout
> this documentation. You can download it here and run the commands described
> thereafter in the same directory.
>
> <div class="center-block text-center py-2"><a class="nounderline btn-outline-secondary btn-lg border px-4 py-2" role="button" href="http://www.mediafire.com/file/01ahkjrsjjzt56l/all_simulated_data.zip">Download .zip</a></div>



### Transport maps ###

```sh
wot optimal_transport --matrix matrix.txt \
 --cell_days days.txt --out tmaps
```

This command will create a file `tmaps_{A}_{B}.loom` for each pair `{A}`, `{B}`
of consecutive timepoints. These maps can then be translated to any format you find
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
      <td>--batch_size</td>
      <td>Number of iterations performed between two duality_gap checks<br/>default : 50</td>
    </tr>
    <tr>
      <td>--tolerance</td>
      <td>Threshold on the ratio between duality gap and primal objective<br/>default : 0.01</td>
    </tr>
    <tr>
      <td>--local_pca</td>
      <td>Number of dimensions to use when doing PCA<br/>default : 30</td>
    </tr>
    <tr>
      <td>--config</td>
      <td>Configuration file for fine-graing control over OT parameters at each timepoint. See <a href="#ot-configuration-file">formats</a></td>
    </tr>
  </tbody>
</table>

##### Convergence #####

The convergence of the calculation of the transport maps is checked
through the use of the duality gap for the optimization problem.
When this gap reaches 0, the transport map is guaranteed to be optimal.
However, close-to-optimal transport maps are often acceptable and are
a lot faster to compute. **wot** lets you choose the tolerance you
want on the duality gap through the `--tolerance` parameter.

Checking the duality gap is computationally expensive, that is why
it is not done on every iteration, but rather after each batch of
iterations. You can change the batch size with the `--batch_size`
parameter. It is difficult to estimate the best batch size, as this
depends heavily on the data being considered. Nonetheless, batch sizes
too small (&lt; 10) or too big (&gt; 1000) are not recommended.


##### Local PCA #####

The default transport cost uses Principal Component Analysis to reduce the
dimension of the data before computing distances between cells.
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
 --cell_set cell_sets.gmt --time 10 --out traj.txt
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
      <td><b>--time</b></td>
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
 --out census --time 10
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
      <td><b>--time</b></td>
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


### Local regulatory model via differential expression ###

The local enrichment command finds the genes that are differentially expressed between two sets of cells.

The input matrices must have timepoints on rows, and genes on columns. This is the format created
by the [trajectory trends command](#trajectory-trends).

```sh
wot local_enrichment --score t_test \
 --matrix1 trends_set1.mean.txt --variance1 trends_set1.variance.txt \
 --matrix2 trends_set2.mean.txt --variance2 trends_set2.variance.txt
```

This will create files `<timepoint>.rnk` for each timepoint, containing each gene's score.

<table class="table table-hover" style="display: table">
  <thead class="thead-light">
    <tr>
      <th>Option</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>--matrix1</b></td>
      <td>A matrix with timepoints on rows and features, such as genes or pathways on columns See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--variance1</b></td>
      <td>A matrix with timepoints on rows and features on columns with the variance of each gene, associated with matrix1. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td>--matrix2</td>
      <td>A matrix with timepoints on rows and features, such as genes or pathways on columns See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td>--variance2</td>
      <td>A matrix with timepoints on rows and features on columns with the variance of each gene, associated with matrix2. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--score</b></td>
      <td>Method to compute differential expression score.<br/>
      Available:
      <ul>
      <li>signal to noise (s2n)</li>
      <li>mean difference (mean_difference)</li>
      <li>t-test (t_test)</li>
      <li>fold change (fold_change)</li>
      </ul></td>
    </tr>
    <tr>
      <td>--comparisons</td>
      <td>The timepoints to compare. See detailled description below</td>
    </tr>
  </tbody>
</table>

This commands accepts two types of input configurations. The first one is as presented in the example with two matrices.
It will compare the two matrices for all matching timepoints, and output a `<timepoint>.rnk` file for each of those.

Alternatively, you can specify a single matrix, and the default behavior will be to compare entries of the matrix
for all consecutives timepoints. This would create files names `<timepoint1>_<timepoint2>.rnk` instead.

For more control over which comparisons are performed, you can specify a file with a tab-separated pair
of timepoints on each line with the `--comparisons` parameter. The first of the two will refer to the
timepoint of the selected entry in the first matrix, and the second to the timpepoints of the selected entry in the
second matrix, which is identical to the first one if only one was specified. You will then get a file named
`<timepoint1>_<timepoint2>.rnk` for each pair of timepoints in the comparisons file.

### Validation ###

You can easily validate the transport maps that have been computed above.

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
    <tr>
      <td>--day_triplets</td>
      <td>Three column file without a header containing start time (t0), interpolation time, and end time (t1). Default is to compute transport maps from time t to time t+2
      and interpolate at time t+1</td>
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


### Force-directed Layout Embedding ###

In order to visualize data in two dimensions, **wot** uses
[Force-directed Layout Embedding (FLE)](https://en.wikipedia.org/wiki/Force-directed_graph_drawing). Intuitively, the FLE searches for a 2D layout of a nearest neighbor graph of the data. 
We compute this nearest neighbor graph in diffusion component space.

```sh
wot force_layout --matrix matrix.txt --out fdlayout
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
      <td>--neighbors</td>
      <td>Number of nearest neighbors to consider<br/>default : 100</td>
    </tr>
    <tr>
      <td>--neighbors_diff</td>
      <td>Number of nearest neighbors to use in diffusion component space<br/>default : 20</td>
    </tr>
    <tr>
      <td>--n_comps</td>
      <td>Number of diffusion components<br/>default : 20</td>
    </tr>
    <tr>
      <td>--n_steps</td>
      <td>Force-directed layout iteration count<br/>default : 1000</td>
    </tr>
    <tr>
      <td>--out</td>
      <td>Output filename prefix<br/>default : 'wot'</td>
    </tr>
  </tbody>
</table>

##### Output #####

This command will create two files containing the projection of the data
in two dimensions: `<prefix>.csv` and `<prefix>.h5ad`

The **h5ad** is just a regular 2D matrix that can be converted to any other
format with the *convert_matrix* tool. The **csv** is provided for convenience
and can be given directly to the *wot_server* interactive tool.


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

### <a name="geneset_file">Gene sets</a> ###

Gene sets can be in **gmx** (Gene MatriX), **gmt** (Gene Matrix Transposed), or **grp** format.

The **gmt** format is convenient to store large databases of gene sets.
However, for a handful of sets, the **gmx** format might offer better
excel-editablity.

More information on the gene set formats can be found
in the [Broad Institute Software Documentation](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats)

##### GMT #####

The **gmt** format consists of one gene set per line. Each line is a
tab-separated list composed as follows :

- The gene set name (can contain spaces)
- A commentary / description of the gene set (may be empty or contain spaces)
- A tab-separated list of genes

Example:

<table class="table" style="display: table">
<tr><td>Tip1</td><td>The first tip</td><td>gene_2</td><td>gene_1</td></tr>
<tr><td>Tip2</td><td>The second tip</td><td>gene_3</td></tr>
<tr><td>Tip3</td><td>The third tip</td><td>gene_4</td><td>gene_1</td></tr>
</table>

##### GMX #####

The **gmx** format is the transposed of the **gmx** format.
Each column represents a gene set. It is also tab-separated.

Example:

<table class="table" style="display: table">
<tr><td>Tip1</td><td>Tip2</td><td>Tip3</td></tr>
<tr><td>The first tip</td><td>The second tip</td><td>The third tip</td></tr>
<tr><td>gene_2</td><td>gene_3</td><td>gene_4</td></tr>
<tr><td>gene_1</td><td></td><td>gene_1</td></tr>
</table>


##### GRP #####

The **grp** format contains a single gene set in a simple newline-delimited text format. 

Example:

<table class="table" style="display: table">
<tr><td>gene_1</td></tr>
<tr><td>gene_2</td></tr>
<tr><td>gene_3</td></tr>
</table>

### <a name="cellset_file">Cell sets</a> ###

Cell sets can be described using the same formats as gene sets.
Simply list the ids of the cell in a set where you would have listed
the name of the genes in a gene set.

##### <a name="cells_by_gene_set">Cell selecting tool</a> #####

If you want to select a cell sets corresponding to a list of gene sets,
you may use the **cells_by_gene_set** command-line tool provided byt **wot**.

```sh
wot cells_by_gene_set --matrix matrix.txt --gene_sets gene_sets.gmt \
 --out cell_sets.gmt --format gmt --quantile 0.99
```

You can select which proportion of the cells having each gene to select
with the `--quantile` option. The default value is 0.99, which would
select the top 1% of each gene. Choosing 0.5 for instance would
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
This file must be a tab-separated plain text file, with two header fields: "id" and "covariate".

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>covariate</td></tr>
<tr><td>cell_1</td><td>0</td></tr>
<tr><td>cell_2</td><td>1</td></tr>
</table>

[pip-install]: https://pip.pypa.io/en/stable/installing/

