---
noheader: true
title: Command Line
permalink: cli_documentation/
layout: documentation
location: Documentation
---


## Command Line ##
-----------

`wot` provides command line interface that offers fine-grained control over calculations. Each tool can be used with the syntax `wot <tool>`.

Help is available for each tool with `wot <tool> -h`. For example, `wot optimal_transport -h`


Type  `wot` to list all available tools.

<hr />


### wot gene_set_scores ###
Score each cell according to its expression of input gene signatures
<hr />

### wot optimal_transport ###

This command will create a file *tmaps_{A}_{B}.h5ad* for each pair *{A}, {B}*
of consecutive timepoints.


##### Local PCA #####

The default transport cost uses Principal Component Analysis to reduce the
dimension of the data before computing distances pairs of timepoints.
By default, optimal transport uses 30 dimensions.

Dimensionality reduction can be disabled with `--local_pca 0`

If you specify more dimensions for the PCA than your dataset has genes, wot will skip PCA and print a warning.
<hr />

### wot trajectory ###
Generate trajectories for cell sets generated at the given time. 
<hr />

### wot trajectory_trends ###
Generate mean expression profiles for ancestors and descendants of each
trajectory. Please use the trajectory tool to compute trajectories.

The command creates a file trends_<trajectory name>.txt with the mean expression
profile among ancestors/descendants.
<hr />


### wot fates ###
Generate fates for cell sets generated at the given time.
<hr />

### wot trajectory_divergence ###
Computes the distance between trajectories across time

### wot census ###
Computes the amount of mass of an
ancestor distribution falls into each cell set.

This command creates several census files named <prefix>_<cellset>.txt,
for instance census_tip1.txt. See <a href="#census_file">formats</a>
for more information.
<hr />

### wot diff_exp ###
Finds the genes that are differentially expressed between two sets of cells. You input one or more fates created using
the fates tool, an expression matrix, and the tool outputs a table with statistics about differentially expressed genes.
<hr />

### wot optimal_transport_validation ###
You can validate the transport maps that have been computed above.

Say, for instance, that you have data at time points 0, 1, and 2.
You could compute the transport map from 0 to 2, and then interpolate with
optimal transport to see what distribution is expected at time 1 by this model.
It is then easy to compare with the actual distributions at those time points,
and thus check that the transport maps properly predict the cell's trajectory.


##### Covariate #####

To measure the quality of interpolation, we compare the distance between
batches of cells at the same time point.

The covariate values may be specified in a tab-separated text file.
It must have exactly two headers : "id" and "covariate".
Each subsequent line must consist of a cell name, a tab, and a covariate value.

##### Null hypothesises #####

wot tests the results of the Optimal Transport predictions
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
<hr/>

### wot neighborhood_graph ###
One method of visualizing data in two or three dimensions is
[Force-directed Layout Embedding (FLE)](https://en.wikipedia.org/wiki/Force-directed_graph_drawing).
We compute this nearest neighbor graph in diffusion component space.
We recommend using [forceatlas2](https://github.com/klarman-cell-observatory/forceatlas2) or [Gephi](https://gephi.org/) to run the force directed layout.
You can generate a graph to use as input to the force directed layout using the command neighborhood_graph
<hr/>

### wot convert_matrix ###
Convert matrix data formats
<hr/>

## File formats ##
----------------------------

### <a name="matrix_file">Gene expression matrices</a> ###

The *matrix* file specifies the gene expression matrix to use.

The following formats are accepted by all tools: *mtx*, *txt*, *h5ad*, *loom*, and *gct*. Please note that wot expects 
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


### <a name="days_file">Cell Days</a> ###

The timestamp associated with each cell of the matrix file is specified in the *days* file.
This file must be a tab or comma separated plain text file, with two header fields: "id" and "day".

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>day</td></tr>
<tr><td>cell_1</td><td>1</td></tr>
<tr><td>cell_2</td><td>2.5</td></tr>
</table>

### <a name="geneset_file">Gene/Cell sets</a> ###

Gene or cell sets can be in **gmx** (Gene MatriX), **gmt** (Gene Matrix Transposed), or **grp** format.

The **gmt** format is convenient to store large databases of sets.
However, for a handful of sets, the **gmx** format might offer better
excel-editablity.

More information on these formats can be found [here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats)

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

### <a name="covariate_file">Covariate file</a> ###

The batch associated with each cell of the matrix file is specified in the *covariate* file.
This file must be a tab or comma separated plain text file, with two header fields: "id" and "covariate".

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>covariate</td></tr>
<tr><td>cell_1</td><td>0</td></tr>
<tr><td>cell_2</td><td>1</td></tr>
</table>


### OT Configuration file ###

There are several options to specify Optimal Transport parameters in wot.

The easiest is to just use constant parameters and specify them when
computing transport maps with the `--epsilon` or `--lambda1` options.

If you want more control over what parameters are used, you can use a
detailed configuration file. There are two kinds of configuration files
accepted by wot.

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




