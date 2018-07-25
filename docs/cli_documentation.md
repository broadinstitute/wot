---
layout: documentation
location: Documentation
---

Waddington Optimal Transport uses time-course data to infer how the
probability distribution of cells in gene-expression space evolves
over time, by using the mathematical approach of Optimal Transport (OT).


## Installation ##
------------------

### Dependencies ###

This packages depends on [Python 3](https://www.python.org/downloads/).

Several other python packages are required, but they can easily be installed through [pip][pip-install]

### Install from pip ###

The recommended and easiest way to install **wot** is via [pip][pip-install]

```sh
pip install --user wot
```

**wot** is then installed and ready to use.


### Install from source ###

Alternatively, you can install **wot** from source, and eventually modify it :

```sh
git clone https://github.com/broadinstitute/wot
cd wot && python setup.py sdist
```

The package gets built to `dist/wot-VERSION.tar.gz`, for instance `dist/wot-0.1.2.tar.gz`.

> <button class="btn-info rounded border-0 px-3 py-1" disabled>Dependencies</button>
>
> Once built, if you don't already satisfy all dependecies, you must install required packages :
>
> ```sh
> pip install --user h5py docutils msgpack-python --no-cache-dir
> pip install --user cython --no-cache-dir
> echo "$PATH" | grep -q "$HOME/.local/bin" || \
>   echo -e "\n# Cython path\nPATH=\"\$PATH:\$HOME/.local/bin\"" \
>   >> ~/.bash_profile
> ```
>
> *h5py*, *docutils*, and *msgpack-python* have to be installed separately and
> before cython because of [a known issue](https://github.com/h5py/h5py/issues/535)
> during the automatic installation through pip.

Then install the built package from the *dist/* directory :

```sh
pip install --user dist/wot-*.tar.gz
```

And then **wot** is installed and ready to use.

<hr />

> <button class="btn-info rounded border-0 px-3 py-1" disabled>Note</button>
>
> To improve random numbers generation performance, you might want to install [gslrandom](https://github.com/slinderman/gslrandom).
>
> ```sh
> pip install --user gslrandom
> ```


## Usage ##
-----------

**wot** is decomposed into several tools. Each tool can be used with the syntax `wot <tool>`.

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
> <div class="center-block text-center py-2"><a class="nounderline btn-outline-secondary btn-lg border px-4 py-2" role="button" href="#">Download .zip</a></div>

<br />

> <button class="btn-success rounded border-0 px-3 py-1" disabled>Interactive</button>
>
> Alternatively, **wot** features an interactive web interface to visualize
> your data and perform all the tasks described below.
>
> <div class="center-block text-center py-2">
>   <a class="nounderline btn-outline-secondary btn-md rounded border px-3 py-2"
>      role="button" href="{{site.baseurl}}/interactive_documentation">
>      Learn more &raquo;
>   </a>
> </div>
>

### Transport maps ###

```sh
wot optimal_transport --matrix matrix.txt \
 --cell_days days.txt --out tmaps --local_pca -1
```

This command will create a file `tmaps_{F}_{T}.loom` for each pair `{F}`, `{T}`
in the day pairs file. These maps can then be translated to any format you find
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
      <td>--day_pairs</td>
      <td>Target day pairs. See <a href="#day_pairs">formats</a></td>
    </tr>
    <tr>
      <td>--scaling_iter</td>
      <td>Number of iterations performed while scaling &epsilon;<br/>default : 3000</td>
    </tr>
    <tr>
      <td>--local_pca</td>
      <td>Number of dimensions to use when doing PCA<br/>default : 30</td>
    </tr>
  </tbody>
</table>

##### Scaling iterations #####

The convergence of the calculation of the transport maps can be
accelerated by gradually increasing the entropy regularization
parameter &epsilon;. Iteration count is thus split between those
that use a scaled &epsilon; (scaling iterations) and those that
use the final value of &epsilon; (extra iterations).

By default, **wot** performs 3000 scaling iterations and 1000 extra iterations.
While this is acceptable to ensure convergence in all cases when computing
couplings, it does require a lot of computing power. You may want to perform
fewer iterations to get an approximate result faster.

##### Local PCA #####

Dimensionality reduction is used when computing distances between cells.
The algorithm used for this purpose is Principal Component Analysis.
While using more dimensions for this purpose will make it more precise,
it will also slow the algorithm down. **wot** chooses by default
to use 30 dimensions.

Dimensionality reduction can be disabled with `--local_pca -1`


### Trajectories ###

Ancestors and descendants in **wot** are computed through the use of trajectories.

You can select a **cell set** by specifying a [cell set file](#cellset_file).
You can either manually edit this type of file, or generate it from a gene set file
using the [cells_by_gene_set](#cells_by_gene_set) tool.

```sh
wot trajectory --tmap . --cell_days days.txt \
 --cell_set cell_sets.gmt --out traj --progress
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
      <td>Directory containing the transport maps</td>
    </tr>
    <tr>
      <td><b>--cell_days</b></td>
      <td>Timestamps for each cell. See <a href="#days_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--cell_set</b></td>
      <td>Target cell set. See <a href="#cellset_file">formats</a></td>
    </tr>
    <tr>
      <td>--day_pairs</td>
      <td>Target day pairs. See <a href="#day_pairs">formats</a></td>
    </tr>
    <tr>
      <td>--out</td>
      <td>Output filenames prefix</td>
    </tr>
    <tr>
      <td>--progress</td>
      <td>Display a progress bar while performing the calculation</td>
    </tr>
  </tbody>
</table>


<a class="btn-info rounded border-0 px-3 py-1 btn-example nounderline"
 href="{{site.baseurl}}/examples/ancestor_census">See example code</a>
### Ancestor census ###

The census command lets you find out in which cell sets the ancestors
of a given cell set were located.


```sh
wot census --tmap . --cell_days days.txt \
 --cell_set cell_sets.gmt --matrix matrix.txt --progress
```
![Ancestor census plot]({{site.baseurl}}/images/ancestor_census.png)

This would create several census files named `<prefix>_<cellset>_<timepoint>.txt`,
for instance `census_tip1_100.0.txt`. See <a href="#census_file">formats</a>
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
      <td>Directory containing the transport maps</td>
    </tr>
    <tr>
      <td><b>--cell_days</b></td>
      <td>Timestamps for each cell. See <a href="#days_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--cell_set</b></td>
      <td>Target cell sets. See <a href="#cellset_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--matrix</b></td>
      <td>Normalized gene expression matrix. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td>--out</td>
      <td>Output filenames prefix.<br/>default : 'census'</td>
    </tr>
    <tr>
      <td>--progress</td>
      <td>Display a progress bar while performing the calculation</td>
    </tr>
  </tbody>
</table>

### Trajectory trends ###

Given **cell sets**, the mean value of different tips' ancestors at each time point will be calculated through trajectory trends.

You can select a **cell set** by specifying a [cell set file](#cellset_file).
You can either manually edit this type of file, or generate it from a gene set file
using the [cells_by_gene_set](#cells_by_gene_set) tool.

```sh
wot trajectory_trends --tmap . --cell_days days.txt --cell_set cell_sets.gmt --matrix matrix.txt
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
      <td>Directory containing the transport maps</td>
    </tr>
    <tr>
      <td><b>--cell_days</b></td>
      <td>Timestamps for each cell. See <a href="#days_file">formats</a></td>
    </tr>
    <tr>
      <td><b>--cell_set</b></td>
      <td>Target cell set. See <a href="#cellset_file">formats</a></td>
    </tr>
      <tr>
      <td><b>--matrix</b></td>
      <td>Normalized gene expression matrix. See <a href="#matrix_file">formats</a></td>
    </tr>
    <tr>
      <td>--matrix_transform</td>
      <td>Transform matrix values into certain type: choices={expm1,log1p,rank}
       </td>
    </tr>
  </tbody>
</table>


### Shared ancestry ###

### Trajectory differential expression ###
you can compare two ancestor distributions through local enrichment. 

The ancestor distributions can be two tips' or one tip's but at different time point. Now we have different ways to give the score that measures the difference between two distributions. Besides, the `matrix.txt` and `matrix1.txt`  should be the form of the result of trajectory trends.
```sh
wot optimal_local_enrichment --matrix1 matrix.txt \
(--matrix2 matrix2.txt) --score t_test \
--comparisons comapre.txt (--gsea C1)
```
When we run the cmd, we can get the file like `timepoint.rnk` or `timepoint1_timepoint2.rnk` including each gene 's score.

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
      <td>A matrix with cells on rows and features, such as
                        genes or pathways on columns See <a href="#matrix_file">formats</a></td>
    </tr> 
    <tr>
      <td><b>--score</b></td>
      <td>Method to compute differential gene expression score.
                        Choices are signal to noise, mean difference, t-test,
                        and fold change.{s2n,mean_difference,fold_change,t_test}
    </tr>
           <tr>
      <td>--matrix2</td>
      <td>A matrix with cells on rows and features, such as
                        genes or pathways on columns See <a href="#matrix_file">formats</a></td>
    </tr>
          <td>--gsea</td>
      <td>Run (<a href="http://s
                        oftware.broadinstitute">GSEA</a> on the specified MSigDB collections.<br/>
                         H (hallmark gene sets), C1 (positional gene sets),
                        C2 (curated gene sets), C3 (motif gene sets), C4
                        (computational gene sets), C5 (GO gene sets), C6
                        (oncogenic signatures), C7 (immunologic signatures)
      </td>
    </tr>
       </tr>
           <tr>
      <td>--comparisons</td>
      <td>Comparisons to generate ranked lists for. By default,
                        for one matrix signatures are created for all
                        consecutive timepoints. For two matrices for all
                        matching timepoints.</td>
    </tr>

  </tbody>
</table>

### Local regulatory model ###

### Global regulatory model ###

### Validation ###

You can easily validate the transport maps that have been computed above.

Say, for instance, that you have data at time points 0, 1, and 2.
You could compute the transport map from 0 to 2, and then interpolate with
optimal transport to see what distribution is expected at time 1 by this model.
It is then easy to compare with the actual distributions at those time points,
and thus check that the transport maps properly predict the cell's trajectory.

In order to do that, you would want a [day pairs file](#day_pairs_file)
containing the pair "0 2", tab-separated. You can then ask for interpolation
at time fraction 0.5, which would result in time point 1 here.

```sh
wot optimal_transport_validation --matrix matrix.txt \
 --cell_days days.txt --day_pairs day_pairs.txt \
 --out val_tmaps --t_interpolate .5 \
 --save_interpolated --local_pca -1 \
 --covariate covariate.txt --progress
```

This would create two files : `val_tmaps_I_1.0.txt`
and `val_tmaps_random_1.0.txt`.

They contain the coordinates of respectively the *interpolated* and the
*randomly generated* cells at that time point.

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
      <td><b>--t_interpolate</b></td>
      <td>Interpolation fraction between two timepoints</td>
    </tr>
    <tr>
      <td>--covariate</td>
      <td>Covariate value for each cell. See format below</td>
    </tr>
    <tr>
      <td>--save_interpolated</td>
      <td>Save interpolated population and the summary<br/>By default : save summary, discard populations</td>
    </tr>
    <tr>
      <td>--progress</td>
      <td>Display a progress bar while performing the calculation</td>
    </tr>
  </tbody>
</table>

The validation tool also accepts all options of the
[optimal_transport tool](#transport-maps). Mandatory options are
reproduced here for convenience.

##### Covariate #####

To enhance the results obtained when performing validation with few timepoints,
a covariate value may be assigned to each cell from each timepoint.

The interpolation can then occur for each pair of batches, and the distance
between the interpolated cells and the real population can be compared
to the distance between real batches at the interpolated time point, to
obtain more meaningful results.

The covariate values may be specified in a tab-separated text file.
It must have exactly two headers : "id" and "covariate".
Each subsequent line must consist of a cell name, a tab, and a covariate value.

> <button class="btn-warning rounded border-0 px-3 py-1" disabled>Warning</button>
>
> While a higher number of batches will greatly enhance the quality of the
> validation performed, the computation time required will increase as the
> square of the number of batches per timepoint.

#### Validation summary ####

Additionally, a validation summary is generated, as a text file.

Each line contains information about the relation between two cell sets :
 - **interval_start** and **interval_end** indicate which day pair is being considered.
 - **pair0** and **pair1** indicate which sets are being considered:
   - **P0** is the population of the dataset at the first timepoint of the day pair
   - **P0.5** is the population of the dataset at the interpolation timepoint
   - **I0.5** is the interpolated population at that timepoint
   - **R0.5** is the randomly generated population at that timepoint
   - **P1** is the population of the dataset at the second timepoint of the day pair
 - **distance** is the Wasserstein distance between the two sets considered


### Force-directed Layout Embedding ###

In order to visualize data in two dimensions, **wot** uses
[Force-directed Layout Embedding](https://en.wikipedia.org/wiki/Force-directed_graph_drawing).


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

The following formats are accepted by all tools: *mtx*, *hdf5*, *h5*, *h5ad*, *loom*, and *gct*.

##### Plain text gene expression matrix #####

Additionally, a plain text format is accepted. It must consist of tab-separated columns.

The first row, the header, must consist of an "id" field, and then the list of genes to be considered.

Each subsequent row will give the expression level of each gene for a given cell.

The first field must be a unique identifier for the cell, and then the tab-separated list
of expression level for each gene/feature.

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>gene_1</td><td>gene_2</td><td>gene_3</td></tr>
<tr><td>cell_1</td><td>1.2</td><td>12.2</td><td>5.4</td></tr>
<tr><td>cell_2</td><td>2.3</td><td>4.1</td><td>5.0</td></tr>
</table>

---

> <button class="btn-info rounded border-0 px-3 py-1" disabled>Note</button>
>
> You can convert input matrices from one format to another :
>
> ```sh
> wot convert_matrix file.mtx --format loom
> ```
>
> This command would create a file named *file.loom* containing the same
> matrix as *file.mtx*, in the same directory.
>
> Supported input formats are *mtx*, *hdf5*, *h5*, *h5ad*, *loom*, and *gct*.
>
> Supported output formats are *txt*, *txt.gz*, *loom*, and *gct*. (default: *loom*)



### <a name="days_file">Cell timestamps (day file)</a> ###

The timestamp associated with each cell of the matrix file is specified in the *days* file.
This file must be a tab-separated plain text file, with two header fields: "id" and "day".

Example:

<table class="table" style="display: table">
<tr><td>id</td><td>day</td></tr>
<tr><td>cell_1</td><td>1</td></tr>
<tr><td>cell_2</td><td>2.5</td></tr>
</table>

### <a name="day_pairs_file">Day pairs file</a> ###

The transitions to considered are specified in the *day\_pairs* file.
This file must be a tab-separated plain text file, without any header.

It indicates which pairs of days should be considered when computing transport maps.

Example:

<table class="table" style="display: table">
<tr><td>0</td><td>2</td></tr>
<tr><td>2</td><td>4</td></tr>
<tr><td>4</td><td>6</td></tr>
</table>

### <a name="geneset_file">Gene sets</a> ###

Gene sets can be in **gmx** (Gene MatriX), or **gmt** (Gene Matrix Transposed) format.

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


## More documentation ##
------------------------------

This document and the [examples]({{site.baseurl}}/examples) section should be more than enough to use **wot**.
However, if you feel the need for a more in-depth documentation about each of the python
functions in this package, it is available and can be generated from the sources of
the package with the [Sphinx](http://www.sphinx-doc.org/en/master/) tool :

```sh
pip install --user sphinx
cd sdocs/
make
```


[pip-install]: https://pip.pypa.io/en/stable/installing/
