---
layout: documentation
---

Waddington Optimal Transport uses time-course data to infer how the
probability distribution of cells in gene-expression space evolves
over time, by using the mathematical approach of Optimal Transport (OT).


## Dependencies ##
------------------

This packages depends on [Python 3](https://www.python.org/downloads/).

Several other python packages are required, but they can easily be installed through *pip*

## Installation ##
------------------

Clone the repo and switch to develop branch :

```sh
git clone https://github.com/broadinstitute/wot
cd wot
git checkout develop
```

Build the package for pip :

```sh
python setup.py sdist
```

It gets built to `dist/wot-VERSION.tar.gz`, for instance `dist/wot-0.1.2.tar.gz`.

Install with :

```sh
pip install --user h5py docutils msgpack-python --no-cache-dir
pip install --user cython --no-cache-dir
echo -e "\n# Required by cython\nPATH=\"\$PATH:\$HOME/.local/bin\"" >> ~/.bashrc
pip install --user dist/wot-*.tar.gz
```

*h5py*, *docutils*, and *msgpack-python* have to be installed separately and
before cython because of [a known issue](https://github.com/h5py/h5py/issues/535)
during the automatic installation through pip.

---

*Note :*
This package depends on recent versions of several other packages.
You might want to upgrade all your python packages first, with :

```sh
pip freeze --local | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U --user
```

*Note :*
To improve random numbers generation performance, you might want to install [gslrandom](https://github.com/slinderman/gslrandom).

```sh
pip install --user gslrandom
```


## Usage ##
-----------

WOT includes several tools. Each tool can be used with the syntax `wot <tool>`.

Help is available for each tool with `wot <tool> -h`. For instance `wot optimal_transport -h`.

### Converting matrices ###

You can convert input matrices from one format to another :

```sh
wot convert_matrix file.mtx --format loom
```

This will create a file named *file.loom* containing the same matrix as *file.mtx*, in the same directory.

Supported input formats are *mtx*, *hdf5*, *h5*, *h5ad*, *loom*, and *gct*.

Supported output formats are *txt*, *txt.gz*, *loom*, and *gct*. (default: *loom*)


You can compute the gene set scores for each of the cells in a gene expression matrix :

```sh
wot gene_set_scores --matrix <matrixfile>
```

For instance :

```sh
wot gene_set_scores --matrix matrix.txt --out matrix_gss.txt
```

If not explicitely specified, the output file will default to the basename of the input
file with "_gene_set_scores" appended to it, and will be stored in your current directory.

### Computing transport maps ###

**TODO:** Add the generation tools to get the simulated matrices

To compute the transport maps, you need to specify the [matrix file](#matrix_file) to use, and the  [days file](#days_file) :

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt --out tmaps --local_pca -1
```

**TODO:** Change options description to a more readable table

#### Restricting to certain day pairs ####

By default, transport maps will be computed for all consecutive days.

You may want to restrict this with a [day pairs file](#day_pairs_file).
You can then use the `--day_pairs <file>` option to point to your file.

#### Scaling iterations ####

By default, WOT performs 3000 scaling iterations when calculating transport maps
to ensure convergence when calculating the couplings. While this is an acceptable
value to ensure convergence in all kinds of situations, it does require a lot of
computing power. You may want to perform fewer iterations to get an approximate
result faster.

You can specify a custom value with `--scaling_iter <n>`, for instance `--scaling_iter 150`.

#### PCA reduction ####

Principal component analysis is used to reduce the dimensionality of the expression matrix locally for consecutive days.

*Note :*
PCA can be disabled as a whole by passing a negative argument with `--local-pca`.
When plotting, bear in mind that this will change the space in which data is located.


### Validating transport maps ###

You can easily validate the transport maps that have been computing above.

Say, for instance, that you have data at time points 0, 1, 2, 3, and 4.
You could compute the maps from 0 to 2 and 2 to 4, and then interpolate with
optimal transport to see what distributions are expected at time 1 and 3.
It is then easy to compare with the actual distributions at those time points,
and thus check that the transport maps properly predict the cell's trajectory.

In order to do that, you would want a [day pairs file](#day_pairs_file) containing
the two pairs "0 2" and "2 4", tab-separated, each on a different line. You can
then ask for interpolation at time 0.5, which would result in the time points 1 and 3.

```sh
wot optimal_transport_validation --matrix matrix.txt --cell_days days.txt --out val_tmaps --t_interpolate .5 --save_interpolated --day_pairs day_pairs.txt --local_pca -1
```

This would create several files `val_tmaps_I_{t}.txt` and `val_tmaps_random_{t}.txt`
where `{t}` takes values 1 and 3.

They contain the coordinates of respectively the *interpolated* and the *randomly generated* points at that time.

#### Validation summary ####

Additionally, a validation summary is generated.

Each line contains information about the relation between two cell sets :
 - *interval_start* and *interval_end* indicate which day pair is being considered.
 - *pair0* and *pair1* indicate which sets are being considered:
   - *P0* is the population of the dataset at the first timepoint of the day pair
   - *P0.5* is the population of the dataset at the interpolation timepoint
   - *I0.5* is the interpolated population at that timepoint
   - *R0.5* is the randomly generated population at that timepoint
   - *P1* is the population of the dataset at the second timepoint of the day pair
 - *distance* is the Wasserstein distance between the two sets considered


### Computing transition tables ###

**TODO:** requires transport maps from one cell set to another. Builds transition tables from these transport maps


## Generating documentation ##
------------------------------

The documentation for each python function can be generated with the [Sphinx](http://www.sphinx-doc.org/en/master/) tool :

```sh
pip install --user sphinx
cd docs/
make
```

## Supported file formats ##
----------------------------

### <a name="matrix_file">Gene expression matrices</a> ###

The *matrix* file specifies the gene expression matrix to use.

The following formats are accepted by all tools: *mtx*, *hdf5*, *h5*, *h5ad*, *loom*, and *gct*.

#### Plain text gene expression matrix ####

Additionally, a plain text format is accepted. It must consist of tab-separated columns.

The first row, the header, must consist of an "id" field, and then the list of genes to be considered.

Each subsequent row will give the expression level of each gene for a given cell.

The first field must be a unique identifier for the cell, and then the tab-separated list
of expression level for each gene/feature.

Example:

<table>
<tr><td>id</td><td>gene_1</td><td>gene_2</td><td>gene_3</td></tr>
<tr><td>cell_1</td><td>1.2</td><td>12.2</td><td>5.4</td></tr>
<tr><td>cell_2</td><td>2.3</td><td>4.1</td><td>5.0</td></tr>
</table>

### <a name="days_file">Cell timestamps</a> ###

The timestamp associated with each cell of the matrix file is specified in the *days* file.
This file must be a tab-separated plain text file, with two header fields: "id" and "day".

Example:

<table>
<tr><td>id</td><td>day</td></tr>
<tr><td>cell_1</td><td>1</td></tr>
<tr><td>cell_2</td><td>2.5</td></tr>
</table>

### <a name="day_pairs_file">Target transitions</a> ###

The transitions to considered are specified in the *day\_pairs* file.
This file must be a tab-separated plain text file, without any header.

It indicates which pairs of days should be considered when computing transport maps.

Example:

<table>
<tr><td>0</td><td>2</td></tr>
<tr><td>2</td><td>4</td></tr>
<tr><td>4</td><td>6</td></tr>
</table>
