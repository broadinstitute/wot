# WOT: Waddington-OT #

Waddington-OT uses time-course data to infer how the probability distribution
of cells in gene-expression space evolves over time, by using the mathematical
approach of Optimal Transport (OT).

## Install ##

Waddington-OT depends on [Python 3](https://www.python.org/downloads/).


### Dependencies ###

You can install dependencies for **wot** with [conda](https://conda.io/docs/) :
```sh
conda install numpy pandas h5py cython scikit-learn scipy matplotlib
conda install -c conda-forge pot
```


### Install the **wot** package ###

```sh
pip install --user wot
```

## Usage ##

### Initializing an OT Model ###

**wot** uses an OTModel as its interface for computing transport maps. 

You can initialize an OT Model in python with :

```python
ot_model = wot.ot.initialize_ot_model('matrix.txt', 'days.txt')
```

All Optimal Transport parameters can be customized when initializing the model.
For instance, you could explicitely specify the defaults :

```python
ot_model = wot.ot.initialize_ot_model('matrix.txt', 'days.txt', tmap_prefix='tmaps',
    epsilon=.05, lambda1=10, lambda2=50, batch_size=50, tolerance=1e-2)
```

You can compute all transport maps with :

```python
ot_model.compute_all_transport_maps()
```

### Loading Transport Maps ###

Once the transport maps have been created, you can operate on the transport maps using the TransportMapModel interface :

```python
tmap_model = wot.tmap.TransportMapModel.from_directory('.')
```


All previously computed transport maps will be available.

### Changing parameters ###


If you want to keep the previously computed transport maps, simply initialize
a new model with a different prefix. Any model will only affect files that use
its `tmap_prefix`, there is no interaction between models with different prefixes.

### Using wot.commands ###

All data-processing functions are located in the `wot.commands` subpackage.
These include :

- ancestor census
- gene set scores
- gene regulatory networks (grn)
- local enrichment
- optimal transport validation
- trajectory
- trajectory trends
- convert matrix
- force layout
- wot server (interactive version of **wot**)

All of these are documented on [wot's github pages website](http://broadinstitute.github.io/wot), with examples using simulated data to show how to use and plot the results of these commands.

## Documentation ##

The full documentation for **wot** is available on Github Pages : <http://broadinstitute.github.io/wot>

For more advanced usage, you may also browse the source code to read each
function's documentation. Most of **wot**'s internal functions have docstrings
with a description of their parameters, output and examples on how to use them.

## Developer Notes ##

For more information about the internal functionning of **wot**, please refer
to the [Developer Notes](developer_notes.md)

[pip-install]: https://pip.pypa.io/en/stable/installing/
