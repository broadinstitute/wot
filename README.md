# WOT: Waddington-OT #

Waddington-OT uses time-course data to infer how the probability distribution
of cells in gene-expression space evolves over time, by using the mathematical
approach of Optimal Transport (OT).

## Install ##

Waddington-OT depends on [Python 3](https://www.python.org/downloads/).

Several other python packages are required, but they can easily be installed through [pip][pip-install]

### Dependencies ###

You can install all dependencies for **wot** with [conda](https://conda.io/docs/) :
```sh
conda install cython h5py flask gunicorn numexpr numpy pandas scikit-learn scipy simplejson psutil
conda install -c conda-forge pot
```

Or with [pip][pip-install] :
```sh
pip install --user h5py docutils msgpack-python
pip install --user cython
pip install --user flask gunicorn numexpr numpy pandas scikit-learn scipy psutil pot
```

### Install the **wot** package ###

```sh
pip install --user wot
```

## Usage ##

### Initializing an OT Model ###

**wot** uses an OTModel as its user interface. This is meant to provide black-boxed
functions to compute ancestors and descendants, while still using efficient caching
of transport maps under the hood, as these take a long time to compute. It also
allows **wot** to compute only the transport maps that are needed at any given time,
speeding up short pull-backs where only a few of those transport maps are needed.

You can initialize an OT Model in python with :

```python
ot_model = wot.initialize_ot_model('matrix.txt', 'days.txt')
```

All Optimal Transport parameters can be customized when initializing the model.
For instance, you could explicitely specify the defaults :

```python
ot_model = wot.initialize_ot_model('matrix.txt', 'days.txt', tmap_prefix='tmaps',
    epsilon=.05, lambda1=10, lambda2=50, batch_size=50, tolerance=1e-2)
```

This initialization steps store the configuration in a `tmaps.yml` file.
The basename of this file can be customized with the `tmap_prefix` option used above.

Each `tmap_prefix` references a different OT Model. If you plan to use several
models at the same time, you should a different prefix for each one of them.

### Loading an OT Model ###

Once the configuration file has been created, you can reload the previous
OT Model in a different script from the YAML configuration file :

```python
ot_model = wot.load_ot_model('matrix.txt', 'days.txt', 'tmaps')
```

Note that it is not necessary to use the '.yml' suffix when loading an OT
Model, the prefix is enough.

All previously computed transport maps will be available, even across calls
to the script.

### Changing parameters ###

If you want to change the OT parameters, you can simply re-initialize a model
with the new parameters. **wot** will check for compatibility between the cached
transport maps and the current configuration, and recompute the transport maps
that do not match the current configuration, so you need not care about these issues.

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
