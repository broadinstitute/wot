---
layout: documentation
location: Documentation
---

Waddington Optimal Transport uses time-course data to infer how the
probability distribution of cells in gene-expression space evolves
over time, by using the mathematical approach of Optimal Transport (OT).

## Installation ##
------------------

The recommended and easiest way to install **wot** is via [pip][pip-install]

```sh
pip install --user wot
```

## Setup ##
-----------

Before the interactive web interface for **wot** is usable, some data has to
be precomputed.

> <button class="btn-info rounded border-0 px-3 py-1" disabled>Note</button>
>
> To help you get started with **wot**, we provide an archive containing all
> the simulated data that was used to compute all the pictures shown throughout
> this documentation. You can download it here and run the commands described
> thereafter in the same directory.
>
> <div class="center-block text-center py-2"><a class="nounderline btn-outline-secondary btn-lg border px-4 py-2" role="button" href="#">Download .zip</a></div>

##### Force-directed layout embedding #####

```sh
wot force_layout --matrix matrix.txt --out fdlayout
```

This will create two files: `fdlayout.csv` and `fdlayout.h5ad`,
containing the data required to plot the datasets in two dimensions.

##### Transport maps #####

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt \
 --out tmaps --local_pca -1
```

This will create several files describin the transport maps
between consecutive time points, that can later be used
by the visualization tool.


##### Interactive web server #####

```sh
wot wot_server --matrix matrix.txt --tmaps . \
 --cell_meta fdlayout.csv --cell_meta days.txt
```

This will launch the server of port 8080 of your localhost by default.
To use the web interface, you can then open a web browser and navigate to :

> <http://127.0.0.1:8080/web/index.html>

[pip-install]: https://pip.pypa.io/en/stable/installing/
