---
noheader: true
permalink: notebook/
layout: documentation
location: Examples
---

# Using wot from Python
-----------------------

This notebook provides a running example with simulated data to walk
you through the most important tools of **wot**.

You can read the code and the detailed comments below, and simply
copy-paste it into a python interpreter, or download the whole
python script for a specific section by clicking the download
button next to the section title.

The examples below will show you how to generate your own simulated
data using **wot**'s simulate module, and compute transport maps
for it to generate your own plots.

Alternatively, we also provide simulated data with all transport maps
precomputed, and an archive with all of the scripts presented in this
notebook. Just click the button below to download the whole archive.

<div class="center-block text-center py-2">
  <a class="nounderline btn-outline-secondary btn-lg border px-3 py-2 mx-3"
     role="button" href="#">Download all data</a>
  <a class="nounderline btn-outline-secondary btn-lg border px-3 py-2 mx-3"
     role="button" href="#">Download all scripts</a>
</div>

<br />

<a class="btn-info rounded border-0 px-3 py-1 btn-example nounderline"
 href="01_generating_data.py">Download python code</a>
## Simulating data ##
---------------------

```python
{% include_relative notebook/01_generating_data.py %}
```

<a class="btn-info rounded border-0 px-3 py-1 btn-example nounderline"
 href="02_plotting_two_features.py">Download python code</a>
## Plotting two features ##
---------------------------

```python
{% include_relative notebook/02_plotting_two_features.py %}
```

![Generated data](images/notebook_generated_data.png)


<a class="btn-info rounded border-0 px-3 py-1 btn-example nounderline"
 href="03_plotting_cell_sets.py">Download python code</a>
## Cell sets ##
---------------

```python
{% include_relative notebook/03_plotting_cell_sets.py %}
```

![Cell sets plots](images/notebook_cell_sets.png)


<a class="btn-info rounded border-0 px-3 py-1 btn-example nounderline"
 href="04_plotting_ancestors.py">Download python code</a>
## Ancestors of a cell set ##
-----------------------------

```python
{% include_relative notebook/04_plotting_ancestors.py %}
```

![Ancestors plot](images/notebook_ancestors.png)
