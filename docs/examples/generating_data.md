---
noheader: true
layout: documentation
location: Examples
---

# Generating simulated data
---------------------------

## Gaussians around a path ##

One way to generate simulated data is to pick a few "paths" in gene-expression
space, defined by piecewise linear interpolation between a set of fixed points,
and then generate cells for those paths by simply adding Gaussian noise to the
paths.

You can generate this kind of simulation along with cell days, cell batches, and a 2-d embedding of cells with the following script :

```python
{% include_relative code/00_generating_data.py %}

```



