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

You can generate this kind of simulation, along with days and covariates data
for future use, with the following script :

```python
{% include_relative code/00_generating_data.py %}

```

 
## Plotting two features ##

One of the best ways to visualize your dataset is to use Force-directed Layout Embedding
to project your data into 2 dimensions.

However, you might also want to just plot 2 specific features of you dataset.
You can easily add colors to this kind of plot providing an array of colors to plot. For instance, here is how to color according to time :

```python
{% include_relative code/02_plotting_two_features.py %}

```

## Result ##

![generated data plot]({{site.baseurl}}/images/generated_pop.png)

