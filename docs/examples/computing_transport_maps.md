---
noheader: true
layout: documentation
location: Examples
---

# Computing Transport maps
--------------------------

## Usage example ##

**wot** uses the Optimal Transport approach to infer relationships between cells.
It thus relies on the computation of transport maps between pairs of timepoints.

These transport maps can be computed from the command line with :

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt
```

Several options are available for this command, to choose the Optimal Transport parameters
that best suit your needs. Please refer to [the CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps)
for more information on those. You can also generate transport maps within python:

```python
{% include_relative code/01_create_tmaps.py %}

```

