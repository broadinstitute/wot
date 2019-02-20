---
noheader: true
layout: documentation
location: Examples
---

# Plotting a validation summary
-------------------------------

## Computing validation summary ##

You can compute the validation summary for a set of constant parameters
from the command line :

```sh
wot optimal_transport_validation --matrix matrix.txt --cell_days days.txt \
 --covariate covariate.txt --out validation_summary.txt
```


```python
```python
{% include_relative code/08_plotting_validation_summary.py %}

```



## Result ##

![Validation summary plot]({{site.baseurl}}/images/validation_summary.png)
