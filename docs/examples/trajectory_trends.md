---
noheader: true
layout: documentation
location: Examples
---

# Plot Trajectory Trends
------------------------

If we don't have real data and cell sets, we can generate the simulate data and cell sets following the steps in [generate data](generate_data) and [plotting_cell_sets](plotting_cell_sets).

## Calculate the Transport Maps

We can get the transport maps through [optimal_transport]({{site.baseurl}}/cli_documentation#transport-maps).

```sh
wot optimal_transport --matrix matrix.txt \
 --cell_days days.txt --out tmaps 
```

## Calculate the Trajectories
We can get the trajectories through
[trajectory]({{site.baseurl}}/cli_documentation#trajectory).
```sh
wot trajectory --tmap . --cell_set cell_sets.gmt --day 7
```

## Calculate Trajectory Trends
Now we can get the trajectory trends through
[trajectory_trends]({{site.baseurl}}/cli_documentation#trajectory-trends).
```sh
wot trajectory_trends --trajectory trajectory.txt --cell_days days.txt --matrix matrix.txt
```

## Plot Trajectory Trends
The result of trajectory trends are mean value of one tip's ancestors and descendants. Specifically, given time point, and the tips in this time point, we will find each tip's ancestors and descendants, and calculate their mean value of each gene expression rate.



```python
{% include_relative code/06_plotting_trajectory_trends.py %}

```



## The Plot





![Trajectory Trends]({{site.baseurl}}/images/trajectory_trends.png)

