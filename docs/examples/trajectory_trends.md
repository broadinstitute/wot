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
 --cell_days days.txt --out tmaps --local_pca -1
```

## Calculate the Trajecotory Trends
Now we can get the trajectory trends through
[trajectory_trends]({{site.baseurl}}/cli_documentation#trajectory-trends).
```sh
wot trajectory_trends --tmap . --cell_days days.txt --cell_set cell_sets.gmt --matrix matrix.txt
```

## Plot Trajectory Trends
The result of trajectory trends are mean value of one tip's ancestors and descendants. Specifically, given time point, and the tips in this time point, we will find each tip's ancestors and descendants, and calculate their mean value of each gene expression rate.

So, we can plot the data, here's an example to plot the 
Tip1's gene expression rate (cells of Tip1 from day 100) for gene 2,3,4 at different time point.

```python
import wot.io
import matplotlib
import numpy as np
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


data_tip1 = wot.io.read_dataset("matrix_Tip1_100.0_trajectory_trends.loom")

a=data_tip1.x
b=data_tip1.row_meta
print(a)

x=np.linspace(0,100,101)

fig,ax=plt.subplots()

data = data_tip1.x
plt.plot(data[:,2],c='#E00000',label='Gene 2',linewidth=1)
plt.plot(data[:,3],c='#00E000',label='Gene 3',linewidth=1)
plt.plot(data[:,4],c='#0000E0',label='Gene 4',linewidth=1)
plt.xticks(fontsize=6)
plt.yticks(fontsize=6)
plt.xlabel("Time point",fontsize=10)
plt.ylabel("Gene expression",fontsize=10)
plt.title(" Trajectory trends from cells of tip1 at day 100",fontsize=10)
plt.legend(loc='best',edgecolor='black',fontsize=10)
plt.savefig("Trjectory_trends_Tip1_100.png")
plt.show()
plt.close()
```



## The Plot

Here are some examples of plot we get.



![Trajectory Trends for Tip 1]({{site.baseurl}}/images/trajectory_trends_1.png)
![Trajectory Trends for Tip 2]({{site.baseurl}}/images/trajectory_trends_2.png)
![Trajectory Trends for Tip 3]({{site.baseurl}}/images/trajectory_trends_3.png)
