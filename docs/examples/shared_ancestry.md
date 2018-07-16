# Shared Ancestry

We can get a plot of different cells' ancestors through shared ancestry.  Before we get the ancestors, we need cell sets only on the day 100. Assume we already have cell cets, let's filter cells not on day 100.



## Filter Cells

We filter cells not on day 100 firstly.

```python
import wot.io

cell_set = wot.io.read_cell_sets('cell_sets.gmt')
# print(cell_set)
value_tip1=cell_set['Tip1']
value_tip2=cell_set['Tip2']
value_tip3=cell_set['Tip3']

value_tip1_new=[]
value_tip2_new=[]
value_tip3_new=[]

for val in value_tip1:
    print(val[6])
    if val[6]=='1':
        value_tip1_new.append(val)

for val in value_tip2:
    print(val[6])
    if val[6]=='1':
        value_tip2_new.append(val)

for val in value_tip3:
    print(val[6])
    if val[6]=='1':
        value_tip3_new.append(val)

cell_set['Tip1']=value_tip1_new
cell_set['Tip2']=value_tip2_new
cell_set['Tip3']=value_tip3_new

wot.io.write_gene_sets(cell_set,'cell_sets_new.gmt',format='gmt')
```



## Find the Trajectory

We use the cmd command to find new cell sets' trajectory.

```
wot trajectory --tmap . --cell_set cell_sets_new.gmt --cell_days days.txt --out traj --progress
```

p.s. We must have the transport maps already.



## Plot the Shared Ancestry

Here we use functions from **wot** to get the plot.

```python
import numpy as np
import wot.io
import wot.graphics
from matplotlib import pyplot

ds = wot.io.read_dataset('matrix.txt')
bgcolor = "#80808020"
wot.io.incorporate_days_information_in_dataset(ds, 'days.txt')
wot.set_cell_metadata(ds, 'color', bgcolor)

mix_data = wot.io.read_dataset('traj.txt').x




for i in range(101):
    for j in range(1000):
        weight = mix_data[1000*(100-i)+j]*10000
        if np.sum(weight)==0:
            continue
        weight = weight/np.sum(weight)
        color  = (weight[0],weight[2],weight[1],1)
        color  = wot.graphics.hexstring_of_rgba(color)
        indice = 1000*i+j
        wot.set_cell_metadata(ds, 'color',
                              color, indices=indice)



pyplot.axis('off')
wot.graphics.plot_2d_dataset(pyplot, ds, x=0, y=1)
pyplot.show()
```



## Results

Here is one of the plot.

![Shared Ancestry]({{site.baseurl}}/images/shared_ancestry.png)

