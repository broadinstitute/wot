# Transport Maps of different $\epsilon$

Here we give an example to show how the size of $\epsilon$ influence transport maps.



## Regenerate 1-D Data

We need the matrix.txt from Generate Data to generate one dimensional data. In fact, we just choose the first column of the original data.

~~~python
# ------ Configuration variables -------
matrix_file = 'matrix.txt'
new_matrix_file = 'matrix1.txt'
# --------------------------------------

import wot.io
import pandas as pd

# Get the first column of matrix
ds = wot.io.read_dataset(matrix_file)
df = pd.DataFrame({'x':ds.x[:,0]},index=ds.row_meta.index.values)
df = df.sort_values(by=['x'])
dg = pd.DataFrame([],index=list(df.index))
dk = wot.Dataset(df.values,row_meta=dg,col_meta=pd.DataFrame([],index=['x']))

# Save data
wot.io.write_dataset(dk,new_matrix_file)
~~~



# Plot the Map 

Because of one dimensional data, we can get good plots. Here is an example with $\epsilon=0.05$

~~~python
# ------ Configuration variables -------
matrix_file = 'matrix1.txt'
days_file = 'days.txt'
bg_color = "#80808080"
gene_x_plot = 0
gene_y_plot = 1
cell_sets_file = 'cell_sets.gmt'
target_cell_set = "Red blood cells"
target_timepoint = 50
destination_file = "ancestors.png"
# --------------------------------------

import numpy as np
import wot.io
import matplotlib.pyplot as plt
from matplotlib import cm

# Calculate the Transport Map
ds = wot.io.read_dataset(matrix_file)

core = wot.initialize_core(matrix_file, days_file,
                           scaling_iter=300, epsilon=.05, lambda1=50)

tmap = core.transport_map(46,47)
data = tmap.x
data=data[0:4900,0:4900]
# print(data)

# Amplify the Signal
for i in range(490):
    for j in range(490):
        data[i * 10:(i + 1) *10, j * 10:(j + 1) * 10]\
            =np.average(data[i*10:(i+1)*10,j*10:(j+1)*10])

# Give The plot
fig = plt.figure()
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
data=data/np.sum(data)
data=np.power(data,0.7)
line1 = ax.imshow(data,cmap=cm.BuGn)
ax.set_title('Transport Maps from Day 46 to 47 (epislon=0.05)')
ax.set_xticks([1000,2500,4000])
ax.set_yticks([1000,2500,4000])
ax.set_xticklabels(['Tip 1','Tip 2',' Tip 3'])
ax.set_yticklabels(['Tip 1','Tip 2',' Tip 3'])
plt.colorbar(line1)
plt.savefig('tmaps_46_47_05.png')
plt.show()
~~~



## Results 

Here are some results.

![Trajectory Trends for Tip 1]({{site.baseurl}}/images/tmaps_46_47_1.png)

![Trajectory Trends for Tip 1]({{site.baseurl}}/images/tmaps_46_47_05.png)

![Trajectory Trends for Tip 1]({{site.baseurl}}/images/tmaps_46_47_01.png)

![Trajectory Trends for Tip 1]({{site.baseurl}}/images/tmaps_46_47_001.png)