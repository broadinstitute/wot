---
noheader: true
layout: documentation
location: Examples
---

# Computing Transport maps
--------------------------

## Usage example ##

**wot** uses the Optimal Transport approach to infer relationships between cells.
It thus relies on the computation of transport maps between consecutive timepoints.

These transport maps can be computed from the command line with :

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt
```

Several options are available for this command, to choose the Optimal Transport parameters
that best suit your needs. Please refer to [the CLI documentation]({{site.baseurl}}/cli_documentation#transport-maps)
for more information on those.

Here is another example, specifying the default values explicitely :

```sh
wot optimal_transport --matrix matrix.txt --cell_days days.txt \
 --out tmaps --local_pca 30 --max_iter 1000000 \
 --epsilon 0.05 --lambda1 1 --lambda2 50 \
 --batch_size 50 --tolerance 0.01
```

## Influence of epsilon ##

Here is an example to show the influence of epsilon on the computed
transport map. The lower the epsilon, the closer the transport map is
to a one-on-one mapping of cells.

Cells have been sorted by increasing value of the first coordinate for these plots.

#### Computing the transport map ####

```python
import wot.ot

matrix_file = 'matrix.txt'
days_file = 'days.txt'

ds = wot.io.read_dataset(matrix_file)
ot_model = wot.ot.initialize_ot_model(matrix_file, days_file,
    epsilon=.05, lambda1=50)

data = ot_model.compute_transport_map(46, 47).X[:4900, :4900]
```

#### Amplifying the signal ####

We extracted 4900 points, but the points are still hard to see on the plot,
so we use a running average to amplify the importance of each non-null
coefficient of the transport map.

```python
import numpy as np

for i in range(490):
    for j in range(490):
        data[i * 10:(i + 1) *10, j * 10:(j + 1) * 10] \
            = np.average(data[i*10:(i+1)*10, j*10:(j+1)*10])
```

#### Plotting ####

We can then repeat this for different values of epsilon, and plot each one :

```python
from matplotlib import cm
import matplotlib.pyplot as plt

fig = plt.figure()
ax  = fig.add_axes([0.1,0.1,0.8,0.8])
data=data/np.sum(data)
line1 = ax.imshow(data,cmap=cm.BuGn)
ax.set_title('Transport Maps from Day 46 to 47 (epsilon=0.05)')
ax.set_xticks([1000, 2500, 4000])
ax.set_yticks([1000, 2500, 4000])
ax.set_xticklabels(['Tip 1', 'Tip 2', 'Tip 3'])
ax.set_yticklabels(['Tip 1', 'Tip 2', 'Tip 3'])
plt.colorbar(line1)
plt.savefig('tmaps_46_47_05.png')
plt.show()
```

#### Results  ####

![Transport map for epsilon .1  ]({{site.baseurl}}/images/tmaps_46_47_1.png)
![Transport map for epsilon .05 ]({{site.baseurl}}/images/tmaps_46_47_05.png)
![Transport map for epsilon .01 ]({{site.baseurl}}/images/tmaps_46_47_01.png)
![Transport map for epsilon .001]({{site.baseurl}}/images/tmaps_46_47_001.png)
