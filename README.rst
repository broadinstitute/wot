============
WOT
============


Uses time-course data to infer how the probability distribution of cells in gene-expression space evolves over time,
by using the mathematical approach of Optimal Transport (OT)

========================
Getting Started
========================
.. code-block:: bash

    git clone https://github.com/broadinstitute/wot.git
    cd wot
    git checkout develop
    pip install -e .

========================
Command Line Usage
========================

Compute transport maps

.. code-block:: bash

    python bin/ot.py \
       --matrix my_expression_matrix.txt \
       --cell_growth_rates my_growth_rates.txt \
       --cell_days my_cell_days.txt \
       --day_pairs my_day_pairs.txt \
       --prefix my_transport_maps \

Summarize by cluster

.. code-block:: bash

    python bin/summarize_by_cluster.py \
       --dir my_transport_maps \
       --clusters my_clusters.txt \
       --prefix my_summary

Compute trajectories

.. code-block:: bash

    python bin/trajectory.py \
       --dir my_transport_maps_dir \
       --time 9 \
       --prefix my_trajectory \
       --id my_cell_ids.txt

========================
Help
========================
Pass --help to any of the commands for a full description of all command line arguments. For example

.. code-block:: bash

    python bin/ot.py --help

========================
File Formats
========================
Expression matrix - tab delimitted text file with genes on the columns and cells on the rows, loom file, or mtx file.

Example:
id TAB gene_1 TAB gene_2 TAB gene_n...

cell_1 TAB 41.2 TAB 12.2 TAB 3

cell_2 TAB 15 TAB 2 TAB 3.0

cell_n TAB 12.2 TAB 2 TAB 3



Cell growth rates - Two column tab delimited file without header with cell ids and growth rates per day.

Example:

cell_1 TAB 3.1

cell_2 TAB 2.1

cell_n TAB 1.2


Cell days - Two column tab delimited file without header with cell ids and days.

Example:

cell_1 TAB 1

cell_2 TAB 1

cell_n TAB 2


Day pairs - Two column tab delimited file without header with pairs of days to compute transport maps for.

Example:

0 TAB 2

2 TAB 4

4 TAB 6


Clusters - Two column tab delimited file without header with cell id and cluster id. Used to summarize transport maps by cluster.

Example:

cell_1 TAB 1

cell_2 TAB 1

cell_n TAB 2
