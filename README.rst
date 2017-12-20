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
    git checkout develop
    cd wot
    pip install -e .
    export PYTHONPATH=$PYTHONPATH:.

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
Expression matrix - tab delimitted text file with genes on the columns and cells on the rows.

Example:
+----+--------+--------+---------+
| id | gene_1 | gene_2 | gene_n...
+----+--------+--------+---------+
| cell_1 | 41.2 | 12.2 | 3
+----+--------+--------+---------+
| cell_2 | 15 | 2 | 3.0
+----+--------+--------+---------+
| cell_n | 12.2 | 2 | 3
+----+--------+--------+---------+

Cell growth rates - Two column tab delimited file without header with cell ids and growth rates per day.

Example:
+----+--------+
| cell_1 | 3.1 |
+----+--------+
| cell_2 | 2.1 |
+----+--------+
| cell_n | 1.2 |
+----+--------+

Cell days - Two column tab delimited file without header with cell ids and days.

Example:
+----+--------+
| cell_1 | 1 |
+----+--------+
| cell_2 | 1 |
+----+--------+
| cell_n | 2 |
+----+--------+

Clusters - Two column tab delimited file without header with cell id and cluster id. Used to summarize transport maps by cluster.

Example:
+----+--------+
| cell_1 | 1 |
+----+--------+
| cell_2 | 1 |
+----+--------+
| cell_n | 2 |
+----+--------+

Day pairs - Two column tab delimited file without header with pairs of days to compute transport maps for.

Example:
+----+--------+
| 0 | 2 |
+----+--------+
| 2 | 4 |
+----+--------+
| 4 | 6 |
+----+--------+

