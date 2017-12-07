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

