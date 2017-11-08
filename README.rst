============
WOT
============


Uses time-course data to infer how the probability distribution of cells in gene-expression space evolves over time,
by using the mathematical approach of Optimal Transport (OT)

========================
Command Line Usage
========================

Compute transport maps

.. code-block:: bash

    python bin/ot
       --expression_file paper/serum_free_dmap_20.txt
       --growth_file paper/growth_scores.txt
       --days_file paper/days.txt
       --prefix my_output_file_name_prefix
       --min_transport_fraction 0.05
       --max_transport_fraction 0.4
       --min_growth_fit 0.9
       --l0_max 100
       --lambd1 1
       --lambda2 1
       --epsilon 0.1
       --growth_ratio 2.5
       --scaling_iter 250

Compute trajectories

.. code-block:: bash

    python bin/trajectory
       --dir paper/transport_maps/2i/
       --id day-9-c1-2i_6 --id day-9-c1-2i_11
       --time 9
       --prefix my_trajectory

============
API Usage
============

.. code-block:: python

    import numpy as np
    import scipy.stats
    import sklearn.metrics
    import pandas
    import wot

    gene_expression_file = "paper/2i_dmap_20.txt"
    growth_scores_file = "paper/growth_scores.txt"
    days_file = "paper/days.txt"

    # gene_expression has cells on rows, features (e.g. diffusion components) on columns
    gene_expression = pandas.read_table(gene_expression_file, index_col=0)
    # growth scores file has 2 columns: cell ids and scores
    growth_scores = pandas.read_table(growth_scores_file, index_col=0)
    growth_score_field_name = growth_scores.columns[0]
    # days file has 2 columns: cell ids and days
    days = pandas.read_table(days_file, index_col=0)
    day_field_name = days.columns[0]

    gene_expression = gene_expression.join(growth_scores).join(days)
    group_by_day = gene_expression.groupby(day_field_name)
    timepoints = group_by_day.groups.keys()
    timepoints.sort()
    max_transport_fraction = 0.4
    min_transport_fraction = 0.05
    min_growth_fit = 0.9
    l0_max = 100
    lambda1 = 1
    lambda2 = 1
    epsilon = 0.1
    growth_ratio = 2.5
    scaling_iter = 250

    # compute transport maps
    transport_maps = list()
    for i in range(len(timepoints) - 1):
        m1 = group_by_day.get_group(i)
        m2 = group_by_day.get_group(i+1)
        delta_t = timepoints[i + 1] - timepoints[i]
        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
            m1.drop([day_field_name, growth_score_field_name], axis=1),
            Y=m2.drop([day_field_name, growth_score_field_name], axis=1),
            metric="sqeuclidean")
        cost_matrix = cost_matrix / np.median(cost_matrix)
        growth_rate = m1.growth_score.values
        result = wot.optimal_transport(cost_matrix, growth_rate,
                                       delta_days=delta_t,
                                       max_transport_fraction=max_transport_fraction,
                                       min_transport_fraction=min_transport_fraction,
                                       min_growth_fit=min_growth_fit,
                                       l0_max=l0_max, lambda1=lambda1,
                                       lambda2=lambda2, epsilon=epsilon,
                                       growth_ratio=growth_ratio,
                                       scaling_iter=scaling_iter)
        transport = pandas.DataFrame(result["transport"], index=m1.index,
                                     columns=m2.index)
        transport_maps.append(
                {"transport_map": transport,
                 "t_minus_1": timepoints[i], "t": timepoints[i + 1]})
        transport.to_csv("transport" + str(timepoints[i]) + "_" + str(
            timepoints[i + 1]) + ".txt", index_label="id", sep="\t")


    # compute trajectories:
    trajectory = wot.trajectory(["day-9-c1-2i_6", "day-9-c1-2i_11"], transport_maps, 9)
    trajectory["ancestors"].to_csv(prefix + "_ancestors.txt", index_label="id",
                           sep="\t")
    trajectory["descendants"].to_csv(prefix + "_descendants.txt", index_label="id",
                             sep="\t")
