============
WOT
============


Uses time-course data to infer how the probability distribution of cells in gene-expression space evolves over time,
by using the mathematical approach of Optimal Transport (OT)


============
Example Usage
============

.. code-block:: python

    import numpy as np
    import scipy.stats
    import sklearn.metrics
    import pandas
    import wot

    diffusion_matrix = pandas.read_csv("tests/data/diffusion_map.csv.gz",
                                               index_col=0)  # cells on rows, diffusion components on columns
    growth_scores = pandas.read_csv("tests/data/growth_scores.csv.gz",
                                    index_col=0)
    days = pandas.read_csv("tests/data/days.csv.gz", index_col=0)
    diffusion_matrix = diffusion_matrix.join(growth_scores).join(days)
    group_by_day = diffusion_matrix.groupby('day')
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

    for i in range(len(timepoints) - 1):
        m1 = group_by_day.get_group(i)
        m2 = group_by_day.get_group(i+1)
        delta_t = timepoints[i + 1] - timepoints[i]
        cost_matrix = sklearn.metrics.pairwise.pairwise_distances(
            m1.drop(["day", "growth_score"], axis=1),
            Y=m2.drop(["day", "growth_score"], axis=1),
            metric='sqeuclidean')
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
        transport.to_csv("transport" + str(timepoints[i]) + "_" + str(
            timepoints[i + 1]) + ".csv", index_label="id")
