# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import scipy


def transport_map_distance(transport_map_1, transport_map_2, column_weights):
    """
        Compute the distance between two transport maps

        Args:
          transport_map_1 (DataFrame): A transport map
          transport_map_2 (DataFrame): A transport map
          column_weights (ndarray): An array of cluster weights by column
        Returns:
          float: The distance between the two transport maps.
    """
    column_distances = np.zeros(transport_map_1.shape[1])
    for j in range(transport_map_1.shape[1]):
        # Kullback-Leibler divergence
        val = scipy.stats.entropy(
            transport_map_1.iloc[:, j],
            qk=transport_map_2.iloc[:, j])
        if not np.isinf(val):
            column_distances[j] = val

    return np.average(column_distances,
                      weights=column_weights)


def transport_maps_by_time(cluster_transport_maps,
                           cluster_weights_by_time):
    """
          Convert many cluster by cluster transport maps across timepoints,
          to a single cluster by cluster transport map. Weight for each
          cluster at timepoint t is determined by fraction of total cells in
          cluster c at timepoint t

          Args:
              cluster_transport_maps (list): A list of cluster by
              cluster transport maps
              cluster_weights_by_time (list): A list of cluster weights by
              time. For
              example cluster_weights_by_time[1][3] gives the weight for the 2nd
              timepoint and the 4th cluster. Weights must sum to 1 for each
              cluster across time.
          Returns:
              DataFrame: A cluster by cluster transport map
    """
    cluster_ids = cluster_transport_maps[0].index
    combined_cluster_map = pd.DataFrame(index=cluster_ids,
                                        columns=cluster_ids,
                                        data=0)

    for cluster_index in range(combined_cluster_map.shape[0]):
        per_cluster_total_weight = 0
        for transport_index in range(len(cluster_transport_maps)):
            transport_map = cluster_transport_maps[transport_index]
            weight = cluster_weights_by_time[transport_index][cluster_index]
            per_cluster_total_weight += weight
            c = transport_map.iloc[:, cluster_index]
            combined_cluster_map.iloc[:, cluster_index] += c * weight

        if per_cluster_total_weight == 0:
            combined_cluster_map.iloc[:, cluster_index] = 0
        else:
            combined_cluster_map.iloc[:, cluster_index] /= per_cluster_total_weight

    return combined_cluster_map


def transport_map_by_cluster(transport_map, grouped_by_cluster, cluster_ids):
    """
      Convert a cell by cell transport map, to a cluster by cluster transport
      map

      Args:
          transport_map (DataFrame): A transport map
          grouped_by_cluster (GroupBy): A GroupBy object
              returned from pd.DataFrame.groupby(). The object maps
              cluster ids to cell ids.
          cluster_ids (list): A list of unique cluster ids
      Returns:
          DataFrame: A cluster by cluster transport map
    """
    result = pd.DataFrame(index=cluster_ids, columns=cluster_ids,
                          data=0)

    for cluster_index_1 in range(len(cluster_ids)):
        cluster_group_1 = grouped_by_cluster.get_group(
            cluster_ids[cluster_index_1])
        # subset rows

        row_subset = transport_map.loc[
            transport_map.index.intersection(
                cluster_group_1.index)]
        if row_subset.shape[0] == 0:
            continue
        for cluster_index_2 in range(len(cluster_ids)):
            cluster_group_2 = cluster_group_1 if \
                cluster_index_1 == cluster_index_2 else \
                grouped_by_cluster.get_group(
                    cluster_ids[cluster_index_2])

            # subset columns
            row_and_column_subset = row_subset[
                row_subset.columns.intersection(
                    cluster_group_2.index)]
            result.iloc[cluster_index_1, cluster_index_2] = \
                row_and_column_subset.values.sum()

    return result


def get_column_weights(all_cell_ids, grouped_by_cluster, cluster_ids):
    total_cluster_size = np.zeros(len(cluster_ids))
    for cluster_index in range(len(cluster_ids)):
        cluster_group = grouped_by_cluster.get_group(cluster_ids[cluster_index])
        total_cluster_size[cluster_index] = cluster_group.loc[cluster_group.index.intersection(all_cell_ids)].shape[0]
    return total_cluster_size


def get_weights(all_cell_ids, column_cell_ids_by_time, grouped_by_cluster,
                cluster_ids):
    """
          Compute the distance between transport maps.

          Args:
              all_cell_ids (set) : Set of all cell ids contained in
              uncollapsed transport maps.
              column_cell_ids_by_time (list): A list of uncollapsed
              cell ids in transport maps.
              grouped_by_cluster (GroupBy): A GroupBy object
              returned from pd.DataFrame.groupby(). The object maps
              cell ids to cluster ids.
              cluster_ids (list): A list of unique cluster ids
          Returns:
              list: A list for each timepoint, containing an array of
              weights for each cluster
    """
    cluster_weights_by_time = []  # each time, then each cluster
    total_cluster_size = get_column_weights(all_cell_ids, grouped_by_cluster, cluster_ids)

    for time_index in range(len(column_cell_ids_by_time)):
        weights = np.zeros(len(cluster_ids))
        cluster_weights_by_time.append(weights)
        for cluster_index in range(len(cluster_ids)):
            cluster_group = grouped_by_cluster.get_group(
                cluster_ids[cluster_index])
            timepoint_cluster_size = cluster_group.loc[
                cluster_group.index.intersection(
                    column_cell_ids_by_time[time_index])].shape[0]

            # join ids in transport map to ids in cluster
            fraction = timepoint_cluster_size / total_cluster_size[cluster_index] if total_cluster_size[
                                                                                         cluster_index] > 0 else 0
            weights[cluster_index] = fraction

    return {'cluster_size': total_cluster_size,
            'cluster_weights_by_time': cluster_weights_by_time}


def cluster_distance(transport_maps_1, transport_maps_2, grouped_by_cluster,
                     cluster_ids):
    """
       Compute the distance between transport maps.

       Args:
           transport_maps_1 (list): A list of transport maps.
           transport_maps_2 (list): A list of transport maps
           grouped_by_cluster (GroupBy): A GroupBy object
             returned from pd.DataFrame.groupby(). The object maps
             cell ids to cluster ids .
           cluster_ids (list): A list of unique cluster ids
       Returns:
           float: Distance between the transport maps
    """
    if len(transport_maps_1) != len(transport_maps_2):
        raise ValueError('Length of transport maps are not equal')

    # column weights are used to compare two transport maps, weight equals size
    #  of clusters
    weights = get_weights(transport_maps_1, grouped_by_cluster, cluster_ids)
    cluster_weights_by_time = weights['cluster_weights_by_time']
    cluster_size = weights['cluster_size']
    cluster_transport_maps_1 = []
    cluster_transport_maps_2 = []
    for i in range(len(transport_maps_1)):
        # compute total mass in each cluster to create a cluster by cluster
        # matrix
        cluster_transport_maps_1.append(transport_map_by_cluster(
            transport_maps_1[i], grouped_by_cluster, cluster_ids))
        cluster_transport_maps_2.append(transport_map_by_cluster(
            transport_maps_2[i], grouped_by_cluster, cluster_ids))
    # combine each set of cluster transport maps into one transport map
    # weight for each cluster at timepoint t is determined by fraction of
    # cells in cluster c at timepoint t

    combined_cluster_map_1 = transport_maps_by_time(cluster_transport_maps_1,
                                                    cluster_weights_by_time)
    combined_cluster_map_2 = transport_maps_by_time(cluster_transport_maps_2,
                                                    cluster_weights_by_time)
    return transport_map_distance(combined_cluster_map_1,
                                  combined_cluster_map_2, cluster_size)
