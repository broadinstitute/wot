import scipy
import numpy as np
import wot


class TrajectoryTrends:

    @staticmethod
    def __weighted_average(t, weights, cell_set_name, datasets=None, dataset_names=None, value_transform=None):
        if weights is not None:
            weights = weights / np.sum(weights)
        results = []
        for ds_index in range(len(datasets)):
            ds = datasets[ds_index]
            values = ds.x

            for feature_index in range(ds.x.shape[1]):
                array = values[:, feature_index]
                array = array.toarray().flatten() if scipy.sparse.isspmatrix(array) else array
                if value_transform is not None:
                    array = value_transform(array)
                mean = np.average(array, weights=weights)
                variance = np.average((array - mean) ** 2, weights=weights)
                trace = {
                    "dataset": dataset_names[ds_index],
                    "cell_set_name": cell_set_name,
                    "feature": str(ds.col_meta.index.values[feature_index]),
                    "name": cell_set_name + '_' + str(ds.col_meta.index.values[feature_index]),
                    "x": t,
                    "y": mean,
                    "variance": variance
                }
                results.append(trace)
        return results

    @staticmethod
    def compute(trajectory_results, unaligned_datasets, dataset_names, value_transform=None):
        """

        Args:
            trajectory_results (list): Results from trajectory_for_cell_sets_at_time_t
            unaligned_datasets (list): List of datasets that can be unaligned with transport map matrices to compute the mean and variance
            dataset_names (list): List of dataset names
            value_transform (function) A function to transform the values in the dataset that takes a numpy array and returns a numpy array of the same shape.
        Return:
            Dict that maps dataset name to traces
        """

        traces = []
        for trajectory_result in trajectory_results:
            p = trajectory_result['p']
            cell_ids = trajectory_result['cell_ids']
            cell_set_name = trajectory_result['cell_set']
            time = trajectory_result['t']
            aligned_datasets = []
            for ds_index in range(len(unaligned_datasets)):
                unaligned_ds = unaligned_datasets[ds_index]
                ds_order = unaligned_ds.row_meta.index.get_indexer_for(cell_ids)
                ds_order = ds_order[ds_order != -1]
                aligned_datasets.append(
                    wot.Dataset(unaligned_ds.x[ds_order], unaligned_ds.row_meta.iloc[ds_order], unaligned_ds.col_meta))
                traces += TrajectoryTrends.__weighted_average(t=time, weights=p,
                                                              datasets=aligned_datasets,
                                                              cell_set_name=cell_set_name,
                                                              dataset_names=dataset_names,
                                                              value_transform=value_transform)
            trace_fields_to_concat = ['x', 'y', 'variance']
            dataset_name_to_all_traces = {}
            for trace in traces:
                dataset_traces = dataset_name_to_all_traces.get(trace['dataset'])
                if dataset_traces is None:
                    dataset_traces = []
                    dataset_name_to_all_traces[trace['dataset']] = dataset_traces
                dataset_traces.append(trace)
            dataset_name_to_line_traces = {}

            for dataset_name in dataset_name_to_all_traces:
                all_traces = dataset_name_to_all_traces[dataset_name]
                trace_name_to_line_trace = {}
                dataset_name_to_line_traces[dataset_name] = trace_name_to_line_trace
                for trace in all_traces:
                    line_trace = trace_name_to_line_trace.get(trace['name'])
                    if line_trace is None:
                        line_trace = {'name': trace['name'],
                                      "mode": "lines+markers",
                                      "showlegend": True,
                                      "type": 'scatter',
                                      "dataset": trace['dataset'],
                                      "cell_set_name": trace['cell_set_name'],
                                      'feature': trace['feature']
                                      }
                        for field in trace_fields_to_concat:
                            line_trace[field] = [trace[field]]
                        trace_name_to_line_trace[line_trace['name']] = line_trace
                    else:
                        for field in trace_fields_to_concat:
                            line_trace[field] += [trace[field]]
            dataset_name_to_traces = {}
            for dataset_name in dataset_name_to_line_traces:
                trace_name_to_line_trace = dataset_name_to_line_traces[dataset_name]
                traces = []
                dataset_name_to_traces[dataset_name] = traces
                for trace_name in trace_name_to_line_trace:
                    trace = trace_name_to_line_trace[trace_name]
                    traces.append(trace)
                    # sort by time
                    for field in trace_fields_to_concat:
                        trace[field] = np.array(trace[field])
                    sort_order = np.argsort(trace['x'])
                    for field in trace_fields_to_concat:
                        trace[field] = trace[field][sort_order]

                    trace['ncells'] = len(p)
                    # if smooth:
                    #     x = trace['x']
                    #     xsmooth, ysmooth = wot.ot.TrajectoryUtil.kernel_smooth(x, trace['y'], stop=x[len(x) - 1])
                    #     trace['x'] = xsmooth
                    #     trace['y'] = ysmooth

        return dataset_name_to_traces
