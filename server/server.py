# import matplotlib
#
# matplotlib.use('Agg')
import flask
import numpy as np
import pandas as pd
import os
import wot.ot
import re
import h5py
from flask_compress import Compress
import json

# run locally:
# export FLASK_APP=server.py
# export FLASK_DEBUG=1
# flask run

# run with gunicorn:
# gunicorn -w 4 -b 127.0.0.1:4000 server:app


app = flask.Flask(__name__)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
compress = Compress()
compress.init_app(app)

with open('config.json') as json_data:
    server_config = json.load(json_data)

name_to_transport_maps = {}
name_to_cell_scores = {}
name_to_cell_sets = {}
signatures_config = server_config['cell_scores']
for key in signatures_config:
    name_to_cell_scores[key] = signatures_config[key]

transport_maps_config = server_config['transport_maps']
for key in transport_maps_config:
    name_to_transport_maps[key] = wot.io.list_transport_maps(transport_maps_config[key])

gene_matrix_path = server_config['gene_matrix']
cell_sets_config = server_config['cell_sets']

for key in cell_sets_config:
    name_to_cell_sets[key] = cell_sets_config[key]

coords = pd.read_csv(server_config['cell_coords'], index_col=0)
Nx = 400
Ny = 400
xmin = np.min(coords['x'])
xmax = np.max(coords['x'])
ymin = np.min(coords['y'])
ymax = np.max(coords['y'])
coords['px'] = np.floor(np.interp(coords['x'].values, [xmin, xmax], [0, Nx])).astype(int)
coords['py'] = np.floor(np.interp(coords['y'].values, [ymin, ymax], [0, Ny])).astype(int)
coords = coords.drop(['x', 'y'], axis=1)

if __name__ == "__main__":
    app.run()


@app.route("/list_genes/", methods=['GET'])
def list_genes():
    f = h5py.File(gene_matrix_path, 'r')
    ids = f['col_attrs/id'][()]
    if ids.dtype.kind == 'S':
        ids = ids.astype(str)
    f.close()
    return flask.jsonify(list(np.char.upper(ids)))


@app.route("/list_signatures/", methods=['GET'])
def list_signatures():
    results = {}
    for key in name_to_cell_scores:
        f = h5py.File(name_to_cell_scores[key], 'r')
        ids = f['col_attrs/id'][()]
        # ids.sort()
        ids = ids.astype(str)
        f.close()

        results[key] = list(ids)
    return flask.jsonify(results)


@app.route("/signatures/", methods=['GET'])
def modules():
    getter = flask.request.args
    s = getter.get('s')  # key;signature_name
    tokens = s.split(';')
    path = name_to_cell_scores[tokens[0]]
    ds = wot.io.read_dataset(path, col_filter={'id': lambda x: x == tokens[1]})
    return flask.jsonify({'ids': ds.row_meta.index.values.tolist(), 'values': ds.x[:, 0].tolist()})


@app.route("/list_cell_sets/", methods=['GET'])
def list_cell_sets():
    result = {}
    for key in name_to_cell_sets:
        with open(name_to_cell_sets[key]) as fp:
            names = []
            for line in fp:
                tokens = line.split('\t')
                names.append(tokens[0])
            result[key] = names
    return flask.jsonify(result)


@app.route("/trajectory/", methods=['GET', 'POST'])
def trajectory():
    getter = flask.request.args
    if flask.request.method == 'POST':
        getter = flask.request.form
    prefined_cell_sets = getter.getlist('cell_set[]')
    ncustom = int(getter.get('ncustom', '0'))
    ncells = int(getter.get('ncells', '1000'))
    genes = getter.getlist('gene[]')

    transport_map_name = getter.get('transport_map')
    dt = None
    try:
        dt = float(getter.get('dt', ''))
    except ValueError:
        pass
    cell_scores = getter.getlist('score[]')  # key;id
    cell_sets = []
    for cell_set in prefined_cell_sets:
        tokens = cell_set.split(';')  # key;set_name;time
        gmt_path = name_to_cell_sets[tokens[0]]
        cell_set_gmt = os.path.abspath(gmt_path)
        cell_set_name = tokens[1]
        cell_set_ds = wot.io.read_gene_sets(cell_set_gmt)
        cell_set_ds = wot.ot.filter_sets(cell_set_ds, cell_set_name)
        cell_sets.append({'ds': cell_set_ds, 'name': cell_set_name,
                          'time': float(re.compile('.*Day ([0-9.]+)').match(cell_set_name).group(1))})  # FIXME

    if ncustom > 0:
        for i in range(ncustom):
            ids = getter.getlist('id' + str(i) + '[]')
            name = getter.get('name' + str(i), '')
            cell_set_ds = wot.Dataset(np.ones(shape=(len(ids), 1)), pd.DataFrame(index=ids), pd.DataFrame(index=[name]))
            m = re.compile('.* Day ([0-9.]+)').match(name)  # FIXME
            if m is not None:
                time = float(m.group(1))
                cell_sets.append({'ds': cell_set_ds, 'name': name, 'time': time})
            else:
                print(name)
    results = []
    cell_score_key_to_ids = {}
    for s in cell_scores:
        tokens = s.split(';')
        ids = cell_score_key_to_ids.get(tokens[0])
        if ids is None:
            ids = set()
            cell_score_key_to_ids[tokens[0]] = ids
        ids.add(tokens[1])
    # gene_set_filter = '|'.join(gene_sets)
    datasets = []
    summaries = []
    if len(genes) > 0:
        genes = set(np.char.lower(np.array(genes)))
        # if 'oskm' in genes:
        #     datasets.append(wot.io.read_dataset(oskm_matrix_path))
        #     summaries.append('mean')
        #     genes.remove('oskm')

        if len(genes) > 0:
            datasets.append(wot.io.read_dataset(gene_matrix_path, col_filter={'id': lambda x: x.lower() in genes},
                                                row_filter=None))
            summaries.append('mean')

    for key in cell_score_key_to_ids:
        path = name_to_cell_scores[key]
        ids = cell_score_key_to_ids[key]
        # gene_set_filter_regex = re.compile(gene_set_filter)
        datasets.append(
            wot.io.read_dataset(path, row_filter=None, col_filter={'id': lambda x: x in ids}))
        summaries.append(None)

    transport_maps = name_to_transport_maps[transport_map_name]
    for cell_set_dict in cell_sets:
        cell_set_ds = cell_set_dict['ds']
        result = wot.ot.TrajectorySampler.compute(cell_set_ds=cell_set_ds,
                                                  transport_maps=transport_maps,
                                                  time=cell_set_dict['time'], unaligned_datasets=datasets,
                                                  summaries=summaries,
                                                  verbose=True, sampling_loader=None, ncells=ncells,
                                                  end=time + dt if dt is not None else None,
                                                  start=time - dt if dt is not None else None)
        results.append(result)

    trace_data = []
    for result_dict in results:  # each cell set
        for trace in result_dict['traces']:
            trace_data.append(trace)

    cell_set_name_to_traces = {}

    for result_dict in results:  # each cell set

        for p in result_dict['pvecs']:
            cell_set = p['cell_set']
            traces = cell_set_name_to_traces.get(cell_set)
            if traces is None:
                traces = []
                cell_set_name_to_traces[cell_set] = traces

            t = p['t']
            v = p['v']
            cell_ids = p['cell_ids']
            joined = coords.join(pd.DataFrame(index=cell_ids, data={'v': v}), how='right')
            df_sum = joined.groupby(['px', 'py']).sum()
            traces.append({'t': t, 'x': df_sum.index.get_level_values(0).tolist(),
                           'y': df_sum.index.get_level_values(1).tolist(),
                           'marker': {'color': df_sum['v'].values.tolist()}})

    for name in cell_set_name_to_traces:
        traces = cell_set_name_to_traces[name]
        traces.sort(key=lambda x: x['t'])

    # group scatter plots by name
    gene_name_to_trace = {}
    violin_name_to_traces = {}
    scatter_trace_fields = ['x', 'y', 'size', 'text']
    for trace in trace_data:
        if trace.get('type') == 'scatter':
            scatter = gene_name_to_trace.get(trace['name'])
            if scatter is None:
                for field in scatter_trace_fields:
                    trace[field] = np.array([trace[field]])

                gene_name_to_trace[trace['name']] = trace
            else:
                for field in scatter_trace_fields:
                    scatter[field] = np.concatenate((scatter[field], [trace[field]]))

        else:
            violin_traces = violin_name_to_traces.get(trace['set_name'])
            if violin_traces is None:
                violin_traces = []
                violin_name_to_traces[trace['set_name']] = violin_traces
            trace['showlegend'] = False
            trace['bandwidth'] = ncells ** (-1. / (1 + 4))
            violin_traces.append(trace)

    cell_scores_line_traces = []
    for key in violin_name_to_traces:
        traces = violin_name_to_traces[key]
        traces.sort(key=lambda x: x['name'])
        x = []
        y = []
        for trace in traces:
            x.append(trace['name'])
            y.append(trace['median'])
        xsmooth, ysmooth = wot.ot.TrajectorySampler.kernel_smooth(x, y, stop=x[len(x) - 1])
        cell_scores_line_traces.append({'name': key, 'x': xsmooth.tolist(), 'y': ysmooth.tolist(), 'mode': 'lines'})

    for key in gene_name_to_trace:
        trace = gene_name_to_trace[key]
        sort_order = np.argsort(trace['x'])
        # max_size = max(max_size, np.max(trace['size']))
        for field in scatter_trace_fields:
            trace[field] = trace[field][sort_order]
        x = trace['x']
        xsmooth, ysmooth = wot.ot.TrajectorySampler.kernel_smooth(x, trace['y'], stop=x[len(x) - 1])
        trace['x'] = xsmooth.tolist()
        trace['y'] = ysmooth.tolist()
        trace['mode'] = 'lines'
        del trace['size']
        del trace['text']
        # trace['size'] = trace['size'].tolist()
        # trace['text'] = trace['text'].tolist()
        trace['showlegend'] = True
        # trace['sizemode'] = 'area'
        # trace['sizemin'] = 4
        # trace['marker'] = {'size': trace['size'], 'sizeref': (2 * 100) / (4 * 4), 'size_min': 4}


    return flask.jsonify(
        {'line_traces': cell_scores_line_traces, 'force': cell_set_name_to_traces, 'violins': violin_name_to_traces,
         'scatters': list(gene_name_to_trace.values())})
