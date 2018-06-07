var $trajectoryEl = $('#trajectory_vis');
var cellInfo; // id, x, y, t
var cellInfoHeaderNames;
var cellIdToIndex = {};
var featureIds;
var $cellSet = $('#cell_sets');
var $features = $('#features');
var $force_layout_group_by = $('#force_layout_group_by');
var forceLayoutColorScale = ['rgb(217,217,217)', 'rgb(255,0,0)'];

var interpolate = function (x, xi, yi, sigma) {
    var n = xi.length;
    var diff = new Float64Array(n);

    for (var i = 0; i < n; i++) {
        var val = x - xi[i];
        val *= -val;
        diff[i] = val;
    }
    var sigma2 = 2 * Math.pow(sigma, 2);
    var wi = new Float64Array(n);
    var wiSum = 0;
    for (var i = 0; i < n; i++) {
        var val = Math.exp(diff[i] / sigma2);
        wi[i] = val;
        wiSum += val;
    }
    var fx = 0;
    for (var i = 0; i < n; i++) {
        fx += yi[i] * wi[i];
    }
    return fx / wiSum;
};

var kernelSmooth = function (xi, yi, stop, start, steps, sigma) {
    var xlist = new Float64Array(steps);
    var stepSize = (stop - start) / (steps - 1);
    for (var i = 0; i < steps; i++) {
        xlist[i] = start;
        start += stepSize;
    }
    var fhat = new Float64Array(xlist.length);
    for (var i = 0, length = xlist.length; i < length; i++) {
        fhat[i] = interpolate(xlist[i], xi, yi, sigma);
    }
    return [xlist, fhat]

};


var createPlotAnimation = function (backgroundTrace, traces, elem, layout) {
    var $controls = $('<div style="display: inline;"><button class="btn btn-default btn-sm" name="play">Play</button>  <select style="width:auto;" class="form-control input-sm" data-name="time"></select></div>');
    var index = 0;
    var times = [];
    times.push('All');
    for (var i = 0; i < traces.length; i++) {
        times.push(traces[i].t);
    }
    var $time = $controls.find('[data-name=time]');
    $time.html(times.map(function (value, timeIndex) {
        return '<option value="' + timeIndex + '">' + value + '</option>';
    }).join(''));

    function showFrame() {
        if (index > traces.length) {
            index = 0; // reset
        }
        var concatTraces;
        var t;
        if (index === 0) {
            traces[0].marker.showscale = true;
            for (var i = 1; i < traces.length; i++) {
                traces[i].marker.showscale = false;
            }
            concatTraces = traces; // all
        } else {
            traces[index - 1].marker.showscale = true;
            concatTraces = [traces[index - 1]]
        }
        $time.val(index);

        Plotly.newPlot(elem, {
            data: [backgroundTrace].concat(concatTraces),
            layout: layout
        });
    }

    function nextTick() {
        if ($playBtn.text() === 'Pause') {
            index++;
            showFrame();
            setTimeout(nextTick, 500);
        }
    }

    var $playBtn = $controls.find('[name=play]');
    $playBtn.on('click', function () {
        if ($playBtn.text() === 'Play') {
            $playBtn.text('Pause');
            nextTick();
        } else {
            $playBtn.text('Play');
        }

    });


    $time.val('All');
    $time.on('change', function (e) {
        index = parseInt($time.val());
        $playBtn.text('Play');
        showFrame();
    });
    showFrame();
    return $controls;
};
var createForceLayoutPlotObject = function (showLegend) {
    var layout =
        {
            xaxis: {
                autorange: true,
                showgrid: false,
                zeroline: false,
                showline: false,
                autotick: false,
                ticks: '',
                showticklabels: false
            },
            yaxis: {
                autorange: true,
                showgrid: false,
                zeroline: false,
                showline: false,
                autotick: false,
                ticks: '',
                showticklabels: false
            },
            title: '',
            width: showLegend ? 1100 : 840, // leave room for legend
            height: 840,
            margin: {
                l: 0,
                b: 0,
                r: showLegend ? 300 : 0,
                t: 15,
                pad: 0
            },
            autosize: true
        };

    var backgroundTrace = {
        hoverinfo: 'none',
        showlegend: false,
        marker: {size: 2, color: 'rgb(217,217,217)', opacity: 0.5},
        mode: 'markers',
        type: 'scattergl',
        name: 'All Cells',
        x: cellInfo.x,
        y: cellInfo.y
    };
    return {layout: layout, backgroundTrace: backgroundTrace}

};
var createForceLayoutTrajectory = function (force, key) {

    var traces = force[key];
    var $div = $('<li style="list-style: none;"><h4>' + key +
        ' Trajectory</h4><div class="plot"></div><div data-name="controls"></div></li>'
    );
    traces.forEach(function (trace) {
        trace.mode = 'markers';
        trace.type = 'scattergl';
        trace.hoverinfo = 'none';
        trace.showlegend = false;
        trace.marker = {
            cmax: trace.marker.cmax,
            cmin: trace.marker.cmin,
            showscale: true,
            colorscale: forceLayoutColorScale,
            size: 2,
            color: trace.marker.color
        };
    });


    var elem = $div.find('.plot')[0];
    var forceLayoutInfo = createForceLayoutPlotObject(false);
    var backgroundTrace = forceLayoutInfo.backgroundTrace;
    var $controls = createPlotAnimation(backgroundTrace, traces, elem, forceLayoutInfo.layout);
    $controls.appendTo($div.find('[data-name=controls]'));

};
var cellForceLayoutInfo = null;
var featureForceLayoutInfo = null;
$.ajax('/cell_info/').done(function (result) {
    cellInfoHeaderNames = [];
    for (var key in result) {
        if (key !== 'id' && key !== 'x' && key !== 'y') {
            cellInfoHeaderNames.push(key);
        }
    }
    $force_layout_group_by.html(cellInfoHeaderNames.map(function (value) {
        return '<option value="' + value + '">' + value + '</option>'
    }).join(''));
    $force_layout_group_by.selectpicker('refresh');
    $force_layout_group_by.selectpicker('render');
    cellInfo = result;
    for (var i = 0, length = cellInfo.id.length; i < length; i++) {
        cellIdToIndex[cellInfo.id[i]] = i;
    }
    cellForceLayoutInfo = createForceLayoutPlotObject(true);

    featureForceLayoutInfo = createForceLayoutPlotObject(false);
    Plotly.newPlot('trajectory_set_vis', {
        data: [cellForceLayoutInfo.backgroundTrace],
        layout: cellForceLayoutInfo.layout
    });


});
$.ajax('/list_features/').done(function (results) {
    featureIds = results;

});
$('body').find('.selectpicker').selectpicker({
    iconBase: 'fa',
    tickIcon: 'fa-check',
    style: 'btn-default btn-sm'
});
var cellSetHtml;
$.ajax('/list_cell_sets/').done(function (result) {
    var options = [];
    for (var key in result) {
        options.push('<optgroup label="' + key + '">');
        var names = result[key];
        for (var i = 0; i < names.length; i++) {
            options.push('<option value="' + names[i] + '">');
            options.push(names[i]);
            options.push('</option>');
        }
        options.push('</optgroup>');
    }
    cellSetHtml = options.join('');
    $cellSet.html(cellSetHtml);
    $cellSet.selectpicker('refresh');
    $cellSet.selectpicker('render');
});

$cellSet.on('change', function () {
    var selectedSets = $cellSet.val();
    if (selectedSets == null || selectedSets.length === 0) {
        Plotly.newPlot('trajectory_set_vis', {
            data: [cellForceLayoutInfo.backgroundTrace],
            layout: cellForceLayoutInfo.layout
        });
        return;
    }
    $.ajax({url: '/cell_set_members/', data: {cell_set: selectedSets}}).done(function (results) {
        var traces = [];
        for (var setIndex = 0, nsets = results.length; setIndex < nsets; setIndex++) {
            var name = results[setIndex].name;
            var cellIds = results[setIndex].ids;
            var forceLayoutX = [];
            var forceLayoutY = [];
            for (var i = 0, ncells = cellIds.length; i < ncells; i++) {
                var index = cellIdToIndex[cellIds[i]];
                if (index != null) {
                    forceLayoutX.push(cellInfo.x[index]);
                    forceLayoutY.push(cellInfo.y[index]);
                }
            }
            traces.push({
                mode: 'markers',
                name: name,
                type: 'scattergl',
                hoverinfo: 'none',
                showlegend: true,
                marker: {
                    size: 2
                },
                x: forceLayoutX,
                y: forceLayoutY
            });


        }

        Plotly.newPlot('trajectory_set_vis', {
            data: [cellForceLayoutInfo.backgroundTrace].concat(traces),
            layout: cellForceLayoutInfo.layout
        });
    });

});

function split(val) {
    return val.split(/,\s*/);
}

function extractLast(term) {
    return split(term).pop();
}

function autocompleteFilter(term) {
    term = term.toUpperCase();
    var filteredResults = [];
    for (var i = 0, length = featureIds.length; i < length; i++) {
        if (featureIds[i].toUpperCase().startsWith(term)) {
            filteredResults.push(featureIds[i]);
            if (filteredResults.length === 10) {
                return filteredResults;
            }
        }
    }

    return filteredResults;

}

$features
    .on('keydown', function (event) {
        if (event.keyCode === $.ui.keyCode.TAB &&
            $(this).autocomplete('instance').menu.active) {
            event.preventDefault();
        }
    }).autocomplete({
    minLength: 1,
    source: function (request, response) {
        // delegate back to autocomplete, but extract the last term
        response(autocompleteFilter(extractLast(request.term)));
    },
    focus: function () {
        // prevent value inserted on focus
        return false;
    },
    select: function (event, ui) {
        var terms = split(this.value);
        // remove the current input
        terms.pop();
        // add the selected item
        terms.push(ui.item.value);
        // add placeholder to get the comma-and-space at the end
        terms.push('');
        this.value = terms.join(', ');
        return false;
    }
});

var showTrajectoryPlots = function (result) {
    var ancestryDivergenceTraces = result.ancestry_divergence_traces;
    var force = result.force;
    var datasetNameToTraces = result.dataset_name_to_traces;
    if (ancestryDivergenceTraces && ancestryDivergenceTraces.length > 0) {
        var $div = $('<li style="list-style: none;"><h4>Ancestry Divergence</h4><div class="plot"></div></li>');
        $div.appendTo($trajectoryEl);

        Plotly.newPlot($div.find('.plot')[0], ancestryDivergenceTraces, {
                title: '',
                showlegend: true,
                autosize: true,
                margin: {t: 15},
                yaxis: {
                    range: [-0.05, 1.05], autorange: false, 'zeroline': false,
                    title: 'Divergence'
                },
                xaxis: {
                    title: 'Time'
                }
            }
        );
    }
    if (datasetNameToTraces) {
        for (var key in datasetNameToTraces) {
            var $div = $('<li style="list-style: none;"><h4>Trajectory Trends</h4><div class="plot"></div></li>');
            $div.appendTo($trajectoryEl);
            var traces = datasetNameToTraces[key];
            traces.forEach(function (trace) {
                var smoothed = kernelSmooth(trace.x, trace.y, trace.x[trace.x.length - 1], 0, 1000, 0.7);
                trace.x = smoothed[0];
                trace.y = smoothed[1];
            });

            Plotly.newPlot($div.find('.plot')[0], traces,
                {
                    title: '',
                    autosize: true,
                    xaxis: {title: 'Time'},
                    yaxis: {title: 'Value', autorange: true, 'zeroline': false},
                    showlegend: true,
                    margin: {t: 15}
                });
        }

    }
    if (cellInfo != null) {
        for (var key in force) {
            createForceLayoutTrajectory(force, key);
        }
    }
    $trajectoryEl.sortable({handle: 'h4'});
};

var fetchTrajectoryData = function () {

    var selectedCellSets = $cellSet.val();
    var selectedFeatures = $features.val().split(',');
    var _selectedFeatures = [];
    selectedFeatures.forEach(function (value) {
        value = value.trim();
        if (value !== '') {
            _selectedFeatures.push(value);
        }
    });
    if (selectedCellSets.length > 0) {
        $trajectoryEl.empty();
        $('#trajectory_loading').show();
        var predefinedCellSets = [];
        var ncustom_cell_sets = 0;
        var data = {feature: _selectedFeatures};
        selectedCellSets.forEach(function (name) {
            var ids = customCellSetNameToIds[name];
            if (ids != null) {
                data['cell_set_name' + ncustom_cell_sets] = name;
                data['cell_set_ids' + ncustom_cell_sets] = ids;
                ncustom_cell_sets++;
            } else {
                predefinedCellSets.push(name);
            }
        });
        data.ncustom_cell_sets = ncustom_cell_sets;
        data.cell_set = selectedCellSets;
        $.ajax({url: '/trajectory/', data: data, method: 'POST'}).done(function (results) {
            showTrajectoryPlots(results);
            $('#trajectory_loading').hide();
        });
    }
};


$('#trajectory_form').on('submit', function (e) {
    e.preventDefault();
    fetchTrajectoryData();
});
$('#trajectory_submit').on('click', function (e) {
    e.preventDefault();
    fetchTrajectoryData();
});


var setSelectedFeature = null;
var $setFeature = $('#set_features');
var featureResult;
var userFilterValue = null;
var excludeZeros = false;
var timeToCellIds = {}; // time -> {ids:[], ntotal:Number}
var customCellSetNameToIds = {};
var $filterOp = $('#filter_op');
var filterOp = $filterOp.val();
var enterQuantile = true; // enter a quantile or value
var createSets = function () {
    var selectedTimes = [];
    $('.time').filter(':checked').each(function () {
        selectedTimes.push($(this).prop('name'))
    });

    selectedTimes.forEach(function (t) {
        var name = setSelectedFeature + '_' + filterOp + '_' + userFilterValue + '_' + t;
        customCellSetNameToIds[name] = timeToCellIds[t].ids;
    });
    var options = [];
    options.push('<optgroup label="Custom Sets">');
    for (var name in customCellSetNameToIds) {
        options.push('<option value="' + name + '">');
        options.push(name);
        options.push('</option>');
    }
    options.push('</optgroup>');
    var html = $cellSet.html();
    var val = $cellSet.val();
    $cellSet.html(options.join('') + cellSetHtml);
    $cellSet.val(val);
    $cellSet.selectpicker('refresh');
    $cellSet.selectpicker('render');
};

var featureForceLayoutData = null;
var showFeature = function () {
    if (featureResult == null) {
        return;
    }

    var sortedValues = featureResult.values.slice(0).sort(function (a, b) {
        return (a === b ? 0 : (a < b ? -1 : 1));
    });
    if (excludeZeros) {
        sortedValues = sortedValues.filter(function (value) {
            return value !== 0;
        });
    }
    var filterValue = null;
    var f = function () {
        return true;
    };
    if (userFilterValue != null && !isNaN(userFilterValue)) {

        filterValue = enterQuantile ? d3.quantile(sortedValues, userFilterValue / 100.0) : userFilterValue;
        if ($('#filter_op').val() == 'gt') {
            f = function (d) {
                return d > filterValue;
            };
        } else {
            f = function (d) {
                return d < filterValue;
            };
        }
    }
    timeToCellIds = {};

    var forceLayoutX = [];
    var forceLayoutY = [];

    var values = [];
    var times = [];

    for (var idIndex = 0, length = featureResult.ids.length; idIndex < length; idIndex++) {
        var index = cellIdToIndex[featureResult.ids[idIndex]];
        if (index != null && index !== -1) {
            var t = cellInfo.t[index];
            var obj = timeToCellIds[t];
            if (obj == null) {
                times.push(t);
                obj = {ids: [], ntotal: 0};
                timeToCellIds[t] = obj;
            }
            if (f(featureResult.values[idIndex])) {
                obj.ids.push(cellInfo.id[index]);
            }
            obj.ntotal++;
            if (featureForceLayoutData == null) {
                forceLayoutX.push(cellInfo.x[index]);
                forceLayoutY.push(cellInfo.y[index]);
            }
        } else {
            console.log(featureResult.ids[idIndex] + ' missing')
        }
    }
    times.sort(function (a, b) {
        return (a === b ? 0 : (a < b ? -1 : 1));
    });
    var html = [];
    if (filterValue != null) {
        var percentFormatter = d3.format('.1f');
        var groupedThousands = d3.format(',');
        html.push('<h4>Cell Selection Summary <small> ' + ($('#filter_op').val() === 'gt' ? 'Greater then' : 'Less then') + ' ' + d3.format('.2f')(filterValue) + '</small></h4><table class="table table-condensed table-bordered"><tr><th><input name="select_all" type="checkbox" checked></th><th>Time</th><th># Cells Selected</th><th>% Cells Selected</th></tr>');
        var totalPass = 0;
        var total = 0;
        for (var i = 0; i < times.length; i++) {
            var t = times[i];
            html.push('<tr>');
            var obj = timeToCellIds[t];
            html.push('<td><input class="time" name="' + t + '" type="checkbox" ' + (obj.ids.length === 0 ? 'disabled' : 'checked') + '>');
            html.push('<td>');
            html.push(t);
            html.push('</td>');
            html.push('<td>');
            html.push(groupedThousands(obj.ids.length) + '/' + groupedThousands(obj.ntotal));
            html.push('</td>');
            html.push('<td>');
            html.push(percentFormatter(100 * (obj.ids.length / obj.ntotal)));
            html.push('</td>');
            html.push('</tr>');
            totalPass += obj.ids.length;
            total += obj.ntotal;
        }
        html.push('<td></td>');
        html.push('<td>All</td>');
        html.push('<td>');
        html.push(groupedThousands(totalPass) + '/' + groupedThousands(total));
        html.push('</td>');
        html.push('<td>');
        html.push(percentFormatter(100 * (totalPass / total)));
        html.push('</td>');
        html.push('</tr>');
        html.push('</table>');


    }
    $('#table_vis').html(html.join(''));

    Plotly.react('histogram_vis', [
        {
            name: 'All Cells',
            x: featureResult.values,
            type: 'histogram',
            marker: {color: 'LightGrey', opacity: 0.7}
        }], {
        shapes: filterValue == null ? null : [{
            type: 'line',
            xref: 'x',
            yref: 'paper',
            x0: filterValue,
            x1: filterValue,
            y0: 0,
            y1: 1,
            line: {
                color: 'rgb(0, 0, 0)',
                width: 2
            }
        }]
    });
    if (featureForceLayoutData == null) {
        featureForceLayoutData = {
            data: [featureForceLayoutInfo.backgroundTrace, {
                mode: 'markers',
                type: 'scattergl',
                hoverinfo: 'none',
                showlegend: false,
                marker: {
                    cmin: sortedValues[0],
                    cmax: sortedValues[sortedValues.length - 1],
                    color: featureResult.values,
                    showscale: true,
                    colorscale: forceLayoutColorScale,
                    size: 2
                },
                x: forceLayoutX,
                y: forceLayoutY,
                ids: featureResult.ids
            }],
            layout: featureForceLayoutInfo.layout
        };
        var $controls = $('#force_layout_vis_controls');
        var groupBy = ['t']; //$force_layout_group_by.val();
        if (groupBy != null && groupBy.length > 0) {
            featureForceLayoutData.groupedTraces = groupTraces(featureForceLayoutData, groupBy);
            $controls.show().html(createPlotAnimation(featureForceLayoutData.data[0], featureForceLayoutData.groupedTraces, 'force_layout_vis', featureForceLayoutInfo.layout));
        } else {
            $controls.hide();
            Plotly.newPlot('force_layout_vis', featureForceLayoutData);
        }

    }
};

var groupTraces = function (plotData, groupingFields) {
    var traceToTile = plotData.data[1];
    var ids = traceToTile.ids;
    if (ids.length !== traceToTile.x.length) {
        throw new Error('Length of ids not equal to coordinates');
    }
    var nfields = groupingFields.length;
    var traceNameToTrace = {};
    for (var i = 0, n = ids.length; i < n; i++) {
        var index = cellIdToIndex[ids[i]];
        var keyArray = [];
        for (var fieldIndex = 0; fieldIndex < nfields; fieldIndex++) {
            var field = groupingFields[fieldIndex];
            var value = cellInfo[field][index];
            keyArray.push(value);
        }
        var key = keyArray.join(', ');
        var trace = traceNameToTrace[key];
        if (trace == null) {
            trace = {
                x: [],
                y: [],
                keyArray: keyArray,
                t: key,
                mode: 'markers',
                type: 'scattergl',
                hoverinfo: 'none',
                showlegend: false,
                marker: {
                    cmin: traceToTile.marker.cmin,
                    cmax: traceToTile.marker.cmax,
                    color: [],
                    showscale: true,
                    colorscale: forceLayoutColorScale,
                    size: 2
                }
            };
            traceNameToTrace[key] = trace;
        }
        trace.x.push(traceToTile.x[i]);
        trace.y.push(traceToTile.y[i]);
        trace.marker.color.push(traceToTile.marker.color[i]);
    }
    var traces = [];
    for (var key in traceNameToTrace) {
        traces.push(traceNameToTrace[key]);
    }
    traces.sort(function (t1, t2) {
        for (var i = 0; i < t1.keyArray.length; i++) {
            var a = t1.keyArray[i];
            var b = t2.keyArray[i];
            var val = (a === b ? 0 : (a < b ? -1 : 1));
            if (val !== 0) {
                return val;
            }
        }
        return 0;

    });
    return traces;
};

var fetchFeatureData = function () {
    var text = $setFeature.val().trim();
    if (text !== '') {
        if (text !== setSelectedFeature) {
            setSelectedFeature = text;
            $('#set_loading').show();
            $.ajax({url: '/feature_value/', data: {feature: text}}).done(function (result) {
                featureForceLayoutData = null;
                featureResult = result;
                showFeature();
                $('#set_loading').hide();
            });
        }
    }
};


$setFeature.on('keydown', function (event) {
    if (event.keyCode === $.ui.keyCode.TAB &&
        $(this).autocomplete('instance').menu.active) {
        event.preventDefault();
    }
}).autocomplete({
    minLength: 1,
    focus: function () {
        // prevent value inserted on focus
        return false;
    },
    source: function (request, response) {
        // delegate back to autocomplete, but extract the last term
        response(autocompleteFilter(request.term));
    }
});
var $filterValue = $('#filter_quantile');

$('#enter-quantile').on('click', function (e) {
    e.preventDefault();
    enterQuantile = true;
    $filterValue.prop('placeholder', 'Enter a quantile');
    $('#enter-value').show();
    $(this).hide();
    showFeature();
});
$('#enter-value').on('click', function (e) {
    e.preventDefault();
    enterQuantile = false;
    $filterValue.prop('placeholder', 'Enter a number');
    $('#enter-quantile').show();
    $(this).hide();
    showFeature();
});

$('#set_form').on('submit', function (e) {
    e.preventDefault();
    fetchFeatureData();
});
$('#set_submit').on('click', function (e) {
    e.preventDefault();
    fetchFeatureData();
});

$('#create_set').on('click', function (e) {
    createSets();
});

function debounce(func, wait, immediate) {
    var timeout;
    return function () {
        var context = this, args = arguments;
        var later = function () {
            timeout = null;
            if (!immediate) func.apply(context, args);
        };
        var callNow = immediate && !timeout;
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
        if (callNow) func.apply(context, args);
    };
};
$filterValue.on('keyup', debounce(function (e) {
    if (e.which === 13) {
        e.preventDefault();
    }
    userFilterValue = parseFloat($(this).val().trim());
    showFeature();
}, 500));

$filterOp.on('change', function (e) {
    filterOp = $(this).val();
    showFeature();
});

$('#table_vis').on('click', 'input[name=select_all]', function (e) {
    var selected = $(this).prop('checked');
    $('#table_vis').find('.time').prop('checked', selected);
});

$('#export_set').on('click', function () {
    var text = [];
    for (var name in customCellSetNameToIds) {
        var ids = customCellSetNameToIds[name];
        text.push(name);
        text.push('\t');
        text.push('');
        text.push(ids.join('\t'));
        text.push('\n');
    }
    if (text.length > 0) {
        var blob = new Blob([text.join('')], {type: "text/plain;charset=utf-8"});
        saveAs(blob, 'custom_cell_sets.gmt');
    }
});
