var createModuleOptions = function (options) {
  ['XU', 'W'].forEach(function (name) {
    options.push('<optgroup label="' + name + '">');
    ['6h', '48h', '96h'].forEach(function (hour) {
      for (var i = 0; i < 50; i++) {
        options.push('<option value="' + hour + ';' + name + ';' + i + '">');
        options.push(hour + ' ' + name + ' ' + i);
        options.push('</option>');
      }
    });

    options.push('</optgroup>');
  });
};
$('#cell_set_creator').on('click', function () {
  $('#customCellSet').
    dialog(
      {
        width: 600, height: Math.min(1000, $(window).height() - 50), position: {
          my: 'center', at: 'top', of: window,
          collision: 'none'
        }
      });
});
var $cellSet = $('#cell_set');
var $geneSet = $('#gene_set');
var $customGeneSet = $('#custom_gene_set');

$.ajax('/list_signatures/').done(function (results) {
  var geneSetOptions = [];
  for (var key in results) {
    geneSetOptions.push('<optgroup label="' + key + '">');
    results[key].
      forEach(function (name) {
        geneSetOptions.push('<option value="' + key + ';' + name + '">');
        geneSetOptions.push(name);
        geneSetOptions.push('</option>');
      });
    geneSetOptions.push('</optgroup>');
  }

  $customGeneSet.html(geneSetOptions.join(''));
  $customGeneSet.selectpicker('refresh');
  $customGeneSet.selectpicker('render');
  $geneSet.html(geneSetOptions.join(''));
  $geneSet.selectpicker('refresh');
  $geneSet.selectpicker('render');
});

var $gene = $('#gene');
$.ajax('/list_cell_sets/').done(function (json) {
  var options = [];
  for (var key in json) {
    options.push('<optgroup label="' + key + '">');
    var results = json[key];
    for (var i = 0; i < results.length; i++) {
      options.push('<option value="' + key + ';' + results[i] + '">');
      options.push(key + ' ' + results[i]);
      options.push('</option>');
    }
    options.push('</optgroup>');
  }
  $cellSet.html(options.join(''));
  $cellSet.selectpicker('refresh');
  $cellSet.selectpicker('render');
});

function split (val) {
  return val.split(/,\s*/);
}

function extractLast (term) {
  return split(term).pop();
}

function autocompleteFilter (results, term) {
  term = term.toUpperCase();
  var filteredResults = [];
  for (var i = 0, length = results.length; i < length; i++) {
    if (results[i].startsWith(term)) {
      filteredResults.push(results[i]);
      if (filteredResults.length === 10) {
        return filteredResults;
      }
    }
  }
  return filteredResults;

}

var $customTime = $('#custom_time');
$.ajax('/list_times/').done(function (results) {
  var html = [];
  html.push('<option>All</option>');
  results.forEach(function (t) {
    html.push('<option>');
    html.push(t);
    html.push('</option>');
  });
  $customTime.html(html.join(''));
  $customTime.selectpicker('refresh');
  $customTime.selectpicker('render');
});
var $customType = $('[name=custom_type]');
var $quantileValue = $('#quantile-value');
var $customCellSetSelector = $('#custom_cell_set');
var $globalSset = $('#global_set');
var $create = $('#create');
customCellSets = {};
var splitTime = 8;
var q = '';
var dayToCountSerum;
var dayToCount2i;
var splitByTime = function () {
  // split by time
  var dayToIds = {};

  for (var i = 0; i < customCellSetIds.length; i++) {
    var id = customCellSetIds[i];
    var day = id.split('_')[0];

    var ids = dayToIds[day];
    if (ids == null) {
      ids = [];
      dayToIds[day] = ids;
    }
    ids.push(id);
  }
  return dayToIds;
};
var createCustomCellSet = function () {
  var geneSetName = $customGeneSet.val();
  var time = parseFloat($customTime.val());
  var isCustomSetGlobal = isNaN(time);
  var type = $customType.filter(':checked').val();

  if (!isCustomSetGlobal) {
    var name = geneSetName + ' Day ' + time + ' q ' + q;
    if (time > splitTime) {
      name += '_' + type;
    }
    customCellSets[name] = customCellSetIds;
  } else {
    var formatter = d3.format('.2f');
    var dayToIds = splitByTime();
    var dayToCount = type.toLowerCase() === 'serum' ? dayToCountSerum : dayToCount2i;
    for (var day in dayToIds) {
      var ids = dayToIds[day];
      name = geneSetName + ' Day ' + day.substring(1) + ' q ' + q + '_' + type + '_ncells_' + ids.length +
        '_percent_' + formatter(100 * (ids.length / dayToCount[day]));
      customCellSets[name] = ids;
    }
  }

  var selection = $customCellSetSelector.val();
  var options = [];
  for (var key in customCellSets) {
    options.push('<option>' + key + '</option>');
  }

  $customCellSetSelector.html(options.join(''));
  $customCellSetSelector.val(selection);
  $customCellSetSelector.selectpicker('refresh').selectpicker('render');
};

$create.on('click', createCustomCellSet);

var savedQuantileInfo = {};
var updateCustomCellSetsPlot = function () {

  var geneSetName = $customGeneSet.val();
  var time = parseFloat($customTime.val());
  var isCustomSetGlobal = isNaN(time);
  $customType.attr('disabled', time <= splitTime);
  var type = $customType.filter(':checked').val();
  if (!isCustomSetGlobal && time <= splitTime) {
    type = '';
  }

  if (geneSetName == null) {
    return;
  }
  var d = $.Deferred();
  if (savedQuantileInfo[geneSetName] != null) {
    d.resolve(savedQuantileInfo[geneSetName]);
  } else {
    $.ajax({url: '/signatures/?s=' + encodeURIComponent(geneSetName), context: geneSetName}).done(function (result) {
      savedQuantileInfo[this] = result;
      d.resolve(result);
    });
  }

  var title = geneSetName;
  if (title.indexOf(';') !== -1) {
    title = title.split(';');
    title = title[0] + ' ' + title[1] + ' ' + title[2];
  }
  $('#histogram').html('Loading...');
  $('#force').empty();
  $('#filter_summary').empty();
  q = $quantileValue.val();
  var quantile = parseFloat(q) / 100.0;
  var filterFunction;
  if (!isCustomSetGlobal) {
    var dayFilter = 'D' + time;
    var typeRegex;
    if (time > splitTime) {
      typeRegex = new RegExp('.*' + type + '.*', 'i');
    } else {
      typeRegex = new RegExp('.*dox.*', 'i');
    }
    filterFunction = function (x) {
      return x.startsWith(dayFilter) && x.match(typeRegex);
    };
  } else {
    // DiPSC
    var typeRegex = new RegExp('.*dox.*|.*' + type + '.*', 'i');
    filterFunction = function (x) {
      return x.match(typeRegex);
    };
  }

  d.done(function (result) {
    var allValues = result.values;
    var arrayCellIds = result.ids;
    var contains = type;

    var ids = [];
    var values = [];
    for (var i = 0; i < allValues.length; i++) {
      var id = arrayCellIds[i];
      if (filterFunction(id)) {
        values.push(allValues[i]);
        ids.push(id);
      }
    }
    if (isCustomSetGlobal) {
      var dayToCount = type.toLowerCase() === 'serum' ? dayToCountSerum : dayToCount2i;
      if (dayToCount == null) {
        dayToCount = {};
        for (var i = 0; i < ids.length; i++) {
          var id = ids[i];
          var day = id.split('_')[0];

          var count = dayToCount[day];
          if (count == null) {
            dayToCount[day] = 0;
          }
          dayToCount[day]++;
        }
        if (type.toLowerCase() === 'serum') {
          dayToCountSerum = dayToCount;
        } else {
          dayToCount2i = dayToCount;
        }
      }
    }
    var quantileValue = d3.quantile(values.slice(0).sort(function (a, b) {
      return (a === b ? 0 : (a < b ? -1 : 1));
    }), quantile);

    var selectedValues = [];
    customCellSetIds = [];
    var selectedForceLayoutX = [];
    var selectedForceLayoutY = [];

    for (var i = 0; i < ids.length; i++) {
      var id = ids[i];
      if (values[i] > quantileValue) {
        customCellSetIds.push(id);
        selectedValues.push(values[i]);
        var xy = cellIdToXY[id];
        selectedForceLayoutX.push(xy[0]);
        selectedForceLayoutY.push(xy[1]);
      }
    }

    if (isCustomSetGlobal) {
      var formatter = d3.format('.2f');
      var groupedThousands = d3.format(',');
      var dayToIds = splitByTime();
      var dayToCount = type.toLowerCase() === 'serum' ? dayToCountSerum : dayToCount2i;
      var html = [];
      html.push('<table class="table"><tr><th>Day</th><th># Cells pass</th><th>% Cells pass</th></tr>');
      var totalPass = 0;
      var total = 0;
      for (var day in dayToIds) {
        html.push('<tr>');
        var ids = dayToIds[day];
        html.push('<td>');
        html.push(day.substring(1));
        html.push('</td>');
        html.push('<td>');
        html.push(groupedThousands(ids.length) + '/' + groupedThousands(dayToCount[day]));
        html.push('</td>');
        html.push('<td>');
        html.push(formatter(100 * (ids.length / dayToCount[day])));
        html.push('</td>');
        html.push('</tr>');
        totalPass += ids.length;
        total += dayToCount[day];
      }
      html.push('<td>All</td>');
      html.push('<td>');
      html.push(groupedThousands(totalPass) + '/' + groupedThousands(total));
      html.push('</td>');
      html.push('<td>');
      html.push(formatter(100 * (totalPass / total)));
      html.push('</td>');
      html.push('</tr>');
      html.push('</table>');
      $('#filter_summary').html(html.join(''));
    }
    var traces = [
      {
        name: title,
        x: allValues,
        type: 'histogram',
        marker: {color: 'LightGrey', opacity: 0.7}
      }];
    if (!isCustomSetGlobal) {
      traces.push({
        name: title + ' Day ' + time + (time > splitTime ? ' ' + type : ''),
        x: values,
        type: 'histogram',
        marker: {color: '#d95f02'}
      });
    }
    traces.push({
      name: title + (!isCustomSetGlobal ? ' Day ' + time + (time > splitTime ? ' ' + type : '') : '') + ' q ' + q +
      '<br />' +
      customCellSetIds.length + ' cells',
      x: selectedValues,
      type: 'histogram',
      marker: {color: '#1b9e77'}
    });
    $('#histogram').empty();
    Plotly.newPlot('histogram', traces, {margin: {t: 0, b: 15}});
    var forceLayoutTraces = [
      {
        hoverinfo: 'none',
        showlegend: false,
        marker: {symbol: 'square', size: 2, color: 'rgb(217,217,217)', opacity: 0.5},
        mode: 'markers',
        type: 'scattergl',
        name: 'All Cells',
        x: forceLayoutX,
        y: forceLayoutY
      },
      {
        hoverinfo: 'none', showlegend: false, marker: {size: 3, color: '#1b9e77', symbol: 'square'},
        mode: 'markers', type: 'scattergl', name: 'Selected Cells', x: selectedForceLayoutX, y: selectedForceLayoutY
      }];
    $('#force').empty();
    Plotly.newPlot('force', forceLayoutTraces, {
      xaxis: {
        autorange: true,
        showgrid: false,
        zeroline: false,
        showline: false,
        autotick: true,
        ticks: '',
        showticklabels: false
      },
      yaxis: {
        autorange: true,
        showgrid: false,
        zeroline: false,
        showline: false,
        autotick: true,
        ticks: '',
        showticklabels: false
      },
      margin: {t: 0}
    });

  });

};
$customTime.on('change', updateCustomCellSetsPlot);
$customGeneSet.on('change', updateCustomCellSetsPlot);
$customType.on('change', updateCustomCellSetsPlot);
$quantileValue.on('keyup', function (e) {
  if (e.which === 13) {
    updateCustomCellSetsPlot();
  }
});
forceLayoutIds = [];
forceLayoutCounts = [];
forceLayoutX = [];
forceLayoutY = [];
var opacityScale = d3.scaleLinear().range([0.2, 0.9]).domain([10, 80]).clamp(true);
$.ajax({url: '/static/FLE.overview.txt', dataType: 'text'}).done(function (txt) {
  var lines = txt.split('\n');
  var tab = /\t/;
  for (var i = 1; i < lines.length; i++) {
    var tokens = lines[i].split(tab);

    forceLayoutX.push(tokens[0]);
    forceLayoutY.push(tokens[1]);
    //  forceLayoutCounts.push(opacityScale(parseInt(tokens[2])));
  }
});
cellIdToXY = {};
$.ajax({url: '/static/FLE.cells.txt', dataType: 'text'}).done(function (txt) {
  var lines = txt.split('\n');
  var tab = /\t/;
  for (var i = 1; i < lines.length; i++) {
    var tokens = lines[i].split(tab);
    cellIdToXY[tokens[0]] = [tokens[1], tokens[2]];
  }
});

$.ajax('/list_genes/').done(function (results) {

  $gene
  // don't navigate away from the field on tab when selecting an item
    .on('keydown', function (event) {
      if (event.keyCode === $.ui.keyCode.TAB &&
        $(this).autocomplete('instance').menu.active) {
        event.preventDefault();
      }
    }).autocomplete({
    minLength: 2,
    source: function (request, response) {
      // delegate back to autocomplete, but extract the last term
      response(autocompleteFilter(
        results, extractLast(request.term)));
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
});

var $violin = $('#violin');
var $loading = $('#loading');
var $ncells = $('#ncells');

var createForceLayoutTrajectory = function (force, key, isSerum) {

  var traces = force[key];
  var $div = $('<li><h4>' + key + ' ' + (isSerum ? 'Serum' : '2i') +
    '</h4><div class="plot"></div><div data-name="controls"></div></li>'
  );
  var max = -Number.MAX_VALUE;
  var min = Number.MAX_VALUE;

  // traces.forEach(function (trace) {
  //
  //   trace.marker.color.forEach(function (v) {
  //     max = Math.max(max, v);
  //     min = Math.min(min, v);
  //   });
  // });

  traces.forEach(function (trace) {
    trace.mode = 'markers';
    trace.type = 'scattergl';
    trace.hoverinfo = 'none';
    trace.showlegend = false;
    var max = -Number.MAX_VALUE;
    var min = Number.MAX_VALUE;

    trace.marker.color.forEach(function (v) {
      max = Math.max(max, v);
      min = Math.min(min, v);
    });
    trace.marker = {
      symbol: 'dot',
      cmax: max,
      cmin: min,
      showscale: true,
      colorscale: 'Hot',
      size: 2,
      color: trace.marker.color
    };
  });
  var $controls = $div.find('[data-name=controls]');
  $('<button class="btn btn-default" name="play">Play</button>  <label data-name="time"></label>').
    appendTo($controls);
  $controls.after($div.find('.plot'));
  $div.appendTo($violin);

  var layout =
    {
      xaxis: {
        autorange: true,
        showgrid: false,
        zeroline: false,
        showline: false,
        autotick: true,
        ticks: '',
        showticklabels: false
      },
      yaxis: {
        autorange: true,
        showgrid: false,
        zeroline: false,
        showline: false,
        autotick: true,
        ticks: '',
        showticklabels: false
      },
      width: 500,
      height: 400,
      margin: {
        l: 0,
        r: 0,
        b: 0,
        t: 0,
        pad: 0
      },
      autosize: false
      // updatemenus: [
      //   {
      //     x: 0,
      //     y: 0,
      //     yanchor: 'top',
      //     xanchor: 'left',
      //     showactive: false,
      //     direction: 'left',
      //     type: 'buttons',
      //     pad: {t: 87, r: 10},
      //     buttons: [
      //       {
      //         method: 'animate',
      //         args: [
      //           null, {
      //             mode: 'immediate',
      //             fromcurrent: true,
      //             transition: {duration: 300},
      //             frame: {duration: 500, redraw: false}
      //           }],
      //         label: 'Play'
      //       }, {
      //         method: 'animate',
      //         args: [
      //           [null], {
      //             mode: 'immediate',
      //             transition: {duration: 0},
      //             frame: {duration: 0, redraw: false}
      //           }],
      //         label: 'Pause'
      //       }]
      //   }],
      // // Finally, add the slider and use `pad` to position it
      // // nicely next to the buttons.
      // sliders: [
      //   {
      //     pad: {l: 130, t: 55},
      //     currentvalue: {
      //       visible: true,
      //       prefix: 'Time: ',
      //       xanchor: 'right',
      //       font: {size: 20, color: '#666'}
      //     },
      //     steps: sliderSteps
      //   }]
    };
  var elem = $div.find('.plot')[0];

  var backgroundTrace = {
    hoverinfo: 'none', showlegend: false, marker: {symbol: 'square', size: 2, color: 'rgb(217,217,217)', opacity: 0.5},
    mode: 'markers', type: 'scattergl', name: 'All Cells', x: forceLayoutX, y: forceLayoutY
  };

  var allTraces = [backgroundTrace].concat([traces[0]]);
  Plotly.react(elem, {
    data: allTraces,
    layout: layout
  });

  var index = 0;

  function update () {
    if ($playBtn.text() === 'Pause') {
      index++;
      if (index >= traces.length) {
        index = 0;
      }
      $time.text('Time: ' + traces[index].t);
      Plotly.newPlot(elem, {
        data: [backgroundTrace].concat([traces[index]]),
        layout: layout
      });

      setTimeout(update, 1000);
    }
  }

  var $playBtn = $controls.find('[name=play]');
  $playBtn.on('click', function () {
    if ($playBtn.text() === 'Play') {
      $playBtn.text('Pause');
      update();
    } else {
      $playBtn.text('Play');
    }

  });
  var $time = $controls.find('[data-name=time]');
  $time.text('Time: ' + traces[0].t);

};
var fetchData = function (cellSets, genes, geneSets, isSerum) {

  var url = {};
  var predefinedSets = [];
  var customCounter = 0;

  var ncells = parseInt($ncells.val());
  cellSets.forEach(function (cellSet) {
    if (cellSet.custom) {
      var cellIds = customCellSets[cellSet.name];
      url['id' + customCounter] = cellIds;
      url['name' + customCounter] = cellSet.name;
      customCounter++;
    } else {
      predefinedSets.push(cellSet.name);
    }
  });
  url.ncells = '' + ncells;
  url.ncustom = '' + customCounter;
  url.cell_set = predefinedSets;
  if (genes != null) {
    url.gene = genes;
  }
  if (geneSets != null) {
    url.score = geneSets;
  }
  url.transport_map = isSerum ? 'serum' : '2i';

  return $.ajax({url: '/trajectory/', method: 'POST', data: url}).done(function (result) {

    var violins = result.violins;
    var scatters = result.gene_traces;
    var ancestry_divergence_traces = result.ancestry_divergence_traces;
    var force = result.force;
    // also separate median
    for (var key in violins) {
      var traces = violins[key];
      traces.forEach(function (trace) {
        trace.bandwidth = Math.pow(ncells, (-1. / (1 + 4)));
      });
      var $div = $('<li><h4>' + key + ' ' + (isSerum ? 'Serum' : '2i') +
        '</h4><div class="plot"></div></li>'
      );
      $div.appendTo($violin);

      Plotly.newPlot($div.find('.plot')[0], traces, {
        title: '', yaxis: {autorange: true, 'zeroline': false}
      }, {
        modeBarButtonsToRemove: ['toImage'],
        modeBarButtonsToAdd: [
          {
            name: 'Save Image',
            icon: Plotly.Icons.camera,
            click: function (gd) {
              Plotly.downloadImage(gd, {format: 'svg', width: $(gd).width(), height: $(gd).height()});
            }
          }]
      });
    }
    if (ancestry_divergence_traces && ancestry_divergence_traces.length > 0) {
      var $div = $('<li><h4>' + (isSerum ? 'Serum' : '2i') +
        '</h4><div class="plot"></div></li>'
      );
      $div.appendTo($violin);

      Plotly.newPlot($div.find('.plot')[0], ancestry_divergence_traces, {
          title: '', yaxis: {
            range: [0, 1], autorange: false, 'zeroline': false,
            title: 'Divergence'
          },
          showlegend: true
        }, {xaxis: {title: 'Time'}}, {
          modeBarButtonsToRemove: ['toImage'],
          modeBarButtonsToAdd: [
            {
              name: 'Save Image',
              icon: Plotly.Icons.camera,
              click: function (gd) {
                Plotly.downloadImage(gd, {format: 'svg', width: $(gd).width(), height: $(gd).height()});
              }
            }]
        }
      );
    }
    if (result.line_traces && result.line_traces.length > 0) {
      var $div = $('<li><h4>' + (isSerum ? 'Serum' : '2i') +
        '</h4><div class="plot"></div></li>'
      );
      $div.appendTo($violin);

      Plotly.newPlot($div.find('.plot')[0], result.line_traces, {
          title: '', yaxis: {
            autorange: true, 'zeroline': false
          },
          showlegend: true
        }, {
          modeBarButtonsToRemove: ['toImage'],
          modeBarButtonsToAdd: [
            {
              name: 'Save Image',
              icon: Plotly.Icons.camera,
              click: function (gd) {
                Plotly.downloadImage(gd, {format: 'svg', width: $(gd).width(), height: $(gd).height()});
              }
            }]
        }
      );
    }
    if (scatters.length > 0) {
      var $div = $('<li><h4>' + (isSerum ? 'Serum' : '2i') +
        '</h4><div class="plot"></div></li>'
      );
      $div.appendTo($violin);

      Plotly.newPlot($div.find('.plot')[0], scatters, {
        title: '', yaxis: {autorange: true, 'zeroline': false}, showlegend: true
      }, {
        modeBarButtonsToRemove: ['toImage'],
        modeBarButtonsToAdd: [
          {
            name: 'Save Image',
            icon: Plotly.Icons.camera,
            click: function (gd) {
              Plotly.downloadImage(gd, {format: 'svg', width: $(gd).width(), height: $(gd).height()});
            }
          }]
      });
    }
    for (var key in force) {
      createForceLayoutTrajectory(force, key, isSerum);
    }
    $violin.sortable({handle: 'h4'});
  });
};

var updateViolinPlots = function () {
  var predefinedSets = $cellSet.val();
  var customCellSets = $customCellSetSelector.val();

  var cellSets = [];
  cellSets = cellSets.concat(predefinedSets.map(function (s) {
    return {name: s, custom: false};
  }));
  cellSets = cellSets.concat(customCellSets.map(function (s) {
    return {name: s, custom: true};
  }));
  if (cellSets.length === 0) {
    return;
  }
  var geneSets = $geneSet.val();
  var genes = $gene.val().split(',');
  var _genes = [];
  genes.forEach(function (gene) {
    var val = gene.trim();
    if (val !== '') {
      _genes.push(val);
    }
  });
  genes = _genes;
  $violin.html('');
  var serumCellSets = [];
  var twoISets = [];
  cellSets.forEach(function (cellSet) {
    var calls = [];
    if (cellSet.name.toLowerCase().indexOf('2i') !== -1) {
      twoISets.push(cellSet);
    } else if (cellSet.name.toLowerCase().indexOf('serum') !== -1) {
      serumCellSets.push(cellSet);
    } else {
      twoISets.push(cellSet);
      serumCellSets.push(cellSet);
    }
  });
  $violin.empty();
  $loading.show();
  var promises = [];
  if (serumCellSets.length > 0) {
    promises.push(fetchData(serumCellSets, genes, geneSets, true));
  }
  if (twoISets.length > 0) {
    promises.push(fetchData(twoISets, genes, geneSets, false));
  }
  $.when.apply($, promises).done(function () {
    $loading.hide();
  });

};
$('#go').on('click', function (e) {
  e.preventDefault();
  updateViolinPlots();
});
$('#ancestor_form').on('submit', function (e) {
  e.preventDefault();
  updateViolinPlots();
});



































