// Move the IPR colorbar on the spectrum plot to the left side so it does not overlap the legend.
(function(){
  if (!window.Plotly || window.__spectrumColorbarPatch) return;
  window.__spectrumColorbarPatch = true;
  const oldNewPlot = window.Plotly.newPlot.bind(window.Plotly);
  window.Plotly.newPlot = function(div, data, layout, config){
    const id = (typeof div === 'string') ? div : (div && div.id);
    if (id === 'spectrumPlot') {
      if (Array.isArray(data)) {
        data.forEach(trace => {
          if (trace && trace.marker && trace.marker.showscale) {
            trace.marker.colorbar = Object.assign({}, trace.marker.colorbar || {}, {
              title: trace.marker.colorbar && trace.marker.colorbar.title ? trace.marker.colorbar.title : 'IPR',
              x: -0.14,
              xanchor: 'right',
              y: 0.52,
              len: 0.72,
              thickness: 16
            });
          }
        });
      }
      layout = Object.assign({}, layout || {});
      layout.margin = Object.assign({t: 38, l: 95, r: 135, b: 38}, layout.margin || {}, {l: 95, r: 135});
      layout.legend = Object.assign({x: 1.03, y: 1.0, xanchor: 'left', yanchor: 'top'}, layout.legend || {});
    }
    return oldNewPlot(div, data, layout, config);
  };
})();
