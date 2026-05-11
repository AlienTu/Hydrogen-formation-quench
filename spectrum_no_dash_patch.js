// Keep both dphi1 and dphi2 visible as ordinary solid curves in the selected mode plot.
// This is a display-only patch; it does not modify the eigensolver or the mode vectors.
(function(){
  if (!window.Plotly || !Plotly.newPlot) return;
  const originalNewPlot = Plotly.newPlot.bind(Plotly);
  Plotly.newPlot = function(id, data, layout, config){
    const target = (typeof id === 'string') ? id : (id && id.id);
    if (target === 'modePlot' && Array.isArray(data)) {
      data.forEach(trace => {
        if (trace && trace.line && Object.prototype.hasOwnProperty.call(trace.line, 'dash')) {
          delete trace.line.dash;
        }
      });
    }
    return originalNewPlot(id, data, layout, config);
  };
})();
