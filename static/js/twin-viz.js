// Cumulative N(H) plot for the H test (twin-check.js), in the style of the
// N(z) plot (e-statistics.js). For twin fraction alpha, H = |I1-I2|/(I1+I2)
// of twin-related acentric pairs is uniform on [0, 1 - 2 alpha], so
// N(H) = H / (1 - 2 alpha): straight lines, steeper for larger alpha. The
// observed curve of each twin law is drawn on top; one that bends away from
// the straight lines points to pseudo-symmetry rather than a clean twin.
let _lastNHPlotArgs = null;

function renderNHPlot(container, laws) {
  _lastNHPlotArgs = [container, laws];
  const canvas = document.createElement('canvas');
  container.innerHTML = '';
  container.appendChild(canvas);
  drawNHPlot(canvas, laws);
  attachPlotZoom(canvas, '<strong>N(H) — H test</strong>', 480, (cv, zoom) => drawNHPlot(cv, laws, zoom));
}

function drawNHPlot(canvas, laws, zoom = 1) {
  const colors = getPlotColors();
  const cssW = 480, cssH = 320;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);
  const margin = { left: 50, right: 15, top: 15, bottom: 40 };
  const plotW = cssW - margin.left - margin.right;
  const plotH = cssH - margin.top - margin.bottom;
  const HMAX = 1;
  const px = h => margin.left + (h / HMAX) * plotW;
  const py = y => margin.top + plotH - y * plotH;

  ctx.strokeStyle = colors.grid; ctx.fillStyle = colors.text; ctx.font = '10px sans-serif'; ctx.lineWidth = 1;
  for (let i = 0; i <= 5; i++) {
    const x = px(i / 5), y = py(i / 5);
    ctx.beginPath(); ctx.moveTo(x, margin.top); ctx.lineTo(x, margin.top + plotH); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText((i / 5).toFixed(1), x, cssH - margin.bottom + 14);
    ctx.beginPath(); ctx.moveTo(margin.left, y); ctx.lineTo(margin.left + plotW, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText((i / 5).toFixed(1), margin.left - 6, y + 3);
  }

  // reference lines for alpha = 0 ... 0.4, labelled where they reach N = 1
  ctx.font = '9px sans-serif';
  for (const alpha of [0, 0.1, 0.2, 0.3, 0.4]) {
    const end = 1 - 2 * alpha;
    ctx.strokeStyle = colors.axis; ctx.lineWidth = 1; ctx.setLineDash(alpha === 0 ? [] : [4, 3]);
    ctx.beginPath(); ctx.moveTo(px(0), py(0)); ctx.lineTo(px(end), py(1)); ctx.lineTo(px(HMAX), py(1)); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = colors.text; ctx.textAlign = 'center';
    ctx.fillText(`α=${alpha.toFixed(1)}`, px(end), py(1) - 4);
  }

  // observed curves, one per twin law
  const lineColors = [colors.point, colors.fitLine, colors.present, colors.symFilled || colors.axis];
  const shown = laws.filter(L => L.hValues && L.hValues.length >= 20).slice(0, lineColors.length);
  shown.forEach((L, i) => {
    const v = [...L.hValues].sort((a, b) => a - b);
    ctx.strokeStyle = lineColors[i]; ctx.lineWidth = 2.5;
    ctx.beginPath();
    const steps = 200;
    for (let g = 0; g <= steps; g++) {
      const h = HMAX * g / steps;
      // fraction of values <= h (binary search)
      let lo = 0, hi = v.length;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (v[mid] <= h) lo = mid + 1; else hi = mid; }
      const n = lo / v.length;
      g === 0 ? ctx.moveTo(px(h), py(n)) : ctx.lineTo(px(h), py(n));
    }
    ctx.stroke();
  });

  // legend, bottom right
  ctx.font = '10px sans-serif'; ctx.textAlign = 'left';
  shown.forEach((L, i) => {
    const y = margin.top + plotH - 12 - (shown.length - 1 - i) * 14;
    const x0 = margin.left + plotW - 205;
    ctx.strokeStyle = lineColors[i]; ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.moveTo(x0, y - 3); ctx.lineTo(x0 + 22, y - 3); ctx.stroke();
    ctx.fillStyle = colors.text;
    ctx.fillText(`TWIN ${L.shelxMatrix.flat().join(' ')}`, x0 + 28, y);
  });

  ctx.fillStyle = colors.text; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('H = |I₁ − I₂| / (I₁ + I₂)', margin.left + plotW / 2, cssH - 6);
  ctx.save(); ctx.translate(12, margin.top + plotH / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText('N(H) = fraction of pairs with H′ ≤ H', 0, 0); ctx.restore();
}

onPlotThemeChange(() => {
  if (_lastNHPlotArgs) renderNHPlot(..._lastNHPlotArgs);
});
