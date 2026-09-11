// absences-viz.js
// Graphical companion to absences.js: renders the systematic-absence
// evidence as reciprocal-lattice pictures instead of just tables, so
// patterns (checkerboards from glide planes, missing-every-other-spot
// from screw axes) are visible at a glance, precession-photo style.
//
// Deliberately reuses the exact same selectors/index functions that
// absences.js already uses for its statistical scan (ZONE_SELECTORS,
// AXIAL_CLASSES, isWeak) — the picture and the table are guaranteed to
// agree because they're built from the same filters, not a second
// hand-copied set of conditions.
//
// Unmerged data means several observations can land on the same
// integer grid point (duplicate/symmetry-equivalent measurements).
// Rather than plotting every raw point (slow, and visually just a
// blur of overlapping dots), reflections are aggregated per grid cell
// first: a cell counts as "present" if ANY observation there is
// significant, "absent" if every observation there is weak. That's
// the same all-or-nothing logic a person reads off a real precession
// photograph.

// Colors are resolved per-render via getPlotColors() (plot-theme.js) so
// each plot follows the OS light/dark setting instead of a fixed palette.

function setupCanvasForDPR(canvas, cssWidth, cssHeight, colors) {
  const dpr = window.devicePixelRatio || 1;
  canvas.width = cssWidth * dpr;
  canvas.height = cssHeight * dpr;
  canvas.style.width = cssWidth + 'px';
  canvas.style.height = cssHeight + 'px';
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssWidth, cssHeight);
  ctx.fillStyle = colors.bg;
  ctx.fillRect(0, 0, cssWidth, cssHeight);
  return ctx;
}

// Picks a "nice" integer tick step (1, 2, 5, 10, 20, 50, ...) so axis
// labels stay readable regardless of how wide the index range is.
function niceIntStep(range, targetTicks = 6) {
  if (range <= 0) return 1;
  const rough = range / targetTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(rough)));
  const norm = rough / mag;
  const step = norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10;
  return Math.max(1, Math.round(step * mag));
}

function integerTicks(min, max, step) {
  const ticks = [];
  const start = Math.ceil(min / step) * step;
  for (let t = start; t <= max; t += step) ticks.push(t);
  return ticks;
}

// Aggregates reflections belonging to one zone/axis onto an integer
// grid, keyed by the coordinate(s) that vary in that zone.
function aggregateByKey(relevant, keyOf) {
  const grid = new Map();
  for (const r of relevant) {
    const key = keyOf(r);
    let cell = grid.get(key.id);
    if (!cell) { cell = { ...key, strong: 0, weak: 0 }; grid.set(key.id, cell); }
    if (isWeak(r.f2, r.sigma)) cell.weak++; else cell.strong++;
  }
  return grid;
}

// --- 2D zone plots (h0l, 0kl, hk0) ---------------------------------------

const ZONE_PLOT_CONFIG = {
  h0l: { select: ZONE_SELECTORS.h0l, xOf: r => r.h, yOf: r => r.l, xLabel: 'h', yLabel: 'l', fixed: 'k = 0' },
  '0kl': { select: ZONE_SELECTORS['0kl'], xOf: r => r.k, yOf: r => r.l, xLabel: 'k', yLabel: 'l', fixed: 'h = 0' },
  hk0: { select: ZONE_SELECTORS.hk0, xOf: r => r.h, yOf: r => r.k, xLabel: 'h', yLabel: 'k', fixed: 'l = 0' },
};

function renderZonePlot(canvas, reflections, zoneName) {
  const colors = getPlotColors();
  const cfg = ZONE_PLOT_CONFIG[zoneName];
  const relevant = reflections.filter(cfg.select);
  const cssW = 300, cssH = 300;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors);

  if (relevant.length === 0) {
    ctx.fillStyle = colors.text;
    ctx.font = '12px sans-serif';
    ctx.fillText(`No reflections in zone ${zoneName} (${cfg.fixed}).`, 10, cssH / 2);
    return;
  }

  const grid = aggregateByKey(relevant, r => ({ id: `${cfg.xOf(r)},${cfg.yOf(r)}`, x: cfg.xOf(r), y: cfg.yOf(r) }));
  const cells = [...grid.values()];
  const xs = cells.map(c => c.x), ys = cells.map(c => c.y);
  let xMin = Math.min(...xs) - 1, xMax = Math.max(...xs) + 1;
  let yMin = Math.min(...ys) - 1, yMax = Math.max(...ys) + 1;
  if (xMin === xMax) { xMin -= 1; xMax += 1; }
  if (yMin === yMax) { yMin -= 1; yMax += 1; }

  const margin = { left: 34, right: 10, top: 10, bottom: 28 };
  const plotW = cssW - margin.left - margin.right;
  const plotH = cssH - margin.top - margin.bottom;
  const px = x => margin.left + ((x - xMin) / (xMax - xMin)) * plotW;
  const py = y => margin.top + plotH - ((y - yMin) / (yMax - yMin)) * plotH;

  // gridlines + tick labels
  const xStep = niceIntStep(xMax - xMin);
  const yStep = niceIntStep(yMax - yMin);
  ctx.strokeStyle = colors.grid;
  ctx.fillStyle = colors.text;
  ctx.font = '9px sans-serif';
  ctx.lineWidth = 1;
  for (const t of integerTicks(xMin, xMax, xStep)) {
    const x = px(t);
    ctx.beginPath(); ctx.moveTo(x, margin.top); ctx.lineTo(x, margin.top + plotH); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText(String(t), x, cssH - margin.bottom + 12);
  }
  for (const t of integerTicks(yMin, yMax, yStep)) {
    const y = py(t);
    ctx.beginPath(); ctx.moveTo(margin.left, y); ctx.lineTo(margin.left + plotW, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(String(t), margin.left - 6, y + 3);
  }

  // zero axes, emphasized
  ctx.strokeStyle = colors.axis;
  ctx.lineWidth = 1.2;
  if (xMin <= 0 && xMax >= 0) {
    const x0 = px(0);
    ctx.beginPath(); ctx.moveTo(x0, margin.top); ctx.lineTo(x0, margin.top + plotH); ctx.stroke();
  }
  if (yMin <= 0 && yMax >= 0) {
    const y0 = py(0);
    ctx.beginPath(); ctx.moveTo(margin.left, y0); ctx.lineTo(margin.left + plotW, y0); ctx.stroke();
  }

  // one square per aggregated grid cell
  const cellSize = Math.max(2, Math.min(plotW / (xMax - xMin), plotH / (yMax - yMin)) * 0.75);
  for (const c of cells) {
    ctx.fillStyle = c.strong > 0 ? colors.present : colors.absent;
    ctx.fillRect(px(c.x) - cellSize / 2, py(c.y) - cellSize / 2, cellSize, cellSize);
  }

  // axis labels
  ctx.fillStyle = colors.text;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText(cfg.xLabel, margin.left + plotW / 2, cssH - 4);
  ctx.save();
  ctx.translate(10, margin.top + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(cfg.yLabel, 0, 0);
  ctx.restore();
}

// --- 1D axial strips (h00, 0k0, 00l) -------------------------------------

const AXIAL_PLOT_CONFIG = {
  h00: { select: AXIAL_CLASSES.h00.select, indexOf: AXIAL_CLASSES.h00.indexOf, label: 'h' },
  '0k0': { select: AXIAL_CLASSES['0k0'].select, indexOf: AXIAL_CLASSES['0k0'].indexOf, label: 'k' },
  '00l': { select: AXIAL_CLASSES['00l'].select, indexOf: AXIAL_CLASSES['00l'].indexOf, label: 'l' },
};

function renderAxialStrip(canvas, reflections, axisName) {
  const colors = getPlotColors();
  const cfg = AXIAL_PLOT_CONFIG[axisName];
  const relevant = reflections.filter(cfg.select);
  const cssW = 300, cssH = 70;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors);

  if (relevant.length === 0) {
    ctx.fillStyle = colors.text;
    ctx.font = '12px sans-serif';
    ctx.fillText(`No reflections on the ${axisName} axis.`, 10, cssH / 2);
    return;
  }

  const grid = aggregateByKey(relevant, r => ({ id: String(cfg.indexOf(r)), x: cfg.indexOf(r) }));
  const cells = [...grid.values()];
  const xs = cells.map(c => c.x);
  let xMin = Math.min(...xs) - 1, xMax = Math.max(...xs) + 1;
  if (xMin === xMax) { xMin -= 1; xMax += 1; }

  const margin = { left: 20, right: 20, top: 8, bottom: 22 };
  const plotW = cssW - margin.left - margin.right;
  const midY = margin.top + (cssH - margin.top - margin.bottom) / 2;
  const px = x => margin.left + ((x - xMin) / (xMax - xMin)) * plotW;

  ctx.strokeStyle = colors.axis;
  ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(margin.left, midY); ctx.lineTo(margin.left + plotW, midY); ctx.stroke();

  const step = niceIntStep(xMax - xMin);
  ctx.fillStyle = colors.text;
  ctx.font = '9px sans-serif';
  ctx.textAlign = 'center';
  for (const t of integerTicks(xMin, xMax, step)) {
    const x = px(t);
    ctx.beginPath(); ctx.moveTo(x, midY - 3); ctx.lineTo(x, midY + 3); ctx.stroke();
    ctx.fillText(String(t), x, cssH - 6);
  }

  const dotRadius = Math.max(2, Math.min(6, (plotW / (xMax - xMin)) * 0.35));
  for (const c of cells) {
    ctx.fillStyle = c.strong > 0 ? colors.present : colors.absent;
    ctx.beginPath();
    ctx.arc(px(c.x), midY, dotRadius, 0, 2 * Math.PI);
    ctx.fill();
  }

  ctx.fillStyle = colors.text;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText(cfg.label, margin.left, margin.top + 8);
}

// --- Top-level entry point ------------------------------------------------
// Builds the canvases into the given container elements and draws all
// six plots. Called after a fresh absences scan, using the same
// `reflections` array the scan itself used.
let _lastAbsenceVizArgs = null;

function renderAbsenceVisuals(zoneContainerEl, axialContainerEl, reflections) {
  _lastAbsenceVizArgs = [zoneContainerEl, axialContainerEl, reflections];
  zoneContainerEl.innerHTML = '';
  for (const [zoneName, cfg] of Object.entries(ZONE_PLOT_CONFIG)) {
    const wrap = document.createElement('div');
    wrap.className = 'zonePlotWrap';
    const caption = document.createElement('div');
    caption.className = 'plotCaption';
    caption.innerHTML = `<strong>${zoneName}</strong> <small>(${cfg.fixed})</small>`;
    const canvas = document.createElement('canvas');
    wrap.appendChild(caption);
    wrap.appendChild(canvas);
    zoneContainerEl.appendChild(wrap);
    renderZonePlot(canvas, reflections, zoneName);
  }

  axialContainerEl.innerHTML = '';
  for (const axisName of Object.keys(AXIAL_PLOT_CONFIG)) {
    const wrap = document.createElement('div');
    wrap.className = 'axialStripWrap';
    const caption = document.createElement('div');
    caption.className = 'plotCaption';
    caption.innerHTML = `<strong>${axisName}</strong>`;
    const canvas = document.createElement('canvas');
    wrap.appendChild(caption);
    wrap.appendChild(canvas);
    axialContainerEl.appendChild(wrap);
    renderAxialStrip(canvas, reflections, axisName);
  }
}

onPlotThemeChange(() => {
  if (_lastAbsenceVizArgs) renderAbsenceVisuals(..._lastAbsenceVizArgs);
});
