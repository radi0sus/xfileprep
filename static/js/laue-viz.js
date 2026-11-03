// laue-viz.js
// "Laue" tab: intensity-weighted reciprocal-lattice pictures of the three
// principal zones (0kl, h0l, hk0), drawn with the true reciprocal metric.
// The symmetry of the intensity pattern is what defines the Laue class,
// so this is the qualitative companion to the Rint table in tab 1.
//
// Intensities are averaged ONLY over exact duplicate indices — never
// over symmetry equivalents — otherwise the very symmetry being judged
// would be averaged in. Reuses ZONE_SELECTORS (absences.js) and
// ZONE_PLOT_CONFIG (absences-viz.js).

const LAUE_ZONE_ORDER = ['0kl', 'h0l', 'hk0'];
let _lastLaueArgs = null;

// Marker radius for mean intensity I relative to the reference intensity.
// 'area': diameter ∝ sqrt(I) (area ∝ I); 'linear': diameter ∝ I; 'log': compressed.
function laueRadius(I, Iref, rmax, mode) {
  if (!(I > 0)) return 0;
  const x = Math.min(I / Iref, 1);
  let r;
  if (mode === 'linear') r = rmax * x;
  else if (mode === 'log') r = rmax * Math.log10(1 + 999 * x) / 3;
  else r = rmax * Math.sqrt(x);
  return Math.max(r, 0.8);
}

function aggregateZoneIntensity(reflections, cfg) {
  const map = new Map();
  for (const r of reflections) {
    if (!cfg.select(r)) continue;
    const u = cfg.xOf(r), v = cfg.yOf(r), id = u + ',' + v;
    let e = map.get(id);
    if (!e) { e = { u, v, sum: 0, n: 0 }; map.set(id, e); }
    e.sum += r.f2; e.n++;
  }
  return [...map.values()].map(e => ({ u: e.u, v: e.v, I: e.sum / e.n, n: e.n }));
}

function drawLaueZone(canvas, points, cfg, cell, Iref, mode, zoom = 1) {
  const colors = getPlotColors();
  const cssW = 340, cssH = 340;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);
  if (zoom === 1) { canvas.style.width = '100%'; canvas.style.height = 'auto'; } // scale with the grid cell (aspect ratio kept)
  if (points.length === 0) {
    ctx.fillStyle = colors.text;
    ctx.font = '12px sans-serif';
    ctx.fillText(`No reflections in zone (${cfg.fixed}).`, 10, cssH / 2);
    return;
  }
  const us = points.map(p => p.u), vs = points.map(p => p.v);
  const r = { uMin: Math.min(...us) - 1, uMax: Math.max(...us) + 1, vMin: Math.min(...vs) - 1, vMax: Math.max(...vs) + 1 };
  const basis = reciprocalPlaneBasis(cell, cfg.uIdx, cfg.vIdx);
  const frame = makeLatticeFrame(basis, r.uMin, r.uMax, r.vMin, r.vMax, { x: 34, y: 12, w: cssW - 46, h: cssH - 40 });
  drawLatticeGrid(ctx, colors, frame, r, { x: cfg.xLabel, y: cfg.yLabel },
    (lo, hi) => integerTicks(lo, hi, niceIntStep(hi - lo)), cssW, cssH);

  const rmax = Math.max(2, frame.spacing * 0.5);
  ctx.fillStyle = colors.point;
  ctx.strokeStyle = colors.text;
  ctx.lineWidth = 1;
  for (const p of points) {
    const x = frame.px(p.u, p.v), y = frame.py(p.u, p.v);
    const rad = laueRadius(p.I, Iref, rmax, mode);
    ctx.beginPath();
    if (rad > 0) { ctx.arc(x, y, rad, 0, 2 * Math.PI); ctx.fill(); }
    else { ctx.arc(x, y, 1.2, 0, 2 * Math.PI); ctx.fill(); } // I <= 0: just a small dot
  }
}

function renderLaueOverview(containerEl, noteEl, reflections, cell, laueClass) {
  _lastLaueArgs = [containerEl, noteEl, reflections, cell, laueClass];
  const modeEl = document.getElementById('laueScaleSelect');
  const mode = modeEl ? modeEl.value : 'area';

  // redundancy of the raw data: tells whether symmetry equivalents / repeats are present at all
  const seen = new Set();
  for (const r of reflections) seen.add(r.h + ',' + r.k + ',' + r.l);
  const redundancy = reflections.length / Math.max(1, seen.size);
  let note = laueClass
    ? `<p>Laue class from merging: <strong>${laueClass}</strong></p>`
    : `<p><em>Laue class not determined from merging.</em></p>`;
  if (redundancy < 1.05) {
    note += `<p style="color:#a60;"><small>⚠ Almost no repeated observations (mean redundancy ${redundancy.toFixed(2)}) — the data look already merged to a unique set, ` +
      `so symmetry-equivalent reflections are missing and the pattern will not show the full Laue symmetry.</small></p>`;
  }
  noteEl.innerHTML = note;

  const zones = LAUE_ZONE_ORDER.map(z => ({ name: z, cfg: ZONE_PLOT_CONFIG[z], pts: aggregateZoneIntensity(reflections, ZONE_PLOT_CONFIG[z]) }));
  // one shared reference intensity (99th percentile of positive means) so the three zones are comparable
  const pos = zones.flatMap(z => z.pts.map(p => p.I)).filter(I => I > 0).sort((a, b) => a - b);
  const Iref = pos.length ? pos[Math.min(pos.length - 1, Math.floor(pos.length * 0.99))] : 1;

  containerEl.innerHTML = '';
  for (const z of zones) {
    const wrap = document.createElement('div');
    wrap.className = 'zonePlotWrap';
    const caption = document.createElement('div');
    caption.className = 'plotCaption';
    caption.innerHTML = `<strong>${z.name}</strong> <small>(${z.cfg.fixed}, ${z.pts.length} points)</small>`;
    const canvas = document.createElement('canvas');
    wrap.appendChild(caption);
    wrap.appendChild(canvas);
    containerEl.appendChild(wrap);
    drawLaueZone(canvas, z.pts, z.cfg, cell, Iref, mode);
    attachPlotZoom(wrap, caption.innerHTML, 340,
      (cv, zoom) => drawLaueZone(cv, z.pts, z.cfg, cell, Iref, mode, zoom));
  }
}

function rerenderLaueOverview() {
  if (_lastLaueArgs) renderLaueOverview(..._lastLaueArgs);
}
(function () {
  const sel = document.getElementById('laueScaleSelect');
  if (sel) sel.addEventListener('change', rerenderLaueOverview);
})();
onPlotThemeChange(rerenderLaueOverview);
