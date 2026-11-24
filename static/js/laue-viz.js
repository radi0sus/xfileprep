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

function laueFillSymmetryEnabled() {
  const el = document.getElementById('laueFillSym');
  return !!(el && el.checked);
}

// Optional display aid for MERGED data (one reflection per symmetry class,
// so most of each zone is empty): returns the positions that were NOT
// measured themselves but have a measured symmetry equivalent under the
// Laue group, each with the mean intensity of the measured members of its
// class, flagged `generated: true`. Measured positions are never touched —
// they keep their own (repeat-averaged) intensity, so the measured picture
// stays exactly as it was and the filled-in points are drawn in a separate
// colour. `group`: the Laue group as 3x3 matrices (generateGroup), applied
// as h' = h M (the convention of canonicalHKL). Off by default: switched
// on, the plots show the assumed symmetry, not the measured one.
const _symFillCache = { reflections: null, group: null, result: null };
function symmetryFilledReflections(reflections, group) {
  if (_symFillCache.reflections === reflections && _symFillCache.group === group) return _symFillCache.result;
  const measured = new Set();
  const classes = new Map(); // canonical key -> { sum, n, members: Set }
  for (const r of reflections) {
    measured.add(r.h + ',' + r.k + ',' + r.l);
    const c = canonicalHKL(r.h, r.k, r.l, group);
    const key = c.join(',');
    let e = classes.get(key);
    if (!e) { e = { sum: 0, n: 0, h: r.h, k: r.k, l: r.l }; classes.set(key, e); }
    e.sum += r.f2; e.n++;
  }
  const out = [];
  for (const e of classes.values()) {
    const I = e.sum / e.n;
    const seen = new Set();
    for (const M of group) {
      const h = e.h * M[0][0] + e.k * M[1][0] + e.l * M[2][0];
      const k = e.h * M[0][1] + e.k * M[1][1] + e.l * M[2][1];
      const l = e.h * M[0][2] + e.k * M[1][2] + e.l * M[2][2];
      const id = h + ',' + k + ',' + l;
      if (seen.has(id) || measured.has(id)) continue;
      seen.add(id);
      out.push({ h, k, l, f2: I, sigma: 0, generated: true });
    }
  }
  _symFillCache.reflections = reflections; _symFillCache.group = group; _symFillCache.result = out;
  return out;
}

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

// Filled-in positions of one zone (see symmetryFilledReflections); never
// overlap a measured position by construction.
function aggregateZoneGenerated(generated, cfg) {
  return generated.filter(cfg.select).map(r => ({ u: cfg.xOf(r), v: cfg.yOf(r), I: r.f2, generated: true }));
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
    ctx.fillStyle = p.generated ? colors.symFilled : colors.point;
    ctx.beginPath();
    if (rad > 0) { ctx.arc(x, y, rad, 0, 2 * Math.PI); ctx.fill(); }
    else { ctx.arc(x, y, 1.2, 0, 2 * Math.PI); ctx.fill(); } // I <= 0: just a small dot
  }
}

function renderLaueOverview(containerEl, noteEl, reflections, cell, laueClass, laueSource = 'merging', laueGroup = null) {
  _lastLaueArgs = [containerEl, noteEl, reflections, cell, laueClass, laueSource, laueGroup];
  const fill = laueFillSymmetryEnabled() && !!laueGroup;
  const generated = fill ? symmetryFilledReflections(reflections, laueGroup) : [];
  const modeEl = document.getElementById('laueScaleSelect');
  const mode = modeEl ? modeEl.value : 'area';

  // redundancy of the raw data: tells whether symmetry equivalents / repeats are present at all
  const seen = new Set();
  for (const r of reflections) seen.add(r.h + ',' + r.k + ',' + r.l);
  const redundancy = reflections.length / Math.max(1, seen.size);
  let note = laueClass
    ? `<p>Laue class ${laueSource === 'cif' ? 'from the CIF space group' : 'from merging'}: <strong>${formatPG(laueClass)}</strong></p>`
    : `<p><em>Laue class not determined from merging.</em></p>`;
  if (redundancy < 1.05) {
    note += `<p style="color:#a60;"><small>⚠ Almost no repeated observations (mean redundancy ${redundancy.toFixed(2)}) — the data look already merged to a unique set, ` +
      `so symmetry-equivalent reflections are missing and the pattern will not show the full Laue symmetry.` +
      (fill ? '' : ` Tick <em>fill in symmetry equivalents</em> above to complete the picture from the Laue class.`) + `</small></p>`;
  }
  if (fill) {
    note += `<p><small><span class="symFilledKey">●</span> orange = not measured itself, filled in from its measured symmetry equivalents under ${laueClass ? formatPG(laueClass) : 'the assumed Laue class'} ` +
      `(${generated.length} positions overall). This shows the <em>assumed</em> symmetry, not a test of it — switch it off to judge the Laue class.</small></p>`;
  }
  noteEl.innerHTML = note;

  const zones = LAUE_ZONE_ORDER.map(z => {
    const measuredPts = aggregateZoneIntensity(reflections, ZONE_PLOT_CONFIG[z]);
    const genPts = fill ? aggregateZoneGenerated(generated, ZONE_PLOT_CONFIG[z]) : [];
    return { name: z, cfg: ZONE_PLOT_CONFIG[z], pts: measuredPts.concat(genPts), nMeasured: measuredPts.length, nGen: genPts.length };
  });
  // one shared reference intensity (99th percentile of positive means, measured points only — the
  // scale must not change when the fill-in is switched on) so the three zones are comparable
  const pos = zones.flatMap(z => z.pts.filter(p => !p.generated).map(p => p.I)).filter(I => I > 0).sort((a, b) => a - b);
  const Iref = pos.length ? pos[Math.min(pos.length - 1, Math.floor(pos.length * 0.99))] : 1;

  containerEl.innerHTML = '';
  for (const z of zones) {
    const wrap = document.createElement('div');
    wrap.className = 'zonePlotWrap';
    const caption = document.createElement('div');
    caption.className = 'plotCaption';
    caption.innerHTML = `<strong>${z.name}</strong> <small>(${z.cfg.fixed}, ${z.nMeasured} points${z.nGen ? ` + <span class="symFilledKey">${z.nGen} filled in</span>` : ''})</small>`;
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
  const fill = document.getElementById('laueFillSym');
  if (fill) fill.addEventListener('change', rerenderLaueOverview);
})();
onPlotThemeChange(rerenderLaueOverview);
