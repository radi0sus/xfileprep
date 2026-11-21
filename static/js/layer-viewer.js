// layer-viewer.js
// Layer viewer ("Laue" tab): page through the reciprocal lattice layer by
// layer (hk0, hk1, hk2 … / h0l, h1l … / 0kl, 1kl …) with a slider, drawn in
// the true reciprocal metric like the three zone plots — a bit like
// browsing precession images.
//
//  • Circles = mean intensity per lattice point (averaged ONLY over exact
//    repeats, never over symmetry equivalents; same scale as the zone plots).
//  • Red rings = GAPS: positions the lattice centering allows, inside the
//    measured resolution range, that were not measured — and neither was any
//    symmetry equivalent of them under the Laue group. Positions that the
//    absence scan found to be systematically forbidden are not marked.
//    (A reflection that is missing in the file but whose symmetry equivalent
//    was measured is not a gap: it is filled in by symmetry.)
//
// Context from app.js: { reflections, cell, group, latticeTest, scan, laue, assumed }.
// Needs: reciprocalPlaneBasis/makeLatticeFrame/drawLatticeGrid (recip-geometry.js),
// canonicalHKL (laue-merge.js), isForbiddenByHeldConditions (absences.js),
// sinThetaOverLambda (d-spacing.js), laueRadius (laue-viz.js).

const LAYER_AXES = {
  l: { name: 'hkL', fixed: 'l', u: 'h', v: 'k', uIdx: 0, vIdx: 1, layerOf: r => r.l, uOf: r => r.h, vOf: r => r.k, hkl: (u, v, n) => [u, v, n] },
  k: { name: 'hKl', fixed: 'k', u: 'h', v: 'l', uIdx: 0, vIdx: 2, layerOf: r => r.k, uOf: r => r.h, vOf: r => r.l, hkl: (u, v, n) => [u, n, v] },
  h: { name: 'Hkl', fixed: 'h', u: 'k', v: 'l', uIdx: 1, vIdx: 2, layerOf: r => r.h, uOf: r => r.k, vOf: r => r.l, hkl: (u, v, n) => [n, u, v] },
};

const _layer = { ctx: null, data: null, dataMeasured: null, dataFilled: null, axis: 'l', n: 0, showGaps: true, timer: null, canvas: null, els: null };

// Optional fill-in (shared checkbox "fill in symmetry equivalents", laue-viz.js):
// the layer data are then built from measured + generated reflections. The
// gaps don't change (they are computed from symmetry classes either way) and
// the intensity scale stays the measured one.
function currentLayerData() {
  if (!laueFillSymmetryEnabled()) return _layer.dataMeasured;
  if (!_layer.dataFilled) {
    const generated = symmetryFilledReflections(_layer.ctx.reflections, _layer.ctx.group);
    _layer.dataFilled = prepareLayerData({ ..._layer.ctx, reflections: _layer.ctx.reflections.concat(generated) });
  }
  return _layer.dataFilled;
}

// Everything that does not depend on the current layer, computed once.
function prepareLayerData(ctx) {
  const recip = reciprocalCellFromDirect(ctx.cell);
  const classes = new Set();
  let sMax = 0;
  const pos = [];
  const axes = {};
  for (const key of Object.keys(LAYER_AXES)) axes[key] = { buckets: new Map(), nMin: Infinity, nMax: -Infinity, uMin: Infinity, uMax: -Infinity, vMin: Infinity, vMax: -Infinity };
  for (const r of ctx.reflections) {
    if (!r.generated) { // filled-in positions change neither the gaps, the resolution limit nor the scale
      const c = canonicalHKL(r.h, r.k, r.l, ctx.group);
      classes.add(c[0] + ',' + c[1] + ',' + c[2]);
      const s = sinThetaOverLambda(r.h, r.k, r.l, recip);
      if (s > sMax) sMax = s;
      if (r.f2 > 0) pos.push(r.f2);
    }
    for (const [key, cfg] of Object.entries(LAYER_AXES)) {
      const a = axes[key], n = cfg.layerOf(r), u = cfg.uOf(r), v = cfg.vOf(r);
      let b = a.buckets.get(n);
      if (!b) { b = []; a.buckets.set(n, b); }
      b.push(r);
      if (n < a.nMin) a.nMin = n; if (n > a.nMax) a.nMax = n;
      if (u < a.uMin) a.uMin = u; if (u > a.uMax) a.uMax = u;
      if (v < a.vMin) a.vMin = v; if (v > a.vMax) a.vMax = v;
    }
  }
  pos.sort((x, y) => x - y);
  const Iref = pos.length ? pos[Math.min(pos.length - 1, Math.floor(pos.length * 0.99))] : 1;
  return { recip, classes, sMax, axes, Iref };
}

// Points and gaps of one layer. Pure (no drawing), so it can be tested.
function computeLayerInfo(ctx, data, axisKey, n) {
  const cfg = LAYER_AXES[axisKey], a = data.axes[axisKey];
  const bucket = a.buckets.get(n) || [];
  const map = new Map();
  for (const r of bucket) {
    const u = cfg.uOf(r), v = cfg.vOf(r), id = u + ',' + v;
    let e = map.get(id);
    if (!e) { e = { u, v, sum: 0, n: 0, generated: !!r.generated }; map.set(id, e); }
    e.sum += r.f2; e.n++;
  }
  const points = [...map.values()].map(e => ({ u: e.u, v: e.v, I: e.sum / e.n, generated: e.generated }));
  const gaps = [];
  let possible = 0;
  const lim = data.sMax + 1e-9;
  for (let u = a.uMin - 1; u <= a.uMax + 1; u++) {
    for (let v = a.vMin - 1; v <= a.vMax + 1; v++) {
      const [h, k, l] = cfg.hkl(u, v, n);
      if (h === 0 && k === 0 && l === 0) continue;
      if (sinThetaOverLambda(h, k, l, data.recip) > lim) continue;
      if (ctx.latticeTest && !ctx.latticeTest({ h, k, l })) continue;
      if (ctx.scan && isForbiddenByHeldConditions(ctx.scan, h, k, l)) continue;
      possible++;
      const c = canonicalHKL(h, k, l, ctx.group);
      if (!data.classes.has(c[0] + ',' + c[1] + ',' + c[2])) gaps.push({ u, v });
    }
  }
  return { points, gaps, possible };
}

function drawLayer(canvas, zoom = 1) {
  const { ctx, data, axis, n, showGaps } = _layer;
  const colors = getPlotColors();
  const cfg = LAYER_AXES[axis], a = data.axes[axis];
  const cssW = 420, cssH = 420;
  const g = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);
  if (zoom === 1) { canvas.style.width = '100%'; canvas.style.maxWidth = cssW + 'px'; canvas.style.height = 'auto'; }
  const info = computeLayerInfo(ctx, data, axis, n);
  const r = { uMin: a.uMin - 1, uMax: a.uMax + 1, vMin: a.vMin - 1, vMax: a.vMax + 1 };
  const basis = reciprocalPlaneBasis(ctx.cell, cfg.uIdx, cfg.vIdx);
  const frame = makeLatticeFrame(basis, r.uMin, r.uMax, r.vMin, r.vMax, { x: 34, y: 12, w: cssW - 46, h: cssH - 40 });
  drawLatticeGrid(g, colors, frame, r, { x: cfg.u, y: cfg.v },
    (lo, hi) => integerTicks(lo, hi, niceIntStep(hi - lo)), cssW, cssH);
  const modeEl = document.getElementById('laueScaleSelect');
  const mode = modeEl ? modeEl.value : 'area';
  const rmax = Math.max(2, frame.spacing * 0.5);

  if (showGaps) {
    g.strokeStyle = '#d33'; g.lineWidth = 1.4;
    for (const p of info.gaps) {
      g.beginPath(); g.arc(frame.px(p.u, p.v), frame.py(p.u, p.v), Math.max(2.5, rmax * 0.6), 0, 2 * Math.PI); g.stroke();
    }
  }
  for (const p of info.points) {
    g.fillStyle = p.generated ? colors.symFilled : colors.point;
    const x = frame.px(p.u, p.v), y = frame.py(p.u, p.v);
    const rad = laueRadius(p.I, data.Iref, rmax, mode);
    g.beginPath(); g.arc(x, y, rad > 0 ? rad : 1.2, 0, 2 * Math.PI); g.fill();
  }
  g.fillStyle = colors.text; g.font = '13px sans-serif'; g.textAlign = 'right';
  g.fillText(`${cfg.fixed} = ${n}`, cssW - 8, 20);
  return info;
}

function refreshLayer() {
  if (!_layer.canvas) return;
  const info = drawLayer(_layer.canvas);
  const { els, axis, n } = _layer;
  const cfg = LAYER_AXES[axis];
  els.label.textContent = `${cfg.fixed} = ${n}`;
  const total = info.points.filter(p => !p.generated).length;
  const filled = info.points.length - total;
  const complete = info.possible ? Math.min(100, (100 * (info.possible - info.gaps.length)) / info.possible) : null;
  els.status.innerHTML = total === 0 && info.gaps.length === 0
    ? '<small>No reflections in this layer.</small>'
    : `<small>${total} positions measured${filled ? ` · <span class="symFilledKey">${filled} filled in by symmetry</span>` : ''} · ${info.gaps.length} gap${info.gaps.length === 1 ? '' : 's'} (red rings)` +
      `${complete === null ? '' : ` · ${complete.toFixed(1)} % of the ${info.possible} possible unique-reflection positions in this layer covered (symmetry equivalents count)`}</small>`;
}

function stopLayerPlay() {
  if (_layer.timer) { clearInterval(_layer.timer); _layer.timer = null; }
  if (_layer.els) _layer.els.play.textContent = '▶ Play';
}

function setLayerAxis(axis) {
  _layer.axis = axis;
  const a = _layer.data.axes[axis];
  _layer.n = Math.max(a.nMin, Math.min(a.nMax, 0));
  const sl = _layer.els.slider;
  sl.min = a.nMin; sl.max = a.nMax; sl.value = _layer.n;
}

function initLayerViewer(container, ctx) {
  stopLayerPlay();
  _layer.ctx = ctx;
  _layer.dataMeasured = prepareLayerData(ctx);
  _layer.dataFilled = null;
  _layer.data = currentLayerData();
  container.innerHTML = `
    <p class="layerControls">
      <label>Layers:
        <select id="layerAxisSel">
          <option value="l">hk<i>L</i> (l fixed)</option>
          <option value="k">h<i>K</i>l (k fixed)</option>
          <option value="h"><i>H</i>kl (h fixed)</option>
        </select></label>
      <button type="button" id="layerPrev" title="previous layer">◀</button>
      <input type="range" id="layerSlider" step="1" style="vertical-align:middle;width:min(260px,40vw);">
      <button type="button" id="layerNext" title="next layer">▶</button>
      <button type="button" id="layerPlay">▶ Play</button>
      <strong id="layerLabel"></strong>
      <label><input type="checkbox" id="layerGaps" checked> mark gaps</label>
    </p>
    <canvas id="layerCanvas"></canvas>
    <div id="layerStatus"></div>`;
  const q = id => container.querySelector('#' + id);
  _layer.els = { slider: q('layerSlider'), label: q('layerLabel'), status: q('layerStatus'), play: q('layerPlay') };
  _layer.canvas = q('layerCanvas');
  _layer.showGaps = true;
  setLayerAxis(_layer.axis in LAYER_AXES ? _layer.axis : 'l');
  q('layerAxisSel').value = _layer.axis;

  const step = d => {
    const a = _layer.data.axes[_layer.axis];
    _layer.n += d;
    if (_layer.n > a.nMax) _layer.n = a.nMin;
    if (_layer.n < a.nMin) _layer.n = a.nMax;
    _layer.els.slider.value = _layer.n;
    refreshLayer();
  };
  q('layerAxisSel').addEventListener('change', e => { stopLayerPlay(); setLayerAxis(e.target.value); refreshLayer(); });
  _layer.els.slider.addEventListener('input', e => { _layer.n = parseInt(e.target.value, 10); refreshLayer(); });
  q('layerPrev').addEventListener('click', () => { stopLayerPlay(); step(-1); });
  q('layerNext').addEventListener('click', () => { stopLayerPlay(); step(1); });
  _layer.els.play.addEventListener('click', () => {
    if (_layer.timer) { stopLayerPlay(); return; }
    _layer.els.play.textContent = '⏸ Pause';
    _layer.timer = setInterval(() => step(1), 450);
  });
  q('layerGaps').addEventListener('change', e => { _layer.showGaps = e.target.checked; refreshLayer(); });
  attachPlotZoom(_layer.canvas, `<strong>Layer viewer</strong>`, 420, (cv, zoom) => drawLayer(cv, zoom));
  refreshLayer();
}

(function () {
  const sel = document.getElementById('laueScaleSelect');
  if (sel) sel.addEventListener('change', refreshLayer);
  const fill = document.getElementById('laueFillSym');
  if (fill) fill.addEventListener('change', () => {
    if (!_layer.ctx || !_layer.els) return;
    _layer.data = currentLayerData();
    // the layer range can grow (e.g. negative h layers that exist only by symmetry); keep the current layer
    const keep = _layer.n;
    setLayerAxis(_layer.axis);
    const a = _layer.data.axes[_layer.axis];
    _layer.n = Math.max(a.nMin, Math.min(a.nMax, keep));
    _layer.els.slider.value = _layer.n;
    refreshLayer();
  });
})();
onPlotThemeChange(refreshLayer);
