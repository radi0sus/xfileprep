// bijvoet-zones.js
// Bijvoet differences in one reciprocal-lattice layer, drawn like the Laue
// layer viewer ("Twins & Abs. Structure" tab, inversion-twin section):
// circle area ∝ |I(h) − I(−h)|/σ (full size at BZ_ZMAX), filled where
// I(h) > I(−h), open where I(h) < I(−h), small grey dots below the chosen
// threshold (1, 2 or 3σ; 3σ by default, so noise pairs of a centrosymmetric
// structure stay dots). Every pair is drawn at all symmetry equivalents of h
// and (opposite sign) of −h. As in the Laue viewer, indices present in the data
// get the normal colour, positions filled in from an equivalent the "generated"
// (orange) one. The status line counts pairs, not positions.
// Centric layers have no pairs, so the plot starts at the layer nearest to 0
// that has some (otherwise e.g. hk0 in P212121 would open empty).
// Entry: initBijvoetZones(container, cell, check), check = checkAnomalousSignal(...).

const BZ_ZMAX = 8;
const BZ_THRESHOLDS = [1, 2, 3];
const BZ_AXES = {
  l: { fixed: 'l', u: 'h', v: 'k', uIdx: 0, vIdx: 1, fIdx: 2 },
  k: { fixed: 'k', u: 'h', v: 'l', uIdx: 0, vIdx: 2, fIdx: 1 },
  h: { fixed: 'h', u: 'k', v: 'l', uIdx: 1, vIdx: 2, fIdx: 0 },
};
const _bz = { sig: null, cell: null, pg: null, pts: [], axis: 'l', n: 0, thr: 3, els: null, hit: [], themeHooked: false };

// all equivalents h·R of h (+z) and of −h (−z)
// id = pair number; generated = index not in the data (measured = set of "h,k,l")
function bzExpand(pairs, pg, measured) {
  const out = new Map();
  pairs.forEach((p, id) => {
    for (const s of [1, -1]) {
      const h = p.h.map(x => s * x);
      for (const R of pg) {
        const e = [0, 1, 2].map(j => h[0] * R[0][j] + h[1] * R[1][j] + h[2] * R[2][j]);
        const key = e.join(',');
        if (!out.has(key)) out.set(key, { h: e, z: s * p.z, id, generated: !measured.has(key) });
      }
    }
  });
  return [...out.values()];
}

function bzLayers(axis) {
  const f = BZ_AXES[axis].fIdx;
  return [...new Set(_bz.pts.map(p => p.h[f]))].sort((a, b) => a - b);
}

function bzNearestLayer(axis, n) {
  const layers = bzLayers(axis);
  if (!layers.length || layers.includes(n)) return n;
  return layers.reduce((best, x) => (Math.abs(x - n) < Math.abs(best - n) || (Math.abs(x - n) === Math.abs(best - n) && x > best) ? x : best));
}

function bzIsCentric(axis, n) { // is −h equivalent to h for a general position of the layer?
  const c = BZ_AXES[axis], h = [0, 0, 0];
  h[c.uIdx] = 3; h[c.vIdx] = 5; h[c.fIdx] = n;
  const neg = h.map(x => -x).join(',');
  return _bz.pg.some(R => [0, 1, 2].map(j => h[0] * R[0][j] + h[1] * R[1][j] + h[2] * R[2][j]).join(',') === neg);
}

function drawBijvoetZone(canvas) {
  const colors = getPlotColors();
  const cfg = BZ_AXES[_bz.axis];
  const cssW = 420, cssH = 440;
  const g = setupCanvasForDPR(canvas, cssW, cssH, colors, 1);
  canvas.style.width = '100%'; canvas.style.maxWidth = cssW + 'px'; canvas.style.height = 'auto';
  // frame from all pairs (not just this layer) so it stays put while stepping
  const range = i => [Math.min(..._bz.pts.map(p => p.h[i])) - 1, Math.max(..._bz.pts.map(p => p.h[i])) + 1];
  const [uMin, uMax] = range(cfg.uIdx), [vMin, vMax] = range(cfg.vIdx);
  const r = { uMin, uMax, vMin, vMax };
  const frame = makeLatticeFrame(reciprocalPlaneBasis(_bz.cell, cfg.uIdx, cfg.vIdx), uMin, uMax, vMin, vMax, { x: 34, y: 12, w: cssW - 46, h: cssH - 60 });
  drawLatticeGrid(g, colors, frame, r, { x: cfg.u, y: cfg.v }, (lo, hi) => integerTicks(lo, hi, niceIntStep(hi - lo)), cssW, cssH - 20);
  const rmax = Math.max(2, frame.spacing * 0.5);
  const hit = [];
  const layer = _bz.pts.filter(p => p.h[cfg.fIdx] === _bz.n).sort((a, b) => b.generated - a.generated); // measured on top
  const pairs = new Set(layer.map(p => p.id)), sig = new Set();
  for (const p of layer) {
    const x = frame.px(p.h[cfg.uIdx], p.h[cfg.vIdx]), y = frame.py(p.h[cfg.uIdx], p.h[cfg.vIdx]);
    const az = Math.abs(p.z);
    let rad;
    if (az < _bz.thr) {
      rad = 1.2; g.fillStyle = colors.axis;
      g.beginPath(); g.arc(x, y, rad, 0, 2 * Math.PI); g.fill();
    } else {
      rad = Math.max(1.5, rmax * Math.sqrt(Math.min(az, BZ_ZMAX) / BZ_ZMAX));
      const col = p.generated ? colors.symFilled : colors.point;
      g.beginPath(); g.arc(x, y, rad, 0, 2 * Math.PI);
      if (p.z > 0) { g.fillStyle = col; g.fill(); }
      else { g.strokeStyle = col; g.lineWidth = 1.3; g.stroke(); }
      sig.add(p.id);
    }
    hit.push({ x, y, r: rad, p });
  }
  _bz.hit = hit;
  g.fillStyle = colors.text; g.font = '13px sans-serif'; g.textAlign = 'right';
  g.fillText(`${cfg.fixed} = ${_bz.n}`, cssW - 8, 20);
  g.font = '11px sans-serif'; g.textAlign = 'left';
  g.fillText(`● I(h) > I(−h)   ○ I(h) < I(−h)   area ∝ |ΔI|/σ (max ${BZ_ZMAX}σ)`, 34, cssH - 6);
  return { n: layer.length, pairs: pairs.size, sig: sig.size };
}

function bzRefresh() {
  const els = _bz.els;
  const info = drawBijvoetZone(els.canvas);
  els.label.textContent = `${BZ_AXES[_bz.axis].fixed} = ${_bz.n}`;
  els.status.innerHTML = '<small>' + (info.n
    ? `${info.pairs} pair${info.pairs === 1 ? '' : 's'} in this layer, ${info.sig} with |ΔI| ≥ ${_bz.thr}σ · drawn at ${info.n} positions, orange = not in the data, filled in from an equivalent`
    : bzIsCentric(_bz.axis, _bz.n) ? 'Centric layer: −h is equivalent to h, so there are no Bijvoet pairs here.' : 'No Bijvoet pairs in this layer.') + '</small>';
}

function bzSetAxis(axis) {
  _bz.axis = axis;
  const layers = bzLayers(axis);
  _bz.n = bzNearestLayer(axis, _bz.n);
  _bz.els.slider.min = layers[0]; _bz.els.slider.max = layers[layers.length - 1]; _bz.els.slider.value = _bz.n;
}

function initBijvoetZones(container, cell, check) {
  if (!check || !check.enough || !check.bijvoetPairs || !check.bijvoetPairs.length || !cell) { container.innerHTML = ''; _bz.els = null; return; }
  const pairs = check.bijvoetPairs;
  // same data as before (the tab re-renders often)? then keep axis and layer
  const sig = `${pairs.length}|${pairs[0].h}|${pairs[0].z}|${cell.a},${cell.b},${cell.c}`;
  if (sig !== _bz.sig) { _bz.sig = sig; _bz.axis = 'l'; _bz.n = 0; }
  _bz.cell = cell; _bz.pg = check.pointGroup; _bz.pts = bzExpand(pairs, check.pointGroup, check.measured || new Set());
  container.innerHTML = `
    <div class="layerControls">
      <span class="ctlTitle">Bijvoet zone</span>
      <div class="ctlRow">
        <select class="bzAxis" title="which index is fixed per layer">
          <option value="l">hkL (l fixed)</option><option value="k">hKl (k fixed)</option><option value="h">Hkl (h fixed)</option>
        </select>
        <span class="ctlGroup">
          <button type="button" class="bzPrev" title="previous layer">◀</button>
          <input type="range" class="bzSlider" step="1">
          <button type="button" class="bzNext" title="next layer">▶</button>
        </span>
        <strong class="layerLabel bzLabel"></strong>
        <label title="pairs with a smaller |I(h) − I(−h)|/σ are drawn as grey dots">|ΔI|/σ ≥
          <select class="bzThr">${BZ_THRESHOLDS.map(t => `<option value="${t}"${t === _bz.thr ? ' selected' : ''}>${t}</option>`).join('')}</select></label>
      </div>
    </div>
    <div class="layer2dWrap"><canvas class="bzCanvas"></canvas></div>
    <div class="bzStatus"></div>`;
  const q = s => container.querySelector(s);
  _bz.els = { canvas: q('.bzCanvas'), slider: q('.bzSlider'), label: q('.bzLabel'), status: q('.bzStatus') };
  const els = _bz.els;
  q('.bzAxis').value = _bz.axis;
  bzSetAxis(_bz.axis);
  q('.bzAxis').addEventListener('change', e => { bzSetAxis(e.target.value); bzRefresh(); });
  q('.bzThr').addEventListener('change', e => { _bz.thr = Number(e.target.value); bzRefresh(); });
  els.slider.addEventListener('input', () => { _bz.n = Number(els.slider.value); bzRefresh(); });
  const step = d => { _bz.n = Math.max(Number(els.slider.min), Math.min(Number(els.slider.max), _bz.n + d)); els.slider.value = _bz.n; bzRefresh(); };
  q('.bzPrev').addEventListener('click', () => step(-1));
  q('.bzNext').addEventListener('click', () => step(1));
  const bzHit = e => {
    const b = els.canvas.getBoundingClientRect(), sc = 420 / b.width;
    const x = (e.clientX - b.left) * sc, y = (e.clientY - b.top) * sc;
    const slop = typeof touchInput === 'function' && touchInput() ? touchSlop(5) * sc : 5;
    let best = null, bestD = Infinity;
    for (const t of _bz.hit) { const d = Math.hypot(t.x - x, t.y - y); if (d <= Math.max(t.r + 2, slop) && d < bestD) { best = t; bestD = d; } }
    return best;
  };
  const bzTipHtml = it => `<strong>${it.p.h.join(' ')}</strong>${it.p.generated ? ' <small>(not in the data, from an equivalent)</small>' : ''}<br>(I(h) − I(−h))/σ = <strong>${it.p.z.toFixed(2)}</strong>`;
  els.canvas.addEventListener('mousemove', e => {
    if (typeof touchInput === 'function' && touchInput()) return;
    const it = bzHit(e);
    if (it) showReflectionTip(bzTipHtml(it), e.clientX, e.clientY, false);
    else hideReflectionTip(false);
  });
  const wasSwipe = typeof attachSwipeStep === 'function' ? attachSwipeStep(els.canvas, d => { step(d); hideInfoBar('refl'); }) : null;
  els.canvas.addEventListener('click', e => {
    if (wasSwipe && wasSwipe()) return;
    const it = bzHit(e);
    if (it) showReflectionTip(bzTipHtml(it), e.clientX, e.clientY, true);
    else if (typeof touchInput === 'function' && touchInput()) hideInfoBar('refl');
  });
  els.canvas.addEventListener('mouseleave', () => hideReflectionTip(false));
  if (!_bz.themeHooked) { _bz.themeHooked = true; onPlotThemeChange(() => { if (_bz.els && _bz.els.canvas.isConnected) bzRefresh(); }); }
  bzRefresh();
}
