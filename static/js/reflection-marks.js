// reflection-marks.js
// Per-position statistics and "this doesn't fit" marks for the layer viewer
// (2D, layer-viewer.js) and its 3D view (laue-3d.js). Both views read the
// same object, so they always mark the same reflections.
//
//  • OUTLIER: one observation disagrees with the OTHER observations of its
//    symmetry class under the Laue group (exact repeats count as members):
//        z = (I − ⟨I⟩others) / sqrt(σ² + σ⟨I⟩²)
//    with σ⟨I⟩ the larger of the propagated σ and the scatter of the others
//    (σ values of real data are often too small). Needs at least two other
//    observations, i.e. unmerged data. A position is an outlier when its
//    worst observation exceeds the threshold. Typical causes: a spot of a
//    second domain/twin component on top of it (too strong), beamstop or
//    detector edge (too weak), a bad frame (one batch).
//  • FORBIDDEN BUT OBSERVED: the position is systematically absent under the
//    observed lattice and the absence conditions that held (ctx.latticeTest,
//    the same test the gaps use), but its repeat-averaged intensity is
//    significant (not isWeak, absences.js: I/σ ≥ 3).
//
// Generated (symmetry-filled) reflections are never marked.
// Needs: canonicalHKL (laue-merge.js), isWeak (absences.js).

const REFL_MARK_THRESHOLDS = [4, 6, 10];
const REFL_MARK_DEFAULT_THRESHOLD = 6;
const _reflMarkCache = { key: null, value: null };

function reflKey(h, k, l) { return h + ',' + k + ',' + l; }

// Everything that does not depend on the threshold: positions (exact
// repeats merged) and the sums per symmetry class.
function _reflMarkBase(ctx) {
  const positions = new Map(); // "h,k,l" -> { h, k, l, cls, obs: [reflection] }
  const classes = new Map();   // canonical key -> { n, S, S2, Ss2, members: Set of position keys }
  for (const r of ctx.reflections) {
    if (r.generated) continue;
    const pk = reflKey(r.h, r.k, r.l);
    let p = positions.get(pk);
    if (!p) {
      const c = canonicalHKL(r.h, r.k, r.l, ctx.group);
      p = { h: r.h, k: r.k, l: r.l, cls: reflKey(c[0], c[1], c[2]), obs: [] };
      positions.set(pk, p);
    }
    p.obs.push(r);
    let c = classes.get(p.cls);
    if (!c) { c = { n: 0, S: 0, S2: 0, Ss2: 0, members: new Set() }; classes.set(p.cls, c); }
    c.n++; c.S += r.f2; c.S2 += r.f2 * r.f2; c.Ss2 += r.sigma * r.sigma; c.members.add(pk);
  }
  let unmerged = false;
  for (const c of classes.values()) if (c.n >= 3) { unmerged = true; break; }
  return { positions, classes, unmerged };
}

// Leave-one-out comparison of one observation with the rest of its class.
function reflOutlierZ(r, c) {
  const n = c.n - 1;
  if (n < 2) return null;
  const mean = (c.S - r.f2) / n;
  const scatter = Math.max(0, (c.S2 - r.f2 * r.f2 - n * mean * mean) / (n - 1));
  const s2mean = Math.max((c.Ss2 - r.sigma * r.sigma) / (n * n), scatter / n);
  const denom = Math.sqrt(r.sigma * r.sigma + s2mean);
  if (!(denom > 0)) return null;
  return { z: (r.f2 - mean) / denom, mean, nOthers: n };
}

// Returns { positions, classes, unmerged, outliers: Map, forbidden: Map, threshold }.
// outliers: position key -> { z, obs (the worst observation), mean, nOthers }
// forbidden: position key -> { I, sigma }
function getReflectionMarks(ctx, threshold) {
  if (_reflMarkCache.key && _reflMarkCache.key.ctx === ctx && _reflMarkCache.key.threshold === threshold) return _reflMarkCache.value;
  const base = (_reflMarkCache.value && _reflMarkCache.key.ctx === ctx) ? _reflMarkCache.value : _reflMarkBase(ctx);
  const outliers = new Map(), forbidden = new Map();
  for (const [pk, p] of base.positions) {
    const c = base.classes.get(p.cls);
    if (base.unmerged) {
      let worst = null;
      for (const r of p.obs) {
        const o = reflOutlierZ(r, c);
        if (o && Math.abs(o.z) > threshold && (!worst || Math.abs(o.z) > Math.abs(worst.z))) worst = { ...o, obs: r };
      }
      if (worst) outliers.set(pk, worst);
    }
    if (ctx.latticeTest && !ctx.latticeTest(p)) {
      const st = reflPositionStats(p);
      if (!isWeak(st.I, st.sigma)) forbidden.set(pk, st);
    }
  }
  const value = { positions: base.positions, classes: base.classes, unmerged: base.unmerged, outliers, forbidden, threshold };
  _reflMarkCache.key = { ctx, threshold };
  _reflMarkCache.value = value;
  return value;
}

// Mean of the exact repeats at one position and its propagated σ.
function reflPositionStats(p) {
  let S = 0, Ss2 = 0;
  for (const r of p.obs) { S += r.f2; Ss2 += r.sigma * r.sigma; }
  const n = p.obs.length;
  return { I: S / n, sigma: Math.sqrt(Ss2) / n, n };
}

// ---------- hover / click info ----------

const _reflTip = { el: null, pinned: false, keyHandler: null };

const _fmtI = x => (Math.abs(x) >= 1000 ? x.toFixed(0) : Math.abs(x) >= 10 ? x.toFixed(1) : x.toFixed(2));

// item: { h, k, l, kind: 'measured' | 'generated' | 'gap' }; marks: getReflectionMarks result
// (or null when no mark is switched on: the statistics are then built on demand).
function reflectionInfoHtml(ctx, item, marks) {
  const m = marks || getReflectionMarks(ctx, REFL_MARK_DEFAULT_THRESHOLD);
  const pk = reflKey(item.h, item.k, item.l);
  const head = `<strong>${item.h} ${item.k} ${item.l}</strong>`;
  if (item.kind === 'gap') return `${head}<br>gap: not measured, and no symmetry equivalent either`;
  const p = m.positions.get(pk);
  if (!p) {
    if (item.kind === 'generated') return `${head}<br><span class="symFilledKey">filled in by symmetry</span> (not measured itself)`;
    return head;
  }
  const st = reflPositionStats(p);
  const c = m.classes.get(p.cls);
  const lines = [head,
    `I = ${_fmtI(st.I)}, σ = ${_fmtI(st.sigma)}, I/σ = ${st.sigma > 0 ? (st.I / st.sigma).toFixed(1) : '—'}` + (st.n > 1 ? ` <small>(mean of ${st.n} repeats)</small>` : ''),
  ];
  if (c && c.n > st.n) lines.push(`symmetry class: ${c.n} observations at ${c.members.size} positions, ⟨I⟩ = ${_fmtI(c.S / c.n)}`);
  const batches = [...new Set(p.obs.map(r => r.batch).filter(b => b !== null && b !== undefined))];
  if (batches.length) lines.push(`batch${batches.length > 1 ? 'es' : ''}: ${batches.slice(0, 8).join(', ')}${batches.length > 8 ? ' …' : ''}`);
  const o = marks && marks.outliers.get(pk);
  if (o) {
    lines.push(`<span class="reflMarkOutlier">outlier</span>: one observation (I = ${_fmtI(o.obs.f2)}${o.obs.batch !== null && o.obs.batch !== undefined ? `, batch ${o.obs.batch}` : ''})` +
      ` is ${o.z > 0 ? 'stronger' : 'weaker'} than the ${o.nOthers} other observations of its class (⟨I⟩ = ${_fmtI(o.mean)}), z = ${o.z.toFixed(1)}`);
  }
  if (marks && marks.forbidden.has(pk)) lines.push('<span class="reflMarkForbidden">forbidden but observed</span>: systematically absent under the observed lattice / absence conditions, yet I/σ ≥ 3');
  return lines.join('<br>');
}

function showReflectionTip(html, clientX, clientY, pin) {
  if (_reflTip.pinned && !pin) return;
  if (!_reflTip.el) {
    _reflTip.el = document.createElement('div');
    _reflTip.el.className = 'reflTip';
    _reflTip.el.setAttribute('role', 'tooltip');
    document.body.appendChild(_reflTip.el);
  }
  const el = _reflTip.el;
  el.innerHTML = html + (pin ? '<div class="reflTipHint">click elsewhere or Esc to close</div>' : '');
  el.style.display = 'block';
  el.classList.toggle('pinned', !!pin);
  // keep inside the viewport
  const pad = 12, w = el.offsetWidth, h = el.offsetHeight;
  let x = clientX + 14, y = clientY + 14;
  if (x + w > window.innerWidth - pad) x = Math.max(pad, clientX - w - 14);
  if (y + h > window.innerHeight - pad) y = Math.max(pad, clientY - h - 14);
  el.style.left = x + 'px'; el.style.top = y + 'px';
  if (pin) {
    _reflTip.pinned = true;
    if (!_reflTip.keyHandler) {
      _reflTip.keyHandler = e => { if (e.key === 'Escape') { e.stopPropagation(); hideReflectionTip(true); } }; // closes only the tip, not an enlarged view behind it
      document.addEventListener('keydown', _reflTip.keyHandler, true);
      setTimeout(() => document.addEventListener('pointerdown', _reflTipOutside, true), 0);
    }
  }
}

function _reflTipOutside(e) {
  if (_reflTip.el && _reflTip.el.contains(e.target)) return;
  hideReflectionTip(true);
}

function hideReflectionTip(force) {
  if (!_reflTip.el || (_reflTip.pinned && !force)) return;
  _reflTip.el.style.display = 'none';
  if (force) {
    _reflTip.pinned = false;
    if (_reflTip.keyHandler) document.removeEventListener('keydown', _reflTip.keyHandler, true);
    _reflTip.keyHandler = null;
    document.removeEventListener('pointerdown', _reflTipOutside, true);
  }
}
