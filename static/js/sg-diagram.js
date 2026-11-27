// sg-diagram.js
// Combined space-group diagram for the selected candidate: the unit cell in
// projection with its symmetry elements AND the general positions in one
// drawing, plus the list of coordinates — in the graphical language of
// International Tables Vol. A (lens = 2-fold axis, hooks = 2₁, half arrow =
// in-plane 2₁, dashed/dotted/dash-dot = glide planes, thick line = mirror,
// "+ −  ½+" and comma for the general positions).
//
// Everything is derived at runtime from the Seitz operators of the Hall
// symbol (hall-symbol-parser.js) and drawn as inline SVG in the app's own
// style (CSS variables, light/dark). Scope: triclinic, monoclinic and
// orthorhombic groups in any axis setting — i.e. all groups whose rotation
// parts are diagonal (±1). Other groups get a short note instead.
//
// Geometry: schematic cell (both in-plane edges the same length) with the
// REAL angle between them (β for monoclinic b-unique etc.); projection along
// the unique axis (monoclinic) or along c (triclinic, orthorhombic).
// In-plane axes: horizontal = next axis after the viewing axis, vertical =
// the one after that (b-unique: c →, a ↑; along c: a →, b ↑).

const _SGD_FRACS = { '1/2': '½', '1/4': '¼', '3/4': '¾', '1/3': '⅓', '2/3': '⅔', '1/6': '⅙', '5/6': '⅚', '1/8': '⅛', '3/8': '⅜', '5/8': '⅝', '7/8': '⅞' };
const _sgdMod1 = v => { const r = ((v % 1) + 1) % 1; return r > 1 - 1e-9 ? 0 : r; };
function _sgdFrac(v) {
  const f = _sgdMod1(v);
  if (f < 1e-9) return '0';
  for (const den of [2, 3, 4, 6, 8, 12]) {
    const num = Math.round(f * den);
    if (Math.abs(num / den - f) < 1e-6) {
      const g = (a, b) => (b ? g(b, a % b) : a);
      const d = g(num, den);
      const k = `${num / d}/${den / d}`;
      return _SGD_FRACS[k] || k;
    }
  }
  return f.toFixed(3);
}

const _sgdDet = R => R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1])
  - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0])
  + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
const _sgdIsDiagonal = R => R.every((row, i) => row.every((v, j) => (i === j ? Math.abs(v) === 1 : v === 0)));

// Quarter turn about c: R = (A ⊕ ±1) with A a 90° rotation of the ab plane.
// R₃₃ = +1 is a proper 4-fold rotation, R₃₃ = −1 a 4̄ rotoinversion.
// plus: 4⁺ is x,y,z → −y,x,z (IT sense), 4̄⁺ = −1·4⁺ is x,y,z → y,−x,−z.
function _sgdQuarterTurn(op) {
  const R = op.R;
  if (R[0][2] || R[1][2] || R[2][0] || R[2][1] || R[0][0] || R[1][1]) return null;
  if (Math.abs(R[0][1]) !== 1 || R[0][1] !== -R[1][0] || Math.abs(R[2][2]) !== 1) return null;
  const proper = R[2][2] === 1;
  return { proper, plus: proper ? R[1][0] === 1 : R[1][0] === -1 };
}
const _sgdIsTetragonal = ops => ops.some(op => _sgdQuarterTurn(op));

// Diagonal operation of a tetragonal group: R swaps x and y,
// R = ([[0, σ], [σ, 0]] ⊕ ±1). With R₃₃ = +1 it is a mirror/glide plane
// perpendicular to the drawing (σ = +1: plane x − y = const, normal [1̄10];
// σ = −1: plane x + y = const, normal [110]); with R₃₃ = −1 a 2-fold axis in
// the drawing plane (σ = +1: along [110]; σ = −1: along [1̄10]).
function _sgdSwap(op) {
  const R = op.R;
  if (R[0][2] || R[1][2] || R[2][0] || R[2][1] || R[0][0] || R[1][1]) return null;
  if (Math.abs(R[0][1]) !== 1 || R[0][1] !== R[1][0] || Math.abs(R[2][2]) !== 1) return null;
  return { sigma: R[0][1], mirror: R[2][2] === 1 };
}

// Drawn: all groups whose rotation parts are diagonal (triclinic,
// monoclinic, orthorhombic in any axis setting) and all tetragonal groups
// (Nos. 75–142), viewed along c.
function sgdIsSupported(ops) {
  if (_sgdIsTetragonal(ops)) return ops.every(op => _sgdIsDiagonal(op.R) || _sgdQuarterTurn(op) || _sgdSwap(op));
  return ops.every(op => _sgdIsDiagonal(op.R));
}

// Viewing axis: the unique axis for monoclinic, c otherwise.
function sgdViewAxis(ops) {
  if (_sgdIsTetragonal(ops)) return 2;
  const special = new Set();
  for (const op of ops) {
    const d = [0, 1, 2].map(i => op.R[i][i]);
    const nNeg = d.filter(v => v < 0).length;
    if (nNeg === 0 || nNeg === 3) continue; // identity / inversion
    // 2-fold: one +1 (its axis); mirror: one −1 (its normal)
    const odd = nNeg === 2 ? d.indexOf(1) : d.indexOf(-1);
    special.add(odd);
  }
  return special.size === 1 ? [...special][0] : 2;
}

function sgdCrystalSystem(ops) {
  if (_sgdIsTetragonal(ops)) return 'tetragonal';
  const special = new Set();
  for (const op of ops) {
    const d = [0, 1, 2].map(i => op.R[i][i]);
    const nNeg = d.filter(v => v < 0).length;
    if (nNeg === 1 || nNeg === 2) special.add(nNeg === 2 ? d.indexOf(1) : d.indexOf(-1));
  }
  return special.size === 0 ? 'triclinic' : special.size === 1 ? 'monoclinic' : 'orthorhombic';
}

// Point-group symbol from the operations (for the header).
function sgdPointGroup(ops) {
  const sys = sgdCrystalSystem(ops);
  const hasInv = ops.some(op => [0, 1, 2].every(i => op.R[i][i] === -1));
  const nRot = new Set(ops.filter(op => [0, 1, 2].filter(i => op.R[i][i] === -1).length === 2).map(op => [0, 1, 2].findIndex(i => op.R[i][i] === 1))).size;
  if (sys === 'tetragonal') {
    // 4, 4̄, 4/m; with elements in the drawing plane or perpendicular to it
    // along a/b or the diagonals: 422, 4mm, 4̄2m / 4̄m2, 4/mmm
    const has4 = ops.some(op => { const q = _sgdQuarterTurn(op); return q && q.proper; });
    const side = ops.filter(op => _sgdSwap(op) || (_sgdIsDiagonal(op.R) && op.R[0][0] !== op.R[1][1]));
    if (!side.length) return formatPG(!has4 ? '-4' : hasInv ? '4/m' : '4');
    if (hasInv) return formatPG('4/mmm');
    if (has4) return formatPG(side.some(op => _sgdDet(op.R) < 0) ? '4mm' : '422');
    // 4̄2m: 2-fold axes along a and b (mirrors on the diagonals); 4̄m2: the reverse
    return formatPG(side.some(op => _sgdIsDiagonal(op.R) && _sgdDet(op.R) > 0) ? '-42m' : '-4m2');
  }
  if (sys === 'triclinic') return formatPG(hasInv ? '-1' : '1');
  if (sys === 'monoclinic') return formatPG(hasInv ? '2/m' : nRot ? '2' : 'm');
  return formatPG(hasInv ? 'mmm' : nRot === 3 ? '222' : 'mm2');
}

// Short H-M symbol from the (full-ish) table entry: drop the "1" place
// holders of monoclinic settings and the origin-choice suffix; subscripts.
function sgdShortSymbol(hm) {
  const base = hm.replace(/:.*$/, '').trim();
  const parts = base.split(/\s+/);
  const lattice = parts.shift();
  const rest = parts.length === 3 && parts.filter(p => p === '1').length === 2 ? parts.filter(p => p !== '1') : parts;
  return formatSG(lattice + rest.join('')); // sg-format.js
}

// --- symmetry elements ------------------------------------------------------
// Rotation parts are diagonal here, so each operation is easy to read:
//   inversion (−1 −1 −1): centre at (t + L)/2
//   2-fold along i (+1 at i): line through x_j = (t_j + L_j)/2, screw part t_i
//   mirror ⊥ i (−1 at i): plane x_i = (t_i + L_i)/2, glide vector t_j (j ≠ i)
// L runs over {0,1}³ (lattice translations change the location by ½).
function sgdElements(ops, view) {
  const H = (view + 1) % 3, V = (view + 2) % 3;
  const t12 = op => op.t.map(v => v / 12);
  const halves = c => [_sgdMod1(c), _sgdMod1(c + 0.5)];
  const inv = new Map(), axesPerp = new Map(), axesH = new Map(), axesV = new Map(), planesH = new Map(), planesV = new Map(), planesPar = new Map();
  const add = (map, key, obj, h) => { let e = map.get(key); if (!e) { e = { ...obj, heights: new Set() }; map.set(key, e); } if (h !== undefined) e.heights.add(_sgdFrac(h)); };
  const axes4 = new Map(), rotoinv4 = new Map(), diagPlanes = new Map(), diagAxes = new Map();
  const pq = (p, q) => `${p.toFixed(4)}|${q.toFixed(4)}`;
  for (const op of ops) {
    // quarter turns about c (view = c): location from the fixed point of the
    // in-plane part, (I − A)·x = t + L for the four lattice shifts L ∈ {0,1}²
    // (each a 4-fold axis of the same kind), screw part t_z for 4, inversion
    // point at z = t_z/2 (and +½) for 4̄
    const qt = _sgdQuarterTurn(op);
    if (qt) {
      const t = t12(op), R = op.R;
      const m00 = 1 - R[0][0], m01 = -R[0][1], m10 = -R[1][0], m11 = 1 - R[1][1];
      const det2 = m00 * m11 - m01 * m10;
      for (const L0 of [0, 1]) for (const L1 of [0, 1]) {
        const b0 = t[0] + L0, b1 = t[1] + L1;
        const p = _sgdMod1((m11 * b0 - m01 * b1) / det2), q = _sgdMod1((-m10 * b0 + m00 * b1) / det2);
        if (qt.proper) {
          const s = _sgdMod1(t[2]);
          // 4⁺ with screw s lies on a 4_{4s} axis, 4⁻ with screw s on a 4_{4(1−s)} axis
          const n = s < 1e-9 ? 0 : Math.round(4 * (qt.plus ? s : 1 - s)) % 4;
          const e = axes4.get(pq(p, q)) || { p, q, kinds: new Set() };
          e.kinds.add(n);
          axes4.set(pq(p, q), e);
        } else {
          for (const h of halves(t[2] / 2)) add(rotoinv4, pq(p, q), { p, q }, h);
        }
      }
      continue;
    }
    // diagonal planes and axes (tetragonal, view = c). With w = (R + I)·t/2
    // the intrinsic (glide/screw) part, the element lies at (I − R)·x = t − w:
    //   σ = +1: line x − y = (t₀ − t₁)/2, in-plane glide/screw (g, g) with g = (t₀ + t₁)/2
    //   σ = −1: line x + y = (t₀ + t₁)/2, in-plane glide/screw (g, −g) with g = (t₀ − t₁)/2
    // A lattice shift L = (1,0) moves the line by ½ and changes g by ½, so
    // mirrors alternate with glide planes and 2 with 2₁ axes. Vertical part:
    // glide t_z for planes, height t_z/2 (and +½) for axes.
    const sw = _sgdSwap(op);
    if (sw) {
      const t = t12(op);
      for (const L0 of [0, 1]) for (const L1 of [0, 1]) {
        const a = t[0] + L0, b = t[1] + L1;
        const c = _sgdMod1(sw.sigma > 0 ? (a - b) / 2 : (a + b) / 2);
        const g = _sgdMod1(sw.sigma > 0 ? (a + b) / 2 : (a - b) / 2);
        if (sw.mirror) {
          const gz = _sgdMod1(t[2]);
          const q = v => Math.abs(v - 0.25) < 1e-9 || Math.abs(v - 0.75) < 1e-9;
          const kind = q(g) || q(gz) ? 'd' : g > 1e-9 && gz > 1e-9 ? 'n' : g > 1e-9 ? 'dash' : gz > 1e-9 ? 'dot' : 'm';
          add(diagPlanes, `${sw.sigma}|${c.toFixed(4)}|${kind}`, { sigma: sw.sigma, c, kind, g, gz });
        } else {
          const screw = Math.abs(g - 0.5) < 1e-9;
          for (const h of halves(t[2] / 2)) add(diagAxes, `${sw.sigma}|${c.toFixed(4)}|${screw}`, { sigma: sw.sigma, c, screw }, h);
        }
      }
      continue;
    }
    const d = [0, 1, 2].map(i => op.R[i][i]);
    const t = t12(op);
    const nNeg = d.filter(v => v < 0).length;
    if (nNeg === 3) {
      for (const p of halves(t[H] / 2)) for (const q of halves(t[V] / 2)) for (const h of halves(t[view] / 2)) add(inv, `${p}|${q}`, { p, q }, h);
    } else if (nNeg === 2) {
      const ax = d.indexOf(1);
      const screw = Math.abs(_sgdMod1(t[ax]) - 0.5) < 1e-9;
      if (ax === view) {
        for (const p of halves(t[H] / 2)) for (const q of halves(t[V] / 2)) add(axesPerp, `${p}|${q}|${screw}`, { p, q, screw });
      } else {
        // in-plane axis: along H (fixed q, height) or along V (fixed p, height)
        const other = ax === H ? V : H;
        for (const c of halves(t[other] / 2)) for (const h of halves(t[view] / 2)) {
          add(ax === H ? axesH : axesV, `${c}|${screw}`, { c, screw }, h);
        }
      }
    } else if (nNeg === 1) {
      const n = d.indexOf(-1);
      const g = [0, 1, 2].map(i => (i === n ? 0 : _sgdMod1(t[i])));
      const gQuarter = g.some(v => Math.abs(v - 0.25) < 1e-9 || Math.abs(v - 0.75) < 1e-9);
      for (const c of halves(t[n] / 2)) {
        if (n === view) {
          const gh = g[H] > 1e-9, gv = g[V] > 1e-9;
          const kind = gQuarter ? 'd' : gh && gv ? 'n' : gh ? 'h' : gv ? 'v' : 'm';
          add(planesPar, `${kind}`, { kind }, c);
        } else {
          const inPlane = n === H ? g[V] : g[H];
          const along = g[view];
          const kind = gQuarter ? 'd' : inPlane > 1e-9 && along > 1e-9 ? 'n' : inPlane > 1e-9 ? 'dash' : along > 1e-9 ? 'dot' : 'm';
          add(n === H ? planesH : planesV, `${c}|${kind}`, { c, kind, gIn: inPlane, along });
        }
      }
    }
  }
  return {
    H, V,
    inversions: [...inv.values()], axesPerp: [...axesPerp.values()],
    axesH: [...axesH.values()], axesV: [...axesV.values()],
    planesH: [...planesH.values()], planesV: [...planesV.values()], planesPar: [...planesPar.values()],
    axes4: [...axes4.values()], rotoinv4: [...rotoinv4.values()],
    diagPlanes: [...diagPlanes.values()], diagAxes: [...diagAxes.values()],
  };
}

// Height label: nothing if 0 is among the heights, else the smallest
// (¼ stands for ¼ and ¾, as in IT).
function _sgdHeightLabel(set) {
  if (!set || set.size === 0 || set.has('0')) return '';
  const order = ['⅛', '⅙', '¼', '⅓', '⅜', '½', '⅝', '⅔', '¾', '⅚', '⅞'];
  return [...set].sort((a, b) => order.indexOf(a) - order.indexOf(b))[0];
}

// General positions: images of a generic point; label = height along the
// viewing axis as (offset)± z, comma for operations with det < 0.
function sgdPositions(ops, view, start = { p: 0.14, q: 0.135, h: 0.17 }, margin = 0.24) {
  const H = (view + 1) % 3, V = (view + 2) % 3;
  const x0 = [0, 0, 0]; x0[H] = start.p; x0[V] = start.q; x0[view] = start.h;
  const pts = new Map();
  for (const op of ops) {
    const x = [0, 1, 2].map(i => op.R[i][0] * x0[0] + op.R[i][1] * x0[1] + op.R[i][2] * x0[2] + op.t[i] / 12);
    const s = op.R[view][view]; // ±1: the viewing axis is mapped onto itself
    const off = _sgdMod1(x[view] - s * start.h);
    const label = (off > 1e-9 ? _sgdFrac(off) : '') + (s > 0 ? '+' : '−');
    const comma = _sgdDet(op.R) < 0;
    const p = _sgdMod1(x[H]), q = _sgdMod1(x[V]);
    for (const dp of [-1, 0, 1]) for (const dq of [-1, 0, 1]) {
      const P = p + dp, Q = q + dq;
      if (P < -margin || P > 1 + margin || Q < -margin || Q > 1 + margin) continue;
      const key = `${P.toFixed(4)}|${Q.toFixed(4)}`;
      let e = pts.get(key);
      if (!e) { e = { p: P, q: Q, list: [] }; pts.set(key, e); }
      if (!e.list.some(o => o.label === label && o.comma === comma)) e.list.push({ label, comma, op });
    }
  }
  return [...pts.values()];
}

// Coordinate list: one operation per coset of the centring translations,
// in IT order (1, 2-folds z/y/x, −1, mirrors z/y/x; tetragonal see below),
// plus the centring vectors.
function sgdCoordinateList(ops) {
  const vars = ['x', 'y', 'z'];
  const centring = ops.filter(op => op.R.every((row, i) => row.every((v, j) => v === (i === j ? 1 : 0)))).map(op => op.t);
  const byR = new Map();
  for (const op of ops) {
    const key = op.R.map(r => r.join(',')).join(';');
    const sum = op.t.reduce((a, b) => a + b, 0);
    const cur = byR.get(key);
    if (!cur || sum < cur.t.reduce((a, b) => a + b, 0)) byR.set(key, op);
  }
  // IT order: 1, 2 ∥ c, 4⁺, 4⁻, 2 ∥ b, 2 ∥ a, 2 ∥ [110], 2 ∥ [1̄10], 1̄,
  // m ⊥ c, 4̄⁺, 4̄⁻, m ⊥ b, m ⊥ a, m ⊥ [110], m ⊥ [1̄10]; without a proper
  // 4-fold axis (4̄, 4̄2m, 4̄m2) the 4̄ operations take the place of 4⁺, 4⁻
  const has4 = ops.some(op => { const q = _sgdQuarterTurn(op); return q && q.proper; });
  const rank = op => {
    const qt = _sgdQuarterTurn(op);
    if (qt) return (qt.proper || !has4 ? 10.3 : 30.3) + (qt.plus ? 0 : 0.1);
    const sw = _sgdSwap(op);
    if (sw) return sw.mirror ? (sw.sigma < 0 ? 33 : 34) : (sw.sigma > 0 ? 12.3 : 12.4);
    const d = [0, 1, 2].map(i => op.R[i][i]);
    const nNeg = d.filter(v => v < 0).length;
    if (nNeg === 0) return 0;
    if (nNeg === 2) return 10 + (2 - d.indexOf(1));
    if (nNeg === 3) return 20;
    return 30 + (2 - d.indexOf(-1));
  };
  const reps = [...byR.values()].sort((a, b) => rank(a) - rank(b));
  // each row of the rotation part has exactly one ±1 here (diagonal or
  // quarter turn about c): ±(that coordinate) + translation
  const fmt = op => [0, 1, 2].map(i => {
    const j = op.R[i].findIndex(v => v !== 0);
    const s = op.R[i][j];
    const t = op.t[i] / 12;
    const v = vars[j];
    if (Math.abs(t) < 1e-9) return s > 0 ? v : `${v}̄`;
    return `${_sgdFrac(t)} ${s > 0 ? '+' : '−'} ${v}`;
  }).join(', ');
  const cv = centring.filter(t => t.some(v => v !== 0)).map(t => `(${t.map(v => _sgdFrac(v / 12)).join(', ')})`);
  const rKey = op => op.R.map(r => r.join(',')).join(';');
  const indexOf = op => reps.findIndex(r => rKey(r) === rKey(op)) + 1;
  return { lines: reps.map(fmt), symbols: reps.map(sgdOperationSymbol), centring: cv, format: fmt, indexOf };
}

// Short name of the symmetry operation of a (diagonal) Seitz operator,
// IT style: 1, 1̄, 2 / 2₁ with its axis, m / a b c n d with the plane normal.
function sgdOperationSymbol(op) {
  const names = ['a', 'b', 'c'];
  const t = op.t.map(v => _sgdMod1(v / 12));
  const qt = _sgdQuarterTurn(op);
  if (qt) {
    const sense = qt.plus ? '⁺' : '⁻';
    if (!qt.proper) return `<span class="sgdBar">4</span>${sense} ∥ c`;
    return `4${sense}${t[2] > 1e-9 ? `(0,0,${_sgdFrac(t[2])})` : ''} ∥ c`;
  }
  const sw = _sgdSwap(op);
  if (sw) {
    // intrinsic part w = (R + I)·t/2: (g, σg, t_z) for a plane, (g, σg, 0) for an axis
    const g = _sgdMod1(sw.sigma > 0 ? (t[0] + t[1]) / 2 : (t[0] - t[1]) / 2);
    if (!sw.mirror) return `${Math.abs(g - 0.5) < 1e-9 ? '2₁' : '2'} ∥ [${sw.sigma > 0 ? '110' : '1̄10'}]`;
    const gz = t[2];
    const q = v => Math.abs(v - 0.25) < 1e-9 || Math.abs(v - 0.75) < 1e-9;
    const vec = `(${_sgdFrac(g)},${sw.sigma > 0 ? _sgdFrac(g) : _sgdFrac(1 - g)},${_sgdFrac(gz)})`;
    const letter = q(g) || q(gz) ? `d${vec}` : g > 1e-9 && gz > 1e-9 ? 'n' : g > 1e-9 ? `g${vec}` : gz > 1e-9 ? 'c' : 'm';
    return `${letter} ⊥ [${sw.sigma > 0 ? '1̄10' : '110'}]`;
  }
  const d = [0, 1, 2].map(i => op.R[i][i]);
  const nNeg = d.filter(v => v < 0).length;
  if (nNeg === 0) return t.some(v => v > 1e-9) ? `t(${t.map(_sgdFrac).join(', ')})` : '1';
  if (nNeg === 3) return '<span class="sgdBar">1</span>';
  if (nNeg === 2) {
    const i = d.indexOf(1);
    return `${Math.abs(t[i] - 0.5) < 1e-9 ? '2₁' : '2'} ∥ ${names[i]}`;
  }
  const n = d.indexOf(-1);
  const g = [0, 1, 2].filter(j => j !== n && t[j] > 1e-9);
  const quarter = g.some(j => Math.abs(t[j] - 0.25) < 1e-9 || Math.abs(t[j] - 0.75) < 1e-9);
  const letter = quarter ? 'd' : g.length === 0 ? 'm' : g.length === 1 ? names[g[0]] : 'n';
  return `${letter} ⊥ ${names[n]}`;
}

// --- drawing ----------------------------------------------------------------

// SVG arc sweep flag for an arc from a through m to b (screen coordinates)
function _sgdSweep(a, m, b) {
  const cross = (m[0] - a[0]) * (b[1] - a[1]) - (m[1] - a[1]) * (b[0] - a[0]);
  return cross > 0 ? 1 : 0;
}

// comma inside a general-position circle: dot with a short tail
function _sgdComma(x, y, k = 1) {
  const r = 2.1 * k;
  return `<circle cx="${(x).toFixed(1)}" cy="${(y - 0.6 * k).toFixed(1)}" r="${r.toFixed(1)}" class="sgdCommaDot"/>`
    + `<path d="M ${(x + r * 0.9).toFixed(1)} ${(y - 0.2 * k).toFixed(1)} q ${(0.4 * k).toFixed(1)} ${(3.2 * k).toFixed(1)} ${(-2.6 * k).toFixed(1)} ${(4.6 * k).toFixed(1)}" class="sgdCommaTail"/>`;
}

function renderSpaceGroupDiagram(container, entry, cell) {
  let ops;
  try { ops = generateSpaceGroupOperators(entry.hall); } catch (e) { container.innerHTML = ''; return; }
  const short = sgdShortSymbol(entry.hm);
  const fullHM = formatSG(entry.hm.replace(/:.*$/, ''));
  if (!sgdIsSupported(ops)) {
    container.innerHTML = `<div class="sgdBox"><div class="sgdTitle"><span class="sgdShort">${short}</span> <span class="sgdNo">No. ${entry.number}</span></div>
      <p><small>The diagram is drawn for triclinic, monoclinic, orthorhombic and tetragonal space groups (Nos. 1–142) so far.</small></p></div>`;
    return;
  }
  const view = sgdViewAxis(ops);
  const el = sgdElements(ops, view);
  const { H, V } = el;
  const names = ['a', 'b', 'c'];
  const angleNames = ['α', 'β', 'γ'];
  // Fixed, schematic angle as in International Tables: 90° where symmetry
  // forces it (orthorhombic), a standard oblique 100° for the free angle of
  // monoclinic (unique axis = viewing axis) and triclinic cells. The real
  // angle of the data is only quoted in the note below.
  const sys = sgdCrystalSystem(ops);
  const ang = sys === 'orthorhombic' || sys === 'tetragonal' ? 90 : 100;
  const cellAngles = cell ? [cell.alpha, cell.beta, cell.gamma] : [null, null, null];
  const realAng = Number.isFinite(cellAngles[view]) ? cellAngles[view] : null;
  const rad = ang * Math.PI / 180;

  // screen frame: origin bottom-left; H to the right, V at angle `ang` from H
  const L = 300;
  const hx = L, hy = 0, vx = L * Math.cos(rad), vy = -L * Math.sin(rad);
  // extra room when a 2 and a 2₁ axis share a line in the drawing plane
  // (their arrows are drawn one behind the other, further out)
  const dup = (list, keyOf) => { const k = list.map(keyOf); return k.length !== new Set(k).size; };
  const pad = dup(el.axesH, a => a.c.toFixed(4)) || dup(el.axesV, a => a.c.toFixed(4)) || dup(el.diagAxes, a => `${a.sigma}|${a.c.toFixed(4)}`) ? 26 : 0;
  const left = 112 + 1.6 * pad + Math.max(0, -vx), top = 70 + pad;
  const Ox = left, Oy = top - vy; // vy < 0
  const P = (p, q) => [Ox + p * hx + q * vx, Oy + p * hy + q * vy];
  const width = left + L + Math.max(0, vx) + 96 + pad;
  const height = Oy + 112 + 1.6 * pad;

  const out = [];
  const f = n => n.toFixed(1);
  const line = (a, b, cls, extra = '') => `<line x1="${f(a[0])}" y1="${f(a[1])}" x2="${f(b[0])}" y2="${f(b[1])}" class="${cls}"${extra}/>`;

  // --- tooltip texts for the symmetry elements ---
  const X = ['x', 'y', 'z'];
  const fr = v => (v > 1 - 1e-9 ? '1' : _sgdFrac(v));
  const heightOrder = ['0', '⅛', '⅙', '¼', '⅓', '⅜', '½', '⅝', '⅔', '¾', '⅚', '⅞'];
  const hList = set => [...set].sort((a, b) => heightOrder.indexOf(a) - heightOrder.indexOf(b)).join(' and ');
  const glideText = (letters) => (letters.length === 1 ? `${letters[0]}/2` : `(${letters.join(' + ')})/2`);
  // "x = ¼, z = ¾" in x, y, z order from {axisIndex: text}
  const at = obj => Object.keys(obj).map(Number).sort().map(i => `${X[i]} = ${obj[i]}`).join(', ');
  const planeText = (kind, inPlaneAxis, normal, c) => {
    // kind for planes ⊥ drawing: m / dash (glide along the other in-plane axis) / dot (along view) / n / d
    const where = `at ${X[normal]} = ${fr(c)}`;
    if (kind === 'm') return `mirror plane ⊥ ${names[normal]} ${where}`;
    if (kind === 'd') return `d-glide plane ⊥ ${names[normal]} ${where}`;
    const letters = kind === 'dash' ? [names[inPlaneAxis]] : kind === 'dot' ? [names[view]] : [names[inPlaneAxis], names[view]].sort();
    const letter = letters.length === 1 ? letters[0] : 'n';
    return `${letter}-glide plane ⊥ ${names[normal]} ${where}, glide ${glideText(letters)}`;
  };
  const titled = (title, inner) => `<g class="sgdEl"><title>${title}</title>${inner}</g>`;
  const hitLine = (a, b) => line(a, b, 'sgdHitLine');

  // cell + half lines
  for (const q of [0, 0.5, 1]) out.push(line(P(0, q), P(1, q), q === 0.5 ? 'sgdGrid' : 'sgdCell'));
  for (const p of [0, 0.5, 1]) out.push(line(P(p, 0), P(p, 1), p === 0.5 ? 'sgdGrid' : 'sgdCell'));

  // planes perpendicular to the drawing
  const dash = { m: '', dash: '7 4', dot: '1.5 3.5', n: '8 3 1.5 3', e: '9 3 1.5 2.5 1.5 3', d: '8 3 1.5 3 1.5 3' };
  const planeLine = (a, b, kind) => line(a, b, `sgdPlane ${kind === 'm' ? 'mirror' : ''}`, dash[kind] ? ` stroke-dasharray="${dash[kind]}"` : '');
  const across = c => (c < 1e-9 ? [0, 1] : [c]);
  // planes on the same line: IT shows one symbol, by priority
  // m > e > a, b, c, g (dashed, dotted) > n > d; a dashed and a dotted glide
  // on the same line form a double glide plane e. The tooltip names all.
  const PRIO = { m: 0, e: 1, dash: 2, dot: 2, n: 3, d: 4 };
  const merged = (list, keyOf) => {
    const groups = new Map();
    for (const pl of list) { const k = keyOf(pl); if (!groups.has(k)) groups.set(k, []); groups.get(k).push(pl); }
    return [...groups.values()].map(g => {
      g.sort((a, b) => PRIO[a.kind] - PRIO[b.kind]);
      const kinds = new Set(g.map(p => p.kind));
      const kind = g[0].kind !== 'm' && kinds.has('dash') && kinds.has('dot') ? 'e' : g[0].kind;
      return { members: g, kind, c: g[0].c, sigma: g[0].sigma };
    });
  };
  const mergedText = (grp, textOf) => (grp.kind === 'e' ? 'double glide plane e: ' : '') + grp.members.map(textOf).join('; ');
  // d-glide plane: dash-dot-dot line with small arrowheads along it that
  // point in the direction of the in-plane glide component belonging to the
  // glide vector with +¼ along the viewing axis (as in IT); the centring
  // makes (u, ¼) and (−u, ¾) the same plane, so the sign is taken relative
  // to the component along the viewing axis
  const dArrows = (A, B, sign) => {
    const l = Math.hypot(B[0] - A[0], B[1] - A[1]);
    const u = [sign * (B[0] - A[0]) / l, sign * (B[1] - A[1]) / l], n = [-u[1], u[0]];
    let s = '';
    for (let k = 1; k <= 3; k++) {
      const m = [A[0] + (B[0] - A[0]) * k / 4 + u[0] * 4, A[1] + (B[1] - A[1]) * k / 4 + u[1] * 4];
      const back = [m[0] - u[0] * 8, m[1] - u[1] * 8];
      s += `<polygon points="${f(m[0])},${f(m[1])} ${f(back[0] + n[0] * 3.6)},${f(back[1] + n[1] * 3.6)} ${f(back[0] - n[0] * 3.6)},${f(back[1] - n[1] * 3.6)}" class="sgdArrowHead"/>`;
    }
    return s;
  };
  const signOf = (g, along) => (Math.abs(g - 0.75) < 1e-9 ? -1 : 1) * (Math.abs(along - 0.75) < 1e-9 ? -1 : 1);
  for (const grp of merged(el.planesH, pl => pl.c.toFixed(4))) for (const c of across(grp.c)) { // plane ⊥ H: x_H = c
    const t = mergedText(grp, pl => planeText(pl.kind, V, H, c));
    out.push(titled(t, hitLine(P(c, 0), P(c, 1)) + planeLine(P(c, 0), P(c, 1), grp.kind)
      + (grp.kind === 'd' ? dArrows(P(c, 0), P(c, 1), signOf(grp.members[0].gIn, grp.members[0].along)) : '')));
  }
  for (const grp of merged(el.planesV, pl => pl.c.toFixed(4))) for (const c of across(grp.c)) {
    const t = mergedText(grp, pl => planeText(pl.kind, H, V, c));
    out.push(titled(t, hitLine(P(0, c), P(1, c)) + planeLine(P(0, c), P(1, c), grp.kind)
      + (grp.kind === 'd' ? dArrows(P(0, c), P(1, c), signOf(grp.members[0].gIn, grp.members[0].along)) : '')));
  }

  // diagonal planes and axes (tetragonal): pieces of the line x − y = c
  // (σ = +1, running along [110]) or x + y = c (σ = −1, along [1̄10]) inside
  // the cell; each line comes with its translate by one cell edge
  const diagSegs = (sigma, c) => {
    const segs = [];
    if (sigma > 0) {
      for (const cc of [c, c - 1]) if (cc > -1 + 1e-9) segs.push(cc >= 0 ? [[cc, 0], [1, 1 - cc]] : [[0, -cc], [1 + cc, 1]]);
    } else {
      for (const s of [c, c + 1]) if (s > 1e-9 && s < 2 - 1e-9) segs.push(s <= 1 ? [[s, 0], [0, s]] : [[1, s - 1], [s - 1, 1]]);
    }
    return segs.filter(([a, b]) => Math.hypot(a[0] - b[0], a[1] - b[1]) > 1e-6);
  };
  const diagWhere = (sigma, c) => `${sigma > 0 ? 'x − y' : 'x + y'} = ${fr(c)}`;
  const diagNormal = sigma => (sigma > 0 ? '[1̄10]' : '[110]');
  const diagDir = sigma => (sigma > 0 ? '[110]' : '[1̄10]');
  const diagText = pl => {
    const vec = `(${_sgdFrac(pl.g)}, ${_sgdFrac(pl.sigma > 0 ? pl.g : 1 - pl.g)}, ${_sgdFrac(pl.gz)})`;
    const what = pl.kind === 'm' ? 'mirror plane' : `${{ dash: 'g', dot: 'c', n: 'n', d: 'd' }[pl.kind]}-glide plane`;
    return `${what} ⊥ ${diagNormal(pl.sigma)} at ${diagWhere(pl.sigma, pl.c)}${pl.kind === 'm' ? '' : `, glide ${vec}`}`;
  };
  for (const grp of merged(el.diagPlanes, pl => `${pl.sigma}|${pl.c.toFixed(4)}`)) {
    const t = mergedText(grp, diagText);
    // segments run with +x along [110] (σ = +1) and with −x along [1̄10]
    // (σ = −1); the in-plane glide is (g, σg)
    const sg = grp.kind === 'd' ? signOf(grp.members[0].g, grp.members[0].gz) * (grp.sigma > 0 ? 1 : -1) : 0;
    for (const [a, b] of diagSegs(grp.sigma, grp.c)) out.push(titled(t, hitLine(P(...a), P(...b)) + planeLine(P(...a), P(...b), grp.kind)
      + (sg ? dArrows(P(...a), P(...b), sg) : '')));
  }

  // unit vectors on screen
  const uH = [1, 0], lenV = Math.hypot(vx, vy), uV = [vx / lenV, vy / lenV];
  const nH = [uH[1], -uH[0]]; // perpendicular, for arrow barbs

  // in-plane 2-fold axes: arrows outside the cell at both ends
  const arrow = (base, dir, screw, labelText, gap = 14) => {
    const len = 30;
    const a = [base[0] + dir[0] * gap, base[1] + dir[1] * gap];
    const b = [a[0] + dir[0] * len, a[1] + dir[1] * len];
    const nx = -dir[1], ny = dir[0];
    const hl = 8, hw = 4;
    const back = [b[0] - dir[0] * hl, b[1] - dir[1] * hl];
    const head = screw
      ? `${f(b[0])},${f(b[1])} ${f(back[0] + nx * hw)},${f(back[1] + ny * hw)} ${f(back[0])},${f(back[1])}`
      : `${f(b[0])},${f(b[1])} ${f(back[0] + nx * hw)},${f(back[1] + ny * hw)} ${f(back[0] - nx * hw)},${f(back[1] - ny * hw)}`;
    let s = line(a, back, 'sgdAxisArrow') + `<polygon points="${head}" class="sgdArrowHead"/>`;
    if (labelText) {
      const lx = b[0] + dir[0] * 6 + (Math.abs(dir[0]) > 0.5 ? (dir[0] > 0 ? 2 : -10) : 4);
      const ly = b[1] + dir[1] * 6 + (Math.abs(dir[1]) > 0.5 ? (dir[1] > 0 ? 10 : 0) : 4);
      s += `<text x="${f(lx)}" y="${f(ly)}" class="sgdSmall"${Math.abs(dir[0]) > 0.5 && dir[0] < 0 ? ' text-anchor="end"' : ''}>${labelText}</text>`;
    }
    return s;
  };
  const axisName = (screw, along) => `${screw ? '2₁ screw axis' : '2-fold rotation axis'} ∥ ${names[along]}`;
  // a 2 and a 2₁ axis on the same line (centred cells, 4/mmm): two arrows
  // one behind the other, the rotation axis next to the cell (as in IT)
  const outer = (list, keyOf) => {
    const seen = new Map();
    for (const ax of list) { const k = keyOf(ax); seen.set(k, (seen.get(k) || 0) + 1); }
    return ax => ax.screw && seen.get(keyOf(ax)) > 1;
  };
  const GAP2 = 14 + 30 + 10;
  const isOuterH = outer(el.axesH, ax => ax.c.toFixed(4)), isOuterV = outer(el.axesV, ax => ax.c.toFixed(4));
  for (const ax of el.axesH) {
    const lab = _sgdHeightLabel(ax.heights), gap = isOuterH(ax) ? GAP2 : 14;
    for (const c of across(ax.c)) {
      const t = `${axisName(ax.screw, H)} at ${at({ [V]: fr(c), [view]: hList(ax.heights) })}`;
      out.push(titled(t, arrow(P(1, c), uH, ax.screw, lab, gap)));
      out.push(titled(t, arrow(P(0, c), [-uH[0], -uH[1]], ax.screw, lab, gap)));
    }
  }
  for (const ax of el.axesV) {
    const lab = _sgdHeightLabel(ax.heights), gap = isOuterV(ax) ? GAP2 : 14;
    for (const c of across(ax.c)) {
      const t = `${axisName(ax.screw, V)} at ${at({ [H]: fr(c), [view]: hList(ax.heights) })}`;
      out.push(titled(t, arrow(P(c, 1), uV, ax.screw, lab, gap)));
      out.push(titled(t, arrow(P(c, 0), [-uV[0], -uV[1]], ax.screw, lab, gap)));
    }
  }
  // diagonal 2-fold axes in the drawing plane: arrows outside the cell at
  // both ends of each piece of their line
  const isOuterD = outer(el.diagAxes, ax => `${ax.sigma}|${ax.c.toFixed(4)}`);
  for (const ax of el.diagAxes) {
    const lab = _sgdHeightLabel(ax.heights), gap = isOuterD(ax) ? GAP2 : 14;
    const t = `${ax.screw ? '2₁ screw axis' : '2-fold rotation axis'} ∥ ${diagDir(ax.sigma)} at ${diagWhere(ax.sigma, ax.c)}, z = ${hList(ax.heights)}`;
    for (const [a, b] of diagSegs(ax.sigma, ax.c)) {
      const A = P(...a), B = P(...b), l = Math.hypot(B[0] - A[0], B[1] - A[1]);
      const u = [(B[0] - A[0]) / l, (B[1] - A[1]) / l];
      out.push(titled(t, arrow(B, u, ax.screw, lab, gap)));
      out.push(titled(t, arrow(A, [-u[0], -u[1]], ax.screw, lab, gap)));
    }
  }

  // 2-fold axes perpendicular to the drawing: lens (2) or lens with hooks (2₁),
  // long axis along the vertical cell edge
  // 2-fold axis perpendicular to the drawing: a lens made of two circular
  // arcs, long axis along the vertical cell edge. 2₁: each arc runs on past
  // its tip as a short hook (top one to the right, bottom one to the left).
  const lens = (p, q, screw, k = 1) => {
    const [x, y] = P(p, q);
    const u = uV, w = [-uV[1], uV[0]]; // along the axis / to the right of it
    const Lh = 11 * k, W = 4.6 * k;
    const c = (Lh * Lh - W * W) / (2 * W), rho = c + W, th0 = Math.asin(Lh / rho);
    const pt = (ww, uu) => [x + ww * w[0] + uu * u[0], y + ww * w[1] + uu * u[1]];
    const leftArc = th => pt(c - rho * Math.cos(th), rho * Math.sin(th));   // bulges to the left
    const rightArc = th => pt(-c + rho * Math.cos(th), rho * Math.sin(th)); // bulges to the right
    const top = leftArc(th0), bottom = leftArc(-th0);
    const sweepL = _sgdSweep(leftArc(-th0), leftArc(0), leftArc(th0));
    const sweepR = _sgdSweep(rightArc(th0), rightArc(0), rightArc(-th0));
    let s = `<path d="M ${f(bottom[0])} ${f(bottom[1])} A ${f(rho)} ${f(rho)} 0 0 ${sweepL} ${f(top[0])} ${f(top[1])} A ${f(rho)} ${f(rho)} 0 0 ${sweepR} ${f(bottom[0])} ${f(bottom[1])} Z" class="sgdLens"/>`;
    if (screw) {
      // the hooks are the lens arcs drawn on past the tips: each starts
      // inside the lens on its arc (shifted inwards by half the stroke so it
      // lies flush with the filled edge) and runs on past the tip, so that
      // it grows out of the lens (as in IT A)
      const d = 0.5, back = 0.45, r = rho - 0.75; // radians past / before the tip, flush radius
      const lArc = th => pt(c - r * Math.cos(th), r * Math.sin(th));
      const rArc = th => pt(-c + r * Math.cos(th), r * Math.sin(th));
      const hook = (arc, t0, t1) => {
        const a = arc(t0), b = arc(t1), sw = _sgdSweep(a, arc((t0 + t1) / 2), b);
        return `<path d="M ${f(a[0])} ${f(a[1])} A ${f(r)} ${f(r)} 0 0 ${sw} ${f(b[0])} ${f(b[1])}" class="sgdHook"/>`;
      };
      s += hook(lArc, th0 - back, th0 + d);   // top: left arc runs on to the right
      s += hook(rArc, -th0 + back, -th0 - d); // bottom: right arc runs on to the left
    }
    return s;
  };
  const inRange = (v) => v >= -1e-9 && v <= 1 + 1e-9;
  const spread = (p, q) => {
    const res = [];
    for (const dp of [0, 1]) for (const dq of [0, 1]) { const a = p + dp, b = q + dq; if (inRange(a) && inRange(b)) res.push([a, b]); }
    return res;
  };
  // a 4-fold or 4̄ axis contains the 2-fold axis at the same place (4² = 2,
  // 4₁² = 2₁, 4̄² = 2): only the higher symbol is drawn
  const pqKey = (p, q) => `${_sgdMod1(p).toFixed(4)}|${_sgdMod1(q).toFixed(4)}`;
  const quadSites = new Set([...el.axes4, ...el.rotoinv4].map(a => pqKey(a.p, a.q)));
  for (const ax of el.axesPerp) for (const [p, q] of spread(ax.p, ax.q)) {
    if (quadSites.has(pqKey(p, q))) continue;
    const [x, y] = P(p, q);
    const t = `${axisName(ax.screw, view)} at ${at({ [H]: fr(p), [V]: fr(q) })}`;
    out.push(titled(t, `<circle cx="${f(x)}" cy="${f(y)}" r="12" class="sgdHit"/>` + lens(p, q, ax.screw)));
  }

  // 4-fold axes perpendicular to the drawing, as in IT A: a square standing
  // on its corner (diamond, corners on the cell axes); filled for 4; for
  // 4₁ / 4₃ every edge runs on straight past its end corner as a tail
  // (pinwheel counter-clockwise for 4₁, clockwise for 4₃); 4₂: two tails at
  // opposite corners (top and bottom); 4̄: open diamond with a small filled lens (its 2-fold
  // axis) inside. An inversion centre on the axis is drawn on top as the
  // usual small circle (→ 4/m and 4₂/m symbols).
  const quad = (p, q, kind) => {
    const [x, y] = P(p, q);
    const a = 11; // centre to corner
    const pt = (s, r) => [x + s * uH[0] + r * uV[0], y + s * uH[1] + r * uV[1]];
    const c = [pt(a, 0), pt(0, a), pt(-a, 0), pt(0, -a)]; // right, top, left, bottom: counter-clockwise on screen
    let s = `<path d="M ${c.map(v => `${f(v[0])} ${f(v[1])}`).join(' L ')} Z" class="${kind === -1 ? 'sgdQuadOpen' : 'sgdQuad'}"/>`;
    if (kind === -1) return s + lens(p, q, false, 0.6);
    // the tail is the edge itself drawn on: it starts well inside the edge
    // and runs half an edge past the corner, shifted inwards by half its
    // stroke so that it lies flush with the filled edge (grows out of it)
    const tail = (k, ccw) => {
      const from = c[(k + (ccw ? 3 : 1)) % 4], at = c[k];
      const e = [at[0] - from[0], at[1] - from[1]], el = Math.hypot(e[0], e[1]);
      const u = [e[0] / el, e[1] / el];
      let n = [-u[1], u[0]];
      if ((at[0] - x) * n[0] + (at[1] - y) * n[1] > 0) n = [-n[0], -n[1]]; // inward normal
      const sh = 0.75; // half the stroke width
      const a0 = [at[0] - e[0] * 0.45 + n[0] * sh, at[1] - e[1] * 0.45 + n[1] * sh];
      const a1 = [at[0] + e[0] * 0.5 + n[0] * sh, at[1] + e[1] * 0.5 + n[1] * sh];
      return `<path d="M ${f(a0[0])} ${f(a0[1])} L ${f(a1[0])} ${f(a1[1])}" class="sgdTail"/>`;
    };
    if (kind === 1) for (const k of [0, 1, 2, 3]) s += tail(k, true);
    if (kind === 3) for (const k of [0, 1, 2, 3]) s += tail(k, false);
    if (kind === 2) for (const k of [1, 3]) s += tail(k, false);
    return s;
  };
  const quadName = n => (n === 0 ? '4-fold rotation axis' : `4${['', '₁', '₂', '₃'][n]} screw axis`);
  for (const ax of el.axes4) {
    const kind = Math.min(...ax.kinds);
    for (const [p, q] of spread(ax.p, ax.q)) {
      const [x, y] = P(p, q);
      const t = `${quadName(kind)} ∥ c at ${at({ [H]: fr(p), [V]: fr(q) })}`;
      out.push(titled(t, `<circle cx="${f(x)}" cy="${f(y)}" r="15" class="sgdHit"/>` + quad(p, q, kind)));
    }
  }
  const sites4 = new Set(el.axes4.map(a => pqKey(a.p, a.q)));
  for (const ax of el.rotoinv4) {
    if (sites4.has(pqKey(ax.p, ax.q))) continue; // 4/m: the 4-fold symbol stands for both
    const lab = _sgdHeightLabel(ax.heights);
    for (const [p, q] of spread(ax.p, ax.q)) {
      const [x, y] = P(p, q);
      const t = `4̄ rotoinversion axis ∥ c at ${at({ [H]: fr(p), [V]: fr(q) })}, inversion points at z = ${hList(ax.heights)}`;
      out.push(titled(t, `<circle cx="${f(x)}" cy="${f(y)}" r="15" class="sgdHit"/>` + quad(p, q, -1)
        + (lab ? `<text x="${f(x + 12)}" y="${f(y - 9)}" class="sgdSmall">${lab}</text>` : '')));
    }
  }

  // inversion centres: small open circles, height if not 0
  for (const ic of el.inversions) {
    const lab = _sgdHeightLabel(ic.heights);
    for (const [p, q] of spread(ic.p, ic.q)) {
      const [x, y] = P(p, q);
      const coords = [...ic.heights].sort((a, b) => heightOrder.indexOf(a) - heightOrder.indexOf(b)).map(h => {
        const v = []; v[H] = fr(p); v[V] = fr(q); v[view] = h; return `(${v.join(', ')})`;
      });
      const t = `inversion centre at ${coords.join(' and ')}`;
      out.push(titled(t, `<circle cx="${f(x)}" cy="${f(y)}" r="9" class="sgdHit"/><circle cx="${f(x)}" cy="${f(y)}" r="3.8" class="sgdInv"/>`
        + (lab ? `<text x="${f(x - 7)}" y="${f(y - 6)}" class="sgdSmall" text-anchor="end">${lab}</text>` : '')));
    }
  }

  // planes parallel to the drawing: corner symbols below-left of the origin.
  // Several such planes are nested as in IT: the first (mirror first, then
  // glides a, b, n, d) is the outer "L" with its height above the vertical
  // leg, every further one sits inside it, moved in along the diagonal with
  // its legs shortened so that their ends line up, its height to the right.
  const KORD = { m: 0, h: 1, v: 2, n: 3, d: 4 };
  el.planesPar.sort((a, b) => KORD[a.kind] - KORD[b.kind] || _sgdHeightLabel(a.heights).localeCompare(_sgdHeightLabel(b.heights)));
  const nested = el.planesPar.length > 1;
  el.planesPar.forEach((pl, i) => {
    // the "L" follows the cell edges (H to the right, V along the tilted
    // edge); the glide arrow starts in its corner: along H, along V, or along
    // the diagonal H + V for an n-glide
    const S0 = 30, inset = 10 * i;
    const bx0 = Ox - 84 - pad, by0 = Oy + 84 + pad;
    const bx = bx0 + uH[0] * inset + uV[0] * inset, by = by0 + uH[1] * inset + uV[1] * inset;
    const s = S0 - inset;
    const mirror = pl.kind === 'm';
    const vEnd = [bx + uV[0] * s, by + uV[1] * s];
    {
      const letters = pl.kind === 'h' ? [names[H]] : pl.kind === 'v' ? [names[V]] : pl.kind === 'n' ? [names[H], names[V]].sort() : [];
      const what = mirror ? 'mirror plane' : pl.kind === 'd' ? 'd-glide plane' : `${letters.length === 1 ? letters[0] : 'n'}-glide plane`;
      const t = `${what} ⊥ ${names[view]} (parallel to the drawing) at ${X[view]} = ${hList(pl.heights)}${letters.length ? `, glide ${glideText(letters)}` : ''}`;
      out.push(`<g class="sgdEl"><title>${t}</title><rect x="${f(bx - 4)}" y="${f(by - s - 4)}" width="${f(s + 10)}" height="${f(s + 8)}" class="sgdHit"/>`);
    }
    out.push(`<path d="M ${f(vEnd[0])} ${f(vEnd[1])} L ${f(bx)} ${f(by)} L ${f(bx + s)} ${f(by)}" class="sgdCorner${mirror ? ' mirror' : ''}"/>`);
    if (!mirror) {
      let dir;
      if (pl.kind === 'h') dir = uH;
      else if (pl.kind === 'v') dir = uV;
      else { const d = [uH[0] + uV[0], uH[1] + uV[1]], l = Math.hypot(d[0], d[1]); dir = [d[0] / l, d[1] / l]; }
      const len = pl.kind === 'h' || pl.kind === 'v' ? s + 4 : s * 0.95;
      const from = [bx, by];
      const to = [from[0] + dir[0] * len, from[1] + dir[1] * len];
      const nx = -dir[1], ny = dir[0], back = [to[0] - dir[0] * 7, to[1] - dir[1] * 7];
      out.push(line(from, back, 'sgdCornerArrow'));
      out.push(`<polygon points="${f(to[0])},${f(to[1])} ${f(back[0] + nx * 3.5)},${f(back[1] + ny * 3.5)} ${f(back[0] - nx * 3.5)},${f(back[1] - ny * 3.5)}" class="sgdArrowHead"/>`);
    }
    const lab = _sgdHeightLabel(pl.heights);
    if (lab) {
      if (nested && i === 0) out.push(`<text x="${f(vEnd[0])}" y="${f(vEnd[1] - 5)}" class="sgdSmall" text-anchor="middle">${lab}</text>`);
      else out.push(`<text x="${f(bx + s + 8)}" y="${f(by + 5)}" class="sgdSmall">${lab}</text>`);
    }
    out.push('</g>');
  });

  // general positions (tooltip: which line of the list, and its coordinates)
  const R = 9.5;
  const list = sgdCoordinateList(ops);
  const tip = it => `(${list.indexOf(it.op)}) ${list.format(it.op)}`;
  // tetragonal: start point as in IT A (images cluster around the 4-fold
  // axes), off the diagonals and the lines x, y = ¼, ½ (mirrors in 4/mmm)
  const startPt = sys === 'tetragonal' ? { p: 0.16, q: 0.09, h: 0.17 } : undefined;
  for (const g of sgdPositions(ops, view, startPt, sys === 'tetragonal' ? 0.18 : undefined)) {
    const [x, y] = P(g.p, g.q);
    out.push(`<g class="sgdPosGroup"><title>${g.list.map(tip).join(' and ')}</title>`);
    out.push(`<circle cx="${f(x)}" cy="${f(y)}" r="${R + 3}" class="sgdHit"/>`);
    if (g.list.length === 1) {
      const it = g.list[0];
      out.push(`<circle cx="${f(x)}" cy="${f(y)}" r="${R}" class="sgdPos"/>`);
      if (it.comma) out.push(_sgdComma(x, y));
      out.push(`<text x="${f(x + R + 1.5)}" y="${f(y - R + 4)}" class="sgdPosLabel">${it.label}</text>`);
    } else {
      out.push(`<circle cx="${f(x)}" cy="${f(y)}" r="${R}" class="sgdPos"/>`, line([x, y - R], [x, y + R], 'sgdPosSplit'));
      [...g.list].sort((a, b) => (b.comma ? 1 : 0) - (a.comma ? 1 : 0)).slice(0, 2).forEach((it, i) => {
        const sx = i ? 1 : -1;
        if (it.comma) out.push(_sgdComma(x + sx * 4, y, 0.75));
        out.push(`<text x="${f(x + sx * (R + 1.5))}" y="${f(y - R + 4)}" class="sgdPosLabel"${i ? '' : ' text-anchor="end"'}>${it.label}</text>`);
      });
    }
    out.push('</g>');
  }

  // axis labels and angle
  const smallArrow = (from, dir, len) => {
    const to = [from[0] + dir[0] * len, from[1] + dir[1] * len];
    const back = [to[0] - dir[0] * 6, to[1] - dir[1] * 6], nx = -dir[1], ny = dir[0];
    return line(from, back, 'sgdAxisLabelArrow') + `<polygon points="${f(to[0])},${f(to[1])} ${f(back[0] + nx * 3)},${f(back[1] + ny * 3)} ${f(back[0] - nx * 3)},${f(back[1] - ny * 3)}" class="sgdArrowHead"/>`;
  };
  const mid = P(0, 0.5);
  const aBase = [mid[0] - 100 - 1.6 * pad, mid[1] + 12];
  out.push(smallArrow(aBase, uV, 22));
  out.push(`<text x="${f(aBase[0] + 8)}" y="${f(aBase[1] - 4)}" class="sgdAxisName">${names[V]}</text>`);
  const bm = P(0.5, 0);
  out.push(smallArrow([bm[0] - 12, Oy + 92 + 1.6 * pad], uH, 24));
  out.push(`<text x="${f(bm[0] - 3)}" y="${f(Oy + 86 + 1.6 * pad)}" class="sgdAxisName">${names[H]}</text>`);
  if (ang !== 90) {
    const r = 15;
    const a0 = 0, a1 = -rad; // from H towards V (screen y up = negative)
    const p0 = [Ox + r * Math.cos(a0), Oy + r * Math.sin(a0)], p1 = [Ox + r * Math.cos(a1), Oy + r * Math.sin(a1)];
    out.push(`<path d="M ${f(p0[0])} ${f(p0[1])} A ${r} ${r} 0 0 0 ${f(p1[0])} ${f(p1[1])}" class="sgdAngle"/>`);
    const bis = -rad / 2;
    out.push(`<text x="${f(Ox + 27 * Math.cos(bis) - 5)}" y="${f(Oy + 27 * Math.sin(bis) + 5)}" class="sgdAngleName">${angleNames[view]}</text>`);
  }

  const setting = entry.qualifier ? ` <span class="sgdSetting">${/^[12]$/.test(entry.qualifier) ? 'origin choice' : 'setting'} ${entry.qualifier}</span>` : '';
  container.innerHTML = `
    <div class="sgdBox">
      <div class="sgdTitle"><span class="sgdShort">${short}</span><span class="sgdFull">${fullHM}</span><span class="sgdPg">${sgdPointGroup(ops)}</span><span class="sgdNo">No. ${entry.number}</span>${setting}</div>
      <div class="sgdBody">
        <svg viewBox="0 0 ${f(width)} ${f(height)}" class="sgDiagram" role="img" aria-label="Space-group diagram of ${entry.hm}">${out.join('')}</svg>
        <div class="sgdListBox">
          <div class="sgdListHead">General position <span>(multiplicity ${ops.length})</span></div>
          <div class="sgdListSub">coordinates of the symmetry-equivalent copies of the start point (1) — each red circle is one of them; right: the operation that produces it</div>
          <ol class="sgdList">${list.lines.map((l, i) => `<li><span class="sgdXYZ">${l}</span><span class="sgdOpSym">${list.symbols[i]}</span></li>`).join('')}${list.centring.length ? `<li class="sgdCentring">+ ${list.centring.join(' + ')} <span class="sgdOpSym">(centring)</span></li>` : ''}</ol>
        </div>
      </div>
      <p class="sgdNote"><small>Projection along ${names[view]}; ${names[V]} up, ${names[H]} to the right. Schematic cell${ang !== 90 ? ` (${angleNames[view]} drawn as ${ang}°${realAng ? `, here ${realAng.toFixed(2)}°` : ''})` : ''};
        symmetry elements in the symbols of International Tables Vol. A, general positions in red (height along ${names[view]}; comma = mirror image). Hover a symbol or circle for details. Derived from the operators of ${entry.hall}.</small></p>
    </div>`;
}
