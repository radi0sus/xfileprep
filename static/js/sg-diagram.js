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

// Only groups whose rotation parts are all diagonal (triclinic, monoclinic,
// orthorhombic in any axis setting) are drawn.
function sgdIsSupported(ops) {
  return ops.every(op => op.R.every((row, i) => row.every((v, j) => (i === j ? Math.abs(v) === 1 : v === 0))));
}

// Viewing axis: the unique axis for monoclinic, c otherwise.
function sgdViewAxis(ops) {
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
  for (const op of ops) {
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
          add(n === H ? planesH : planesV, `${c}|${kind}`, { c, kind });
        }
      }
    }
  }
  return {
    H, V,
    inversions: [...inv.values()], axesPerp: [...axesPerp.values()],
    axesH: [...axesH.values()], axesV: [...axesV.values()],
    planesH: [...planesH.values()], planesV: [...planesV.values()], planesPar: [...planesPar.values()],
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
function sgdPositions(ops, view, start = { p: 0.14, q: 0.135, h: 0.17 }) {
  const H = (view + 1) % 3, V = (view + 2) % 3;
  const x0 = [0, 0, 0]; x0[H] = start.p; x0[V] = start.q; x0[view] = start.h;
  const pts = new Map();
  for (const op of ops) {
    const x = [0, 1, 2].map(i => op.R[i][i] * x0[i] + op.t[i] / 12);
    const s = op.R[view][view];
    const off = _sgdMod1(x[view] - s * start.h);
    const label = (off > 1e-9 ? _sgdFrac(off) : '') + (s > 0 ? '+' : '−');
    const comma = op.R[0][0] * op.R[1][1] * op.R[2][2] < 0;
    const p = _sgdMod1(x[H]), q = _sgdMod1(x[V]);
    for (const dp of [-1, 0, 1]) for (const dq of [-1, 0, 1]) {
      const P = p + dp, Q = q + dq;
      if (P < -0.24 || P > 1.24 || Q < -0.24 || Q > 1.24) continue;
      const key = `${P.toFixed(4)}|${Q.toFixed(4)}`;
      let e = pts.get(key);
      if (!e) { e = { p: P, q: Q, list: [] }; pts.set(key, e); }
      if (!e.list.some(o => o.label === label && o.comma === comma)) e.list.push({ label, comma, op });
    }
  }
  return [...pts.values()];
}

// Coordinate list: one operation per coset of the centring translations,
// in IT order (1, 2-folds z/y/x, −1, mirrors z/y/x), plus the centring vectors.
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
  const rank = op => {
    const d = [0, 1, 2].map(i => op.R[i][i]);
    const nNeg = d.filter(v => v < 0).length;
    if (nNeg === 0) return 0;
    if (nNeg === 2) return 10 + (2 - d.indexOf(1));
    if (nNeg === 3) return 20;
    return 30 + (2 - d.indexOf(-1));
  };
  const reps = [...byR.values()].sort((a, b) => rank(a) - rank(b));
  const fmt = op => [0, 1, 2].map(i => {
    const s = op.R[i][i];
    const t = op.t[i] / 12;
    const v = vars[i];
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
  const d = [0, 1, 2].map(i => op.R[i][i]);
  const t = op.t.map(v => _sgdMod1(v / 12));
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
      <p><small>The diagram is drawn for triclinic, monoclinic and orthorhombic space groups only.</small></p></div>`;
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
  const ang = sys === 'orthorhombic' ? 90 : 100;
  const cellAngles = cell ? [cell.alpha, cell.beta, cell.gamma] : [null, null, null];
  const realAng = Number.isFinite(cellAngles[view]) ? cellAngles[view] : null;
  const rad = ang * Math.PI / 180;

  // screen frame: origin bottom-left; H to the right, V at angle `ang` from H
  const L = 300;
  const hx = L, hy = 0, vx = L * Math.cos(rad), vy = -L * Math.sin(rad);
  const left = 112 + Math.max(0, -vx), top = 70;
  const Ox = left, Oy = top - vy; // vy < 0
  const P = (p, q) => [Ox + p * hx + q * vx, Oy + p * hy + q * vy];
  const width = left + L + Math.max(0, vx) + 96;
  const height = Oy + 112;

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
  const dash = { m: '', dash: '7 4', dot: '1.5 3.5', n: '8 3 1.5 3', d: '8 3 1.5 3 1.5 3' };
  const planeLine = (a, b, kind) => line(a, b, `sgdPlane ${kind === 'm' ? 'mirror' : ''}`, dash[kind] ? ` stroke-dasharray="${dash[kind]}"` : '');
  const across = c => (c < 1e-9 ? [0, 1] : [c]);
  for (const pl of el.planesH) for (const c of across(pl.c)) { // plane ⊥ H: x_H = c
    const t = planeText(pl.kind, V, H, c);
    out.push(titled(t, hitLine(P(c, 0), P(c, 1)) + planeLine(P(c, 0), P(c, 1), pl.kind)));
  }
  for (const pl of el.planesV) for (const c of across(pl.c)) {
    const t = planeText(pl.kind, H, V, c);
    out.push(titled(t, hitLine(P(0, c), P(1, c)) + planeLine(P(0, c), P(1, c), pl.kind)));
  }

  // unit vectors on screen
  const uH = [1, 0], lenV = Math.hypot(vx, vy), uV = [vx / lenV, vy / lenV];
  const nH = [uH[1], -uH[0]]; // perpendicular, for arrow barbs

  // in-plane 2-fold axes: arrows outside the cell at both ends
  const arrow = (base, dir, screw, labelText) => {
    const len = 30, gap = 14;
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
  for (const ax of el.axesH) {
    const lab = _sgdHeightLabel(ax.heights);
    for (const c of across(ax.c)) {
      const t = `${axisName(ax.screw, H)} at ${at({ [V]: fr(c), [view]: hList(ax.heights) })}`;
      out.push(titled(t, arrow(P(1, c), uH, ax.screw, lab)));
      out.push(titled(t, arrow(P(0, c), [-uH[0], -uH[1]], ax.screw, lab)));
    }
  }
  for (const ax of el.axesV) {
    const lab = _sgdHeightLabel(ax.heights);
    for (const c of across(ax.c)) {
      const t = `${axisName(ax.screw, V)} at ${at({ [H]: fr(c), [view]: hList(ax.heights) })}`;
      out.push(titled(t, arrow(P(c, 1), uV, ax.screw, lab)));
      out.push(titled(t, arrow(P(c, 0), [-uV[0], -uV[1]], ax.screw, lab)));
    }
  }

  // 2-fold axes perpendicular to the drawing: lens (2) or lens with hooks (2₁),
  // long axis along the vertical cell edge
  // 2-fold axis perpendicular to the drawing: a lens made of two circular
  // arcs, long axis along the vertical cell edge. 2₁: each arc runs on past
  // its tip as a short hook (top one to the right, bottom one to the left).
  const lens = (p, q, screw) => {
    const [x, y] = P(p, q);
    const u = uV, w = [-uV[1], uV[0]]; // along the axis / to the right of it
    const Lh = 11, W = 4.6;
    const c = (Lh * Lh - W * W) / (2 * W), rho = c + W, th0 = Math.asin(Lh / rho);
    const pt = (ww, uu) => [x + ww * w[0] + uu * u[0], y + ww * w[1] + uu * u[1]];
    const leftArc = th => pt(c - rho * Math.cos(th), rho * Math.sin(th));   // bulges to the left
    const rightArc = th => pt(-c + rho * Math.cos(th), rho * Math.sin(th)); // bulges to the right
    const top = leftArc(th0), bottom = leftArc(-th0);
    const sweepL = _sgdSweep(leftArc(-th0), leftArc(0), leftArc(th0));
    const sweepR = _sgdSweep(rightArc(th0), rightArc(0), rightArc(-th0));
    let s = `<path d="M ${f(bottom[0])} ${f(bottom[1])} A ${f(rho)} ${f(rho)} 0 0 ${sweepL} ${f(top[0])} ${f(top[1])} A ${f(rho)} ${f(rho)} 0 0 ${sweepR} ${f(bottom[0])} ${f(bottom[1])} Z" class="sgdLens"/>`;
    if (screw) {
      const d = 0.62; // how far the arcs run on past the tips (radians)
      const h1 = leftArc(th0 + d), h2 = rightArc(-th0 - d);
      const s1 = _sgdSweep(top, leftArc(th0 + d / 2), h1), s2 = _sgdSweep(bottom, rightArc(-th0 - d / 2), h2);
      s += `<path d="M ${f(top[0])} ${f(top[1])} A ${f(rho)} ${f(rho)} 0 0 ${s1} ${f(h1[0])} ${f(h1[1])}" class="sgdHook"/>`;
      s += `<path d="M ${f(bottom[0])} ${f(bottom[1])} A ${f(rho)} ${f(rho)} 0 0 ${s2} ${f(h2[0])} ${f(h2[1])}" class="sgdHook"/>`;
    }
    return s;
  };
  const inRange = (v) => v >= -1e-9 && v <= 1 + 1e-9;
  const spread = (p, q) => {
    const res = [];
    for (const dp of [0, 1]) for (const dq of [0, 1]) { const a = p + dp, b = q + dq; if (inRange(a) && inRange(b)) res.push([a, b]); }
    return res;
  };
  for (const ax of el.axesPerp) for (const [p, q] of spread(ax.p, ax.q)) {
    const [x, y] = P(p, q);
    const t = `${axisName(ax.screw, view)} at ${at({ [H]: fr(p), [V]: fr(q) })}`;
    out.push(titled(t, `<circle cx="${f(x)}" cy="${f(y)}" r="12" class="sgdHit"/>` + lens(p, q, ax.screw)));
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

  // planes parallel to the drawing: corner symbols below-left of the origin
  el.planesPar.sort((a, b) => _sgdHeightLabel(a.heights).localeCompare(_sgdHeightLabel(b.heights)));
  el.planesPar.forEach((pl, i) => {
    // the "L" follows the cell edges (H to the right, V along the tilted
    // edge); the glide arrow starts in its corner: along H, along V, or along
    // the diagonal H + V for an n-glide
    const bx = Ox - 84 + i * 54, by = Oy + 84;
    const s = 26;
    const mirror = pl.kind === 'm';
    const vEnd = [bx + uV[0] * s, by + uV[1] * s];
    {
      const letters = pl.kind === 'h' ? [names[H]] : pl.kind === 'v' ? [names[V]] : pl.kind === 'n' ? [names[H], names[V]].sort() : [];
      const what = mirror ? 'mirror plane' : pl.kind === 'd' ? 'd-glide plane' : `${letters.length === 1 ? letters[0] : 'n'}-glide plane`;
      const t = `${what} ⊥ ${names[view]} (parallel to the drawing) at ${X[view]} = ${hList(pl.heights)}${letters.length ? `, glide ${glideText(letters)}` : ''}`;
      out.push(`<g class="sgdEl"><title>${t}</title><rect x="${bx - 6}" y="${by - s - 6}" width="${s + 26}" height="${s + 12}" class="sgdHit"/>`);
    }
    out.push(`<path d="M ${f(vEnd[0])} ${f(vEnd[1])} L ${bx} ${by} L ${bx + s} ${by}" class="sgdCorner${mirror ? ' mirror' : ''}"/>`);
    if (!mirror) {
      let dir;
      if (pl.kind === 'h') dir = uH;
      else if (pl.kind === 'v') dir = uV;
      else { const d = [uH[0] + uV[0], uH[1] + uV[1]], l = Math.hypot(d[0], d[1]); dir = [d[0] / l, d[1] / l]; }
      const len = pl.kind === 'h' || pl.kind === 'v' ? 34 : 28;
      const from = [bx, by];
      const to = [from[0] + dir[0] * len, from[1] + dir[1] * len];
      const nx = -dir[1], ny = dir[0], back = [to[0] - dir[0] * 7, to[1] - dir[1] * 7];
      out.push(line(from, back, 'sgdCornerArrow'));
      out.push(`<polygon points="${f(to[0])},${f(to[1])} ${f(back[0] + nx * 3.5)},${f(back[1] + ny * 3.5)} ${f(back[0] - nx * 3.5)},${f(back[1] - ny * 3.5)}" class="sgdArrowHead"/>`);
    }
    const lab = _sgdHeightLabel(pl.heights);
    if (lab) out.push(`<text x="${bx + 38}" y="${by + 5}" class="sgdSmall">${lab}</text>`);
    out.push('</g>');
  });

  // general positions (tooltip: which line of the list, and its coordinates)
  const R = 9.5;
  const list = sgdCoordinateList(ops);
  const tip = it => `(${list.indexOf(it.op)}) ${list.format(it.op)}`;
  for (const g of sgdPositions(ops, view)) {
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
      g.list.slice(0, 2).forEach((it, i) => {
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
  const aBase = [mid[0] - 100, mid[1] + 12];
  out.push(smallArrow(aBase, uV, 22));
  out.push(`<text x="${f(aBase[0] + 8)}" y="${f(aBase[1] - 4)}" class="sgdAxisName">${names[V]}</text>`);
  const bm = P(0.5, 0);
  out.push(smallArrow([bm[0] - 12, Oy + 92], uH, 24));
  out.push(`<text x="${f(bm[0] - 3)}" y="${f(Oy + 86)}" class="sgdAxisName">${names[H]}</text>`);
  if (ang !== 90) {
    const r = 15;
    const a0 = 0, a1 = -rad; // from H towards V (screen y up = negative)
    const p0 = [Ox + r * Math.cos(a0), Oy + r * Math.sin(a0)], p1 = [Ox + r * Math.cos(a1), Oy + r * Math.sin(a1)];
    out.push(`<path d="M ${f(p0[0])} ${f(p0[1])} A ${r} ${r} 0 0 0 ${f(p1[0])} ${f(p1[1])}" class="sgdAngle"/>`);
    const bis = -rad / 2;
    out.push(`<text x="${f(Ox + 27 * Math.cos(bis) - 5)}" y="${f(Oy + 27 * Math.sin(bis) + 5)}" class="sgdAngleName">${angleNames[view]}</text>`);
  }

  const setting = entry.qualifier ? ` <span class="sgdSetting">setting ${entry.qualifier}</span>` : '';
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
