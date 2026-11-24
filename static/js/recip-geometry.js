// recip-geometry.js
// Shared geometry for the zone plots (absences-viz.js, laue-viz.js):
// draws a reciprocal-lattice plane with the TRUE metric (a*, b*, c* and
// the starred angles) instead of an index-square grid, so oblique cells
// (monoclinic beta, hexagonal 120°) look the way a precession photo does.
//
// Without a usable cell the basis falls back to a unit square lattice,
// i.e. the old index-grid look.

const _toRad = d => d * Math.PI / 180;

// uIdx/vIdx: which Miller index runs along the plot's x / y lattice
// direction (0 = h, 1 = k, 2 = l). Returns lengths (1/Å) of the two
// reciprocal axes and the angle between them.
function reciprocalPlaneBasis(cell, uIdx, vIdx) {
  const fallback = { lu: 1, lv: 1, cosT: 0, sinT: 1, metric: false };
  if (!cell) return fallback;
  const keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'];
  if (!keys.every(k => Number.isFinite(cell[k]) && cell[k] > 0)) return fallback;

  const [ca, cb, cg] = [cell.alpha, cell.beta, cell.gamma].map(d => Math.cos(_toRad(d)));
  const [sa, sb, sg] = [cell.alpha, cell.beta, cell.gamma].map(d => Math.sin(_toRad(d)));
  const rad = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg;
  if (!(rad > 1e-9) || !(sa > 1e-6 && sb > 1e-6 && sg > 1e-6)) return fallback;
  const V = cell.a * cell.b * cell.c * Math.sqrt(rad);

  const L = [cell.b * cell.c * sa / V, cell.a * cell.c * sb / V, cell.a * cell.b * sg / V];
  // cos of alpha*, beta*, gamma*
  const cosStar = [(cb * cg - ca) / (sb * sg), (ca * cg - cb) / (sa * sg), (ca * cb - cg) / (sa * sb)];
  // angle between axis i* and j* is the starred angle of the THIRD axis
  const cosT = cosStar[3 - uIdx - vIdx];
  return { lu: L[uIdx], lv: L[vIdx], cosT, sinT: Math.sqrt(Math.max(0, 1 - cosT * cosT)), metric: true };
}

// Same, for two arbitrary reciprocal-lattice directions given as index
// vectors, e.g. u = [1, 1, 0] (a* + b*) and v = [0, 0, 1] (c*) for the
// diagonal plane hhl. Lengths and angle from the reciprocal metric tensor
// G* = G⁻¹ (G = direct metric tensor).
function reciprocalBasisForVectors(cell, uVec, vVec) {
  const fallback = { lu: 1, lv: 1, cosT: 0, sinT: 1, metric: false };
  if (!cell) return fallback;
  const keys = ['a', 'b', 'c', 'alpha', 'beta', 'gamma'];
  if (!keys.every(k => Number.isFinite(cell[k]) && cell[k] > 0)) return fallback;
  const { a, b, c } = cell;
  const [ca, cb, cg] = [cell.alpha, cell.beta, cell.gamma].map(d => Math.cos(_toRad(d)));
  const G = [[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]];
  const det = G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1]) - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0]) + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]);
  if (!(det > 1e-9)) return fallback;
  const Gs = [0, 1, 2].map(i => [0, 1, 2].map(j => {
    const r = [0, 1, 2].filter(x => x !== j), cIdx = [0, 1, 2].filter(x => x !== i);
    const m = G[r[0]][cIdx[0]] * G[r[1]][cIdx[1]] - G[r[0]][cIdx[1]] * G[r[1]][cIdx[0]];
    return ((i + j) % 2 ? -m : m) / det;
  }));
  const dot = (x, y) => x.reduce((s, xi, i) => s + xi * y.reduce((t, yj, j) => t + Gs[i][j] * yj, 0), 0);
  const lu = Math.sqrt(dot(uVec, uVec)), lv = Math.sqrt(dot(vVec, vVec));
  const cosT = dot(uVec, vVec) / (lu * lv);
  return { lu, lv, cosT, sinT: Math.sqrt(Math.max(0, 1 - cosT * cosT)), metric: true };
}

// Maps lattice coordinates (u, v) to canvas pixels, equal aspect, fitted
// and centered into `box` {x, y, w, h}. Returns px(u,v), py(u,v) and the
// pixel distance between neighbouring lattice rows (for marker sizing).
function makeLatticeFrame(basis, uMin, uMax, vMin, vMax, box) {
  const pos = (u, v) => [u * basis.lu + v * basis.lv * basis.cosT, v * basis.lv * basis.sinT];
  const xs = [], ys = [];
  for (const u of [uMin, uMax]) for (const v of [vMin, vMax]) { const p = pos(u, v); xs.push(p[0]); ys.push(p[1]); }
  const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
  const scale = Math.min(box.w / (x1 - x0), box.h / (y1 - y0));
  const offX = box.x + (box.w - (x1 - x0) * scale) / 2;
  const offY = box.y + (box.h - (y1 - y0) * scale) / 2;
  return {
    px: (u, v) => offX + (pos(u, v)[0] - x0) * scale,
    py: (u, v) => offY + (y1 - pos(u, v)[1]) * scale,
    spacing: scale * Math.min(basis.lu, basis.lv * basis.sinT),
  };
}

// Lattice lines (u = const, v = const), emphasized zero axes, tick labels
// at the line ends, and axis captions.
function drawLatticeGrid(ctx, colors, frame, r, labels, tickFn, cssW, cssH) {
  const { px, py } = frame;
  ctx.lineWidth = 1;
  ctx.font = '9px sans-serif';
  const uTicks = tickFn(r.uMin, r.uMax), vTicks = tickFn(r.vMin, r.vMax);
  for (const t of uTicks) {
    ctx.strokeStyle = t === 0 ? colors.axis : colors.grid;
    ctx.beginPath(); ctx.moveTo(px(t, r.vMin), py(t, r.vMin)); ctx.lineTo(px(t, r.vMax), py(t, r.vMax)); ctx.stroke();
    ctx.fillStyle = colors.text; ctx.textAlign = 'center';
    ctx.fillText(String(t), px(t, r.vMin), py(t, r.vMin) + 12);
  }
  for (const t of vTicks) {
    ctx.strokeStyle = t === 0 ? colors.axis : colors.grid;
    ctx.beginPath(); ctx.moveTo(px(r.uMin, t), py(r.uMin, t)); ctx.lineTo(px(r.uMax, t), py(r.uMax, t)); ctx.stroke();
    ctx.fillStyle = colors.text; ctx.textAlign = 'right';
    ctx.fillText(String(t), px(r.uMin, t) - 5, py(r.uMin, t) + 3);
  }
  ctx.fillStyle = colors.text;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'right';
  ctx.fillText(labels.x, cssW - 6, cssH - 6);
  ctx.textAlign = 'left';
  ctx.fillText(labels.y, 6, 14);
}
