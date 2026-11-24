// Merohedral / pseudo-merohedral twinning check (H test, Yeates 1988).
//
// When the cell METRIC has higher symmetry than the Laue class found from the
// intensities (e.g. monoclinic with beta = 90.07 in an orthorhombic-looking
// cell, or Laue 4/m in a tetragonal metric that also allows 4/mmm), the
// lattice can be twinned by the operations of the metric holohedry that are
// not in the Laue group. Their cosets give the possible twin laws.
//
// H test: for each pair of acentric reflections h and h' = twin(h), both
// measured, H = |I1 - I2| / (I1 + I2). For an untwinned crystal (independent
// Wilson-distributed I1, I2) <H> = 1/2; for twin fraction alpha
// <H> = 1/2 - alpha, and alpha = 1/2 (perfect twin) gives H = 0.
//
// Pair selection by the SUM: I1 + I2 > 10 sigma(I1 + I2) (3 sigma if that
// leaves fewer than 200 pairs). Given the sum, the
// share I1/(I1+I2) of two independent exponentials is uniform, so this keeps
// <H> unbiased; requiring BOTH reflections strong does not — a simulation of
// untwinned data gave alpha = 0.21 that way (alpha = 0.00 with the sum).
// Known limitation: pseudo-symmetry (e.g. Z' > 1 with near-symmetric
// molecules) makes twin-related intensities similar by itself, so alpha is
// overestimated (yw6, P 21 with beta = 90.07, Z = 8: refined 0.29, H test
// 0.35). Noise works the other way: a crystal with the full higher symmetry
// should give alpha = 0.5, but noisy weak pairs pull it down (0.38 at 3 sigma),
// hence the 10 sigma selection and "near-perfect" from alpha >= 0.40. Use it as a detector and BASF start value, not as a refined number.
//
// Racemic (inversion) twinning cannot be seen in intensities without
// anomalous signal; it is only mentioned for acentric space groups.

const TWIN_MIN_PAIRS = 50;
const TWIN_SIGNIFICANT_ALPHA = 0.05;

function twinMatrixKey(M) {
  return M.map(r => r.join(',')).join(';');
}

function applyRow(h, M) { // h' = h M (the convention of canonicalHKL / generateGroup)
  return [0, 1, 2].map(j => h[0] * M[0][j] + h[1] * M[1][j] + h[2] * M[2][j]);
}

function det3i(M) {
  return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
    - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
    + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
}

function mul3i(A, B) {
  return A.map((r, i) => [0, 1, 2].map(j => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
}

function transpose3i(M) {
  return [0, 1, 2].map(i => [0, 1, 2].map(j => M[j][i]));
}

// Coset representatives of `holohedry` modulo `laueGroup` (right cosets
// G·T: every element of a coset relates the same pairs of reflection
// classes). Prefers a proper 2-fold rotation as representative — the usual
// way to write a twin law (TWIN -1 0 0 0 -1 0 0 0 1).
function twinLawCosets(holohedry, laueGroup) {
  const inG = new Set(laueGroup.map(twinMatrixKey));
  const covered = new Set();
  const laws = [];
  const isTwoFold = M => det3i(M) === 1 && twinMatrixKey(mul3i(M, M)) === '1,0,0;0,1,0;0,0,1';
  for (const T of holohedry) {
    if (inG.has(twinMatrixKey(T)) || covered.has(twinMatrixKey(T))) continue;
    const coset = laueGroup.map(G => mul3i(G, T));
    coset.forEach(M => covered.add(twinMatrixKey(M)));
    const twoFolds = coset.filter(isTwoFold);
    // among 2-folds prefer one along a cell axis (a diagonal matrix), c first
    // (TWIN -1 0 0 0 -1 0 0 0 1, as usually written), then the first
    const diagonal = twoFolds.filter(M => M[0][1] === 0 && M[0][2] === 0 && M[1][0] === 0 && M[1][2] === 0 && M[2][0] === 0 && M[2][1] === 0);
    const rep = diagonal.find(M => M[2][2] === 1) || diagonal[0]
      || twoFolds[0] || coset.find(M => det3i(M) === 1) || coset[0];
    laws.push(rep);
  }
  return laws;
}

// Largest change of a cell parameter under the twin operation: how far
// the lattice is from really having this symmetry (0 = exact merohedry).
function twinMetricDeviation(cell, M) {
  if (typeof transformCellParameters !== 'function') return null;
  let best = null;
  for (const X of [M, transpose3i(M)]) {
    const c2 = transformCellParameters(cell, X);
    const dAng = Math.max(Math.abs(c2.alpha - cell.alpha), Math.abs(c2.beta - cell.beta), Math.abs(c2.gamma - cell.gamma));
    const dLen = Math.max(Math.abs(c2.a - cell.a) / cell.a, Math.abs(c2.b - cell.b) / cell.b, Math.abs(c2.c - cell.c) / cell.c);
    if (!best || dAng + 100 * dLen < best.angle + 100 * best.length) best = { angle: dAng, length: dLen };
  }
  return best;
}

// --- Lattice symmetry independent of the cell choice (Le Page 1982) -------
// A lattice 2-fold axis is a direct lattice vector t (nearly) parallel to a
// reciprocal lattice vector h; the angle between them (delta) is the
// obliquity. Its matrix in the cell basis is R = 2 t h^T / (t.h) - I, which
// must be integral and map the centring vectors onto lattice vectors. The
// cell's own angles miss such pseudo-symmetry whenever the conventional cell
// is not the one that shows it — nm56_25 (C 2/c, beta = 132.6 deg) has a
// 2-fold along [001], perpendicular to (-1 0 1) within 0.012 deg: exactly its
// refined TWIN -1 0 -2 0 -1 0 0 0 1.
const CENTERING_VECTORS = {
  'P (primitive)': [], A: [[0, 0.5, 0.5]], B: [[0.5, 0, 0.5]], C: [[0.5, 0.5, 0]], I: [[0.5, 0.5, 0.5]],
  F: [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
  'R (obverse)': [[2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]], 'R (reverse)': [[1 / 3, 2 / 3, 1 / 3], [2 / 3, 1 / 3, 2 / 3]],
};
const LEPAGE_MAX_DELTA = 3; // degrees, as in Le Page's program and PLATON

function lePageTwoFolds(cell, centering = 'P (primitive)', maxDelta = LEPAGE_MAX_DELTA, N = 2) {
  const d2r = Math.PI / 180;
  const { a, b, c } = cell;
  const al = cell.alpha * d2r, be = cell.beta * d2r, ga = cell.gamma * d2r;
  const cx = c * Math.cos(be), cy = c * (Math.cos(al) - Math.cos(be) * Math.cos(ga)) / Math.sin(ga);
  const A = [[a, 0, 0], [b * Math.cos(ga), b * Math.sin(ga), 0], [cx, cy, Math.sqrt(Math.max(0, c * c - cx * cx - cy * cy))]];
  const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const V = dot(A[0], cross(A[1], A[2]));
  const Rr = [cross(A[1], A[2]), cross(A[2], A[0]), cross(A[0], A[1])].map(v => v.map(x => x / V));
  const cart = (t, B) => [0, 1, 2].map(j => t[0] * B[0][j] + t[1] * B[1][j] + t[2] * B[2][j]);
  const cvecs = CENTERING_VECTORS[centering] || [];
  const isLattice = v => [[0, 0, 0], ...cvecs].some(cv => v.every((x, i) => Math.abs(x - cv[i] - Math.round(x - cv[i])) < 1e-6));
  const ts = [];
  for (let u = -N; u <= N; u++) for (let v = -N; v <= N; v++) for (let w = -N; w <= N; w++) {
    if (!u && !v && !w) continue;
    ts.push([u, v, w]);
    for (const cv of cvecs) ts.push([u + cv[0], v + cv[1], w + cv[2]]);
  }
  const found = new Map();
  for (const t of ts) {
    const tc = cart(t, A);
    const tl = Math.sqrt(dot(tc, tc));
    for (let h = -N; h <= N; h++) for (let k = -N; k <= N; k++) for (let l = -N; l <= N; l++) {
      if (!h && !k && !l) continue;
      const th = t[0] * h + t[1] * k + t[2] * l;
      if (th <= 1e-9) continue;
      const hc = cart([h, k, l], Rr);
      const cos = dot(tc, hc) / (tl * Math.sqrt(dot(hc, hc)));
      const delta = Math.acos(Math.min(1, cos)) / d2r;
      if (delta > maxDelta) continue;
      const hh = [h, k, l];
      const Rf = [0, 1, 2].map(i => [0, 1, 2].map(j => 2 * t[i] * hh[j] / th - (i === j ? 1 : 0)));
      if (!Rf.flat().every(x => Math.abs(x - Math.round(x)) < 1e-6)) continue;
      const R = Rf.map(r => r.map(Math.round));
      if (!cvecs.every(cv => isLattice([0, 1, 2].map(i => R[i][0] * cv[0] + R[i][1] * cv[1] + R[i][2] * cv[2])))) continue;
      const key = twinMatrixKey(R);
      if (!found.has(key) || found.get(key).delta > delta) found.set(key, { R, delta, t, h: hh });
    }
  }
  return [...found.values()].sort((x, y) => x.delta - y.delta);
}

// Closure of a set of integer matrices (h' = h M); null if it exceeds 48
// elements (inconsistent approximate axes).
function closeMatrixGroup(gens) {
  const out = new Map([['1,0,0;0,1,0;0,0,1', [[1, 0, 0], [0, 1, 0], [0, 0, 1]]]]);
  for (const g of gens) out.set(twinMatrixKey(g), g);
  let grown = true;
  while (grown) {
    grown = false;
    for (const x of [...out.values()]) for (const y of [...out.values()]) {
      const z = mul3i(x, y);
      const k = twinMatrixKey(z);
      if (!out.has(k)) { out.set(k, z); grown = true; if (out.size > 48) return null; }
    }
  }
  return [...out.values()];
}

// Lattice holohedry from the Le Page 2-folds (smallest obliquity first, a
// 2-fold is dropped if it would make the group inconsistent), plus the
// inversion and the Laue group itself.
function latticeHolohedry(cell, centering, laueGroup) {
  const folds = lePageTwoFolds(cell, centering);
  // For a 2-fold, x' = R x on coordinates and h' = h R^-1 = h R on indices: R serves directly.
  let gens = [[[-1, 0, 0], [0, -1, 0], [0, 0, -1]], ...laueGroup];
  let H = closeMatrixGroup(gens);
  const used = [];
  for (const f of folds) {
    const trial = closeMatrixGroup([...gens, f.R]);
    if (!trial) continue;
    gens = [...gens, f.R];
    H = trial;
    used.push(f);
  }
  return { H, folds: used };
}

function isDiagonal3(M) {
  return M[0][1] === 0 && M[0][2] === 0 && M[1][0] === 0 && M[1][2] === 0 && M[2][0] === 0 && M[2][1] === 0;
}

// Laue class of a holohedry group from its order and the highest order of a
// proper rotation in it (4/m vs mmm, 6/m vs -3m, m-3 vs 6/mmm).
function laueNameOfGroup(H) {
  const id = '1,0,0;0,1,0;0,0,1';
  const orderOf = M => { let X = M; for (let n = 1; n <= 6; n++) { if (twinMatrixKey(X) === id) return n; X = mul3i(X, M); } return 0; };
  const maxProper = Math.max(...H.filter(M => det3i(M) === 1).map(orderOf));
  const n = H.length;
  if (n === 2) return '-1';
  if (n === 4) return '2/m';
  if (n === 6) return '-3';
  if (n === 8) return maxProper === 4 ? '4/m' : 'mmm';
  if (n === 12) return maxProper === 6 ? '6/m' : '-3m';
  if (n === 16) return '4/mmm';
  if (n === 24) return maxProper === 6 ? '6/mmm' : 'm-3';
  if (n === 48) return 'm-3m';
  return 'order ' + n;
}

// reflections: [{h,k,l,f2,sigma}] (lattice-absent ones may stay — they are
// weak in both members of a pair and drop out of the sum selection).
// centering: observed lattice centring name (absences.js), for the Le Page search.
// Returns { applicable, reason?, laue, holohedry, laws: [...] }.
function checkMerohedralTwinning(reflections, cell, laueClass, monoclinicAxis = 'b', centering = 'P (primitive)') {
  const G = generateGroup(getLaueClassGenerators(laueClass, monoclinicAxis));
  const { H, folds } = latticeHolohedry(cell, centering, G);
  const inLaue = new Set(G.map(twinMatrixKey));
  const allFolds = lePageTwoFolds(cell, centering).map(f => ({ ...f, inLaue: inLaue.has(twinMatrixKey(f.R)), used: folds.includes(f) || folds.some(u => twinMatrixKey(u.R) === twinMatrixKey(f.R)) }));
  const maxFold = folds.length ? Math.max(...folds.map(f => f.delta)) : 0;
  const holoName = laueNameOfGroup(H);
  const holoNote = `Le Page search, obliquity δ ≤ ${LEPAGE_MAX_DELTA}°${folds.length ? `, largest used ${maxFold.toFixed(2)}°` : ''}`;
  if (H.length === G.length) {
    return { applicable: false, laue: laueClass, holohedry: holoName, holohedryNote: holoNote, folds: allFolds,
      reason: `the lattice shows no symmetry beyond the Laue class ${typeof formatPG === 'function' ? formatPG(laueClass) : laueClass} (Le Page search for 2-fold axes, obliquity ≤ ${LEPAGE_MAX_DELTA}°) — no (pseudo-)merohedral twin law` };
  }
  const foldByKey = new Map(folds.map(f => [twinMatrixKey(f.R), f]));

  // merge under the Laue group (weighted mean)
  const merged = new Map();
  for (const r of reflections) {
    if (!(r.sigma > 0)) continue;
    const key = canonicalHKL(r.h, r.k, r.l, G).join(',');
    let e = merged.get(key);
    if (!e) { e = { sw: 0, swi: 0 }; merged.set(key, e); }
    const w = 1 / (r.sigma * r.sigma);
    e.sw += w; e.swi += w * r.f2;
  }
  for (const e of merged.values()) { e.I = e.swi / e.sw; e.s = Math.sqrt(1 / e.sw); }

  // centric (in any group with this Laue class): h g = -h for some g other than -1
  const nonInversion = G.filter(M => twinMatrixKey(M) !== '-1,0,0;0,-1,0;0,0,-1');
  const isCentric = h => nonInversion.some(M => { const v = applyRow(h, M); return v[0] === -h[0] && v[1] === -h[1] && v[2] === -h[2]; });

  // coset representative: the Le Page 2-fold with the smallest obliquity, if the coset has one
  const reps = twinLawCosets(H, G).map(T => {
    const coset = G.map(g => mul3i(g, T));
    // smallest obliquity; on a tie (2-folds related by a Laue operation have the
    // same delta) a 2-fold along c, then along any cell axis — as usually written
    const pref = f => (isDiagonal3(f.R) ? (f.R[2][2] === 1 ? 0 : 1) : 2);
    const best = coset.map(M => foldByKey.get(twinMatrixKey(M))).filter(Boolean)
      .sort((x, y) => (Math.abs(x.delta - y.delta) > 1e-3 ? x.delta - y.delta : pref(x) - pref(y)))[0];
    return best ? { T: best.R, fold: best } : { T, fold: null };
  });
  const laws = reps.map(({ T, fold }) => {
    const Hs = [];
    const seen = new Set();
    for (const [key, e] of merged) {
      const h = key.split(',').map(Number);
      if (isCentric(h)) continue;
      const t = applyRow(h, T);
      const key2 = canonicalHKL(t[0], t[1], t[2], G).join(',');
      if (key2 === key) continue; // twin-invariant class (e.g. on the twin axis): no information
      const pairKey = key < key2 ? key + '|' + key2 : key2 + '|' + key;
      if (seen.has(pairKey)) continue;
      seen.add(pairKey);
      const e2 = merged.get(key2);
      if (!e2) continue;
      const sum = e.I + e2.I;
      const snr = sum / Math.hypot(e.s, e2.s);
      if (!(snr > 3)) continue;
      Hs.push({ H: Math.abs(e.I - e2.I) / sum, strong: snr > 10 });
    }
    const strong = Hs.filter(x => x.strong);
    const used = (strong.length >= 200 ? strong : Hs).map(x => x.H);
    const selection = strong.length >= 200 ? 10 : 3;
    const n = used.length;
    const meanH = n ? used.reduce((a, b) => a + b, 0) / n : null;
    const alpha = meanH === null ? null : Math.max(0, Math.min(0.5, 0.5 - meanH));
    // standard error of alpha (= of <H>)
    const sdH = n > 1 ? Math.sqrt(used.reduce((a, b) => a + (b - meanH) * (b - meanH), 0) / (n - 1)) : null;
    const alphaSe = sdH !== null ? sdH / Math.sqrt(n) : null;
    let verdict;
    if (n < TWIN_MIN_PAIRS) verdict = 'too few pairs';
    else if (alpha >= 0.40) verdict = 'perfect or near-perfect twin — or the true Laue class is higher';
    else if (alpha >= TWIN_SIGNIFICANT_ALPHA && alpha > 3 * alphaSe) verdict = 'twinned';
    else verdict = 'no indication';
    return {
      matrix: T, // index operation h' = h M
      shelxMatrix: transpose3i(T), // SHELX TWIN / HKLF convention: h' = R h (column) — R = M^T
      deviation: fold ? { angle: fold.delta, length: 0, lePage: true } : twinMetricDeviation(cell, T),
      pairs: n, selection, meanH, alpha, alphaSe, verdict,
      hValues: used, // for the cumulative N(H) plot
    };
  });
  return { applicable: true, laue: laueClass, holohedry: holoName, holohedryNote: holoNote, folds: allFolds, laws };
}

// Twin law in the transformed setting (HKLF matrix S: h_new = S h_old):
// R_new = S R S^-1. Used when the .ins changes setting.
function twinMatrixInSetting(R, S) {
  if (!S) return R;
  const Sinv = invert3x3(S);
  return mul3i(mul3i(S, R), Sinv.map(r => r.map(v => Math.round(v))));
}

// --- Non-merohedral twins (partial overlap) ---------------------------------
// A lattice 2-fold with an obliquity of a few degrees (Le Page delta ~1-6 deg)
// does not map the lattice onto itself: rotated reciprocal-lattice points
// land NEAR other points only at low angle, and only there do the spots of
// both domains overlap. Twin law = the exact 180 deg rotation about the
// direct row t or the reciprocal row h (both are tried; non-integral matrix).
// Test: the H test (as above) on the pairs whose twin image lies within
// NONMERO_OVERLAP of a lattice point — overlapping spots share intensity.
// Pseudo-symmetry gives correlation at any distance, a twin only where the
// spots overlap; the H test stays near 0 for untwinned pseudo-symmetric data
// (COD 1561021: alpha -0.11; COD 1543504: -0.17). pf12c_26_o (P -1, PLATON
// TwinRotMat: 2-fold about (1 -1 0)*, 3.40 deg, BASF 0.17): alpha 0.26-0.29 on
// the closest pairs, falling with distance. Refinement needs an HKLF 5 file.
const NONMERO_MIN_DELTA = 1;      // below: the merohedral H test above applies
const NONMERO_MAX_DELTA = 6;      // degrees
const NONMERO_OVERLAP = 0.01;     // 1/Angstrom, distance of the twin image from a lattice point

function checkNonMerohedralTwinning(reflections, cell, laueClass, monoclinicAxis = 'b', centering = 'P (primitive)') {
  const G = generateGroup(getLaueClassGenerators(laueClass, monoclinicAxis));
  const inG = new Set(G.map(twinMatrixKey));
  const folds = lePageTwoFolds(cell, centering, NONMERO_MAX_DELTA)
    .filter(f => f.delta > NONMERO_MIN_DELTA && !inG.has(twinMatrixKey(f.R)));
  if (!folds.length) return { applicable: false, reason: `no lattice 2-fold with an obliquity of ${NONMERO_MIN_DELTA}–${NONMERO_MAX_DELTA}° outside the Laue class` };

  // Cartesian bases (rows): direct A, reciprocal Rr
  const d2r = Math.PI / 180;
  const { a, b, c } = cell;
  const al = cell.alpha * d2r, be = cell.beta * d2r, ga = cell.gamma * d2r;
  const cx = c * Math.cos(be), cy = c * (Math.cos(al) - Math.cos(be) * Math.cos(ga)) / Math.sin(ga);
  const A = [[a, 0, 0], [b * Math.cos(ga), b * Math.sin(ga), 0], [cx, cy, Math.sqrt(Math.max(0, c * c - cx * cx - cy * cy))]];
  const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const V = dot(A[0], cross(A[1], A[2]));
  const Rr = [cross(A[1], A[2]), cross(A[2], A[0]), cross(A[0], A[1])].map(v => v.map(x => x / V));
  const cart = (t, B) => [0, 1, 2].map(j => t[0] * B[0][j] + t[1] * B[1][j] + t[2] * B[2][j]);

  // merge under the Laue group
  const merged = new Map();
  for (const r of reflections) {
    if (!(r.sigma > 0)) continue;
    const h = canonicalHKL(r.h, r.k, r.l, G);
    const key = h.join(',');
    let e = merged.get(key);
    if (!e) { e = { h, sw: 0, swi: 0 }; merged.set(key, e); }
    const w = 1 / (r.sigma * r.sigma);
    e.sw += w; e.swi += w * r.f2;
  }
  for (const e of merged.values()) { e.I = e.swi / e.sw; e.s = Math.sqrt(1 / e.sw); }

  const laws = [];
  for (const f of folds) {
    for (const about of ['reciprocal', 'direct']) {
      const u0 = about === 'direct' ? cart(f.t, A) : cart(f.h, Rr);
      const un = Math.sqrt(dot(u0, u0));
      const u = u0.map(x => x / un);
      // index matrix of the exact rotation: h' = M h (columns = images of unit indices)
      const img = hk => { const hc = cart(hk, Rr); const p = 2 * dot(hc, u); const hp = [0, 1, 2].map(i => p * u[i] - hc[i]); return { hp, fi: [0, 1, 2].map(j => dot(hp, A[j])) }; };
      const cols = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(e => img(e).fi);
      const matrix = [0, 1, 2].map(i => [0, 1, 2].map(j => cols[j][i]));
      let overlaps = 0;
      const Hs = [];
      const seen = new Set();
      for (const [key, e] of merged) {
        const { hp, fi } = img(e.h);
        const g = fi.map(Math.round);
        const gc = cart(g, Rr);
        const dist = Math.hypot(hp[0] - gc[0], hp[1] - gc[1], hp[2] - gc[2]);
        if (dist > NONMERO_OVERLAP) continue;
        overlaps++;
        const key2 = canonicalHKL(g[0], g[1], g[2], G).join(',');
        if (key2 === key) continue;
        const pk = key < key2 ? key + '|' + key2 : key2 + '|' + key;
        if (seen.has(pk)) continue;
        seen.add(pk);
        const e2 = merged.get(key2);
        if (!e2) continue;
        const sum = e.I + e2.I, snr = sum / Math.hypot(e.s, e2.s);
        if (snr > 3) Hs.push({ H: Math.abs(e.I - e2.I) / sum, strong: snr > 10 });
      }
      const strong = Hs.filter(x => x.strong);
      const used = (strong.length >= 200 ? strong : Hs).map(x => x.H);
      const n = used.length;
      const meanH = n ? used.reduce((s, x) => s + x, 0) / n : null;
      const alphaRaw = meanH === null ? null : 0.5 - meanH;
      const sd = n > 1 ? Math.sqrt(used.reduce((s, x) => s + (x - meanH) ** 2, 0) / (n - 1)) : null;
      const alphaSe = sd !== null ? sd / Math.sqrt(n) : null;
      const alpha = alphaRaw === null ? null : Math.max(0, alphaRaw);
      let verdict;
      if (n < 100) verdict = 'too few overlapping pairs';
      else if (alpha >= TWIN_SIGNIFICANT_ALPHA && alpha > 3 * alphaSe) verdict = 'twinned (non-merohedral)';
      else verdict = 'no indication';
      laws.push({ axis: f.t, plane: f.h, about, delta: f.delta, matrix, overlaps, pairs: n, meanH, alpha, alphaSe, verdict, hValues: used });
    }
  }
  laws.sort((x, y) => (y.verdict.startsWith('twinned') - x.verdict.startsWith('twinned')) || ((y.alpha || 0) - (x.alpha || 0)));
  return { applicable: true, laws, overlapTolerance: NONMERO_OVERLAP };
}
