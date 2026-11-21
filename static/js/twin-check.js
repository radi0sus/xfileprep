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

// reflections: [{h,k,l,f2,sigma}] (lattice-absent ones may stay — they are
// weak in both members of a pair and drop out of the sum selection).
// Returns { applicable, reason?, laue, holohedry, laws: [...] }.
function checkMerohedralTwinning(reflections, cell, laueClass, monoclinicAxis = 'b') {
  const classification = classifyCellMetric(cell);
  const holoName = classification.laueClasses[0];
  if (classification.system.includes('rhombohedral')) {
    return { applicable: false, laue: laueClass, holohedry: holoName, reason: 'rhombohedral axes — twin laws are only derived for hexagonal axes' };
  }
  if (holoName === laueClass) {
    return { applicable: false, laue: laueClass, holohedry: holoName, reason: `the Laue class ${laueClass} already has the full symmetry of the cell metric — no (pseudo-)merohedral twin law` };
  }
  const G = generateGroup(getLaueClassGenerators(laueClass, monoclinicAxis));
  const H = generateGroup(getLaueClassGenerators(holoName, monoclinicAxis));
  const inH = new Set(H.map(twinMatrixKey));
  if (!G.every(M => inH.has(twinMatrixKey(M)))) {
    return { applicable: false, laue: laueClass, holohedry: holoName, reason: `${laueClass} (unique axis ${monoclinicAxis}) is not a subgroup of ${holoName} in this orientation` };
  }

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

  const laws = twinLawCosets(H, G).map(T => {
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
      deviation: twinMetricDeviation(cell, T),
      pairs: n, selection, meanH, alpha, alphaSe, verdict,
    };
  });
  return { applicable: true, laue: laueClass, holohedry: holoName, laws };
}

// Twin law in the transformed setting (HKLF matrix S: h_new = S h_old):
// R_new = S R S^-1. Used when the .ins changes setting.
function twinMatrixInSetting(R, S) {
  if (!S) return R;
  const Sinv = invert3x3(S);
  return mul3i(mul3i(S, R), Sinv.map(r => r.map(v => Math.round(v))));
}
