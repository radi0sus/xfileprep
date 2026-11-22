// Anomalous signal from Bijvoet pairs: is the structure acentric at all, and
// is it an inversion (racemic) twin?
//
// Merged under the point group of an ACENTRIC candidate (Friedel mates kept
// apart), every class h is paired with the class of -h. Without anomalous
// scattering both are equal; differences beyond the noise exist only if
// (a) the structure is acentric and contains anomalous scatterers, and
// (b) it is not an inversion twin near 50:50.
//
// 1. Observed signal: chi2 = <(I+ - I-)^2 / (s+^2 + s-^2)>, 1 for noise only
//    (sigmas first rescaled by the scatter of repeated measurements, see below).
//    A clearly significant excess means Friedel's law is broken — which no
//    centrosymmetric structure can do (yw6: E² says "centrosymmetric",
//    Bijvoet chi2 = 17 says acentric; refined in P 21). Absorption errors can
//    fake a small excess, hence the generous threshold.
// 2. Expected signal (needs f'' per element, read from the CIF's
//    _atom_type_scat_dispersion_imag, and the composition): for randomly
//    placed atoms <dI^2>/<I^2> = 4 sum f''^2 / sum f0(s)^2 (acentric, from
//    I+ - I- = -4 sum_jk f0_j f''_k sin(phi_j - phi_k)). The noise-corrected
//    observed/expected ratio estimates |1 - 2x| for an inversion twin with
//    fraction x — no sign, there is no structure model to tell the hand.
//    Checked: jg6 (Cu, Flack -0.006) -> 0.02, nm75_26 (Co, Flack 0.45) ->
//    0.40, COD 1501632 (Sr/Cu, Flack 0.52) -> 0.48. Combined with a
//    (pseudo-)merohedral twin the two effects mix (yw6: 0.27 vs. a refined
//    inversion share of 0.45) — then it is only a rough indication.

// f'' per element from the CIF atom-type loop, or null.
function extractDispersionFromCIF(text) {
  const lines = String(text).split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].trim() !== 'loop_') continue;
    const tags = [];
    let j = i + 1;
    while (j < lines.length && lines[j].trim().startsWith('_')) { tags.push(lines[j].trim().split(/\s+/)[0].toLowerCase()); j++; }
    const iSym = tags.indexOf('_atom_type_symbol');
    const iImag = tags.indexOf('_atom_type_scat_dispersion_imag');
    if (iSym < 0 || iImag < 0) continue;
    // values: tokens (quoted strings and ;-text fields as single tokens) until the next loop_/tag/data_
    const tokens = [];
    for (; j < lines.length; j++) {
      const raw = lines[j];
      const t = raw.trim();
      if (t === 'loop_' || t.startsWith('_') || t.startsWith('data_')) break;
      if (raw.startsWith(';')) { // text field
        let k = j + 1;
        while (k < lines.length && !lines[k].startsWith(';')) k++;
        tokens.push('');
        j = k;
        continue;
      }
      const re = /'([^']*)'|"([^"]*)"|(\S+)/g;
      let m;
      while ((m = re.exec(t))) tokens.push(m[1] ?? m[2] ?? m[3]);
    }
    const out = {};
    for (let r = 0; r + tags.length <= tokens.length; r += tags.length) {
      const el = String(tokens[r + iSym]).replace(/[^A-Za-z]/g, '');
      const v = parseFloat(tokens[r + iImag]);
      if (el && Number.isFinite(v)) out[el.charAt(0).toUpperCase() + el.slice(1).toLowerCase()] = v;
    }
    return Object.keys(out).length ? out : null;
  }
  return null;
}

// reflections [{h,k,l,f2,sigma}], cell, pointGroup: proper+improper rotations
// of an acentric candidate (h' = h R), fpp: {El: f''} or null, composition:
// [{element, count}] per unit cell or null.
function checkAnomalousSignal(reflections, cell, pointGroup, fpp, composition) {
  const canon = (h, k, l) => canonicalHKL(h, k, l, pointGroup).join(',');
  const merged = new Map();
  for (const r of reflections) {
    if (!(r.sigma > 0)) continue;
    const key = canon(r.h, r.k, r.l);
    let e = merged.get(key);
    if (!e) { e = { sw: 0, swi: 0, obs: [] }; merged.set(key, e); }
    const w = 1 / (r.sigma * r.sigma);
    e.sw += w; e.swi += w * r.f2; e.obs.push(r);
  }
  for (const e of merged.values()) { e.I = e.swi / e.sw; e.v = 1 / e.sw; }

  // sigma calibration: within one Bijvoet class all observations measure the
  // same intensity, so sum (I_i - I_mean)^2 / sigma_i^2 should be (n - 1).
  // Data-reduction sigmas are often too small (COD 1561021, centrosymmetric:
  // raw Bijvoet chi2 1.83 -> would look acentric); k rescales the variances.
  // Merged data (one observation per class) leave k = 1 — not calibrated.
  let ssum = 0, dof = 0;
  for (const e of merged.values()) {
    if (e.obs.length < 2) continue;
    for (const r of e.obs) ssum += (r.f2 - e.I) ** 2 / (r.sigma * r.sigma);
    dof += e.obs.length - 1;
  }
  const sigmaScale = dof >= 200 ? Math.max(1, ssum / dof) : 1;
  const calibrated = dof >= 200;
  for (const e of merged.values()) { e.v *= sigmaScale; delete e.obs; }

  const recip = reciprocalCellFromDirect(cell);
  const usable = composition && fpp ? composition.filter(c => c.count > 0) : null;
  const sumFpp2 = usable ? usable.reduce((a, c) => a + c.count * (fpp[c.element] || 0) ** 2, 0) : 0;
  const missingFpp = usable ? usable.filter(c => fpp[c.element] === undefined && c.element !== 'H').map(c => c.element) : [];

  const pairs = [];
  const seen = new Set();
  for (const [key, e] of merged) {
    const h = key.split(',').map(Number);
    const key2 = canon(-h[0], -h[1], -h[2]);
    if (key2 === key || seen.has(key2)) continue; // centric, or pair already taken
    seen.add(key);
    const e2 = merged.get(key2);
    if (!e2) continue;
    const d = e.I - e2.I, v = e.v + e2.v, Ib = (e.I + e2.I) / 2;
    const p = { d, v, Ib, dexp2: null };
    if (sumFpp2 > 0 && Ib > 0) {
      const s = sinThetaOverLambda(h[0], h[1], h[2], recip);
      const sumF02 = usable.reduce((a, c) => a + c.count * (scatteringFactor(c.element, s) || 0) ** 2, 0);
      if (sumF02 > 0) p.dexp2 = 4 * sumFpp2 / sumF02;
    }
    pairs.push(p);
  }
  const n = pairs.length;
  if (n < 50) return { pairs: n, enough: false };

  const chi2 = pairs.reduce((a, p) => a + p.d * p.d / p.v, 0) / n;
  const z = (chi2 - 1) * Math.sqrt(n / 2); // excess in units of its noise spread
  const acentricEvidence = chi2 > 1.5 && z > 10 ? 'strong' : chi2 > 1.2 && z > 5 ? 'weak' : 'none';

  // observed / expected, with a jackknife over 10 interleaved blocks
  let ratio = null, ratioSe = null, expectedChi2 = null;
  const withExp = pairs.filter(p => p.dexp2 !== null);
  if (withExp.length >= 50) {
    const est = list => {
      let num = 0, den = 0;
      for (const p of list) { const w = 1 / p.v; num += w * (p.d * p.d - p.v); den += w * p.Ib * p.Ib * p.dexp2; }
      return den > 0 ? Math.sqrt(Math.max(0, num / den)) : null;
    };
    ratio = est(withExp);
    const blocks = 10;
    const jk = [];
    for (let b = 0; b < blocks; b++) jk.push(est(withExp.filter((_, i) => i % blocks !== b)));
    const mean = jk.reduce((a, x) => a + x, 0) / blocks;
    ratioSe = Math.sqrt((blocks - 1) / blocks * jk.reduce((a, x) => a + (x - mean) ** 2, 0));
    expectedChi2 = 1 + withExp.reduce((a, p) => a + p.Ib * p.Ib * p.dexp2 / p.v, 0) / withExp.length;
  }
  let twin = null;
  if (ratio !== null) {
    const r = Math.min(1, ratio);
    const x = (1 - r) / 2;
    let verdict;
    if (expectedChi2 < 1.5) verdict = 'expected signal too weak to judge';
    else if (x < 0.1) verdict = 'enantiopure — no inversion twinning';
    else if (x < 0.35) verdict = 'partial inversion twin';
    else verdict = 'racemic (inversion) twin near 50:50 — or centrosymmetric';
    twin = { ratio, ratioSe, x, xSe: ratioSe / 2, expectedChi2, verdict };
  }
  return { pairs: n, enough: true, chi2, z, acentricEvidence, twin, missingFpp, hasFpp: sumFpp2 > 0, sigmaScale, calibrated };
}
