// abs-structure.js
// Absolute structure from a REFINED model: Parsons quotients, Hooft y and a
// fixed-model Flack x. Only for a final CIF that carries the SHELXL model
// (_shelx_res_file) together with the unmerged reflections (_shelx_hkl_file).
//
// All three compare the observed Bijvoet pairs I(h), I(-h) with |Fc(h)|²,
// |Fc(-h)|² of the model, so — unlike anomalous-check.js — they have a sign:
// x ≈ 0 correct hand, x ≈ 1 inverted, x ≈ 0.5 racemic twin.
//
//  - Parsons (Parsons, Flack & Wagner 2013, Acta Cryst. B69, 249):
//      Qo = (I+ − I−)/(I+ + I−) = (1 − 2x) Qc, weighted LSQ through the
//      origin; the scale cancels. This is the value SHELXL (2013+) writes as
//      _refine_ls_abs_structure_Flack, so the two should agree closely.
//  - Hooft (Hooft, Straver & Spek 2008, 2009, 2010): Bijvoet differences on
//      the Fc² scale, Do = (1 − 2y) Dc, Student-t likelihood with ν fitted to
//      the probability plot, exactly as PLATON does it (platon.f, see
//      hooftStudentT). P2 (y = 0 vs 1) and P3 (y = 0, ½, 1) from the
//      likelihoods.
//  - Flack, classical: I = k[(1 − x)|Fc(h)|² + x|Fc(−h)|²] fitted with the
//      model held fixed (only k and x free) and SHELXL's WGHT weights. A real
//      refinement frees all parameters, so the su here comes out somewhat
//      smaller and x can differ by a few su.
//
// Supported: SHELXL .res with LATT/SYMM, short or long SFAC, DISP, FVAR free
// variables, riding Uiso, anisotropic U, EXTI, HKLF 4 with a matrix, OMIT hkl,
// SHEL, TWIN only as a pure inversion twin. Not: other twin laws, HKLF 5,
// SWAT, multipoles — those return {applicable: false, reason}.

const SHELX_COMMANDS = new Set(('TITL CELL ZERR LATT SYMM SFAC DISP UNIT LAUE REM MORE TIME END HKLF OMIT SHEL BASF TWIN ' +
  'EXTI SWAT HOPE MERG SPEC RESI MOVE ANIS AFIX HFIX FRAG FEND EXYZ EADP EQIV CONN PART BIND FREE DFIX DANG BUMP SAME ' +
  'SADI CHIV FLAT DELU SIMU DEFS ISOR NCSY SUMP L.S. CGLS BLOC DAMP STIR WGHT FVAR BOND CONF MPLA RTAB HTAB LIST ACTA ' +
  'SIZE TEMP WPDB FMAP GRID PLAN MOLE ABIN ANSC ANSR NEUT XNPD REST RIGU PRIG WIGL BEDT TWST SHEL STIR ESEL EGEN ' +
  'LAUE NOTR FLAT SPIN XPRT FIND PSEE PATT TREF INIT DSUL FVAR BUMP ABIN').split(/\s+/));

const _deg = d => d * Math.PI / 180;

function _normEl(s) {
  const m = String(s).replace(/^\$/, '').match(/^[A-Za-z]{1,2}/);
  if (!m) return null;
  return m[0].charAt(0).toUpperCase() + m[0].slice(1).toLowerCase();
}

// "-X+1/2, Y, -Z" -> {R: 3x3 (x' = R x + t), t}
function parseSymmString(s) {
  const parts = String(s).replace(/\s+/g, '').toLowerCase().split(',');
  if (parts.length !== 3) throw new Error(`bad SYMM "${s}"`);
  const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], t = [0, 0, 0];
  parts.forEach((p, i) => {
    const re = /([+-]?)([0-9.]*(?:\/[0-9.]+)?)\*?([xyz]?)/g;
    let m;
    while ((m = re.exec(p)) !== null) {
      if (m[0] === '') { re.lastIndex++; continue; }
      const sign = m[1] === '-' ? -1 : 1;
      let num = 1;
      if (m[2]) { const [a, b] = m[2].split('/'); num = b ? parseFloat(a) / parseFloat(b) : parseFloat(a); }
      if (m[3]) R[i]['xyz'.indexOf(m[3])] += sign * (m[2] ? num : 1);
      else if (m[2]) t[i] += sign * num;
    }
  });
  return { R, t };
}

// SHELX free-variable coding of an atom parameter
function _decode(v, fv) {
  if (Math.abs(v) < 5) return v;
  const m = Math.round(v / 10), p = v - 10 * m;
  if (Math.abs(m) === 1) return p;
  const f = fv[Math.abs(m) - 1];
  if (f === undefined) throw new Error(`free variable ${Math.abs(m)} not defined`);
  return m > 0 ? p * f : p * (f - 1);
}

function _isNum(t) { return /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/.test(t); }

// Model from a SHELXL .res text (as embedded in the CIF).
function parseShelxResModel(text) {
  const lines = joinContinuations(String(text));
  const m = { cell: null, lambda: null, latt: 1, symm: [], sfac: [], disp: {}, fvar: [], rawAtoms: [],
    hklf: null, twin: null, basf: [], exti: 0, swat: false, wght: [0.1, 0, 0, 0, 0, 1 / 3], omitHKL: [], omit2theta: 180, shel: null, notes: [] };
  for (const raw of lines) {
    const line = raw.replace(/!.*$/, '').trim();
    if (!line) continue;
    const tok = line.split(/\s+/);
    const cmd = tok[0].toUpperCase();
    const nums = tok.slice(1).map(Number);
    if (cmd === 'END') break;
    if (cmd === 'HKLF') {
      const v = nums.filter(Number.isFinite);
      m.hklf = { n: v[0] || 4, s: v[1] || 1, matrix: v.length >= 11 ? [v.slice(2, 5), v.slice(5, 8), v.slice(8, 11)] : null };
      break; // atoms and instructions end here; Q-peaks follow in a .res
    }
    if (cmd === 'REM' || cmd === 'TITL') continue;
    switch (cmd) {
      case 'CELL': m.lambda = nums[0]; m.cell = { a: nums[1], b: nums[2], c: nums[3], alpha: nums[4], beta: nums[5], gamma: nums[6] }; continue;
      case 'LATT': m.latt = nums[0]; continue;
      case 'SYMM': m.symm.push(parseSymmString(line.replace(/^\s*SYMM\s+/i, ''))); continue;
      case 'SFAC': {
        const rest = tok.slice(1);
        if (rest.length > 1 && _isNum(rest[1])) { // long form: El a1 b1 a2 b2 a3 b3 a4 b4 c f' f'' mu r wt
          const v = rest.slice(1).map(Number);
          m.sfac.push({ el: _normEl(rest[0]), label: rest[0], coef: { a: [v[0], v[2], v[4], v[6]], b: [v[1], v[3], v[5], v[7]], c: v[8] }, fp: v[9], fpp: v[10] });
        } else for (const e of rest) m.sfac.push({ el: _normEl(e), label: e });
        continue;
      }
      case 'DISP': { const el = _normEl(tok[1]); if (el) m.disp[el] = { fp: Number(tok[2]), fpp: Number(tok[3]) }; continue; }
      case 'FVAR': m.fvar.push(...nums.filter(Number.isFinite)); continue;
      case 'TWIN': { const v = nums.filter(Number.isFinite); m.twin = { matrix: v.length >= 9 ? [v.slice(0, 3), v.slice(3, 6), v.slice(6, 9)] : [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], n: v.length >= 10 ? v[9] : 2 }; continue; }
      case 'BASF': m.basf.push(...nums.filter(Number.isFinite)); continue;
      case 'EXTI': m.exti = Number.isFinite(nums[0]) ? nums[0] : 0; continue;
      case 'SWAT': m.swat = true; continue;
      case 'WGHT': { const d = [0.1, 0, 0, 0, 0, 1 / 3]; nums.forEach((x, i) => { if (Number.isFinite(x) && i < 6) d[i] = x; }); m.wght = d; continue; }
      case 'OMIT': {
        if (nums.length >= 3 && nums.slice(0, 3).every(Number.isInteger)) m.omitHKL.push(nums.slice(0, 3));
        else if (Number.isFinite(nums[1])) m.omit2theta = nums[1];
        continue;
      }
      case 'SHEL': m.shel = { dmax: nums[0] || 99, dmin: nums[1] || 0 }; continue;
    }
    if (SHELX_COMMANDS.has(cmd) || SHELX_COMMANDS.has(cmd.slice(0, 4))) continue;
    // atom: name sfac x y z sof [Uiso | U11 U22 U33 U23 U13 U12]
    if (tok.length >= 6 && /^\d+$/.test(tok[1]) && tok.slice(2).every(_isNum) && /^[A-Za-z]/.test(tok[0])) {
      const v = tok.slice(1).map(Number);
      if (v.length === 5 || v.length === 6 || v.length === 11) m.rawAtoms.push({ name: tok[0], v });
    }
  }
  if (!m.cell) throw new Error('no CELL in the embedded .res');
  return m;
}

function _metric(cell) {
  const ca = Math.cos(_deg(cell.alpha)), cb = Math.cos(_deg(cell.beta)), cg = Math.cos(_deg(cell.gamma));
  return [[cell.a * cell.a, cell.a * cell.b * cg, cell.a * cell.c * cb],
    [cell.a * cell.b * cg, cell.b * cell.b, cell.b * cell.c * ca],
    [cell.a * cell.c * cb, cell.b * cell.c * ca, cell.c * cell.c]];
}

// atoms with decoded parameters; riding Uiso resolved against the previous atom's Ueq
function decodeAtoms(model, recip) {
  const G = _metric(model.cell), as = [recip.aStar, recip.bStar, recip.cStar];
  const fv = model.fvar;
  const atoms = [];
  let lastUeq = null;
  for (const { name, v } of model.rawAtoms) {
    const sfac = v[0];
    if (sfac < 1 || sfac > model.sfac.length) throw new Error(`atom ${name}: SFAC number ${sfac} out of range`);
    const xyz = [_decode(v[1], fv), _decode(v[2], fv), _decode(v[3], fv)];
    const occ = _decode(v[4], fv);
    const a = { name, type: sfac - 1, xyz, occ };
    if (v.length === 11) {
      const [u11, u22, u33, u23, u13, u12] = v.slice(5).map(x => _decode(x, fv));
      const U = [[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]];
      a.U = U;
      let s = 0;
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) s += U[i][j] * as[i] * as[j] * G[i][j];
      lastUeq = s / 3;
    } else {
      let u = v.length >= 6 ? v[5] : 0.05;
      if (u <= -0.5 && u >= -5) {
        if (lastUeq === null) throw new Error(`atom ${name}: riding Uiso without a preceding atom`);
        u = -u * lastUeq;
      } else {
        u = _decode(u, fv);
        if (_normEl(model.sfac[a.type].el) !== 'H') lastUeq = u;
      }
      a.uiso = u;
    }
    atoms.push(a);
  }
  return atoms;
}

// all symmetry operators (x' = R x + t) from LATT and SYMM
function expandShelxOperators(model) {
  const I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  let ops = [{ R: I, t: [0, 0, 0] }, ...model.symm];
  if (model.latt > 0) ops = ops.concat(ops.map(o => ({ R: o.R.map(r => r.map(x => -x)), t: o.t.map(x => -x) })));
  const CENT = { 1: [], 2: [[0.5, 0.5, 0.5]], 3: [[2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]], 4: [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    5: [[0, 0.5, 0.5]], 6: [[0.5, 0, 0.5]], 7: [[0.5, 0.5, 0]] };
  const cv = CENT[Math.abs(model.latt)];
  if (!cv) throw new Error(`unknown LATT ${model.latt}`);
  const all = [...ops];
  for (const c of cv) for (const o of ops) all.push({ R: o.R, t: o.t.map((x, i) => x + c[i]) });
  return all.map(o => ({ R: o.R, t: o.t.map(x => x - Math.floor(x)) }));
}

// scattering data per SFAC type: f0 coefficients, f', f''
function scatteringTypes(model, cifDisp, lambda) {
  const missing = [], fromTable = [];
  const ita = typeof dispersionForWavelength === 'function' ? dispersionForWavelength(lambda) : null;
  const types = model.sfac.map(s => {
    const coef = s.coef || IT92_SCATTERING_FACTORS[s.el];
    if (!coef) throw new Error(`no scattering factor for ${s.label}`);
    let fp = null, fpp = null;
    if (model.disp[s.el]) { fp = model.disp[s.el].fp; fpp = model.disp[s.el].fpp; }
    else if (cifDisp && cifDisp[s.el]) { fp = cifDisp[s.el].fp; fpp = cifDisp[s.el].fpp; }
    else if (s.coef && Number.isFinite(s.fpp)) { fp = s.fp; fpp = s.fpp; }
    if (fpp === null || !Number.isFinite(fpp)) {
      if (s.el === 'H' || s.el === 'D') { fp = 0; fpp = 0; }
      else if (ita && ita.fpp[s.el] !== undefined) { fpp = ita.fpp[s.el]; fp = 0; fromTable.push(s.el); }
      else { missing.push(s.el); fp = 0; fpp = 0; }
    }
    if (!Number.isFinite(fp)) fp = 0;
    return { el: s.el, coef, fp, fpp };
  });
  return { types, missing, fromTable };
}

// f' and f'' per element from the CIF atom-type loop
function extractDispersionPairsFromCIF(text) {
  const lines = String(text).split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].trim() !== 'loop_') continue;
    const tags = [];
    let j = i + 1;
    while (j < lines.length && lines[j].trim().startsWith('_')) { tags.push(lines[j].trim().split(/\s+/)[0].toLowerCase()); j++; }
    const iSym = tags.indexOf('_atom_type_symbol');
    const iRe = tags.indexOf('_atom_type_scat_dispersion_real');
    const iIm = tags.indexOf('_atom_type_scat_dispersion_imag');
    if (iSym < 0 || iIm < 0) continue;
    const tokens = [];
    for (; j < lines.length; j++) {
      const t = lines[j].trim();
      if (t === 'loop_' || t.startsWith('_') || t.startsWith('data_')) break;
      if (lines[j].startsWith(';')) { let k = j + 1; while (k < lines.length && !lines[k].startsWith(';')) k++; tokens.push(''); j = k; continue; }
      const re = /'([^']*)'|"([^"]*)"|(\S+)/g;
      let mm;
      while ((mm = re.exec(t))) tokens.push(mm[1] ?? mm[2] ?? mm[3]);
    }
    const out = {};
    for (let r = 0; r + tags.length <= tokens.length; r += tags.length) {
      const el = _normEl(tokens[r + iSym]);
      const fpp = parseFloat(tokens[r + iIm]);
      const fp = iRe >= 0 ? parseFloat(tokens[r + iRe]) : 0;
      if (el && Number.isFinite(fpp)) out[el] = { fp: Number.isFinite(fp) ? fp : 0, fpp };
    }
    return Object.keys(out).length ? out : null;
  }
  return null;
}

// |F(h)|² and |F(-h)|² with anomalous dispersion
function makeFcCalculator(atoms, ops, types, recip) {
  const as = [recip.aStar, recip.bStar, recip.cStar];
  const cs = [[1, recip.cosGStar, recip.cosBStar], [recip.cosGStar, 1, recip.cosAStar], [recip.cosBStar, recip.cosAStar, 1]];
  const TWO_PI = 2 * Math.PI, TWO_PI2 = 2 * Math.PI * Math.PI;
  const prep = atoms.map(a => ({ ...a, beta: a.U ? a.U.map((r, i) => r.map((u, j) => TWO_PI2 * u * as[i] * as[j])) : null }));
  const nT = types.length;
  return function (h) {
    const s2 = (h[0] * h[0] * as[0] * as[0] + h[1] * h[1] * as[1] * as[1] + h[2] * h[2] * as[2] * as[2]
      + 2 * h[1] * h[2] * as[1] * as[2] * cs[1][2] + 2 * h[0] * h[2] * as[0] * as[2] * cs[0][2] + 2 * h[0] * h[1] * as[0] * as[1] * cs[0][1]) / 4;
    const A = new Float64Array(nT), B = new Float64Array(nT);
    const hp = ops.map(o => [h[0] * o.R[0][0] + h[1] * o.R[1][0] + h[2] * o.R[2][0],
      h[0] * o.R[0][1] + h[1] * o.R[1][1] + h[2] * o.R[2][1],
      h[0] * o.R[0][2] + h[1] * o.R[1][2] + h[2] * o.R[2][2]]);
    const ht = ops.map(o => h[0] * o.t[0] + h[1] * o.t[1] + h[2] * o.t[2]);
    for (const a of prep) {
      let sa = 0, sb = 0;
      const tIso = a.beta ? 1 : Math.exp(-8 * Math.PI * Math.PI * a.uiso * s2);
      for (let k = 0; k < ops.length; k++) {
        const q = hp[k];
        const phi = TWO_PI * (q[0] * a.xyz[0] + q[1] * a.xyz[1] + q[2] * a.xyz[2] + ht[k]);
        let T = tIso;
        if (a.beta) {
          const b = a.beta;
          T = Math.exp(-(q[0] * q[0] * b[0][0] + q[1] * q[1] * b[1][1] + q[2] * q[2] * b[2][2] + 2 * (q[0] * q[1] * b[0][1] + q[0] * q[2] * b[0][2] + q[1] * q[2] * b[1][2])));
        }
        sa += T * Math.cos(phi); sb += T * Math.sin(phi);
      }
      A[a.type] += a.occ * sa; B[a.type] += a.occ * sb;
    }
    let pr = 0, pi = 0, mr = 0, mi = 0;
    for (let t = 0; t < nT; t++) {
      const ty = types[t];
      const c = ty.coef;
      let f0 = c.c;
      for (let i = 0; i < 4; i++) f0 += c.a[i] * Math.exp(-c.b[i] * s2);
      const fr = f0 + ty.fp, fi = ty.fpp, a = A[t], b = B[t];
      pr += fr * a - fi * b; pi += fr * b + fi * a;   // F(h)
      mr += fr * a + fi * b; mi += fi * a - fr * b;   // F(-h): geometric part conjugated, f unchanged
    }
    return { fp2: pr * pr + pi * pi, fm2: mr * mr + mi * mi, s: Math.sqrt(s2) };
  };
}

// Reported values in the CIF (SHELXL / Olex2)
function reportedAbsStructure(cifText) {
  const val = tag => { const raw = extractCIFTagValue(cifText, tag); if (raw === null || raw === undefined) return null; const r = parseCIFNumberWithEsd(raw); return r && Number.isFinite(r.value) ? { value: r.value, su: r.esd, raw: String(raw).trim() } : null; };
  let details = extractCIFTagValue(cifText, '_refine_ls_abs_structure_details');
  if (details === null || details === undefined || /^\s*(;|$)/.test(details)) details = extractCIFTextField(cifText, '_refine_ls_abs_structure_details');
  details = details ? String(details).replace(/\s+/g, ' ').trim() : '';
  const flack = val('_refine_ls_abs_structure_Flack');
  const nQ = details.match(/(\d+)\s+quotients/i);
  return {
    flack,
    flackIsParsons: /quotient|parsons/i.test(details),
    nQuotients: nQ ? Number(nQ[1]) : null,
    details,
    r1: val('_refine_ls_R_factor_gt'),
  };
}

// Long loops give the browser a moment every ~40 ms, so the page stays
// responsive and a spinner can turn; onProgress(text, fraction) reports the step.
const _now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());
let _lastYield = 0;
async function _maybeYield() {
  if (_now() - _lastYield < 40) return;
  await new Promise(r => setTimeout(r, 0));
  _lastYield = _now();
}

// ---- Student-t machinery for Hooft y (follows platon.f) ----
function _lgamma(x) {
  const g = 7, c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
    12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - _lgamma(1 - x);
  x -= 1;
  let a = c[0];
  const t = x + g + 0.5;
  for (let i = 1; i < 9; i++) a += c[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}
function _betacf(a, b, x) {
  const FPMIN = 1e-300;
  let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= 300; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 3e-14) break;
  }
  return h;
}
function _betai(a, b, x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(_lgamma(a + b) - _lgamma(a) - _lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return x < (a + 1) / (a + b + 2) ? bt * _betacf(a, b, x) / a : 1 - bt * _betacf(b, a, 1 - x) / b;
}
function _tCdf(t, nu) {
  const p = 0.5 * _betai(nu / 2, 0.5, nu / (nu + t * t));
  return t > 0 ? 1 - p : p;
}
// inverse CDF by bracketing and bisection to 1e-6 (GEN140)
function _tInv(p, nu) {
  let lo, hi;
  if (p < 0.5) { hi = 0; lo = -1; while (_tCdf(lo, nu) > p) { hi = lo; lo *= 2; } }
  else { lo = 0; hi = 1; while (_tCdf(hi, nu) < p) { lo = hi; hi *= 2; } }
  while (hi - lo > 1e-6) { const m = (hi + lo) / 2; if (_tCdf(m, nu) < p) lo = m; else hi = m; }
  return (hi + lo) / 2;
}
// quantiles at (i − ½)/n (GEN141); symmetric, so half is computed
async function _tQuantiles(n, nu) {
  const w = new Float64Array(n);
  for (let i = 0; i < Math.ceil(n / 2); i++) { const q = _tInv((i + 0.5) / n, nu); w[i] = q; w[n - 1 - i] = -q; if (i % 64 === 0) await _maybeYield(); }
  if (n % 2) w[(n - 1) / 2] = 0;
  return w;
}
// probability-plot correlation and slope through the origin (GEN119, mode 2)
function _ppFit(y, w) {
  const n = y.length;
  const ybar = y.reduce((a, v) => a + v, 0) / n;
  let s1 = 0, s2 = 0, s3 = 0;
  for (let i = 0; i < n; i++) { s1 += (y[i] - ybar) ** 2; s2 += w[i] * y[i]; s3 += w[i] * w[i]; }
  return { cc: s1 !== 0 ? s2 / Math.sqrt(s3 * s1) : 1, slope: s2 / s3 };
}
// ν with the best probability-plot correlation, golden section on [2, 100] (GEN143)
async function _fitNu(y, onStep) {
  const R = 0.61803399, C = 1 - R, TOL = 0.0005;
  let nStep = 0;
  const f = async nu => { if (onStep) onStep(++nStep); return 1 - _ppFit(y, await _tQuantiles(y.length, nu)).cc; };
  let x0 = 2, x3 = 100;
  const bx = 99;
  let x1, x2;
  if (Math.abs(x3 - bx) > Math.abs(bx - x0)) { x1 = bx; x2 = bx + C * (x3 - bx); } else { x2 = bx; x1 = bx - C * (bx - x0); }
  let f1 = await f(x1), f2 = await f(x2);
  while (Math.abs(x3 - x0) > TOL * (Math.abs(x1) + Math.abs(x2))) {
    if (f2 < f1) { x0 = x1; x1 = x2; x2 = R * x1 + C * x3; f1 = f2; f2 = await f(x2); }
    else { x3 = x2; x2 = x1; x1 = R * x2 + C * x0; f2 = f1; f1 = await f(x1); }
  }
  return f1 < f2 ? x1 : x2;
}
// D: [{Do, Dc, sD}] on the Fc² scale. Grid over g = 1 − 2y from −5 to 5 in
// steps of 0.025; pass 1 with ν = 10, then ν and the slope from the Student-t
// probability plot of (g·Dc − Do)/σ (5 outliers off each end), pass 2 with
// that ν and the log-likelihood divided by slope². y and su from the mean and
// spread of the posterior; P2/P3 from the grid at g = 1, 0, −1.
async function hooftStudentT(D, onProgress) {
  const STEP = 0.025, NSP1 = Math.round(1 / STEP), NSTP = 10 * NSP1 + 1, NSPT = 5 * NSP1;
  const grid = async nu => {
    let c = 0;
    const L = new Float64Array(NSTP);
    for (const d of D) {
      if (++c % 100 === 0) await _maybeYield();
      for (let j = 0; j < NSTP; j++) {
        const g = (j - NSPT) * STEP;
        const z = (g * d.Dc - d.Do) / d.sD;
        L[j] -= Math.log(1 + z * z / nu) * (nu + 1) / 2;
      }
    }
    return L;
  };
  const moments = L => {
    let mx = -Infinity; for (const v of L) if (v > mx) mx = v;
    let s0 = 0, s1 = 0;
    for (let j = 0; j < NSTP; j++) { const e = L[j] - mx; if (e > -500) { const p = Math.exp(e); s0 += p; s1 += (j - NSPT) * STEP * p; } }
    const g = s1 / s0;
    let s2 = 0;
    for (let j = 0; j < NSTP; j++) { const e = L[j] - mx; if (e > -500) s2 += ((j - NSPT) * STEP - g) ** 2 * Math.exp(e); }
    return { g, sd: Math.sqrt(s2 / s0), mx };
  };
  if (onProgress) onProgress('Hooft y: first pass (ν = 10)', 0.8);
  const g0 = moments(await grid(10)).g;
  let r = D.map(d => (g0 * d.Dc - d.Do) / Math.max(d.sD, 1e-12)).sort((a, b) => a - b);
  if (r.length > 30) r = r.slice(5, r.length - 5);
  const nu = await _fitNu(r, n => onProgress && onProgress(`Hooft y: fitting ν to the Student-t probability plot (step ${n})`, 0.85));
  const pp = _ppFit(r, await _tQuantiles(r.length, nu));
  if (onProgress) onProgress('Hooft y: second pass', 0.97);
  const L = await grid(nu);
  for (let j = 0; j < NSTP; j++) L[j] /= pp.slope * pp.slope;
  const m = moments(L);
  const e = j => { const v = L[j] - m.mx; return v > -500 ? Math.exp(v) : 0; };
  const pT = e(NSPT + NSP1), pF = e(NSPT - NSP1), pH = e(NSPT);
  const sum3 = pT + pF + pH;
  const y = (1 - m.g) / 2, su = m.sd / 2;
  // two-value test only for y near 0 or 1 (PLA425)
  const p3twin = sum3 > 0 ? pH / sum3 : null;
  const p2ok = !(Math.abs(Math.abs(y - 0.5) - 0.5) > Math.max(0.1, 3 * su) && p3twin > 0.1);
  return { y, su, n: D.length, nu, slope: pp.slope, cc: pp.cc,
    p2true: p2ok && pT + pF > 0 ? pT / (pT + pF) : null,
    p3true: sum3 > 0 ? pT / sum3 : null, p3twin, p3false: sum3 > 0 ? pF / sum3 : null };
}

// Main entry: cifText with _shelx_res_file + reflections [{h,k,l,f2,sigma}] from its _shelx_hkl_file.
async function computeAbsoluteStructure(cifText, reflections, onProgress) {
  _lastYield = _now();
  const resText = extractCIFTextField(cifText, '_shelx_res_file') || extractCIFTextField(cifText, '_iucr_refine_instructions_details');
  if (!resText) return { applicable: false, missingModel: true, reason: 'the CIF has no refined model (_shelx_res_file)' };
  const model = parseShelxResModel(resText);
  const reported = reportedAbsStructure(cifText);
  if (model.latt > 0) return { applicable: false, reason: 'the space group is centrosymmetric (LATT > 0) — no absolute structure', reported };
  if (model.hklf && model.hklf.n === 5) return { applicable: false, reason: 'HKLF 5 (non-merohedral twin) is not supported', reported };
  if (model.swat) return { applicable: false, reason: 'SWAT (solvent model) is not supported', reported };
  let inversionTwin = false;
  if (model.twin) {
    const isInv = model.twin.matrix.every((r, i) => r.every((x, j) => Math.abs(x - (i === j ? -1 : 0)) < 1e-6)) && Math.abs(model.twin.n) === 2;
    if (!isInv) return { applicable: false, reason: 'refined with a twin law other than inversion (TWIN) — not supported yet', reported };
    inversionTwin = true;
  }
  const recip = reciprocalCellFromDirect(model.cell);
  const atoms = decodeAtoms(model, recip);
  if (!atoms.length) return { applicable: false, reason: 'no atoms found in the embedded .res', reported };
  const ops = expandShelxOperators(model);
  const cifDisp = extractDispersionPairsFromCIF(cifText);
  const sc = scatteringTypes(model, cifDisp, model.lambda);
  const fc = makeFcCalculator(atoms, ops, sc.types, recip);

  // point group (rotations of the non-centring operators) for merging
  const pg = [...new Map(ops.map(o => [o.R.flat().join(','), o.R])).values()];

  // reflections in the refinement setting, OMIT/SHEL applied
  const M = model.hklf && model.hklf.matrix;
  const omitKeys = new Set(model.omitHKL.map(h => canonicalHKL(h[0], h[1], h[2], pg).join(',')));
  const dmin = model.shel ? model.shel.dmin : 0, dmax = model.shel ? model.shel.dmax : 1e9;
  const sinMax = Math.sin(_deg(model.omit2theta / 2)) / model.lambda;
  const merged = new Map();
  let nRaw = 0;
  for (const r of reflections) {
    if (!(r.sigma > 0)) continue;
    let h = [r.h, r.k, r.l];
    if (M) h = [M[0][0] * h[0] + M[0][1] * h[1] + M[0][2] * h[2], M[1][0] * h[0] + M[1][1] * h[1] + M[1][2] * h[2], M[2][0] * h[0] + M[2][1] * h[1] + M[2][2] * h[2]].map(Math.round);
    const s = sinThetaOverLambda(h[0], h[1], h[2], recip);
    const d = 1 / (2 * s);
    if (d < dmin || d > dmax || s > sinMax) continue;
    const key = canonicalHKL(h[0], h[1], h[2], pg).join(',');
    if (omitKeys.has(key)) continue;
    let e = merged.get(key);
    if (!e) { e = { h: key.split(',').map(Number), sw: 0, swi: 0, swii: 0, n: 0 }; merged.set(key, e); }
    const w = 1 / (r.sigma * r.sigma);
    e.sw += w; e.swi += w * r.f2; e.swii += w * r.f2 * r.f2; e.n++; nRaw++;
  }
  const exti = model.exti;
  const lam3 = model.lambda ** 3;
  const extCorr = (f2, s) => {
    if (!(exti > 0)) return f2;
    const th = Math.asin(Math.min(1, model.lambda * s));
    return f2 / Math.sqrt(1 + 0.001 * exti * f2 * lam3 / Math.sin(2 * th));
  };
  // Fc² of the model for every merged class; F(-h) belongs to the Friedel class
  const items = [];
  const seen = new Set();
  let nDone = 0;
  if (onProgress) onProgress('Structure factors of the model', 0);
  for (const [key, e] of merged) {
    if (++nDone % 50 === 0) { if (onProgress) onProgress(`Structure factors of the model (${nDone} / ${merged.size})`, 0.75 * nDone / merged.size); await _maybeYield(); }
    // σ of the mean: the larger of the propagated σ and the scatter of the
    // equivalents. 1/√Σw alone shrinks with every repeat and ignores
    // systematic errors (absorption, scaling) — on jg6 (multiplicity 8.5) it
    // gave σ about 3× too small, G ≈ 5, and Parsons x = −0.017 instead of
    // SHELXL's −0.006; with the scatter included: −0.007.
    const I = e.swi / e.sw;
    let sig = 1 / Math.sqrt(e.sw);
    if (e.n >= 2) sig = Math.max(sig, Math.sqrt(Math.max(0, e.swii - e.sw * I * I) / ((e.n - 1) * e.sw)));
    const keyM = canonicalHKL(-e.h[0], -e.h[1], -e.h[2], pg).join(',');
    const centric = keyM === key;
    const f = fc(e.h);
    const fp2 = extCorr(f.fp2, f.s), fm2 = extCorr(f.fm2, f.s);
    items.push({ key, keyM, I, sig, fp2, fm2, centric, s: f.s });
  }
  const byKey = new Map(items.map(it => [it.key, it]));
  const meanFc2 = items.reduce((a, it) => a + it.fp2, 0) / Math.max(1, items.length);
  const pairs = [];
  for (const it of items) {
    if (it.centric || seen.has(it.key)) continue;
    const mate = byKey.get(it.keyM);
    if (!mate) continue;
    seen.add(it.key); seen.add(mate.key);
    if (it.fp2 + it.fm2 < 1e-6 * meanFc2) continue; // systematically absent
    pairs.push({ Ip: it.I, sp: it.sig, Im: mate.I, sm: mate.sig, cp: it.fp2, cm: it.fm2 });
  }
  if (pairs.length < 20) return { applicable: false, reason: `only ${pairs.length} Bijvoet pairs (both mates measured) — Friedel pairs merged in the data?`, reported };

  // --- scale and model check (all reflections, x = 0) ---
  let sxy = 0, syy = 0;
  for (const it of items) { const w = 1 / (it.sig * it.sig); sxy += w * it.I * it.fp2; syy += w * it.fp2 * it.fp2; }
  let k = sxy / syy;

  // --- classical Flack, model fixed: I = a|F(h)|² + c|F(-h)|², SHELXL weights ---
  const [wa, wb, wc, wd, we, wf] = model.wght;
  const acentric = items.filter(it => !it.centric && it.fp2 + it.fm2 > 1e-6 * meanFc2);
  let flack = null;
  if (acentric.length >= 20) {
    let a = k, c = 0;
    for (let iter = 0; iter < 6; iter++) {
      const kk = a + c;
      let s11 = 0, s12 = 0, s22 = 0, b1 = 0, b2 = 0, chi = 0;
      const rows = acentric.map(it => {
        const Fo2 = it.I / kk, sF = it.sig / kk;
        const Fc2 = ((a * it.fp2 + c * it.fm2) / kk);
        const P = wf * Math.max(Fo2, 0) + (1 - wf) * Fc2;
        let w = 1 / (sF * sF + (wa * P) ** 2 + wb * P + wd + we * (it.s * it.s));
        if (wc > 0) w *= Math.exp(wc * it.s * it.s); else if (wc < 0) w *= 1 - Math.exp(wc * it.s * it.s);
        return { w: w / (kk * kk), it };
      });
      for (const { w, it } of rows) { s11 += w * it.fp2 * it.fp2; s12 += w * it.fp2 * it.fm2; s22 += w * it.fm2 * it.fm2; b1 += w * it.I * it.fp2; b2 += w * it.I * it.fm2; }
      const det = s11 * s22 - s12 * s12;
      if (!(Math.abs(det) > 0)) break;
      a = (b1 * s22 - b2 * s12) / det; c = (s11 * b2 - s12 * b1) / det;
      if (iter === 5) {
        for (const { w, it } of rows) chi += w * (it.I - a * it.fp2 - c * it.fm2) ** 2;
        const gof2 = chi / (rows.length - 2);
        const v11 = s22 / det * gof2, v22 = s11 / det * gof2, v12 = -s12 / det * gof2;
        const S = a + c, x = c / S;
        // x = c/(a+c): dx/da = -c/S², dx/dc = a/S²
        const ga = -c / (S * S), gc = a / (S * S);
        const varx = ga * ga * v11 + gc * gc * v22 + 2 * ga * gc * v12;
        flack = { x, su: Math.sqrt(Math.max(0, varx)), n: rows.length, gof: Math.sqrt(gof2) };
      }
    }
    if (flack) k = a + c;
  }
  // model check: R1 on F for I > 2σ, with the fitted inversion fraction (an inverted model must not look like a bad fit)
  const xm = flack ? Math.min(1, Math.max(0, flack.x)) : 0;
  let rNum = 0, rDen = 0;
  for (const it of items) {
    if (it.I < 2 * it.sig) continue;
    const Fo = Math.sqrt(it.I / k), Fc = Math.sqrt((1 - xm) * it.fp2 + xm * it.fm2);
    rNum += Math.abs(Fo - Fc); rDen += Fo;
  }
  const r1 = rDen > 0 ? rNum / rDen : null;

  // --- Parsons quotients ---
  let sQ = 0, sQQ = 0, nQ = 0;
  const qs = [];
  for (const p of pairs) {
    const S = p.Ip + p.Im;
    if (!(S > 0)) continue;
    const Qo = (p.Ip - p.Im) / S;
    const sQo = 2 * Math.sqrt(p.Im * p.Im * p.sp * p.sp + p.Ip * p.Ip * p.sm * p.sm) / (S * S);
    const Qc = (p.cp - p.cm) / (p.cp + p.cm);
    if (!(sQo > 0)) continue;
    const w = 1 / (sQo * sQo);
    sQ += w * Qo * Qc; sQQ += w * Qc * Qc; nQ++;
    qs.push({ w, Qo, Qc });
  }
  let parsons = null;
  if (nQ >= 20 && sQQ > 0) {
    const b = sQ / sQQ;
    const chi = qs.reduce((acc, q) => acc + q.w * (q.Qo - b * q.Qc) ** 2, 0) / (nQ - 1);
    const suRaw = 0.5 / Math.sqrt(sQQ);
    parsons = { x: (1 - b) / 2, su: suRaw * Math.max(1, Math.sqrt(chi)), suRaw, n: nQ, gof: Math.sqrt(chi) };
  }

  // --- Hooft y, as PLATON does it (platon.f PLA421/424/425, GEN139-143) ---
  const D = pairs.map(p => ({ Do: (p.Ip - p.Im) / k, sD: Math.sqrt(p.sp * p.sp + p.sm * p.sm) / k, Dc: p.cp - p.cm }));
  const hooft = D.length >= 20 ? await hooftStudentT(D, onProgress) : null;

  // Bijvoet signal strength: mean |Dc| relative to the noise
  const meanAbsDc = D.reduce((a, d) => a + Math.abs(d.Dc), 0) / D.length;
  const meanSD = D.reduce((a, d) => a + d.sD, 0) / D.length;

  return {
    applicable: true, reported, inversionTwin, basf: model.basf,
    nAtoms: atoms.length, nOps: ops.length, nRaw, nUnique: items.filter(it => it.fp2 + it.fm2 >= 1e-6 * meanFc2).length, nPairs: pairs.length,
    r1, scale: k, flack, parsons, hooft,
    signal: meanAbsDc / meanSD,
    dispersion: { missing: sc.missing, fromTable: sc.fromTable, fromCif: !!cifDisp },
    exti: model.exti, hklfMatrix: !!M,
  };
}

// Verdict from x and its su (Flack & Bernardinelli 2000). x a few su below 0
// (or above 1) is not a twin: it comes from systematic errors.
function absStructureVerdict(x, su) {
  if (!(su > 0)) return { text: 'no su', color: '#888', level: 'none' };
  if (su > 0.3) return { text: 'too imprecise — absolute structure not determined', color: '#888', level: 'weak' };
  const near0 = x < 3 * su, near1 = x > 1 - 3 * su, nearHalf = Math.abs(x - 0.5) < 3 * su;
  if (near0 && !near1 && su <= 0.1) return { text: 'correct hand', color: '#080', level: 'ok' };
  if (near1 && !near0 && su <= 0.1) return { text: 'inverted — invert the model', color: '#b00', level: 'bad' };
  if (nearHalf) return { text: 'racemic (inversion) twin or wrong space group', color: '#b60', level: 'warn' };
  if (near0 && !near1) return { text: 'probably correct hand (su > 0.1 — enantiopure compound needed)', color: '#080', level: 'ok' };
  if (near1 && !near0) return { text: 'probably inverted', color: '#b00', level: 'bad' };
  if (near0 && near1) return { text: 'too imprecise — absolute structure not determined', color: '#888', level: 'weak' };
  return { text: 'partial inversion twin', color: '#b60', level: 'warn' };
}
