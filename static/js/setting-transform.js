// setting-transform.js
// Standardizes a matched candidate's cell/axis setting to the conventional
// one (e.g. measured I 2 c b -> standard I b a 2).
//
// v2 — DATA-DRIVEN, not table-driven. The first version of this module
// picked the transformation matrix by trusting hall-symbols-data.js's
// `qualifier` string (e.g. "cab") as if it directly named the matrix to
// apply. Checked against a real dataset (space group 45, qualifier "cab",
// I 2 c b), that assumption turned out to be wrong: the qualifier "cab"
// entry is NOT related to the qualifier-null entry by the matrix literally
// named "cab" — a genuine, silent, hard-to-notice error of exactly the
// kind this whole feature was supposed to guard against.
//
// So: don't trust names. DERIVE the matrix from the two entries' actual
// generated symmetry operators (the same operators generateSpaceGroupOperators
// already produces for every candidate, i.e. data the tool has already
// validated by scoring it against the observed reflections) and VERIFY it
// reproduces the standard entry's operator set exactly before using it.
//
// Method: find 3x3 integer matrices M that, applied to the matched
// candidate's operators, exactly reproduce the target entry's operator set
// (as a set, for SOME origin shift). M is always the CELL-AXIS matrix in
// SHELX convention (rows = new axes in terms of old ones, a' = M a), the
// same matrix applied to h,k,l and written on the HKLF line.
//
// v3 — CONVENTION FIX (monoclinic cell choices). Fractional coordinates do
// NOT transform with M but with Q = (M^-1)^T (x' = Q x). v2 conjugated the
// operators with M itself. For signed permutations (orthogonal, Q = M)
// that is harmless, which is why all orthorhombic cases and the monoclinic
// choice 3 -> 1 relabeling were right. For a genuine shear it is not:
// v2's "verified" P 21/n -> P 21/c matrix [[-1,0,-1],[0,-1,0],[-1,0,0]]
// actually re-indexes P 21/n data into P 21/a, and the I -> C matrix kept
// the wrong centering. PLATON's matrix [[1,0,0],[0,-1,0],[-1,0,-1]], which
// v2 rejected, is the correct one (checked against the reflections
// themselves: allowed in the source group <=> allowed in the target group
// after h' = M h).
//
// Monoclinic groups (Nos. 3-15): instead of a fixed pool, ALL integer
// matrices with entries in {-1,0,1} and det = +1 are searched (covers cell
// choices 1/2/3, A/C/I centering and a/b/c-unique axis changes). Several
// valid matrices exist for any pair (a' can be any lattice vector in the
// ac plane that completes the basis); with a cell given, the one giving
// the least oblique cell is chosen (shortest axes, then monoclinic angle
// closest to 90, obtuse preferred) — the IUCr/ITA convention, and what
// XPREP/PLATON do. Without a cell, the sparsest matrix wins.
//
// Other groups: the signed-permutation pool below (from PLATON's TRDAT
// table, platon.f ~line 25612) is tried; every candidate is still checked
// against the real operators before use, a failure returns null.
//
// SHELX convention (SHELX manual, HKLF instruction): new_h = r11*h +
// r12*k + r13*l — the same matrix that transforms the cell edges is
// applied, unmodified, directly to h,k,l.

const CANDIDATE_MATRIX_POOL = [
  [[1, 0, 0], [0, 1, 0], [0, 0, 1]], // identity ("abc")
  [[0, 1, 0], [1, 0, 0], [0, 0, -1]], // "ba-c"
  [[0, 1, 0], [0, 0, 1], [1, 0, 0]], // "cab"
  [[0, 0, 1], [0, 1, 0], [-1, 0, 0]], // "-cba"
  [[0, 0, 1], [1, 0, 0], [0, 1, 0]], // "bca"
  [[1, 0, 0], [0, 0, -1], [0, 1, 0]], // "a-cb"
];

function determinant3(M) {
  return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
       - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
       + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
}

// General 3x3 inverse via the adjugate/cofactor method. Every matrix this
// module works with has determinant +-1, so the result is exact integers —
// round defensively and throw if it turns out not to be, rather than
// silently returning a fractional "transformation matrix".
function invert3x3(M) {
  const det = determinant3(M);
  if (det === 0) throw new Error('setting-transform: matrix is singular, cannot invert');
  const cof = (r, c) => {
    const rows = [0, 1, 2].filter(i => i !== r);
    const cols = [0, 1, 2].filter(j => j !== c);
    const sign = (r + c) % 2 === 0 ? 1 : -1;
    return sign * (M[rows[0]][cols[0]] * M[rows[1]][cols[1]] - M[rows[0]][cols[1]] * M[rows[1]][cols[0]]);
  };
  const adjT = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let r = 0; r < 3; r++)
    for (let c = 0; c < 3; c++)
      adjT[c][r] = cof(r, c); // transpose of the cofactor matrix = adjugate
  const inv = adjT.map(row => row.map(v => v / det));
  const rounded = inv.map(row => row.map(v => Math.round(v)));
  const maxDrift = Math.max(...rounded.map((row, i) => Math.max(...row.map((v, j) => Math.abs(v - inv[i][j])))));
  if (maxDrift > 1e-6) throw new Error('setting-transform: inverse is not integer-valued — unexpected input matrix');
  return rounded;
}

function matMul3(A, B) {
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
  return C;
}

function matVec3(M, v) {
  return [0, 1, 2].map(i => M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2]);
}

function transpose3(M) {
  return M.map((_, i) => M.map(row => row[i]));
}

function mod12(x) {
  return ((x % 12) + 12) % 12;
}

// True when M only ever relabels/negates whole axes (no mixing) — i.e. a
// signed permutation matrix: exactly one nonzero (+-1) entry per row and
// per column. Distinguishes a harmless axis-relabeling from a genuine
// shear that actually changes the cell's shape — used to decide whether
// the original measured cell/esds are worth keeping around in a REM line.
function isPurePermutation(M) {
  for (let i = 0; i < 3; i++) {
    const row = M[i].filter(v => v !== 0);
    if (row.length !== 1 || Math.abs(row[0]) !== 1) return false;
  }
  for (let j = 0; j < 3; j++) {
    const col = [0, 1, 2].map(i => M[i][j]).filter(v => v !== 0);
    if (col.length !== 1 || Math.abs(col[0]) !== 1) return false;
  }
  return true;
}

// --- The actual derivation: search + verify against real operators -------

function rotKey(R) {
  return R.map(r => r.join(',')).join(';');
}

// Fractional-coordinate matrix belonging to cell-axis matrix M: x' = Q x.
function coordinateMatrix(M) {
  return transpose3(invert3x3(M));
}

// Does the basis change M (cell-axis convention, see header) turn candOps
// into exactly stdOps, for SOME origin shift? Operators are conjugated
// with the coordinate matrix Q: R' = Q R Q^-1, t' = Q t. The rotation
// parts are compared first (cheap, rejects almost every wrong M); only
// then is the origin shift searched.
//
// Shift grid: twelfths, and with `fineShifts` also 24ths. Origin choices 1
// and 2 of F d -3 m, F d d d, I 41/a m d, ... differ by 1/8 (e.g. -1/8,-1/8,-1/8
// for F d -3 m), which no multiple of 1/12 hits — these were reported as
// "matrix not found" although M is simply the identity (same hkl indexing).
// Translations are handled in 24ths internally (exact: both operator sets
// are in twelfths). Monoclinic searches don't need it and skip the fine grid
// (it multiplies the shift loop by 8).
function reproducesOperators(M, candOps, stdOps, fineShifts = false) {
  if (candOps.length !== stdOps.length) return false;
  const Q = coordinateMatrix(M);
  const Qinv = invert3x3(Q);
  const m24 = x => ((x % 24) + 24) % 24;
  const key = t => (m24(t[0]) * 24 + m24(t[1])) * 24 + m24(t[2]);
  const stdByRot = new Map();
  for (const op of stdOps) {
    const k = rotKey(op.R);
    if (!stdByRot.has(k)) stdByRot.set(k, new Set());
    stdByRot.get(k).add(key(op.t.map(x => 2 * x)));
  }
  const transformed = [];
  for (const op of candOps) {
    const R = matMul3(Q, matMul3(op.R, Qinv));
    const targets = stdByRot.get(rotKey(R));
    if (!targets) return false;
    transformed.push({ R, t: matVec3(Q, op.t).map(x => 2 * x), targets });
  }
  // identity-rotation ops (pure lattice translations) don't depend on the
  // shift — check them once up front.
  for (const op of transformed) {
    if (isIdentity3(op.R) && !op.targets.has(key(op.t))) return false;
  }
  const moving = transformed.filter(op => !isIdentity3(op.R));
  const step = fineShifts ? 1 : 2; // in 24ths
  for (let s0 = 0; s0 < 24; s0 += step)
    for (let s1 = 0; s1 < 24; s1 += step)
      for (let s2 = 0; s2 < 24; s2 += step) {
        const v = [s0, s1, s2];
        let ok = true;
        for (const op of moving) {
          const R = op.R;
          const t = [0, 1, 2].map(i => op.t[i] + v[i] - (R[i][0] * v[0] + R[i][1] * v[1] + R[i][2] * v[2]));
          if (!op.targets.has(key(t))) { ok = false; break; }
        }
        if (ok) return true;
      }
  return false;
}

function isIdentity3(R) {
  return R[0][0] === 1 && R[1][1] === 1 && R[2][2] === 1
    && R[0][1] === 0 && R[0][2] === 0 && R[1][0] === 0 && R[1][2] === 0 && R[2][0] === 0 && R[2][1] === 0;
}

// All integer 3x3 matrices with entries in {-1,0,1} and determinant +1,
// sparsest first (deterministic order). Built once, lazily.
let SMALL_UNIMODULAR_CACHE = null;
function smallUnimodularMatrices() {
  if (SMALL_UNIMODULAR_CACHE) return SMALL_UNIMODULAR_CACHE;
  const out = [];
  for (let n = 0; n < 19683; n++) {
    let x = n;
    const e = [];
    for (let i = 0; i < 9; i++) { e.push((x % 3) - 1); x = Math.floor(x / 3); }
    const M = [e.slice(0, 3), e.slice(3, 6), e.slice(6, 9)];
    if (determinant3(M) === 1) out.push(M);
  }
  const nonzero = M => M.flat().filter(v => v !== 0).length;
  out.sort((A, B) => nonzero(A) - nonzero(B));
  SMALL_UNIMODULAR_CACHE = out;
  return out;
}

function isMonoclinicNumber(n) {
  return n >= 3 && n <= 15;
}

// Monoclinic angle of a cell = the one furthest from 90 degrees.
function monoclinicAngle(cell) {
  return [cell.alpha, cell.beta, cell.gamma].reduce((w, x) => (Math.abs(x - 90) > Math.abs(w - 90) ? x : w));
}

function cellIsUsable(cell) {
  return !!cell && ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].every(k => Number.isFinite(cell[k]) && cell[k] > 0);
}

// Ranks valid matrices for a given measured cell: least oblique cell
// first (sum of axis lengths, 0.01 A tolerance; then |angle - 90|; then
// obtuse angle preferred, the ITA convention for monoclinic cells), the
// sparsest matrix as final tie-break (input order is already sparsest
// first and Array.prototype.sort is stable).
function pickLeastObliqueMatrix(matrices, cell) {
  if (!cellIsUsable(cell)) return matrices[0];
  const scored = matrices.map(M => {
    const c = transformCellParameters(cell, M);
    const ang = monoclinicAngle(c);
    return { M, len: Math.round((c.a + c.b + c.c) * 100), dev: Math.round(Math.abs(ang - 90) * 100), acute: ang < 90 - 1e-6 ? 1 : 0 };
  });
  scored.sort((x, y) => x.len - y.len || x.dev - y.dev || x.acute - y.acute);
  return scored[0].M;
}

const DERIVE_CACHE = new Map();

// Returns a matrix M (cell-axis / HKLF convention) that transforms the
// candidate setting's cell/hkl into the standard (target) setting's, or
// null if none is found (genuinely unhandled relationship — e.g.
// rhombohedral hex/rhomb axes, tetragonal/cubic origin choices — the
// caller should leave the setting untransformed rather than guess).
// `cell` (optional, as measured) picks the least oblique of several valid
// monoclinic matrices; numbers are only needed to decide whether the
// monoclinic search applies and are derived from the operators otherwise.
function deriveTransformMatrix(candidateHall, standardHall, cell) {
  const cacheKey = candidateHall + '|' + standardHall + '|'
    + (cellIsUsable(cell) ? ['a', 'b', 'c', 'alpha', 'beta', 'gamma'].map(k => cell[k].toFixed(3)).join(',') : '');
  if (DERIVE_CACHE.has(cacheKey)) return DERIVE_CACHE.get(cacheKey);

  const candOps = generateSpaceGroupOperators(candidateHall);
  const stdOps = generateSpaceGroupOperators(standardHall);
  const entry = typeof HALL_SYMBOLS_TABLE !== 'undefined' ? HALL_SYMBOLS_TABLE.find(e => e.hall === standardHall) : null;
  const monoclinic = entry ? isMonoclinicNumber(entry.number) : isMonoclinicOperatorSet(stdOps);

  let result = null;
  if (monoclinic) {
    const valid = smallUnimodularMatrices().filter(M => reproducesOperators(M, candOps, stdOps));
    if (valid.length) result = pickLeastObliqueMatrix(valid, cell);
  } else {
    const pool = [];
    for (const M of CANDIDATE_MATRIX_POOL) {
      pool.push(M);
      pool.push(invert3x3(M));
    }
    const seen = new Set();
    for (const M of pool) {
      const key = M.flat().join(',');
      if (seen.has(key)) continue;
      seen.add(key);
      if (determinant3(M) <= 0) continue; // SHELX requires a positive-determinant HKLF matrix
      if (reproducesOperators(M, candOps, stdOps, true)) { result = M; break; }
    }
  }
  DERIVE_CACHE.set(cacheKey, result);
  return result;
}

// Fallback when the Hall table isn't loaded: a monoclinic point group has
// 2 or 4 distinct rotation parts, one 2-fold/mirror axis.
function isMonoclinicOperatorSet(ops) {
  const rots = new Set(ops.map(op => rotKey(op.R)));
  return rots.size === 2 || rots.size === 4 ? [...rots].some(k => k !== '1,0,0;0,1,0;0,0,1' && k !== '-1,0,0;0,-1,0;0,0,-1') : false;
}

// Transforms cell parameters {a,b,c,alpha,beta,gamma} (degrees) through
// the metric tensor: G' = M * G * M^T, then extracts lengths/angles back
// out of G'. Correct for any integer basis-change matrix, permutation or
// genuine shear alike.
function transformCellParameters(cell, M) {
  const toRad = d => d * Math.PI / 180;
  const toDeg = r => r * 180 / Math.PI;
  const { a, b, c, alpha, beta, gamma } = cell;
  const [ca, cb, cg] = [Math.cos(toRad(alpha)), Math.cos(toRad(beta)), Math.cos(toRad(gamma))];
  const G = [
    [a * a, a * b * cg, a * c * cb],
    [a * b * cg, b * b, b * c * ca],
    [a * c * cb, b * c * ca, c * c],
  ];
  const Gp = matMul3(matMul3(M, G), transpose3(M));
  const na = Math.sqrt(Gp[0][0]);
  const nb = Math.sqrt(Gp[1][1]);
  const nc = Math.sqrt(Gp[2][2]);
  const clamp = v => Math.max(-1, Math.min(1, v));
  return {
    a: na, b: nb, c: nc,
    alpha: toDeg(Math.acos(clamp(Gp[1][2] / (nb * nc)))),
    beta: toDeg(Math.acos(clamp(Gp[0][2] / (na * nc)))),
    gamma: toDeg(Math.acos(clamp(Gp[0][1] / (na * nb)))),
  };
}

// Relabels esds under a PURE PERMUTATION only (see isPurePermutation) —
// exact relabeling/sign-drop of the given length/angle esds. For a
// genuine shear, rigorous esd propagation needs the full cell covariance
// matrix, not available here — callers must not use this for that case.
function permuteEsds(esd, M) {
  const axisSource = i => M[i].findIndex(v => v !== 0);
  const src = [0, 1, 2].map(axisSource);
  const lenEsd = [esd.a, esd.b, esd.c];
  const angEsd = [esd.alpha, esd.beta, esd.gamma];
  const newLen = src.map(j => lenEsd[j]);
  const newAng = [0, 1, 2].map(i => {
    const others = [0, 1, 2].filter(k => k !== i).map(k => src[k]);
    const missingOld = [0, 1, 2].find(k => !others.includes(k));
    return angEsd[missingOld];
  });
  return { a: newLen[0], b: newLen[1], c: newLen[2], alpha: newAng[0], beta: newAng[1], gamma: newAng[2] };
}

// Given a matched candidate entry ({number, qualifier, hm, hall}) and the
// full Hall table, returns the standardization transform, or null if the
// candidate is already standard or no matrix reproduces the standard
// entry's operators from this candidate's (see deriveTransformMatrix).
//
// "Standard" sibling: same `number`, qualifier null/"abc" for orthorhombic-
// style entries; for monoclinic ones b-unique cell choice 1 ("b1", or "b"
// for groups without cell choices) — from any a/b/c-unique, any cell choice,
// including the "-b1".. style entries of Nos. 9 and 15.
//
// Return shape: { matrix, standardEntry, isPurePermutation }
// — matrix transforms THIS candidate's cell/hkl into standardEntry's.
function getStandardizationTransform(entry, hallTable, cell) {
  const q = entry.qualifier;

  if (isMonoclinicNumber(entry.number)) {
    const standardEntry = hallTable.find(e => e.number === entry.number && (e.qualifier === 'b1' || e.qualifier === 'b'));
    if (!standardEntry || standardEntry.hall === entry.hall) return null;
    const matrix = deriveTransformMatrix(entry.hall, standardEntry.hall, cell);
    if (!matrix) return null;
    return { matrix, standardEntry, isPurePermutation: isPurePermutation(matrix) };
  }

  if (q === null || q === 'abc') return null; // already standard (or unqualified, nothing to do)
  const standardEntry = hallTable.find(e => e.number === entry.number && (e.qualifier == null || e.qualifier === 'abc'));
  if (!standardEntry) return null;
  const matrix = deriveTransformMatrix(entry.hall, standardEntry.hall, cell);
  if (!matrix) return null;
  return { matrix, standardEntry, isPurePermutation: isPurePermutation(matrix) };
}

// Formats a 3x3 integer matrix as the 9 numbers on a SHELX HKLF line, in
// the row-major order the manual specifies (new_h = r11 h + r12 k + r13 l).
function formatHKLFMatrix(M) {
  return M.flat().map(v => (Number.isInteger(v) ? v : v.toFixed(4))).join(' ');
}
