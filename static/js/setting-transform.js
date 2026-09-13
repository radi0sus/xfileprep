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
// Method: search a fixed pool of candidate 3x3 integer matrices (signed
// permutations only — swap/negate whole axes, no shear — plus the one
// genuine monoclinic-choice-2 shear) for one that, applied to the matched
// candidate's operators (rotation part conjugated, translation part
// transformed, then optionally re-origin-shifted), exactly reproduces the
// standard entry's operator set as a set. If more than one matrix in the
// pool passes (happens when the space group's own symmetry makes two
// relabelings equivalent — observed for exactly this Iba2 case), the first
// in a fixed, deterministic order is used; any passing matrix is equally
// physically correct, this just makes the choice reproducible.
//
// The pool itself (8 matrices) still comes from PLATON's TRDAT table
// (platon.f ~line 25612) — not because it's trusted blindly anymore, but
// because it's a convenient, ITA-conventional set of *candidates worth
// trying*. Every one of them is now checked against the actual matched
// entries before use, so a wrong guess is caught (returns null) rather
// than silently applied.
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
  [[0, 0, 1], [0, -1, 0], [1, 0, 0]], // monoclinic choice 3 (a-glide) -> choice 1 (c-glide): pure relabeling
  [[-1, 0, -1], [0, -1, 0], [-1, 0, 0]], // monoclinic choice 2 (n-glide) -> choice 1 (c-glide): genuine shear — found by brute-force search over small unimodular matrices and verified against actual operators; an earlier hand-transcription from platon.f ([[1,0,0],[0,-1,0],[-1,0,-1]]) looked plausible but FAILED verification (no origin shift reconciled all 4 operators) — proof this whole derive-and-verify approach earns its keep
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

function seitzKey(op) {
  return op.R.map(r => r.join(',')).join(';') + '|' + op.t.map(mod12).join(',');
}

// Re-origins a transformed operator set by shift v (in twelfths), the same
// way hall-symbol-parser.js's internal shift logic works: t' = t + v - R*v.
function shiftOperators(ops, v) {
  return ops.map(op => {
    const Rv = matVec3(op.R, v);
    return { R: op.R, t: op.t.map((ti, i) => mod12(ti + v[i] - Rv[i])) };
  });
}

// Does applying M to candOps (rotation conjugated: R' = M R M^-1;
// translation: t' = M t) reproduce stdOps exactly, for SOME choice of
// origin? Shifts in quarter-cell steps (0, 1/4, 1/2, 3/4 per axis) are
// tried — sufficient for every case in CANDIDATE_MATRIX_POOL (verified:
// the monoclinic choice-2 shear needs exactly a 1/4-cell shift on two
// axes; the pure permutations need none) — not an exhaustive twelfths
// search, which would be needed for full generality beyond this pool.
function reproducesOperators(M, candOps, stdOps) {
  if (candOps.length !== stdOps.length) return false;
  const Minv = invert3x3(M); // NOT transpose — M isn't orthogonal for the monoclinic shear case
  const transformed = candOps.map(op => ({
    R: matMul3(M, matMul3(op.R, Minv)),
    t: matVec3(M, op.t).map(mod12),
  }));
  const stdKeySet = new Set(stdOps.map(seitzKey));
  // Step by quarters (3 twelfths), not just halves — the monoclinic
  // choice-2 shear needs a 1/4-cell origin shift to line up, found by
  // direct calculation from its translation parts (0/6-only missed it).
  for (let s0 = 0; s0 < 12; s0 += 3)
    for (let s1 = 0; s1 < 12; s1 += 3)
      for (let s2 = 0; s2 < 12; s2 += 3) {
        const keys = shiftOperators(transformed, [s0, s1, s2]).map(seitzKey);
        if (keys.every(k => stdKeySet.has(k))) return true;
      }
  return false;
}

// Tries every matrix in the pool AND its inverse (many pool entries aren't
// self-inverse) against the candidate's actual operators, returns the
// first that reproduces the standard entry's operators exactly, or null
// if the whole pool fails (genuinely unhandled relationship — e.g.
// rhombohedral hex/rhomb axes, tetragonal/trigonal origin choices — the
// caller should leave the setting untransformed rather than guess).
function deriveTransformMatrix(candidateHall, standardHall) {
  const candOps = generateSpaceGroupOperators(candidateHall);
  const stdOps = generateSpaceGroupOperators(standardHall);
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
    if (reproducesOperators(M, candOps, stdOps)) return M;
  }
  return null;
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
// candidate is already standard or no matrix in the pool reproduces the
// standard entry's operators from this candidate's (see deriveTransformMatrix).
//
// "Standard" sibling: same `number`, qualifier === null for orthorhombic-
// style entries; qualifier === "b1" for monoclinic ones (b-unique only —
// a1../c1.. candidates are left untransformed, not implemented).
//
// Return shape: { matrix, standardEntry, isPurePermutation }
// — matrix transforms THIS candidate's cell/hkl into standardEntry's.
function getStandardizationTransform(entry, hallTable) {
  const q = entry.qualifier;

  const monoMatch = typeof q === 'string' && /^([abc])([123])$/.exec(q);
  if (monoMatch) {
    const [, axis, choice] = monoMatch;
    if (axis !== 'b') return null; // a/c-unique cell choice: not implemented
    if (choice === '1') return null; // already standard
    const standardEntry = hallTable.find(e => e.number === entry.number && e.qualifier === 'b1');
    if (!standardEntry) return null;
    const matrix = deriveTransformMatrix(entry.hall, standardEntry.hall);
    if (!matrix) return null;
    return { matrix, standardEntry, isPurePermutation: isPurePermutation(matrix) };
  }

  if (q === null || q === 'abc') return null; // already standard (or unqualified, nothing to do)
  const standardEntry = hallTable.find(e => e.number === entry.number && (e.qualifier == null || e.qualifier === 'abc'));
  if (!standardEntry) return null;
  const matrix = deriveTransformMatrix(entry.hall, standardEntry.hall);
  if (!matrix) return null;
  return { matrix, standardEntry, isPurePermutation: isPurePermutation(matrix) };
}

// Formats a 3x3 integer matrix as the 9 numbers on a SHELX HKLF line, in
// the row-major order the manual specifies (new_h = r11 h + r12 k + r13 l).
function formatHKLFMatrix(M) {
  return M.flat().map(v => (Number.isInteger(v) ? v : v.toFixed(4))).join(' ');
}
