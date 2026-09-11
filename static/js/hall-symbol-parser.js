// hall-symbol-parser.js
// Decodes a Hall symbol (Hall, 1981; International Tables Vol. B, Section
// 1.4) into the full set of space-group symmetry operators (Seitz
// matrices: rotation + translation).
//
// As discussed earlier in this project: NOTHING about space-group
// symmetry is hardcoded as a per-space-group table here. The only "data"
// is the same small set of standard crystallographic generator matrices
// used throughout this project (rotation matrices for principal, face-
// diagonal, and body-diagonal axes; translation fractions; lattice
// centering vectors) — all 230 space groups (530 settings) are produced
// by parsing their Hall symbol and closing the resulting generators under
// group multiplication, the same technique used for the Laue classes in
// laue-merge.js.
//
// Translations are tracked as integers in units of 1/12 (not floats) —
// every standard crystallographic translation fraction (1/2, 1/3, 1/4,
// 1/6, 1/12) divides evenly into twelfths, so this avoids floating-point
// rounding entirely during group closure and deduplication.

function mod12(x) {
  return ((x % 12) + 12) % 12;
}

// --- Table 3: rotation matrices for principal axes ----------------------

const PRINCIPAL_ROTATION_MATRICES = {
  x: {
    1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    2: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
    3: [[1, 0, 0], [0, 0, -1], [0, 1, -1]],
    4: [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    6: [[1, 0, 0], [0, 1, -1], [0, 1, 0]],
  },
  y: {
    1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    2: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
    3: [[-1, 0, 1], [0, 1, 0], [-1, 0, 0]],
    4: [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
    6: [[0, 0, 1], [0, 1, 0], [-1, 0, 1]],
  },
  z: {
    1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    2: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
    3: [[0, -1, 0], [1, -1, 0], [0, 0, 1]],
    4: [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    6: [[1, -1, 0], [1, 0, 0], [0, 0, 1]],
  },
};

// --- Table 4: face-diagonal 2-fold matrices (direction depends on the
// axis of the immediately preceding generator) ---------------------------

const FACE_DIAGONAL_MATRICES = {
  x: {
    "'": [[-1, 0, 0], [0, 0, -1], [0, -1, 0]], // b-c
    '"': [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],   // b+c
  },
  y: {
    "'": [[0, 0, -1], [0, -1, 0], [-1, 0, 0]], // a-c
    '"': [[0, 0, 1], [0, -1, 0], [1, 0, 0]],   // a+c
  },
  z: {
    "'": [[0, -1, 0], [-1, 0, 0], [0, 0, -1]], // a-b
    '"': [[0, 1, 0], [1, 0, 0], [0, 0, -1]],   // a+b
  },
};
const FACE_DIAGONAL_DIRECTION = {
  x: { "'": [0, 1, -1], '"': [0, 1, 1] },
  y: { "'": [1, 0, -1], '"': [1, 0, 1] },
  z: { "'": [1, -1, 0], '"': [1, 1, 0] },
};

// --- Table 5: body-diagonal 3-fold (a+b+c) ------------------------------

const BODY_DIAGONAL_MATRIX = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];

// --- Table 2: translation symbols, in twelfths --------------------------

const FIXED_TRANSLATIONS_12 = {
  a: [6, 0, 0], b: [0, 6, 0], c: [0, 0, 6], n: [6, 6, 6],
  u: [3, 0, 0], v: [0, 3, 0], w: [0, 0, 3], d: [3, 3, 3],
};

// --- Table 1: lattice translation vectors, in twelfths ------------------

const LATTICE_TRANSLATIONS_12 = {
  P: [[0, 0, 0]],
  A: [[0, 0, 0], [0, 6, 6]],
  B: [[0, 0, 0], [6, 0, 6]],
  C: [[0, 0, 0], [6, 6, 0]],
  I: [[0, 0, 0], [6, 6, 6]],
  R: [[0, 0, 0], [8, 4, 4], [4, 8, 8]],
  S: [[0, 0, 0], [4, 4, 8], [8, 8, 4]],
  T: [[0, 0, 0], [4, 8, 4], [8, 4, 8]],
  F: [[0, 0, 0], [0, 6, 6], [6, 0, 6], [6, 6, 0]],
};

// --- Matrix/vector helpers (all translations in twelfths) ---------------

function matMul3(A, B) {
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
  return C;
}
function matVecMul3(M, v) {
  return [0, 1, 2].map(i => M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2]);
}
function negateMat(M) {
  return M.map(row => row.map(x => -x));
}

// --- Parsing a single NAT token ------------------------------------------

function parseTranslationChars(chars, order, axisVector) {
  let t = [0, 0, 0];
  for (const ch of chars) {
    if (FIXED_TRANSLATIONS_12[ch]) {
      t = t.map((v, i) => v + FIXED_TRANSLATIONS_12[ch][i]);
    } else if (/[0-9]/.test(ch)) {
      const frac12 = Math.round((parseInt(ch, 10) / order) * 12);
      t = t.map((v, i) => v + axisVector[i] * frac12);
    }
  }
  return t.map(mod12);
}

// position: 0-based index of this generator among the NAT tokens (used
// for the default-axis rules). previousOrder/previousAxis describe the
// immediately preceding generator, needed to resolve face-diagonal
// direction and the "second rotation" default axis rule.
function parseNatToken(token, position, previousOrder, previousAxis) {
  let i = 0;
  let sign = 1;
  if (token[0] === '-') { sign = -1; i = 1; }
  const order = parseInt(token[i], 10);
  i++;

  let axisChar = null;
  if (i < token.length && /[xyz'"*]/.test(token[i])) {
    axisChar = token[i];
    i++;
  }
  const translationChars = token.slice(i);

  let axisType = axisChar;
  if (!axisType) {
    if (position === 0) axisType = 'z';
    else if (position === 1) {
      if (previousOrder === 2 || previousOrder === 4) axisType = 'x';
      else if (previousOrder === 3 || previousOrder === 6) axisType = "'";
      else axisType = 'z';
    } else {
      axisType = '*';
    }
  }

  let R;
  let axisVector;
  let resolvedPrincipalAxis = null;
  if (axisType === 'x' || axisType === 'y' || axisType === 'z') {
    R = PRINCIPAL_ROTATION_MATRICES[axisType][order];
    axisVector = axisType === 'x' ? [1, 0, 0] : axisType === 'y' ? [0, 1, 0] : [0, 0, 1];
    resolvedPrincipalAxis = axisType;
  } else if (axisType === "'" || axisType === '"') {
    const base = previousAxis || 'z';
    R = FACE_DIAGONAL_MATRICES[base][axisType];
    axisVector = FACE_DIAGONAL_DIRECTION[base][axisType];
  } else if (axisType === '*') {
    R = BODY_DIAGONAL_MATRIX;
    axisVector = [1, 1, 1];
  } else {
    throw new Error(`Unresolvable axis for NAT token "${token}"`);
  }

  if (sign === -1) R = negateMat(R);
  const t = parseTranslationChars(translationChars, order, axisVector);

  return { seitz: { R, t }, order, resolvedPrincipalAxis };
}

// Shifts a Seitz operation's origin by v (in twelfths): S' = T(v) S T(-v),
// giving t' = t + v - R*v.
function shiftOperation(op, v12) {
  const Rv = matVecMul3(op.R, v12);
  const t = op.t.map((ti, i) => mod12(ti + v12[i] - Rv[i]));
  return { R: op.R, t };
}

// Parses a full Hall symbol string into its lattice translations and
// generator Seitz operations (already origin-shifted if a V vector was
// given at the end of the symbol).
function parseHallSymbol(hallSymbol) {
  const rawTokens = hallSymbol.trim().split(/\s+/);

  let originShift = [0, 0, 0];
  let tokens = rawTokens;
  const openIdx = rawTokens.findIndex(t => t.startsWith('('));
  if (openIdx !== -1) {
    const shiftStr = rawTokens.slice(openIdx).join(' ').replace(/[()]/g, '');
    originShift = shiftStr.split(/\s+/).map(Number);
    tokens = rawTokens.slice(0, openIdx);
  }

  let latticeToken = tokens[0];
  let centrosymmetric = false;
  if (latticeToken.startsWith('-')) {
    centrosymmetric = true;
    latticeToken = latticeToken.slice(1);
  }
  const latticeTranslations = LATTICE_TRANSLATIONS_12[latticeToken.toUpperCase()];
  if (!latticeTranslations) throw new Error(`Unknown lattice symbol: ${latticeToken}`);

  const natTokens = tokens.slice(1);
  const generators = [];
  let previousOrder = null;
  let previousAxis = null;

  natTokens.forEach((token, position) => {
    const parsed = parseNatToken(token, position, previousOrder, previousAxis);
    generators.push(parsed.seitz);
    previousOrder = parsed.order;
    if (parsed.resolvedPrincipalAxis) previousAxis = parsed.resolvedPrincipalAxis;
  });

  if (centrosymmetric) {
    generators.push({ R: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], t: [0, 0, 0] });
  }

  const shiftedGenerators = generators.map(g => shiftOperation(g, originShift));

  return { centrosymmetric, latticeTranslations, generators: shiftedGenerators };
}

// --- Group closure (same technique as laue-merge.js, extended to
// include translations) ---------------------------------------------------

function composeOps(a, b) {
  // apply a first, then b
  const R = matMul3(b.R, a.R);
  const Rta = matVecMul3(b.R, a.t);
  const t = Rta.map((v, i) => mod12(v + b.t[i]));
  return { R, t };
}
function seitzKey(op) {
  return op.R.map(r => r.join(',')).join(';') + '|' + op.t.join(',');
}

function generatePointGroupOperators(generators) {
  const identity = { R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] };
  const elements = new Map([[seitzKey(identity), identity]]);
  let frontier = [identity];

  while (frontier.length > 0) {
    const next = [];
    for (const g of frontier) {
      for (const gen of generators) {
        const product = composeOps(g, gen);
        const key = seitzKey(product);
        if (!elements.has(key)) {
          elements.set(key, product);
          next.push(product);
        }
      }
    }
    frontier = next;
  }
  return Array.from(elements.values());
}

// Full pipeline: Hall symbol string -> all general-position Seitz
// operators of the space group (point-group operations x lattice
// centering translations).
function generateSpaceGroupOperators(hallSymbol) {
  const { latticeTranslations, generators } = parseHallSymbol(hallSymbol);
  const pointOps = generatePointGroupOperators(generators);

  const allOps = new Map();
  for (const pOp of pointOps) {
    for (const Lt of latticeTranslations) {
      const t = pOp.t.map((ti, i) => mod12(ti + Lt[i]));
      const op = { R: pOp.R, t };
      allOps.set(seitzKey(op), op);
    }
  }
  return Array.from(allOps.values());
}

// Formats a Seitz operation back to x,y,z notation, mainly for debugging
// and for later SYMM-card generation.
function seitzToXYZ(op) {
  const vars = ['x', 'y', 'z'];
  return [0, 1, 2].map(row => {
    let parts = [];
    for (let col = 0; col < 3; col++) {
      const coeff = op.R[row][col];
      if (coeff === 1) parts.push(`+${vars[col]}`);
      else if (coeff === -1) parts.push(`-${vars[col]}`);
    }
    let s = parts.join('').replace(/^\+/, '');
    const frac = op.t[row];
    if (frac !== 0) {
      s += (frac > 0 ? '+' : '-') + `${Math.abs(frac)}/12`;
    }
    return s;
  }).join(',');
}
