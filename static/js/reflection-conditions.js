// reflection-conditions.js
// Derives systematic reflection conditions directly from a space group's
// symmetry operators (as produced by hall-symbol-parser.js) — the
// theoretical counterpart to the OBSERVED conditions from absences.js.
// Matching the two against each other is how a specific space group gets
// picked out of the candidates compatible with the observed Laue class
// and centering.
//
// Core relation (standard crystallographic result): for a symmetry
// operation (R, t), a reflection h is only affected if h·R = h (the
// operation leaves this reflection's indices unchanged, as opposed to
// mapping it to a different — even if equivalent — reflection). When
// that holds, structure factor invariance requires exp(2*pi*i * h·t) = 1,
// i.e. h·t must be an integer; otherwise the reflection is systematically
// absent. Translations are in twelfths (see hall-symbol-parser.js), so
// this reduces to an integer congruence mod 12 — no floating point.

function gcd(a, b) {
  a = Math.abs(a); b = Math.abs(b);
  while (b) { [a, b] = [b, a % b]; }
  return a;
}
function lcm(a, b) {
  return a / gcd(a, b) * b;
}

// Checks whether operator R leaves the reflection class defined by
// zeroIdx (indices among 0=h,1=k,2=l fixed at zero) exactly invariant —
// not merely mapped to a Friedel-equivalent or symmetry-equivalent
// reflection, which wouldn't impose a phase restriction on THIS index
// combination.
function operatorPreservesClass(R, zeroIdx) {
  const freeIdx = [0, 1, 2].filter(i => !zeroIdx.includes(i));
  for (const i of freeIdx) {
    for (const j of zeroIdx) {
      if (R[i][j] !== 0) return false;
    }
    for (const j of freeIdx) {
      const expected = i === j ? 1 : 0;
      if (R[i][j] !== expected) return false;
    }
  }
  return true;
}

// Collects, for a given class, the class-preserving operators' translation
// components at the free index positions (in twelfths) — the raw material
// for the congruences below.
function collectClassTranslations(operators, zeroIdx) {
  const freeIdx = [0, 1, 2].filter(i => !zeroIdx.includes(i));
  const coeffs = [];
  for (const op of operators) {
    if (operatorPreservesClass(op.R, zeroIdx)) {
      coeffs.push(freeIdx.map(i => op.t[i]));
    }
  }
  return { freeIdx, coeffs };
}

// Single free index (axial classes h00/0k0/00l, or the general/centering
// case reduces trivially): combines all operators' constraints into one
// minimal period. Returns 1 if there's no restriction at all.
function derivePeriodCondition(operators, zeroIdx) {
  const { coeffs } = collectClassTranslations(operators, zeroIdx);
  let period = 1;
  for (const [c] of coeffs) {
    if (c === 0) continue;
    const g = gcd(c, 12);
    period = lcm(period, 12 / g);
  }
  return period;
}

// Two free indices (zonal classes 0kl/h0l/hk0): reports the constraint on
// each individual free index and on their sum, wherever the congruences
// reduce to that simple form. More complex (e.g. d-glide-like) congruences
// are reported as a raw pair rather than force-fit into a named condition.
function deriveZonalConditions(operators, zeroIdx) {
  const { coeffs } = collectClassTranslations(operators, zeroIdx);

  let periodP = 1, periodQ = 1, periodSum = 1;
  const rawCongruences = [];

  for (const [cp, cq] of coeffs) {
    if (cp === 0 && cq === 0) continue;
    if (cq === 0) {
      periodP = lcm(periodP, 12 / gcd(cp, 12));
    } else if (cp === 0) {
      periodQ = lcm(periodQ, 12 / gcd(cq, 12));
    } else if (cp === cq) {
      periodSum = lcm(periodSum, 12 / gcd(cp, 12));
    } else {
      rawCongruences.push({ cp, cq, modulus: 12 / gcd(gcd(cp, cq), 12) });
    }
  }
  return { periodP, periodQ, periodSum, rawCongruences };
}

// General form of the zonal conditions, used for matching: the set of index
// pairs (p, q) — taken mod 12, the grid on which every translation lives —
// that the space group ALLOWS in a zone. A zone is a plane of reflections
// spanned by two parameter vectors (e.g. 0kl: p → (0,1,0), q → (0,0,1); the
// diagonal plane hhl: p → (1,1,0), q → (0,0,1)). An operator constrains the
// zone if it leaves both vectors unchanged (v·R = v); the constraint is then
// c_p·p + c_q·q ≡ 0 (mod 12) with c = v·t. The allowed pairs satisfy every
// constraint — this also covers conditions that combine several operators
// (e.g. a and b glide together, d-glide plus centering), which the one-by-one
// period bookkeeping above cannot express.
function deriveAllowedResidues(operators, paramVectors) {
  const preserves = (R, v) => [0, 1, 2].every(j => v[0] * R[0][j] + v[1] * R[1][j] + v[2] * R[2][j] === v[j]);
  const constraints = [];
  for (const op of operators) {
    if (!paramVectors.every(v => preserves(op.R, v))) continue;
    const [cp, cq] = paramVectors.map(v => v[0] * op.t[0] + v[1] * op.t[1] + v[2] * op.t[2]);
    if (cp % 12 !== 0 || cq % 12 !== 0) constraints.push([cp, cq]);
  }
  const allowed = new Set();
  for (let p = 0; p < 12; p++) {
    for (let q = 0; q < 12; q++) {
      if (constraints.every(([cp, cq]) => (((cp * p + cq * q) % 12) + 12) % 12 === 0)) allowed.add(p + ',' + q);
    }
  }
  return allowed;
}

// Does every allowed pair satisfy  a·p + b·q ≡ 0 (mod m)?  (m divides 12, so
// working mod 12 is exact.) True = the group predicts this condition to hold.
function allowedSatisfies(allowed, [a, b, m]) {
  for (const key of allowed) {
    const [p, q] = key.split(',').map(Number);
    if ((((a * p + b * q) % m) + m) % m !== 0) return false;
  }
  return true;
}

// Named centering types for matching the lattice translation vectors
// straight from the parsed Hall symbol — no need to re-derive this from
// the operator set, since the Hall symbol already states it explicitly.
const CENTERING_NAMES_12 = {
  'P (primitive)': [[0, 0, 0]],
  A: [[0, 0, 0], [0, 6, 6]],
  B: [[0, 0, 0], [6, 0, 6]],
  C: [[0, 0, 0], [6, 6, 0]],
  I: [[0, 0, 0], [6, 6, 6]],
  F: [[0, 0, 0], [0, 6, 6], [6, 0, 6], [6, 6, 0]],
  'R (obverse)': [[0, 0, 0], [8, 4, 4], [4, 8, 8]],
};

function latticeTranslationsToName(translations) {
  const key = t => t.join(',');
  const sortedSet = translations.map(key).sort().join(';');
  for (const [name, ref] of Object.entries(CENTERING_NAMES_12)) {
    if (ref.map(key).sort().join(';') === sortedSet) return name;
  }
  return 'unrecognized';
}

// Full derivation, mirroring the shape of scanAllAbsences() in absences.js
// so the two can be compared directly. Takes the Hall symbol string
// directly (rather than just an operator list) since centering is read
// off the lattice symbol, not re-derived statistically.
function deriveAllConditions(hallSymbol) {
  const { latticeTranslations } = parseHallSymbol(hallSymbol);
  const operators = generateSpaceGroupOperators(hallSymbol);

  return {
    centering: latticeTranslationsToName(latticeTranslations),
    axial: {
      h00: derivePeriodCondition(operators, [1, 2]),
      '0k0': derivePeriodCondition(operators, [0, 2]),
      '00l': derivePeriodCondition(operators, [0, 1]),
    },
    zonal: {
      '0kl': deriveZonalConditions(operators, [0]),
      h0l: deriveZonalConditions(operators, [1]),
      hk0: deriveZonalConditions(operators, [2]),
    },
    // zone → set of allowed (p, q) mod 12; compare with absences.js ZONAL_CONDITIONS[..].test
    zonalAllowed: {
      '0kl': deriveAllowedResidues(operators, [[0, 1, 0], [0, 0, 1]]),
      h0l: deriveAllowedResidues(operators, [[1, 0, 0], [0, 0, 1]]),
      hk0: deriveAllowedResidues(operators, [[1, 0, 0], [0, 1, 0]]),
      hhl: deriveAllowedResidues(operators, [[1, 1, 0], [0, 0, 1]]),
    },
  };
}
