// laue-merge.js
// Determines the true Laue class by merging reflections under each
// candidate point group (from cell-metric.js) and comparing Rint — the
// data-driven step that XPREP itself relies on, rather than trusting cell
// metric alone (a metrically-orthorhombic cell can still be a monoclinic
// structure with beta accidentally close to 90°, for example).
//
// Following the same principle as the space-group symmetry discussion:
// the 11 Laue classes are NOT hardcoded as full element lists. Each is
// generated on the fly from a small set of generator matrices via group
// closure (standard Cayley-graph BFS). Only the generators themselves are
// "data" — a handful of standard crystallographic rotation matrices,
// unavoidable in the same sense a periodic table is unavoidable.
//
// Reflection indices transform under a symmetry operation with the SAME
// integer matrix used for the real-space coordinate operation (applied as
// h' = h * M, row vector times matrix) — this is the standard convention
// used throughout crystallography (International Tables Vol. B) and
// verified here against the well-known hexagonal 3-fold equivalence
// {(h,k,l), (-h-k,h,l), (k,-h-k,l)}.

const IDENTITY = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const INVERSION = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];

// orthogonal-axes two-folds (monoclinic/orthorhombic/tetragonal/cubic)
const TWO_A = [[1, 0, 0], [0, -1, 0], [0, 0, -1]];
const TWO_B = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]];
const TWO_C = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];

// 4-fold along c (tetragonal/cubic): x,y,z -> -y,x,z
const FOUR_C = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];

// hexagonal-axes operations (trigonal/hexagonal): x,y,z -> -y,x-y,z (3-fold)
// and x,y,z -> x-y,x,z (6-fold), both about c
const THREE_HEX = [[0, -1, 0], [1, -1, 0], [0, 0, 1]];
const SIX_HEX = [[1, -1, 0], [1, 0, 0], [0, 0, 1]];
// secondary 2-fold in hexagonal axes, along a: x,y,z -> x-y,-y,-z
const TWO_A_HEX = [[1, -1, 0], [0, -1, 0], [0, 0, -1]];

// cubic 3-fold along the [111] body diagonal: cyclic permutation x,y,z -> z,x,y
const THREE_CUBIC_111 = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];

// Standard 3x3 integer matrix multiplication, A*B.
function matMul(A, B) {
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      let sum = 0;
      for (let k = 0; k < 3; k++) sum += A[i][k] * B[k][j];
      C[i][j] = sum;
    }
  }
  return C;
}

function matKey(M) {
  return M.map(row => row.join(',')).join(';');
}

// Generates the full point group (as a list of matrices) from a small set
// of generators, via breadth-first closure — no group needs to be listed
// element by element by hand, including the 48-element cubic m-3m.
function generateGroup(generators) {
  const elements = new Map([[matKey(IDENTITY), IDENTITY]]);
  let frontier = [IDENTITY];

  while (frontier.length > 0) {
    const next = [];
    for (const g of frontier) {
      for (const gen of generators) {
        const product = matMul(g, gen);
        const key = matKey(product);
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

// Generator sets for each of the 11 Laue classes, in standard settings.
// monoclinicAxis lets the caller pick b-unique (default, IUCr standard) or
// c-unique — every other class is orientation-independent once its
// generators are closed under the group operation.
function getLaueClassGenerators(name, monoclinicAxis = 'b') {
  switch (name) {
    case '-1': return [INVERSION];
    case '2/m': return [monoclinicAxis === 'c' ? TWO_C : TWO_B, INVERSION];
    case 'mmm': return [TWO_A, TWO_B, INVERSION];
    case '4/m': return [FOUR_C, INVERSION];
    case '4/mmm': return [FOUR_C, TWO_A, INVERSION];
    case '-3': return [THREE_HEX, INVERSION];
    case '-3m': return [THREE_HEX, TWO_A_HEX, INVERSION];
    case '6/m': return [SIX_HEX, INVERSION];
    case '6/mmm': return [SIX_HEX, TWO_A_HEX, INVERSION];
    case 'm-3': return [TWO_A, TWO_B, THREE_CUBIC_111, INVERSION];
    case 'm-3m': return [FOUR_C, THREE_CUBIC_111, INVERSION];
    default: throw new Error(`Unknown Laue class: ${name}`);
  }
}

function compareTriples(a, b) {
  for (let i = 0; i < 3; i++) if (a[i] !== b[i]) return a[i] - b[i];
  return 0;
}

// Picks the lexicographically smallest hkl in the symmetry orbit as the
// canonical (merged) index for a reflection.
function canonicalHKL(h, k, l, groupMatrices) {
  let best = null;
  for (const M of groupMatrices) {
    const candidate = [
      h * M[0][0] + k * M[1][0] + l * M[2][0],
      h * M[0][1] + k * M[1][1] + l * M[2][1],
      h * M[0][2] + k * M[1][2] + l * M[2][2],
    ];
    if (best === null || compareTriples(candidate, best) < 0) best = candidate;
  }
  return best;
}

// Groups reflections by their canonical (symmetry-reduced) hkl under the
// given point group — the shared first step for both Rint and an actual
// merged reflection list.
function groupReflectionsByCanonicalHKL(reflections, groupMatrices) {
  const groups = new Map();
  for (const r of reflections) {
    const [ch, ck, cl] = canonicalHKL(r.h, r.k, r.l, groupMatrices);
    const key = `${ch},${ck},${cl}`;
    if (!groups.has(key)) groups.set(key, { h: ch, k: ck, l: cl, members: [] });
    groups.get(key).members.push(r);
  }
  return groups;
}

// Merges reflections under the given point group and computes Rint
// (standard definition: sum|F2 - mean F2| / sum F2, over redundant groups).
function mergeUnderGroup(reflections, groupMatrices) {
  const groups = groupReflectionsByCanonicalHKL(reflections, groupMatrices);

  let rIntNumerator = 0;
  let rIntDenominator = 0;
  for (const { members } of groups.values()) {
    if (members.length < 2) continue; // singletons don't contribute to Rint
    const mean = members.reduce((a, r) => a + r.f2, 0) / members.length;
    for (const r of members) rIntNumerator += Math.abs(r.f2 - mean);
    rIntDenominator += members.reduce((a, r) => a + r.f2, 0);
  }

  return {
    uniqueReflections: groups.size,
    totalReflections: reflections.length,
    redundancy: reflections.length / groups.size,
    rInt: rIntDenominator > 0 ? rIntNumerator / rIntDenominator : null,
  };
}

// Produces the actual merged (symmetry-unique) reflection list — mean F²
// per canonical hkl, with sigma propagated as the standard error of the
// mean — for uses beyond Rint itself. Wilson/E² statistics in particular
// should run on unique data, not raw symmetry-redundant observations,
// since duplicate observations of the same true reflection would
// otherwise be double-counted as if they were independent data.
function mergeReflections(reflections, groupMatrices) {
  const groups = groupReflectionsByCanonicalHKL(reflections, groupMatrices);
  const merged = [];
  for (const { h, k, l, members } of groups.values()) {
    const n = members.length;
    const f2 = members.reduce((a, r) => a + r.f2, 0) / n;
    const sigma = Math.sqrt(members.reduce((a, r) => a + r.sigma * r.sigma, 0)) / n;
    merged.push({ h, k, l, f2, sigma, n });
  }
  return merged;
}

// Evaluates every candidate Laue class (as produced by cell-metric.js) and
// returns merging statistics for each, highest symmetry first — the same
// ranked comparison XPREP shows, left for the user to judge (or for a
// later automatic "jump detection" heuristic to pick from).
function evaluateLaueCandidates(reflections, candidateLaueClasses, options = {}) {
  const monoclinicAxis = options.monoclinicAxis || 'b';
  return candidateLaueClasses.map(name => {
    const generators = getLaueClassGenerators(name, monoclinicAxis);
    const group = generateGroup(generators);
    const stats = mergeUnderGroup(reflections, group);
    return { laueClass: name, groupOrder: group.length, ...stats };
  });
}

// Automatically picks the Laue class from the merging statistics, instead
// of leaving that judgment entirely to the person reading the table.
// Heuristic: starting from the LOWEST symmetry candidate (which merges the
// fewest reflections together and is therefore always a "safe" baseline),
// walk toward higher symmetry as long as Rint doesn't jump — a genuine
// jump means the next candidate's extra symmetry isn't actually present
// in the data. jumpFactor is a judgment call, not a physical law: 1.5
// (50% relative increase tolerated) is a reasonable default, but this is
// exactly the kind of borderline case worth cross-checking against the
// full table rather than trusting blindly.
function determineLaueClass(reflections, candidateLaueClasses, options = {}) {
  const jumpFactor = options.jumpFactor || 1.5;
  const evaluations = evaluateLaueCandidates(reflections, candidateLaueClasses, options);

  // If NONE of the candidates found even one redundant pair to merge,
  // Rint is undefined everywhere and this method has no signal at all to
  // distinguish them — typically because the input is already a final,
  // merged/symmetry-unique reflection list (common for deposited CIF
  // "_refln" loops), not raw multiply-measured data. Walking "up" from
  // the lowest candidate would then silently default to it, which looks
  // exactly like a confident answer but isn't one. Flag this instead of
  // hiding it, so callers can warn the user rather than trust a default
  // that has no real evidence behind it.
  const redundancyInformative = evaluations.some(e => e.rInt !== null);

  // evaluations is ordered highest symmetry first (matches
  // classifyCellMetric's laueClasses order) — walk it in reverse (lowest
  // symmetry first).
  const lowToHigh = evaluations.slice().reverse();
  let chosen = lowToHigh[0];
  let jumpDetectedAt = null;

  if (redundancyInformative) {
    for (let i = 1; i < lowToHigh.length; i++) {
      const candidate = lowToHigh[i];
      if (candidate.rInt === null || chosen.rInt === null) break; // insufficient data to compare
      const ratio = candidate.rInt / chosen.rInt;
      if (ratio > jumpFactor) {
        jumpDetectedAt = candidate.laueClass;
        break;
      }
      chosen = candidate;
    }
  }

  return { chosen: chosen.laueClass, jumpDetectedAt, evaluations, redundancyInformative };
}

